# -*- coding: utf-8 -*-
"""
Created on Wed May 20 10:42:26 2020

@author: Grim.M
"""



import sys
import os
import re
from pathlib import Path
import glob

sys.path.insert(0, os.path.dirname(__file__) + "/../cellular_cycle")
sys.path.insert(0, os.path.dirname(__file__) + "/../fit")

#os.environ["PYOPENCL_CTX"] ="0"
from Cell_Population_Simulation import Cell_Population_Simulation
from retrieve_cellState_dw_2 import retrieve_cellState
from call_cell_events import call_events #
import numpy as np
import pandas as pd
from pysb.simulator import BngSimulator
from pysb.simulator import ScipyOdeSimulator
import fracprolif as fp
from rbf_mmd import mmd_rbf
from SMs_G1s_IMTs_extraction import SMs_G1s_IMTs_extraction
from scipy.stats import multivariate_normal

###############################################################################
# Model run check

# tmp_t=np.linspace(0, 10, 100)
# sim = BngSimulator(model, tspan=tmp_t, verbose=True)
# traj = sim.run()
#
# sim = ScipyOdeSimulator(model, tspan=tmp_t, verbose=True)
# traj = sim.run()

###############################################################################

# omega=200

###############################################################################
# Importing data

DATADIR = os.path.dirname(__file__) + "/../../Data/Exp_Cell_Counts"
#DATADIR = os.getcwd() + "/Documents/git_repos/ECM_cell_cycle_model/Data/Exp_Cell_Counts"
#DATADIR ='C:\\Users\\Grim.M' + "/Documents/git_repos/ECM_cell_cycle_model/Data/Exp_Cell_Counts"
filepaths = os.listdir(DATADIR)
filenames = [f for f in filepaths if re.match(r'^DIP_DS8', f)]
filepaths = [Path(DATADIR, f) for f in filenames]
filepaths.sort()
filepaths

all_response = []
for f in filepaths:
    all_response.append(np.loadtxt(f,delimiter=','))

means=[]
for vals in range(len(all_response)):
    means.append(np.mean(all_response[vals],axis=1))
    
stdevs=[]
for vals in range(len(all_response)):
    stdevs.append(np.std(all_response[vals],axis=1))
  
t = np.arange(0,420/3,1/3)


IMT_DATADIR = Path(os.path.dirname(__file__) + "/../../Data/Exp_Intermitotic_Times")
#IMT_DATADIR = Path(path + "/Data/Exp_Intermitotic_Times")

print(IMT_DATADIR)
file_location = os.path.join(IMT_DATADIR,'*imt_*.csv')
filenames = glob.glob(file_location)
filenames.sort()
#print(filenames)


IMTs = []
for f in filenames:
    tmp_df = pd.read_csv(f)
    tmp_df = tmp_df.replace(to_replace=-1, value=np.nan)
    IMTs.append(tmp_df)

ii=0
for df in IMTs:
    df["IMT"] = df["Division"] - df["Mitotic Frame"]
    df["G1"] = df["Green On"] - df["Mitotic Frame"]
    df["S-M"] = df["Division"] - df["Green On"]
    if ii == 2 or ii == 3 or ii == 6 or ii == 7:
        df["IMT_t"] = df["IMT"] * (32/60)
        df["G1_t"] = df["G1"] * (32/60)
        df["S-M_t"] = df["S-M"] * (32/60)
    else:
        df["IMT_t"] = df["IMT"] * (20/60)
        df["G1_t"] = df["G1"] * (20/60)
        df["S-M_t"] = df["S-M"] * (20/60)
    #df = df[(df["IMT_t"] ==27)].reset_index(drop=True)
    ii+=1
    
stiffs = [0.2, 0.2, 0.6, 0.6, 4, 4, 55, 55]
ii=0
for df in IMTs:
    tmp_stiff = np.ones(len(df))*stiffs[ii]
    df.insert(0, "Stiffness", tmp_stiff)
    if ii % 2 == 0:
        df["Drug"] = "DMSO"
    else:
        df["Drug"] = "erlotinib"    
    ii+=1

all_IMTs = pd.concat(IMTs, ignore_index=True, sort=False)
centers = np.arange(0.5, 100.5, 1)  # histogram bin centers 


all_IMTs = all_IMTs[all_IMTs['IMT_t'].notna()]

all_IMTs = all_IMTs[all_IMTs['S-M_t']>0]
all_IMTs = all_IMTs[all_IMTs['G1_t']>0]

##############################################################################



##############################################################################

def RBS(x, x0, sigma = 5) :
    observations , dimensions= np.shape(x0)
    s=0
    for i in range(observations):
        s+=np.exp(-np.linalg.norm(x-x0[i])**2/(2*sigma**2)) / (sigma*np.sqrt(2*np.pi))
    if observations != 0:
        return s/observations
    else:
        return s

##############################################################################

def log_prior(parameter_vector, mean, variance):
    # Create a multivariate normal distribution object
    mvn = multivariate_normal(mean, np.diag(variance))

    # Calculate the log likelihood for each data point
    log_likelihoods = mvn.logpdf(parameter_vector)
    
    # Calculate the sum of log likelihoods
    log_likelihood_sum = np.sum(log_likelihoods)
    
    return log_likelihood_sum


    

##############################################################################

def log_likelihood_single_condition(init, run_time, inputs, condition, sampled_parameter_names ,chain, iteration, save=True):
    
    stiffness, _, treatment = condition.split(maxsplit=2)  # split the string into two parts
    stiffness = float(stiffness)  # convert the stiffness from string to float
    treatment = treatment.strip()  # remove any leading or trailing white spaces
        
    t_sim = np.arange(0,run_time,1/3)
    observed_cells_time = run_time

    parameters_sim = {}
    cost = {}
    
    # FAK, ERK and AKT inputs

    cond_split = condition.split(' ')
    cond_split_DMSO = condition.split(' ')
    cond_split_DMSO[-1] = 'DMSO'

    
    ###########################################################################
    print(condition)
    (cells, cell_count) = Cell_Population_Simulation(init, run_time, run_time, model, inputs, condition, sampled_parameter_names, return_plot=True);

    ###########################################################################
    print('Cell count')
    print(cell_count)
    
    if type(cell_count)==np.ndarray and save:
        np.save(os.path.dirname(__file__) + "/../fit_data/pso/cell_count/cell_count_particle_" + str(chain) + "_iteration_" + str(iteration) + "_" + '_'.join(cond_split) + ".npy", cell_count)
    
    IMT_eoe_counts = retrieve_cellState(run_time, init, cells, return_all_prot=False, cells_starting_at_mitosis=False) #cell count
    all_events_simulated = call_events(run_time, IMT_eoe_counts);
    all_events_simulated = all_events_simulated[all_events_simulated['event'] != 'Starting cell']
    all_events_simulated = all_events_simulated[np.array(all_events_simulated['Birthtime']<observed_cells_time)]
    mitotic_lifespans = all_events_simulated['Lifespan'][all_events_simulated['event']=='mitosis']
    censored_lifespans = all_events_simulated['Lifespan'][all_events_simulated['event']=='End_of_Expt']
    if save:
        all_events_simulated.to_pickle(os.path.dirname(__file__) + "/../fit_data/pso/distribution/Events_particle_" + str(chain) + "_iteration_" + str(iteration) + "_" + '_'.join(cond_split) + ".pkl") #doesn't save all imts, only the cells that are born before 30 hours


    ###########################################################################
    #Extract useful information from the simulations:

    simulated_IMTs, simulated_SMs, simulated_G1s = SMs_G1s_IMTs_extraction(cells) #G1, SM, IMT

    if save:
        np.save(os.path.dirname(__file__) + "/../fit_data/pso/distribution/G1_distribution_particle_"  + str(chain) + "_iteration_" + str(iteration) + "_" + '_'.join(cond_split) + ".npy", np.array(simulated_G1s))
        np.save(os.path.dirname(__file__) + "/../fit_data/pso/distribution/SM_distribution_particle_" + str(chain) + "_iteration_" + str(iteration) + "_" + '_'.join(cond_split) + ".npy", np.array(simulated_SMs))
        np.save(os.path.dirname(__file__) + "/../fit_data/pso/distribution/IMT_distribution_particle_" + str(chain) + "_iteration_" + str(iteration) + "_" + '_'.join(cond_split) + ".npy", np.array(simulated_IMTs))
        #save density estimation
        # x=np.linspace(0, run_time, 360)
        # y=[RBS(i, np.transpose([mitotic_lifespans.tolist()])) for i in x]
       #  np.save(os.path.dirname(__file__) + "/../fit_data/pso/distribution/IMT_distribution_particle_" + str(chain) + "_iteration_" + str(iteration) + '_' + '_'.join(cond_split) + "_kde.npy", np.array(y))



    ###########################################################################
    #compute frac prolif parameters of the simulation:
    # apoptosis rate
    if len(all_events_simulated) > 0:
        number_of_deaths = len(all_events_simulated) - (all_events_simulated['event']=='apoptosis').value_counts()[0]
        cell_time_total = sum(all_events_simulated['Lifespan'])        
        parameters_sim['a'] = number_of_deaths/cell_time_total
    else:
        parameters_sim['a'] = np.nan

    # mitosis and quiescence rates
    if len(mitotic_lifespans) > 0:
        m = fp.mle(mitotic_lifespans, censored_lifespans)            
        parameters_sim['mu'] = m['x'][0]
        parameters_sim['sigma'] = m['x'][1]
        parameters_sim['lambda'] = m['x'][2]
        parameters_sim.update(fp.q_rates(m))
    else:
        parameters_sim['d'] = np.nan
        parameters_sim['q'] = np.nan
        parameters_sim['mu'] = np.nan
        parameters_sim['sigma'] = np.nan
        parameters_sim['lambda'] = np.nan

    ###########################################################################      
    # Compare distributions

    experimental_G1s = all_IMTs[(all_IMTs['Drug']==treatment) & (all_IMTs['Stiffness']==stiffness)]['G1_t']
    experimental_SMs = all_IMTs[(all_IMTs['Drug']==treatment) & (all_IMTs['Stiffness']==stiffness)]['S-M_t']

    if len(experimental_G1s) > 0 and len(experimental_SMs) > 0:
        cost['cost G1 distributions'] = mmd_rbf(np.array([simulated_G1s]).transpose(), np.array([experimental_G1s]).transpose(), 1) #difference between IMT distribution
        cost['cost S-M distributions'] = mmd_rbf(np.array([simulated_SMs]).transpose(), np.array([experimental_SMs]).transpose(), 1) #difference between IMT distribution
    
    else:
        cost['cost G1 distributions'] = np.inf
        cost['cost S-M distributions'] = np.inf
    
    ###########################################################################
    # Cost computation: rmse dip rates

    # if cell_count is not None:
    if type(cell_count)==np.ndarray:
        cell_count_t = np.zeros(len(t_sim))
        time_tmp_m = 0
        
        for time_tmp_p, cell_count_time_tmp in cell_count:
            cell_count_t[(time_tmp_m <= t_sim) & (t_sim < time_tmp_p)] = cell_count_time_tmp
            time_tmp_m = time_tmp_p
        
        num_condition = ['0.6 kPa DMSO', '0.6 kPa erlotinib', '4 kPa DMSO', '4 kPa erlotinib', '55 kPa DMSO', '55 kPa erlotinib'].index(condition)
        
        DIP_sim = np.log2(cell_count_t/cell_count_t[0])
        cost['cost population dynamics'] = np.sqrt(np.sum((DIP_sim[1:len(t_sim)] - means[num_condition][1:len(t_sim)])**2 / stdevs[num_condition][1:len(t_sim)]**2))/len(t_sim)
    
    else:
        cost['cost population dynamics'] = np.inf
        
    return parameters_sim, cost
    
    
##############################################################################



def log_likelihood(parameter_vector, chain, sampled_parameter_names, inputs, conditions = ['0.6 kPa DMSO', '0.6 kPa erlotinib', '4 kPa DMSO', '4 kPa erlotinib', '55 kPa DMSO', '55 kPa erlotinib'], run_times = [120, 120, 120, 120, 60, 60], init_cells = [4,6,3,3,2,2], save=True, iteration=0):
    
    print("The parameters are ")
    print(parameter_vector)
    # np.save(os.path.dirname(__file__) + "/../fit_data/pso/particle_position/parameter_vector_particle_" + str(chain) + "_iteration_" + str(iteration) + ".npy", parameter_vector)
    
    global model, t, means
    
    print(iteration)
    
    param_dict = {pname: pvalue for pname, pvalue in zip(sampled_parameter_names, parameter_vector)}

    # Change model parameter values to current location in parameter space
    for pname, pvalue in param_dict.items():
        model.parameters[pname].value = (pvalue)
        
    parameters_sim={} 
    cost = {}
###############################################################################
###############################################################################
    
    for init, run_time, condition in zip(init_cells, run_times, conditions):
        print(condition)
        parameters_sim_tmp, cost_tmp = log_likelihood_single_condition(init, run_time, inputs, condition, sampled_parameter_names, iteration, save)
        
        parameters_sim[condition] = parameters_sim_tmp
        cost[condition] = cost_tmp
        
        
    # Save cost summary
    parameters_sim = {param + ' ' + cond_tmp: [val_tmp] for cond_tmp in parameters_sim.keys() for param, val_tmp in parameters_sim[cond_tmp].items()}
    cost = {param + ' ' + cond_tmp: [val_tmp] for cond_tmp in cost.keys() for param, val_tmp in cost[cond_tmp].items()}
    
    frac_prolif_parameters_temp = {"chain":[chain], "iteration":[iteration], **parameters_sim, **cost}

    vect = pd.DataFrame(np.array([parameter_vector]), columns = sampled_parameter_names)
    frac_prolif_parameters_temp = pd.DataFrame(data=frac_prolif_parameters_temp)   
    frac_prolif_parameters_temp = pd.concat([frac_prolif_parameters_temp, vect], axis=1)
    
    
    frac_prolif_parameters_temp.columns = pd.MultiIndex.from_tuples([("pso", "chain"), ("pso", "iteration")] + [("fractional proliferation model", param_name) for param_name in parameters_sim.keys()] + [("cost", cost_tmp_name) for cost_tmp_name in cost.keys()] + [("model parameters", parameter_name) for parameter_name in sampled_parameter_names])
    
          
    
    
    return frac_prolif_parameters_temp['cost'].sum(axis=1).values[0], frac_prolif_parameters_temp


def log_posterior(parameter_vector, chain, sampled_parameter_names, mean, variance, inputs, conditions = ['0.6 kPa DMSO', '0.6 kPa erlotinib', '4 kPa DMSO', '4 kPa erlotinib', '55 kPa DMSO', '55 kPa erlotinib'], run_times = [120, 120, 120, 120, 60, 60], init_cells = [4,6,3,3,2,2], save=False, iteration = 0):
    
    prior = - log_prior(parameter_vector, mean, variance)
    ll, frac_prolif_parameters_temp = log_likelihood(parameter_vector, chain, sampled_parameter_names, inputs, conditions, run_times, init_cells, save, iteration)
    frac_prolif_parameters_temp[('cost', '- log prior')] = prior
    if save:
        frac_prolif_parameters_temp.to_pickle(os.path.dirname(__file__) + "/../fit_data/pso/cost/parameters_particle_" + str(chain) + "_iteration_" + str(iteration) + "_pso.pkl") 
    return prior + ll



###############################################################################
# Test the cost function:

if __name__=="__main__":
    # import model
    from cell_cycle_death_model import model
    import pandas as pd
    import matplotlib.pyplot as plt

    # import pysb sampled parameters
    sampled_parameter_names = np.loadtxt(os.path.dirname(__file__) + '/../model_data/sampled_parameter_names.txt',
                                         dtype=str)

    # Importing the parameters:
    new_param = pd.read_csv(os.path.dirname(__file__) + '/../model_data/parameters.csv')
    del new_param['Stress']
    for key in new_param.keys():
        try:
            if model.parameters[key].value != new_param[key][0]:
                model.parameters[key].value = new_param[key]
        except:
            print(key)
    model.parameters['omega'].value = 200

    # Setting initial concentrations '/../model_data/initial_concentrations.csv'
    initial_concentrations = pd.read_csv(os.path.dirname(__file__) + '/../model_data/initial_concentrations.csv')
    for key in initial_concentrations.keys():
        try:
            model.parameters[key].value = initial_concentrations[key][0] / model.parameters['omega'].value
        except:
            print(key)

    # Import ECM_conditions.csv
    ECM_conditions = pd.read_csv(os.path.dirname(__file__) + '/../model_data/ECM_conditions_2.csv')


    original_parameters = []
    for name in sampled_parameter_names:
        original_parameters += [np.log10(model.parameters[name].value)]
    parameter_vector = 10 ** np.array(original_parameters)
    idx = 'original_parameters'
    mean = 10 ** np.array(original_parameters)

    ECM_conditions = pd.read_csv(os.path.dirname(__file__) + '/../model_data/ECM_conditions_2.csv')

    # Calculate prior
    prior = log_prior(parameter_vector, mean, np.diag(np.ones(len(sampled_parameter_names))))

    # # Calculate likelihood
    ll, frac_prolif_parameters = log_likelihood(parameter_vector, 0, sampled_parameter_names, ECM_conditions, save=False)
    #
    # print("The prior is " + str(prior))
    # print("The likelihood is " + str(ll))
    # print("The posterior is " + str(prior + ll))
