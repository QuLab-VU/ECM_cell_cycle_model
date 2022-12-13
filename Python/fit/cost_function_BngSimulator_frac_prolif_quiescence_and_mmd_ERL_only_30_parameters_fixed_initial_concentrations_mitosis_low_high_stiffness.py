# -*- coding: utf-8 -*-
"""
Created on Wed May 20 10:42:26 2020

@author: Grim.M
"""



import sys
import os

sys.path.insert(0, os.path.dirname(__file__) + "/../cellular_cycle")
sys.path.insert(0, os.path.dirname(__file__) + "/../fit")
#sys.path.insert(0, os.getcwd() + "/cellular_cycle")


from Austin_proliferation_model_concentration_5 import model
model_concentration=model
from Austin_proliferation_model_population_5 import model
model_population=model
#os.environ["PYOPENCL_CTX"] ="0"


#from runAKTERK_2_starting_at_mitosis import runAKTERK
from runAKTERK_2 import runAKTERK
from retrieve_cellState_dw_2 import retrieve_cellState
from call_cell_events import call_events #
import numpy as np
import pandas as pd
from pysb.simulator import BngSimulator
from pysb.simulator import ScipyOdeSimulator
import fracprolif as fp
from rbf_mmd import nmmd2_rbf
import random
from pathlib import Path
from SMs_G1s_IMTs_extraction import SMs_G1s_IMTs_extraction


tmp_t=np.linspace(0, 10, 100)
sim = BngSimulator(model_population, tspan=tmp_t, verbose=True)
traj = sim.run()

sim = ScipyOdeSimulator(model_concentration, tspan=tmp_t, verbose=True)
traj = sim.run()

###############################################################################

omega=200


pysb_sampled_parameter_names = ["V_s2skp","V_Fp","k_M","k_EMY","k_CDF","k_CDS","k_ED","K_S","K_MY","K_DF","d_CD","k_P1", "k_P2", "K_CD", "k_E", "K_M", "K_CE", "V_2E2F", "K_2E2F", "k_RE", "kf1", "kf2", "kb", "kasbaxmbcl2", "kdsbaxmbcl2", "k1sbh3", "k2sbh3", "kdbh3", "kasbh3bcl2", "kdsbh3bcl2"]

###############################################################################
import glob

IMT_DATADIR = Path(os.path.dirname(__file__) + "/../../Data/Exp_Intermitotic_Times")

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
        df["Drug"] = "ERL"    
    ii+=1

all_IMTs = pd.concat(IMTs, ignore_index=True, sort=False)
centers = np.arange(0.5, 100.5, 1)  # histogram bin centers 


all_IMTs = all_IMTs[all_IMTs['IMT_t'].notna()]


##############################################################################
def RBS(x,x0,sigma) :
    observations , dimensions= np.shape(x0)
    s=0
    for i in range(observations):
        s+=np.exp(-np.linalg.norm(x-x0[i])**2/(2*sigma**2)) / (sigma*np.sqrt(2*np.pi))
    if observations != 0:
        return s/observations
    else:
        return s
    
    
##############################################################################


def cost_function_reduced_quiescence_ERL_only_30_parameters(parameter_vector, idx):    
    
    run_time=102.3
    omega=200
    observed_cells_time = 30
    sigma=5
    parameter_names = ['a','d','q']
    
    print("The parameters are ")
    print(parameter_vector)
    #np.save(path + "/fit_data/sobol_seq/cost_function_model_4/particle_position/parameter_vector_" + str(idx) + ".npy", parameter_vector)
    
    global model_concentration, model_population
    
    param_dict = {pname: pvalue for pname, pvalue in zip(pysb_sampled_parameter_names, parameter_vector)}   
    for pname, pvalue in param_dict.items():
        #Change model parameter values to current location in parameter space
        model_population.parameters[pname].value = (pvalue)
        model_concentration.parameters[pname].value = (pvalue)
    

    
###############################################################################
###############################################################################
    #Setting 1st condition
    parameters_sim={}
    cost_mmd_006_ERL=0
    cost_006_ERL=0    
    condition = '0.6 kPa ERL'
    init=16 #40
    parameters_sim[condition] = {}
    
    model_population.parameters['ERK'].value = 0.2 * omega
    model_population.parameters['FAK'].value = 0.4 * omega
    model_population.parameters['AKT'].value = 0.3 * omega
    
    model_concentration.parameters['ERK'].value = 1.0 * 1 ###
    model_concentration.parameters['FAK'].value = 0.4 * 1
    model_concentration.parameters['AKT'].value = 0.6 * 1 ###
    
###############################################################################
    
    try:
        (cells, cell_count) = runAKTERK(init, run_time, observed_cells_time, model_concentration, model_population, idx, condition);
        IMT_eoe_counts = retrieve_cellState(run_time, init, cells, return_all_prot=False, cells_starting_at_mitosis=False) #cell count
        all_events_simulated = call_events(run_time, IMT_eoe_counts); 
        all_events_simulated = all_events_simulated[all_events_simulated['event'] != 'Starting cell']
        all_events_simulated = all_events_simulated[np.array(all_events_simulated['Birthtime']<observed_cells_time)]
        mitotic_lifespans = all_events_simulated['Lifespan'][all_events_simulated['event']=='mitosis']
        censored_lifespans = all_events_simulated['Lifespan'][all_events_simulated['event']=='End_of_Expt']
        
        while len(mitotic_lifespans)<70: #while we have less than 100 timecourses with apoptosis, 
            init = init + 2
            (cells_tmp, cell_count) = runAKTERK(2, run_time, observed_cells_time, model_concentration, model_population, idx, condition);
            
            #merge simulations
            cells = {key: cells[key] + cells_tmp[key] for key in cells.keys()}
        
            IMT_eoe_counts = retrieve_cellState(run_time, init, cells, return_all_prot=False, cells_starting_at_mitosis=False) #cell count
            all_events_simulated = call_events(run_time, IMT_eoe_counts); 
            all_events_simulated = all_events_simulated[all_events_simulated['event'] != 'Starting cell']
            all_events_simulated = all_events_simulated[np.array(all_events_simulated['Birthtime']<observed_cells_time)]
            mitotic_lifespans = all_events_simulated['Lifespan'][all_events_simulated['event']=='mitosis']
            censored_lifespans = all_events_simulated['Lifespan'][all_events_simulated['event']=='End_of_Expt']
            
###############################################################################
        #Extract useful inforation from the simulations:
    
        simulated_IMT, simulated_sm, simulated_g1 = SMs_G1s_IMTs_extraction(cells) #G1, SM, IMT
    
        
        all_events_simulated.to_pickle(os.path.dirname(__file__) + "/../fit_data/sobol_seq/distribution/Events_" + str(idx) + "_006kPa_ERL.pkl") #doesn't save all imts, only the cells that are born before 30 hours
        np.save(os.path.dirname(__file__) + "/../fit_data/sobol_seq/distribution/G1_distribution_" + str(idx) + "_006kPa_ERL.npy",np.array(simulated_g1))
        np.save(os.path.dirname(__file__) + "/../fit_data/sobol_seq/distribution/SM_distribution_" + str(idx) + "_006kPa_ERL.npy",np.array(simulated_sm))
        
        #save density estimation
        x=np.linspace(0,120,360)
        y=[RBS(i, np.transpose([mitotic_lifespans.tolist()]), sigma) for i in x]
        np.save(os.path.dirname(__file__) + "/../fit_data/sobol_seq/distribution/IMT_distribution_" + str(idx) + "_006kPa_ERL_kde.npy",np.array(y))
    
###############################################################################
        #compute frac prolif parameters of the simulation:
    
        number_of_deaths = len(all_events_simulated) - (all_events_simulated['event']=='apoptosis').value_counts()[0]
        cell_time_total = sum(all_events_simulated['Lifespan'])
        m = fp.mle(mitotic_lifespans, censored_lifespans)
        parameters_sim[condition]['a'] = number_of_deaths/cell_time_total
        
        parameters_sim[condition]['mu'] = m['x'][0]
        parameters_sim[condition]['sigma'] = m['x'][1]
        parameters_sim[condition]['lambda'] = m['x'][2]
        parameters_sim[condition].update(fp.q_rates(m))
###############################################################################         
    #Compare G1 and SM istribution with experimental
        
        number_of_imt_to_sample=70
        experimental_IMTs_picked = random.sample(list(all_IMTs[(all_IMTs['Drug']=='ERL') & (all_IMTs['Stiffness']==0.2)]['IMT_t']), number_of_imt_to_sample)
        simulated_IMTs = random.sample(list(mitotic_lifespans), number_of_imt_to_sample)
        
        cost_mmd_006_ERL = nmmd2_rbf(np.array([simulated_IMTs]).transpose(), np.array([experimental_IMTs_picked]).transpose(), 1) #difference between IMT distribution
            
        cost_006_ERL = cost_mmd_006_ERL
    except :
        parameters_sim[condition]['a'] = np.nan
        parameters_sim[condition]['d'] = np.nan
        parameters_sim[condition]['q'] = np.nan
        parameters_sim[condition]['mu'] = np.nan
        parameters_sim[condition]['sigma'] = np.nan
        parameters_sim[condition]['lambda'] = np.nan
        
        cost_mmd_006_ERL = np.inf
        cost_006_ERL = np.inf

###############################################################################
###############################################################################
    #Setting 2nd condition
    
    init = 15
    condition = '55 kPa ERL'
    parameters_sim[condition] = {}
    
    cost_mmd_55_ERL=0
    cost_55_ERL=0
    
    model_population.parameters['ERK'].value = 0.2 * omega
    model_population.parameters['FAK'].value = 1.0 * omega
    model_population.parameters['AKT'].value = 0.6 * omega

    model_concentration.parameters['ERK'].value = 1 * 1
    model_concentration.parameters['FAK'].value = 1 * 1
    model_concentration.parameters['AKT'].value = 1 * 1

###############################################################################
    
    try:
        (cells, cell_count) = runAKTERK(init, run_time, observed_cells_time, model_concentration, model_population, idx, condition);
        IMT_eoe_counts = retrieve_cellState(run_time, init, cells, return_all_prot=False, cells_starting_at_mitosis=False) #cell count
        all_events_simulated = call_events(run_time, IMT_eoe_counts); 
        all_events_simulated = all_events_simulated[all_events_simulated['event'] != 'Starting cell']
        all_events_simulated = all_events_simulated[np.array(all_events_simulated['Birthtime']<observed_cells_time)]
        mitotic_lifespans = all_events_simulated['Lifespan'][all_events_simulated['event']=='mitosis']
        censored_lifespans = all_events_simulated['Lifespan'][all_events_simulated['event']=='End_of_Expt']

        while len(mitotic_lifespans)<70: #while we have less than 100 timecourses with apoptosis, 
            init = init + 2
            (cells_tmp, cell_count) = runAKTERK(2, run_time, observed_cells_time, model_concentration, model_population, idx, condition);
            
            #merge simulations
            cells = {key: cells[key] + cells_tmp[key] for key in cells.keys()}
        
            IMT_eoe_counts = retrieve_cellState(run_time, init, cells, return_all_prot=False, cells_starting_at_mitosis=False) #cell count
            all_events_simulated = call_events(run_time, IMT_eoe_counts); 
            all_events_simulated = all_events_simulated[all_events_simulated['event'] != 'Starting cell']
            all_events_simulated = all_events_simulated[np.array(all_events_simulated['Birthtime']<observed_cells_time)]
            mitotic_lifespans = all_events_simulated['Lifespan'][all_events_simulated['event']=='mitosis']
            censored_lifespans = all_events_simulated['Lifespan'][all_events_simulated['event']=='End_of_Expt']
            
            
###############################################################################
        #Extract useful inforation from the simulations:
    
        simulated_IMT, simulated_sm, simulated_g1 = SMs_G1s_IMTs_extraction(cells) #G1, SM, IMT
    
        #Save distribution
        all_events_simulated.to_pickle(os.path.dirname(__file__) + "/../fit_data/sobol_seq/distribution/Events_" + str(idx) + "_55kPa_ERL.pkl") #doesn't save all imts, only the cells that are born before 30 hours
        np.save(os.path.dirname(__file__) + "/../fit_data/sobol_seq/distribution/G1_distribution_" + str(idx) + "_55kPa_ERL.npy",np.array(simulated_g1))
        np.save(os.path.dirname(__file__) + "/../fit_data/sobol_seq/distribution/SM_distribution_" + str(idx) + "_55kPa_ERL.npy",np.array(simulated_sm))
        
        #save density estimation
        x=np.linspace(0,120,360)
        y=[RBS(i, np.transpose([mitotic_lifespans.tolist()]), sigma) for i in x]
        np.save(os.path.dirname(__file__) + "/../fit_data/sobol_seq/distribution/IMT_distribution_" + str(idx) + "_55kPa_ERL_kde.npy",np.array(y))
        
###############################################################################
        #compute frac prolif parameters of the simulation:
    
        number_of_deaths = len(all_events_simulated) - (all_events_simulated['event']=='apoptosis').value_counts()[0]
        cell_time_total = sum(all_events_simulated['Lifespan'])
        m = fp.mle(mitotic_lifespans, censored_lifespans)
        parameters_sim[condition]['a'] = number_of_deaths/cell_time_total
        
    
        parameters_sim[condition]['mu'] = m['x'][0]
        parameters_sim[condition]['sigma'] = m['x'][1]
        parameters_sim[condition]['lambda'] = m['x'][2]
        parameters_sim[condition].update(fp.q_rates(m))
            
###############################################################################         
    #Compare G1 and SM istribution with experimental
    
        number_of_imt_to_sample=70
        
    #        simulated_G1s = random.sample(simulated_g1, number_of_G1_to_sample)
    #        experimental_G1s_picked = random.sample(experimental_G1s['0.2 kPa ERL'], number_of_G1_to_sample)
    #        simulated_SMs = random.sample(simulated_sm, number_of_SM_to_sample)
    #        experimental_SMs_picked = random.sample(experimental_SMs['0.2 kPa ERL'], number_of_SM_to_sample)
        
        experimental_IMTs_picked = random.sample(list(all_IMTs[(all_IMTs['Drug']=='ERL') & (all_IMTs['Stiffness']==0.2)]['IMT_t']), number_of_imt_to_sample)
        simulated_IMTs = random.sample(list(mitotic_lifespans), number_of_imt_to_sample)
        
        cost_mmd_55_ERL = nmmd2_rbf(np.array([simulated_IMTs]).transpose(), np.array([experimental_IMTs_picked]).transpose(), 1) #difference between IMT distribution
        cost_55_ERL = cost_mmd_55_ERL
        
#            plt.hist(simulated_g1)
#            plt.hist(experimental_G1s_picked)

#        for name in parameter_names:
#            print(name)
#            if parameters_exp[condition][name]<10**-5:
#                parameters_exp[condition][name] = 1e-5
#            if parameters_sim[condition][name]<10**-5:
#                parameters_sim[condition][name] = 1e-5
#            cost_frac_prolif_model_006 += abs(np.log(parameters_sim[condition][name]) - np.log(parameters_exp[condition][name]))
    except :
        parameters_sim[condition]['a'] = np.nan
        parameters_sim[condition]['d'] = np.nan
        parameters_sim[condition]['q'] = np.nan
        parameters_sim[condition]['mu'] = np.nan
        parameters_sim[condition]['sigma'] = np.nan
        parameters_sim[condition]['lambda'] = np.nan
        
        cost_mmd_55_ERL = np.inf
        cost_55_ERL = np.inf
#        
    #save frac prolif parameters and cost + parameter vector
    frac_prolif_parameters_temp = {"index":[idx], "a_006_ERL":[parameters_sim['0.6 kPa ERL']['a']], "d_006_ERL":[parameters_sim['0.6 kPa ERL']['d']], "q_006_ERL":[parameters_sim['0.6 kPa ERL']['q']], "cost_mmd_006_ERL":[cost_mmd_006_ERL], "a_55_ERL":[parameters_sim['55 kPa ERL']['a']], "d_55_ERL":[parameters_sim['55 kPa ERL']['d']], "q_55_ERL":[parameters_sim['55 kPa ERL']['q']], "cost_mmd_55_ERL":[cost_mmd_55_ERL], "cost_total_ERL":[cost_mmd_006_ERL + cost_mmd_55_ERL]}
    vect = pd.DataFrame(np.array([parameter_vector]), columns = pysb_sampled_parameter_names)
    frac_prolif_parameters_temp = pd.DataFrame(data=frac_prolif_parameters_temp)
        
    frac_prolif_parameters_temp = pd.concat([frac_prolif_parameters_temp, vect], axis=1)
        
    frac_prolif_parameters_temp.to_pickle(os.path.dirname(__file__) + "/../fit_data/sobol_seq/cost/parameters_" + str(idx) + "_high_and_low_stiffness_ERL.pkl")        

    return cost_mmd_006_ERL + cost_mmd_55_ERL


