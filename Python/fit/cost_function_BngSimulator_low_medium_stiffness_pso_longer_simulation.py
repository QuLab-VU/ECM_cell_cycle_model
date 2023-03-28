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
import re

tmp_t=np.linspace(0, 10, 100)
sim = BngSimulator(model_population, tspan=tmp_t, verbose=True)
traj = sim.run()

sim = ScipyOdeSimulator(model_concentration, tspan=tmp_t, verbose=True)
traj = sim.run()

###############################################################################

omega=200

pysb_sampled_parameter_names = ["V_s2skp","V_Fp","k_M","k_EMY","k_CDF","k_CDS","k_ED","K_S","K_MY","K_DF","d_CD","k_P1", "k_P2", "K_CD", "k_E", "K_M", "K_CE", "V_2E2F", "K_2E2F", "k_RE", "kf1", "kf2", "kb", "kasbaxmbcl2", "kdsbaxmbcl2", "k1sbh3", "k2sbh3", "kdbh3", "kasbh3bcl2", "kdsbh3bcl2"]

###############################################################################


DATADIR = os.path.dirname(__file__) + "/../../Data/Exp_Cell_Counts"


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


def cost_function_30_parameters_pso(parameter_vector, idx, init_cell_number = [8, 12, 4, 4, 4, 4]):    


#    randomnumber = random.random()
#    print("random number is ")
#    print(randomnumber)
#    
#    print("The parameters are ")
#    print(parameter_vector)
#    np.save(os.path.dirname(__file__) + "/../fit_data/pso/particle_position/parameter_vector_particle_" + str(idx) + "_" + str(randomnumber) + ".npy", parameter_vector)
        
    global model_concentration, model_population, t, means
    
    param_dict = {pname: pvalue for pname, pvalue in zip(pysb_sampled_parameter_names, parameter_vector)}
    for pname, pvalue in param_dict.items():
        #Change model parameter values to current location in parameter space
        model_population.parameters[pname].value = (pvalue)
        model_concentration.parameters[pname].value = (pvalue)
        
    
    conditions = ['0.6 kPa DMSO', '0.6 kPa ERL', '4 kPa DMSO', '4 kPa ERL', '55 kPa DMSO', '55 kPa ERL']
    inputs = {}
    inputs['FAK'] = {'4 kPa DMSO': 0.6, '4 kPa ERL': 0.6, '0.6 kPa DMSO': 0.4, '0.6 kPa ERL': 0.4, '55 kPa DMSO': 1.0, '55 kPa ERL': 1.0}
    inputs['ERK'] = {'4 kPa DMSO': 1.0, '4 kPa ERL': 0.2, '0.6 kPa DMSO': 1.0, '0.6 kPa ERL': 0.2, '55 kPa DMSO': 1.0, '55 kPa ERL': 0.2}
    inputs['AKT'] = {'4 kPa DMSO': 1.0, '4 kPa ERL': 0.6, '0.6 kPa DMSO': 0.6, '0.6 kPa ERL': 0.3, '55 kPa DMSO': 1.0, '55 kPa ERL': 0.6}
    

    parameters_sim={} 
    cost = {}
    
    t_sim = np.arange(0,120,1/3)
    
    run_times = [120, 120, 120, 120, 60, 60] 
    omega = 200
    sigma = 5
    parameter_names = ['a','d','q']
    
#    run_times = [120] 
#    init_cell_number = [4]
#    conditions = ['0.6 kPa ERL']
    
#    condition = '55 kPa DMSO'
#    condition = '0.6 kPa ERL'
#    run_time = 120
#    run_time=60
#    condition = '55 kPa ERL'
###############################################################################
###############################################################################
    
    for init, run_time, condition, num_condition in zip(init_cell_number, run_times, conditions, range(len(conditions))):
        
        t_sim = np.arange(0,run_time,1/3)
        observed_cells_time = run_time
        cond_split = condition.split(' ')
        cond_split_DMSO = condition.split(' ')
        cond_split_DMSO[-1] = 'DMSO'
        parameters_sim[condition] = {}
        
        model_population.parameters['ERK'].value = inputs['ERK'][condition] * omega
        model_population.parameters['FAK'].value = inputs['FAK'][condition] * omega
        model_population.parameters['AKT'].value = inputs['AKT'][condition] * omega
        
        model_concentration.parameters['ERK'].value = inputs['ERK'][' '.join(cond_split_DMSO)] * 1
        model_concentration.parameters['FAK'].value = inputs['FAK'][' '.join(cond_split_DMSO)] * 1
        model_concentration.parameters['AKT'].value = inputs['AKT'][' '.join(cond_split_DMSO)] * 1
        
###############################################################################
            
        #try:
        
        (cells, cell_count) = runAKTERK(init, run_time, run_time, model_concentration, model_population, idx, condition);
        
        print('Cell count')
        print(cell_count)
        
        if type(cell_count)==np.ndarray:
            np.save(os.path.dirname(__file__) + "/../fit_data/pso/longer_sim/cell_count/cell_count_" + str(idx) + "_" + '_'.join(cond_split) + ".npy", cell_count)
        
        try:
            IMT_eoe_counts = retrieve_cellState(run_time, init, cells, return_all_prot=False, cells_starting_at_mitosis=False) #cell count
            all_events_simulated = call_events(run_time, IMT_eoe_counts); 
            all_events_simulated = all_events_simulated[all_events_simulated['event'] != 'Starting cell']
            all_events_simulated = all_events_simulated[np.array(all_events_simulated['Birthtime']<observed_cells_time)]
            mitotic_lifespans = all_events_simulated['Lifespan'][all_events_simulated['event']=='mitosis']
            censored_lifespans = all_events_simulated['Lifespan'][all_events_simulated['event']=='End_of_Expt']
            
            all_events_simulated.to_pickle(os.path.dirname(__file__) + "/../fit_data/pso/longer_sim/distribution/Events_" + str(idx) + "_" + '_'.join(cond_split) + ".pkl") #doesn't save all imts, only the cells that are born before 30 hours
        except:
            pass

        
###############################################################################
        #Extract useful inforation from the simulations:
        
        try:
            simulated_IMT, simulated_sm, simulated_g1 = SMs_G1s_IMTs_extraction(cells) #G1, SM, IMT
            
            np.save(os.path.dirname(__file__) + "/../fit_data/pso/longer_sim/distribution/G1_distribution_" + str(idx) + "_" + '_'.join(cond_split) + ".npy", np.array(simulated_g1))
            np.save(os.path.dirname(__file__) + "/../fit_data/pso/longer_sim/distribution/SM_distribution_" + str(idx) + "_" + '_'.join(cond_split) + ".npy", np.array(simulated_sm))
            
            #save density estimation
            x=np.linspace(0, run_time, 360)
            y=[RBS(i, np.transpose([mitotic_lifespans.tolist()]), sigma) for i in x]
            np.save(os.path.dirname(__file__) + "/../fit_data/pso/longer_sim/distribution/IMT_distribution_" + str(idx) + "_" + '_'.join(cond_split) + "_kde.npy", np.array(y))
        
        except:
            pass
        
###############################################################################
        #compute frac prolif parameters of the simulation:
        try:
            number_of_deaths = len(all_events_simulated) - (all_events_simulated['event']=='apoptosis').value_counts()[0]
            cell_time_total = sum(all_events_simulated['Lifespan'])        
            parameters_sim[condition]['a'] = number_of_deaths/cell_time_total
        except:
            parameters_sim[condition]['a'] = np.nan
        
        try:
            m = fp.mle(mitotic_lifespans, censored_lifespans)            
            parameters_sim[condition]['mu'] = m['x'][0]
            parameters_sim[condition]['sigma'] = m['x'][1]
            parameters_sim[condition]['lambda'] = m['x'][2]
            parameters_sim[condition].update(fp.q_rates(m))
        
        except:
            parameters_sim[condition]['d'] = np.nan
            parameters_sim[condition]['q'] = np.nan
            parameters_sim[condition]['mu'] = np.nan
            parameters_sim[condition]['sigma'] = np.nan
            parameters_sim[condition]['lambda'] = np.nan
###############################################################################         
        #cost computation: rmse dip rates
        
        try:
            cell_count_t = np.zeros(len(t_sim))
            time_tmp_m = 0
            
            for time_tmp_p, cell_count_time_tmp in cell_count:
                cell_count_t[(time_tmp_m <= t_sim) & (t_sim < time_tmp_p)] = cell_count_time_tmp
                time_tmp_m = time_tmp_p
            
            
            DIP_sim = np.log2(cell_count_t/cell_count_t[0])
            
            cost['Cost ' + condition] = np.sqrt(np.sum((DIP_sim[1:len(t_sim)] - means[num_condition][1:len(t_sim)])**2 / stdevs[num_condition][1:len(t_sim)]**2))/len(t_sim)
        
        except:
            cost['Cost ' + condition] = np.inf
        
###############################################################################         
        
        #except :
            #cost['Cost ' + condition] = np.inf



    parameters_sim = {param + ' ' + cond_tmp: [val_tmp] for cond_tmp in parameters_sim.keys() for param, val_tmp in parameters_sim[cond_tmp].items()}
    
    
    frac_prolif_parameters_temp = {"index":[idx], **parameters_sim, **cost}
    
    vect = pd.DataFrame(np.array([parameter_vector]), columns = pysb_sampled_parameter_names)
    
    frac_prolif_parameters_temp = pd.DataFrame(data=frac_prolif_parameters_temp)
        
    frac_prolif_parameters_temp = pd.concat([frac_prolif_parameters_temp, vect], axis=1)
        
    frac_prolif_parameters_temp.to_pickle(os.path.dirname(__file__) + "/../fit_data/pso/longer_sim/cost/parameters_" + str(idx) + "_pso.pkl")        

    return sum(cost.values())


