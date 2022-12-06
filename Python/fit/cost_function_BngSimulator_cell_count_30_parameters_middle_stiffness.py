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


def cost_function_4kPa_30_parameters(parameter_vector, idx):    
    
    run_time=102.3
    omega=200
    observed_cells_time = 102.3
    sigma=5
    parameter_names = ['a','d','q']
    
    
    print("The parameters are ")
    print(parameter_vector)
    #np.save(path + "/fit_data/sobol_seq/cost_function_model_4/particle_position/parameter_vector_" + str(idx) + ".npy", parameter_vector)
    
    global model_concentration, model_population, t, means
    
    param_dict = {pname: pvalue for pname, pvalue in zip(pysb_sampled_parameter_names, parameter_vector)}   
    for pname, pvalue in param_dict.items():
        #Change model parameter values to current location in parameter space
        model_population.parameters[pname].value = (pvalue)
        model_concentration.parameters[pname].value = (pvalue)
    
###############################################################################
###############################################################################
    #Setting 1st condition
    parameters_sim={} 
    condition = '4 kPa DMSO'
    init=3 #40
    parameters_sim[condition] = {}
    
    model_population.parameters['ERK'].value = 1.0 * omega
    model_population.parameters['FAK'].value = 0.6 * omega
    model_population.parameters['AKT'].value = 1.0 * omega
    
    model_concentration.parameters['ERK'].value = 1.0 * 1
    model_concentration.parameters['FAK'].value = 0.6 * 1
    model_concentration.parameters['AKT'].value = 1.0 * 1
    
###############################################################################
    
    try:
        (cells, cell_count) = runAKTERK(init, t[-1], t[-1], model_concentration, model_population, idx, condition);
        
        
        np.save(os.path.dirname(__file__) + "/../fit_data/sobol_seq/cell_count/cell_count_" + str(idx) + "_4kPa_DMSO.npy", cell_count)

        IMT_eoe_counts = retrieve_cellState(run_time, init, cells, return_all_prot=False, cells_starting_at_mitosis=False) #cell count
        all_events_simulated = call_events(run_time, IMT_eoe_counts); 
        all_events_simulated = all_events_simulated[all_events_simulated['event'] != 'Starting cell']
        all_events_simulated = all_events_simulated[np.array(all_events_simulated['Birthtime']<observed_cells_time)]
        mitotic_lifespans = all_events_simulated['Lifespan'][all_events_simulated['event']=='mitosis']
        censored_lifespans = all_events_simulated['Lifespan'][all_events_simulated['event']=='End_of_Expt']
        
            
###############################################################################
        #Extract useful inforation from the simulations:
    
        simulated_IMT, simulated_sm, simulated_g1 = SMs_G1s_IMTs_extraction(cells) #G1, SM, IMT
    
        
        all_events_simulated.to_pickle(os.path.dirname(__file__) + "/../fit_data/sobol_seq/distribution/Events_" + str(idx) + "_4kPa_DMSO.pkl") #doesn't save all imts, only the cells that are born before 30 hours
        np.save(os.path.dirname(__file__) + "/../fit_data/sobol_seq/distribution/G1_distribution_" + str(idx) + "_4kPa_DMSO.npy",np.array(simulated_g1))
        np.save(os.path.dirname(__file__) + "/../fit_data/sobol_seq/distribution/SM_distribution_" + str(idx) + "_4kPa_DMSO.npy",np.array(simulated_sm))
        
        #save density estimation
        x=np.linspace(0,120,360)
        y=[RBS(i, np.transpose([mitotic_lifespans.tolist()]), sigma) for i in x]
        np.save(os.path.dirname(__file__) + "/../fit_data/sobol_seq/distribution/IMT_distribution_" + str(idx) + "_4kPa_DMSO_kde.npy",np.array(y))
    
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
        

    except :
        parameters_sim[condition]['a'] = np.nan
        parameters_sim[condition]['d'] = np.nan
        parameters_sim[condition]['q'] = np.nan
        parameters_sim[condition]['mu'] = np.nan
        parameters_sim[condition]['sigma'] = np.nan
        parameters_sim[condition]['lambda'] = np.nan
        


###############################################################################
###############################################################################
    #Setting 2nd condition
    
    init = 8
    condition = '4 kPa ERL'
    parameters_sim[condition] = {}
    
    
    model_population.parameters['ERK'].value = 0.2 * omega
    model_population.parameters['FAK'].value = 0.6 * omega
    model_population.parameters['AKT'].value = 0.6 * omega
    
    model_concentration.parameters['ERK'].value = 1.0 * 1
    model_concentration.parameters['FAK'].value = 0.6 * 1
    model_concentration.parameters['AKT'].value = 1.0 * 1

###############################################################################
    
    try:
        (cells, cell_count) = runAKTERK(init, t[-1], t[-1], model_concentration, model_population, idx, condition);
        
        np.save(os.path.dirname(__file__) + "/../fit_data/sobol_seq/cell_count/cell_count_" + str(idx) + "_4kPa_ERL.npy", cell_count)
        
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
        all_events_simulated.to_pickle(os.path.dirname(__file__) + "/../fit_data/sobol_seq/distribution/Events_" + str(idx) + "_4kPa_ERL.pkl") #doesn't save all imts, only the cells that are born before 30 hours
        np.save(os.path.dirname(__file__) + "/../fit_data/sobol_seq/distribution/G1_distribution_" + str(idx) + "_4kPa_ERL.npy",np.array(simulated_g1))
        np.save(os.path.dirname(__file__) + "/../fit_data/sobol_seq/distribution/SM_distribution_" + str(idx) + "_4kPa_ERL.npy",np.array(simulated_sm))
        
        #save density estimation
        x=np.linspace(0,120,360)
        y=[RBS(i, np.transpose([mitotic_lifespans.tolist()]), sigma) for i in x]
        np.save(os.path.dirname(__file__) + "/../fit_data/sobol_seq/distribution/IMT_distribution_" + str(idx) + "_4kPa_ERL_kde.npy",np.array(y))
        
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

    except :
        parameters_sim[condition]['a'] = np.nan
        parameters_sim[condition]['d'] = np.nan
        parameters_sim[condition]['q'] = np.nan
        parameters_sim[condition]['mu'] = np.nan
        parameters_sim[condition]['sigma'] = np.nan
        parameters_sim[condition]['lambda'] = np.nan
        
#        """
    #save frac prolif parameters and cost + parameter vector
    frac_prolif_parameters_temp = {"index":[idx], "a_4_ERL":[parameters_sim['4 kPa ERL']['a']], "d_4_ERL":[parameters_sim['4 kPa ERL']['d']], "q_4_ERL":[parameters_sim['4 kPa ERL']['q']], "a_4_DMSO":[parameters_sim['4 kPa DMSO']['a']], "d_4_DMSO":[parameters_sim['4 kPa DMSO']['d']], "q_4_DMSO":[parameters_sim['4 kPa DMSO']['q']]}
    vect = pd.DataFrame(np.array([parameter_vector]), columns = pysb_sampled_parameter_names)
    frac_prolif_parameters_temp = pd.DataFrame(data=frac_prolif_parameters_temp)
        
    frac_prolif_parameters_temp = pd.concat([frac_prolif_parameters_temp, vect], axis=1)
        
    frac_prolif_parameters_temp.to_pickle(os.path.dirname(__file__) + "/../fit_data/sobol_seq/cost/parameters_" + str(idx) + "_4kPa_DMSO_ERL.pkl")        

    return vect


