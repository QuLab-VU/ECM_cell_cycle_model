# -*- coding: utf-8 -*-
"""
Created on Wed May 20 10:42:26 2020

@author: Grim.M
"""
#time stopper:

#
#import signal
#
#class TimeoutException(RuntimeError):
#    """ Time out occurred! """
#    pass
#
#def handler(signum, frame):
#    print('forever is over!')
#    raise TimeoutException()
#


# Register the signal function handler
#signal.signal(signal.SIGALRM, handler)



path="/home/jeremie"
#path=r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code"


import sys
sys.path.insert(0, __file__ + "/../cellular_cycle") #add a pathway to the model
#sys.path.insert(0, path + "/cellular_cycle") #add a pathway to the model
#sys.path.insert(0, path + "/fit") #add a pathway to the model

from Austin_proliferation_model_concentration_4 import model
model_concentration=model
from Austin_proliferation_model_population_4 import model
model_population=model


from pysb.simulator import BngSimulator


import random

import os
os.environ["PYOPENCL_CTX"] ="0"

import pylab as pl
import pandas as pd
import math
import numpy as np
from runAKTERK_2_starting_at_mitosis import runAKTERK
from retrieve_cellState_dw import retrieve_cellState
from plot_frac_div_q_ap1 import all_events
import numpy as np
import pandas as pd
from pysb.simulator import BngSimulator
from pysb.simulator import StochKitSimulator
from pysb.simulator import ScipyOdeSimulator

from pysb.simulator import  CupSodaSimulator
import matplotlib.pyplot as plt
import pylab as pl
import fracprolif as fp
from rbf_mmd import nmmd2_rbf
import random

from SMs_G1s_IMTs_extraction import IMTs


tmp_t=pl.linspace(0, 10, 100)
sim = BngSimulator(model_population, tspan=tmp_t, verbose=True)
traj = sim.run()


sim = ScipyOdeSimulator(model_concentration, tspan=tmp_t, verbose=True)
traj = sim.run()

###############################################################################
all_events_exp = {}

run_time=102.3
omega=200

IMT_eoe_counts = pd.DataFrame(columns=['End_of_Expt', 'Death' ,'Birthtime','Lifespan'])

df = pd.read_csv(path + '/cell_fates/ds8_06_DMSO_fates.csv')
lines = {}
for col in df :
    lines[col]=df[col].dropna().tolist()

for (death, mitotic_time, end_time, lifetime) in zip(lines['Death'], lines['Mitotic Time'], lines['End Time'], lines['Lifetime']):
    tmp_counts = {'End_of_Expt' : end_time==102.3, 'Death' : death==1, 'Birthtime' : mitotic_time, 'Lifespan' : lifetime}
    IMT_eoe_counts = IMT_eoe_counts.append(tmp_counts, ignore_index=True)

all_events_exp['0.6 kPa DMSO'] = all_events(run_time, IMT_eoe_counts)



mitotic_lifespans_exp = all_events_exp['0.6 kPa DMSO']['Lifespan'][all_events_exp['0.6 kPa DMSO']['event']=='mitosis']
censored_lifespans_exp = all_events_exp['0.6 kPa DMSO']['Lifespan'][all_events_exp['0.6 kPa DMSO']['event']=='End_of_Expt']

number_of_deaths = len(all_events_exp['0.6 kPa DMSO']) - (all_events_exp['0.6 kPa DMSO']['event']=='apoptosis').value_counts()[0]
cell_time_total = sum(all_events_exp['0.6 kPa DMSO']['Lifespan'])


m = fp.mle(mitotic_lifespans_exp, censored_lifespans_exp)

parameters_exp={}
parameters_exp['0.6 kPa DMSO'] ={}

parameters_exp['0.6 kPa DMSO']['a'] = number_of_deaths/cell_time_total
parameters_exp['0.6 kPa DMSO']['mu'] = m['x'][0]
parameters_exp['0.6 kPa DMSO']['sigma'] = m['x'][1]
parameters_exp['0.6 kPa DMSO']['lambda'] = m['x'][2]
parameters_exp['0.6 kPa DMSO'].update(fp.q_rates(m))


df = pd.read_csv(path + '/cell_fates/ds8_06_DMSO_fates.csv')
lines = {}
for col in df :
    lines[col]=df[col].dropna().tolist()

pysb_sampled_parameter_names = ["V_s2skp","V_Fp","k_M","k_EMY","k_CDF","k_CDS","k_ED","K_S","K_MY","K_DF","d_CD","k_P1", "k_P2", "K_CD", "k_E", "K_M", "K_CE", "V_2E2F", "K_2E2F", "k_RE", "kf1", "kf2", "kb", "kasbaxmbcl2", "kdsbaxmbcl2", "k1sbh3", "k2sbh3", "kdbh3", "kasbh3bcl2", "kdsbh3bcl2"]

######################################################

experimental_SMs={}
experimental_G1s={}
experimental_IMTs={}


df = pd.read_csv(path+'/IMT_times/DS8_SM_0_2_55.csv')

for col in df.columns :
    experimental_SMs[col]=df[col].dropna().tolist()


df = pd.read_csv(path+'/IMT_times/DS8_G1_0_2_55.csv')

for col in df.columns :
    experimental_G1s[col]=df[col].dropna().tolist()

for k in experimental_G1s.keys() :
    minimum=min(len(experimental_G1s[k]),len(experimental_SMs[k]))
    experimental_IMTs[k]=[experimental_G1s[k][i]+experimental_SMs[k][i] for i in range(minimum)]


pysb_sampled_parameter_names = ["V_s2skp","V_Fp","k_M","k_EMY","k_CDF","k_CDS","k_ED","K_S","K_MY","K_DF","d_CD","k_P1", "k_P2", "K_CD", "k_E", "K_M", "K_CE", "V_2E2F", "K_2E2F", "k_RE", "kf1", "kf2", "kb", "kasbaxmbcl2", "kdsbaxmbcl2", "k1sbh3", "k2sbh3", "kdbh3", "kasbh3bcl2", "kdsbh3bcl2"]


##############################################################################
def RBS(x,x0,sigma) :
    observations , dimensions= np.shape(x0)
    s=0
    for i in range(observations):
        s+=math.exp(-np.linalg.norm(x-x0[i])**2/(2*sigma**2)) / (sigma*math.sqrt(2*math.pi))
    if observations != 0:
        return s/observations
    else:
        return s
    
    
##############################################################################

conditions=['0.6 kPa DMSO', '55 kPa DMSO']
parameter_names = ['a','d','q']

observed_cells_time = 30
sigma=5


#parameter_vector = parameter_vector_error
#parameter_vector = parameter_vector_error_over_60_sec
#parameter_vector = parameter_vector_error_runakterk
#parameter_vector = [0.01796364778625122, 0.11362811190427101, 0.0008276599610971477, 0.6451204157405731, 4.515842910184013, 0.016427851822702934, 6.794676345427822, 0.06526961947442374, 0.1934733242502379, 0.34232423258816486, 82.19791287433931, 0.3610094736457807, 116.75734850864045, 2.368276243820679, 5.816445291899858, 1.5990499381060275, 0.04505357247470116, 43.70269268253706, 19.379783033100328, 23.761725463454734, 0.3705507155537046, 60.679678761195824, 0.021925638754823482, 709.7991617788466, 0.03004544205271144, 5.955282169100989e-05, 0.0020243858391033893, 0.30033748241996405, 0.33677088633493957, 0.0028291772439897416]

def cost_function_reduced_quiescence_DMSO_only_30_parameters(parameter_vector, idx):    
    randomnumber = random.random()
    print("random number is ")
    print(randomnumber)
    
    print("The parameters are ")
    print(parameter_vector)
    np.save(path + "/fit_data/pso/high_and_low_stiffness_cost/particle_position/parameter_vector_particle_" + str(idx) + "_" + str(randomnumber) + ".npy", parameter_vector)
    
    global model_concentration, model_population, init, run_time, observed_cells_time, sigma
    
    param_dict = {pname: pvalue for pname, pvalue in zip(pysb_sampled_parameter_names, parameter_vector)}   
    for pname, pvalue in param_dict.items():
        #Change model parameter values to current location in parameter space
        model_population.parameters[pname].value = (pvalue)
        model_concentration.parameters[pname].value = (pvalue)
    
    parameters_sim={}
    cost_frac_prolif_model_006=0
    cost_mmd_006=0
    cost_006=0
    
    
    
    init = 30
    condition = '55 kPa DMSO'
    
    cost_mmd_55=0
    

    model_concentration.parameters['ERK'].value = 1
    model_concentration.parameters['FAK'].value = 1
    model_concentration.parameters['AKT'].value = 1
    
    model_population.parameters['ERK'].value = 1 * omega
    model_population.parameters['FAK'].value = 1 * omega
    model_population.parameters['AKT'].value = 1 * omega

    try:
        try:
            (cells, cell_count) = runAKTERK(init, run_time, observed_cells_time, model_concentration, model_population, randomnumber, condition)
            #signal.alarm(0)
        except :
            cells=0
            np.save(path + "/fit_data/pydream/frac_prolif_quiescence_DMSO_only_30_parameters/parameter_vector_error_runakterk.npy", parameter_vector)
        
        if cells == 0:
            return np.inf
        
        simulated_IMT, simulated_sm, simulated_g1 = IMTs(cells)
        
        
        IMT_eoe_counts = retrieve_cellState(run_time, init, cells)
        all_events_simulated = all_events(run_time, IMT_eoe_counts)
        all_events_simulated = all_events_simulated[np.array(all_events_simulated['Birthtime']<30)]
        all_events_simulated.to_pickle(path + "/fit_data/pso/high_and_low_stiffness_cost/all_events_simulated/all_events_simulated_particle_" + str(idx) + "_" + str(randomnumber) + "_30_initial_cells_55kPa_DMSO.pkl")        
        #finding the parameters for the EMG and fractional proliferation model
        
        if len(all_events_simulated)<=1 : #not enough events to optimize the fractional proliferation model
            return np.inf
        
        
        mitotic_lifespans = all_events_simulated['Lifespan'][all_events_simulated['event']=='mitosis']
        
        censored_lifespans = all_events_simulated['Lifespan'][all_events_simulated['event']=='End_of_Expt']
        
        if len(mitotic_lifespans)<=1 : #not enough events to optimize the fractional proliferation model
            return np.inf
        
        try :
            number_of_deaths = len(all_events_simulated) - (all_events_simulated['event']=='apoptosis').value_counts()[0]
            cell_time_total = sum(all_events_simulated['Lifespan'])
            
            m = fp.mle(mitotic_lifespans, censored_lifespans)
            
            parameters_sim[condition] = {}
            
            parameters_sim[condition]['a'] = number_of_deaths/cell_time_total
            
        except : 
            np.save(path + "/fit_data/pydream/frac_prolif_quiescence_DMSO_only_30_parameters/parameter_vector_error_mle.npy", parameter_vector)
            np.save(path + "/fit_data/pydream/frac_prolif_quiescence_DMSO_only_30_parameters/mle_censored_lifespans.npy", censored_lifespans)
            np.save(path + "/fit_data/pydream/frac_prolif_quiescence_DMSO_only_30_parameters/mle_mitotic_lifespans.npy", mitotic_lifespans)
            return np.inf
        
        print(parameters_sim)
        
        try :
            parameters_sim[condition]['mu'] = m['x'][0]
            parameters_sim[condition]['sigma'] = m['x'][1]
            parameters_sim[condition]['lambda'] = m['x'][2]
            parameters_sim[condition].update(fp.q_rates(m))
        except :
            np.save(path + "/fit_data/pydream/frac_prolif_quiescence_DMSO_only_30_parameters/parameter_vector_error_qrates.npy", parameter_vector)
            np.save(path + "/fit_data/pydream/frac_prolif_quiescence_DMSO_only_30_parameters/mitotic_lifespans.npy", mitotic_lifespans)
            np.save(path + "/fit_data/pydream/frac_prolif_quiescence_DMSO_only_30_parameters/censored_lifespans.npy", censored_lifespans)
            np.save(path + "/fit_data/pydream/frac_prolif_quiescence_DMSO_only_30_parameters/all_events_simulated.npy", all_events_simulated)
            return np.inf
            
        print(parameters_sim)
        #comparing those parameters and mmd:
        
        try:
            #number_of_IMT_to_sample = min (len(mitotic_lifespans), len(mitotic_lifespans_exp))
            number_of_SM_to_sample=160
            number_of_G1_to_sample=160
            
            simulated_G1s = random.sample(simulated_g1, number_of_G1_to_sample)
            experimental_G1s_picked = random.sample(experimental_G1s['55 kPa DMSO'], number_of_G1_to_sample)
            
            
            simulated_SMs = random.sample(simulated_sm, number_of_SM_to_sample)
            experimental_SMs_picked = random.sample(experimental_SMs['55 kPa DMSO'], number_of_SM_to_sample)
            
            cost_mmd_55 += nmmd2_rbf(np.array([simulated_G1s]).transpose(), np.array([experimental_G1s_picked]).transpose(), 1) #difference between IMT distribution
            cost_mmd_55 += nmmd2_rbf(np.array([simulated_SMs]).transpose(), np.array([experimental_SMs_picked]).transpose(), 1) #difference between IMT distribution
            
        except:
            cost_mmd+=160
            
        for name in parameter_names:
            print(name)
#                if parameters_exp[condition][name]<10**-5:
#                    parameters_exp[condition][name] = 1e-5
            if parameters_sim[condition][name]<10**-5:
                parameters_sim[condition][name] = 1e-5
#
#                cost_frac_prolif_model += abs(np.log(parameters_sim[condition][name]) - np.log(parameters_exp[condition][name]))
            
        print(path + "/fit_data/pso/high_and_low_stiffness_cost/cost_and_IMT_distribution/g1_distribution_particle_" + str(idx) + "_" + str(randomnumber) + "_30_initial_cells_55kPa_DMSO.npy")
        
        x=np.linspace(0,120,360)
        print(x)
        print('y')
        y=[RBS(i, np.transpose([mitotic_lifespans.tolist()]), sigma) for i in x]
        print(y)
        
        np.save(path + "/fit_data/pso/high_and_low_stiffness_cost/cost_and_IMT_distribution/IMT_distribution_particle_" + str(idx) + "_" + str(randomnumber) + "_30_initial_cells_55kPa_DMSO.npy",np.array(y))
        np.save(path + "/fit_data/pso/high_and_low_stiffness_cost/cost_and_IMT_distribution/g1_distribution_particle_" + str(idx) + "_" + str(randomnumber) + "_30_initial_cells_55kPa_DMSO.npy",np.array(simulated_g1))
        np.save(path + "/fit_data/pso/high_and_low_stiffness_cost/cost_and_IMT_distribution/SM_distribution_particle_" + str(idx) + "_" + str(randomnumber) + "_30_initial_cells_55kPa_DMSO.npy",np.array(simulated_sm)) 
        
    except :#TimeoutException as exc:
        np.save(path + "/fit_data/pydream/frac_prolif_quiescence_DMSO_only_30_parameters/parameter_vector_error_over_X_sec" + str(random.random()) + ".npy", parameter_vector)
        return np.inf
    
    
    
    

    
    ##############################################################
    #low stiffness
    
    condition = '0.6 kPa DMSO'
    init=50

    model_concentration.parameters['ERK'].value = 1.06970165050038 
    model_concentration.parameters['FAK'].value = 0.337143388171177
    model_concentration.parameters['AKT'].value = 0.616658111482869
    
    model_population.parameters['ERK'].value = 1.06970165050038 * omega
    model_population.parameters['FAK'].value = 0.337143388171177 * omega
    model_population.parameters['AKT'].value = 0.616658111482869 * omega
    
    #running the simulation with the experiment parameters (run_time)
    #signal.alarm(100)
    try:
        try:
            (cells, cell_count) = runAKTERK(init, run_time, observed_cells_time, model_concentration, model_population, randomnumber, condition)
            #signal.alarm(0)
        except :
            cells=0
            np.save(path + "/fit_data/pydream/frac_prolif_quiescence_DMSO_only_30_parameters/parameter_vector_error_runakterk.npy", parameter_vector)
        
        if cells == 0:
            return np.inf
        
        IMT_eoe_counts = retrieve_cellState(run_time, init, cells)
        all_events_simulated = all_events(run_time, IMT_eoe_counts)
        all_events_simulated = all_events_simulated[np.array(all_events_simulated['Birthtime']<30)]
        all_events_simulated.to_pickle(path + "/fit_data/pso/high_and_low_stiffness_cost/all_events_simulated/all_events_simulated_particle_" + str(idx) + "_" + str(randomnumber) + "_50_initial_cells.pkl")
        
        #finding the parameters for the EMG and fractional proliferation model
        
        
        
        if len(all_events_simulated)<=1 : #not enough events to optimize the fractional proliferation model
            return np.inf
        
        
        mitotic_lifespans = all_events_simulated['Lifespan'][all_events_simulated['event']=='mitosis']
        
        censored_lifespans = all_events_simulated['Lifespan'][all_events_simulated['event']=='End_of_Expt']
        
        if len(mitotic_lifespans)<=1 : #not enough events to optimize the fractional proliferation model
            return np.inf
        
        try :
            number_of_deaths = len(all_events_simulated) - (all_events_simulated['event']=='apoptosis').value_counts()[0]
            cell_time_total = sum(all_events_simulated['Lifespan'])
            
            m = fp.mle(mitotic_lifespans, censored_lifespans)
            
            parameters_sim[condition] = {}
            
            parameters_sim[condition]['a'] = number_of_deaths/cell_time_total
            
        except : 
            np.save(path + "/fit_data/pydream/frac_prolif_quiescence_DMSO_only_30_parameters/parameter_vector_error_mle.npy", parameter_vector)
            np.save(path + "/fit_data/pydream/frac_prolif_quiescence_DMSO_only_30_parameters/mle_censored_lifespans.npy", censored_lifespans)
            np.save(path + "/fit_data/pydream/frac_prolif_quiescence_DMSO_only_30_parameters/mle_mitotic_lifespans.npy", mitotic_lifespans)
            return np.inf
        
        try :
            parameters_sim[condition]['mu'] = m['x'][0]
            parameters_sim[condition]['sigma'] = m['x'][1]
            parameters_sim[condition]['lambda'] = m['x'][2]
            parameters_sim[condition].update(fp.q_rates(m))
        except :
            np.save(path + "/fit_data/pydream/frac_prolif_quiescence_DMSO_only_30_parameters/parameter_vector_error_qrates.npy", parameter_vector)
            np.save(path + "/fit_data/pydream/frac_prolif_quiescence_DMSO_only_30_parameters/mitotic_lifespans.npy", mitotic_lifespans)
            np.save(path + "/fit_data/pydream/frac_prolif_quiescence_DMSO_only_30_parameters/censored_lifespans.npy", censored_lifespans)
            np.save(path + "/fit_data/pydream/frac_prolif_quiescence_DMSO_only_30_parameters/all_events_simulated.npy", all_events_simulated)
            return np.inf
            
            
        #comparing those parameters and mmd:
        
        try:
            #number_of_IMT_to_sample = min (len(mitotic_lifespans), len(mitotic_lifespans_exp))
            number_of_IMT_to_sample=72
            
            simulated_IMT = random.sample(mitotic_lifespans.to_list(), number_of_IMT_to_sample)
            experimental_IMT = random.sample(mitotic_lifespans_exp.to_list(), number_of_IMT_to_sample)
            
            cost_mmd_006 += nmmd2_rbf(np.array([simulated_IMT]).transpose(), np.array([experimental_IMT]).transpose(), 1) #difference between IMT distribution
            
        except:
            cost_mmd+=60
            
        for name in parameter_names:
            print(name)
            if parameters_exp[condition][name]<10**-5:
                parameters_exp[condition][name] = 1e-5
            if parameters_sim[condition][name]<10**-5:
                parameters_sim[condition][name] = 1e-5
            
            cost_frac_prolif_model_006 += abs(np.log(parameters_sim[condition][name]) - np.log(parameters_exp[condition][name]))
            
        cost_006 += cost_frac_prolif_model_006 + cost_mmd_006
        
        x=np.linspace(0,120,360)
        y=[RBS(i, np.transpose([mitotic_lifespans.tolist()]), sigma) for i in x]
        np.save(path + "/fit_data/pso/high_and_low_stiffness_cost/cost_and_IMT_distribution/IMT_distribution_particle_" + str(idx) + "_" + str(randomnumber) + "_50_initial_cells.npy",np.array(y))
        
        print(idx)
        print(parameters_sim)
        print(cost_006)
        print(cost_mmd_006)
        print(cost_mmd_55)
        print(cost_frac_prolif_model_006)
        
        frac_prolif_parameters_temp = {"index":[idx], "a_006":[parameters_sim['0.6 kPa DMSO']['a']], "d_006":[parameters_sim['0.6 kPa DMSO']['d']], "q_006":[parameters_sim['0.6 kPa DMSO']['q']], "cost_006":[cost_006], "cost_mmd_006":[cost_mmd_006], "cost_frac_prolif_model_006":[cost_frac_prolif_model_006], "a_55":[parameters_sim['55 kPa DMSO']['a']], "d_55":[parameters_sim['55 kPa DMSO']['d']], "q_55":[parameters_sim['55 kPa DMSO']['q']], "cost_mmd_55":[cost_mmd_55], "cost_total":[cost_006 + cost_mmd_55 * (3/8)]}
        
        print(frac_prolif_parameters_temp)
        
        frac_prolif_parameters_temp = pd.DataFrame(data=frac_prolif_parameters_temp)
        frac_prolif_parameters_temp.to_pickle(path + "/fit_data/pso/high_and_low_stiffness_cost/cost_and_IMT_distribution/parameters_particle_" + str(idx) + "_" + str(randomnumber) + "_high_and_low_stiffness.pkl")
        
        
        
    except :#TimeoutException as exc:
        np.save(path + "/fit_data/pydream/frac_prolif_quiescence_DMSO_only_30_parameters/parameter_vector_error_over_X_sec" + str(random.random()) + ".npy", parameter_vector)
        return np.inf
    
    return cost_006 + cost_mmd_55 * (3/8)





"""

        
        

        


            
"""

"""
    



#
#
#all_events_simulated = np.load(path + "/fit_data/pydream/frac_prolif_quiescence_DMSO_only_30_parameters/error2/all_events_simulated.npy", allow_pickle=True)
#
#all_events = pd.DataFrame(columns=['event' ,'Birthtime','Lifespan'])
#
#for i in range(len(all_events_simulated)):
#    all_events = all_events.append(dict(zip(all_events.columns,all_events_simulated[i].tolist())), ignore_index=True)
#    
#all_events_simulated = all_events
    

#different types of error
    

#over 60 sec running : parameter_vector_error_over_60_sec: should b okay now
#runakterk doesn't work : parameter_vector_error_runakterk : should b okay now
#fractional proliferation model is not fitted by the data : mle_censored_lifespans and mle_mitotic_lifespans : 
    
"""