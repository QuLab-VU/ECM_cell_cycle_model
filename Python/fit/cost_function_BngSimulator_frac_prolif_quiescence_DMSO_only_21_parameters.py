# -*- coding: utf-8 -*-
"""
Created on Wed May 20 10:42:26 2020

@author: Grim.M
"""

path="/home/jeremie"

#path=r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code"

import sys
sys.path.insert(0, path+"/cellular_cycle") #add a pathway to the model
sys.path.insert(0, path+"/fit") #add a pathway to the model

from Austin_proliferation_model_population_3 import model
model_population=model
from Austin_proliferation_model_concentration_3 import model
model_concentration=model


from pysb.simulator import BngSimulator


import random

import os
os.environ["PYOPENCL_CTX"] ="0"

import pylab as pl
import pandas as pd
import math
import numpy as np
from runAKTERK import runAKTERK
from retrieve_cellState_dw import retrieve_cellState
from plot_frac_div_q_ap1 import all_events
import numpy as np
import pandas as pd
from pysb.simulator import BngSimulator
from pysb.simulator import ScipyOdeSimulator
import matplotlib.pyplot as plt
import pylab as pl
import fracprolif as fp


tmp_t=pl.linspace(0, 500, 2500)
sim = BngSimulator(model_population, tspan=tmp_t, verbose=True)
traj = sim.run()

###############################################################################
all_events_exp = {}

run_time=120
omega=200
init=10


IMT_eoe_counts = pd.DataFrame(columns=['End_of_Expt', 'Death' ,'Birthtime','Lifespan'])

df = pd.read_csv(path+'/cell_fates/ds8_06_DMSO_fates.csv')
lines = {}
for col in df :
    lines[col]=df[col].dropna().tolist()

for (death, mitotic_time, end_time, lifetime) in zip(lines['Death'], lines['Mitotic Time'], lines['End Time'], lines['Lifetime']):
    tmp_counts = {'End_of_Expt' : end_time==102.3, 'Death' : death==1, 'Birthtime' : mitotic_time, 'Lifespan' : lifetime}
    IMT_eoe_counts = IMT_eoe_counts.append(tmp_counts, ignore_index=True)
    
all_events_exp['0.6 kPa DMSO'] = all_events(run_time, IMT_eoe_counts)



mitotic_lifespans = all_events_exp['0.6 kPa DMSO']['Lifespan'][all_events_exp['0.6 kPa DMSO']['event']=='mitosis']
censored_lifespans = all_events_exp['0.6 kPa DMSO']['Lifespan'][all_events_exp['0.6 kPa DMSO']['event']=='End_of_Expt']

number_of_deaths = len(all_events_exp['0.6 kPa DMSO']) - (all_events_exp['0.6 kPa DMSO']['event']=='apoptosis').value_counts()[0]
cell_time_total = sum(all_events_exp['0.6 kPa DMSO']['Lifespan'])

m = fp.mle(mitotic_lifespans, censored_lifespans)

parameters_exp={}
parameters_exp['0.6 kPa DMSO'] ={}

parameters_exp['0.6 kPa DMSO']['a'] = number_of_deaths/cell_time_total
parameters_exp['0.6 kPa DMSO']['mu'] = m['x'][0]
parameters_exp['0.6 kPa DMSO']['sigma'] = m['x'][1]
parameters_exp['0.6 kPa DMSO']['lambda'] = m['x'][2]
parameters_exp['0.6 kPa DMSO'].update(fp.q_rates(m))


df = pd.read_csv(path+'/cell_fates/ds8_06_DMSO_fates.csv')
lines = {}
for col in df :
    lines[col]=df[col].dropna().tolist()

pysb_sampled_parameter_names = ["V_s2skp","V_Fp","k_M","k_EMY","k_CDF","k_CDS","k_ED","K_S","K_MY","K_DF","d_CD","k_P1", "k_P2", "K_CD", "k_E", "K_M", "K_CE", "V_2E2F", "K_2E2F", "k_RE", "k_b"]


##############################################################################
conditions=['0.6 kPa DMSO']
parameter_names = ['a','d','q']


def cost_function_reduced_quiescence_DMSO_only_21_parameters(parameter_vector):    
    
    print("The parameters are ")
    print(parameter_vector)
    
    global model_concentration, model_population, init, run_time
    
    param_dict = {pname: pvalue for pname, pvalue in zip(pysb_sampled_parameter_names, parameter_vector)}   
    for pname, pvalue in param_dict.items():
        #Change model parameter values to current location in parameter space
        model_population.parameters[pname].value = (pvalue)
        model_concentration.parameters[pname].value = (pvalue)

    parameters_sim={}
    cost=0
    
    for condition in conditions:
        print(condition)
        
        if 'DMSO' in condition:            
            if '0.6 kPa' in condition :
                model_concentration.parameters['ERK'].value = 1.06970165050038 
                model_concentration.parameters['FAK'].value = 0.337143388171177
                model_concentration.parameters['AKT'].value = 0.616658111482869
                
                model_population.parameters['ERK'].value = 1.06970165050038 * omega
                model_population.parameters['FAK'].value = 0.337143388171177 * omega
                model_population.parameters['AKT'].value = 0.616658111482869 * omega
            elif '4 kPa' in condition :
                model_concentration.parameters['ERK'].value = 0.897485095798471
                model_concentration.parameters['FAK'].value = 0.594150207445117
                model_concentration.parameters['AKT'].value = 1.01975420186633
                
                model_population.parameters['ERK'].value = 0.897485095798471 * omega
                model_population.parameters['FAK'].value = 0.594150207445117 * omega
                model_population.parameters['AKT'].value = 1.01975420186633 * omega
            else : #55 kPa
                model_concentration.parameters['ERK'].value = 1
                model_concentration.parameters['FAK'].value = 1
                model_concentration.parameters['AKT'].value = 1
                
                model_population.parameters['ERK'].value = 1 * omega
                model_population.parameters['FAK'].value = 1 * omega
                model_population.parameters['AKT'].value = 1 * omega
        else : #ERL treated
            if '0.6 kPa' in condition :
                model_concentration.parameters['ERK'].value = 0.186745227418267
                model_concentration.parameters['FAK'].value = 0.397060244114162
                model_concentration.parameters['AKT'].value = 0.245702711121333
                
                model_population.parameters['ERK'].value = 0.186745227418267 * omega
                model_population.parameters['FAK'].value = 0.397060244114162 * omega
                model_population.parameters['AKT'].value = 0.245702711121333 * omega
            elif '4 kPa' in condition :
                model_concentration.parameters['ERK'].value = 0.211453830597783
                model_concentration.parameters['FAK'].value = 0.709172297576571
                model_concentration.parameters['AKT'].value = 0.598518145949655
                
                model_population.parameters['ERK'].value = 0.211453830597783 * omega
                model_population.parameters['FAK'].value = 0.709172297576571 * omega
                model_population.parameters['AKT'].value = 0.598518145949655 * omega
            else : #55 kPa
                model_concentration.parameters['ERK'].value = 0.273618423688386
                model_concentration.parameters['FAK'].value = 1.075820621939
                model_concentration.parameters['AKT'].value = 0.585067267459499
                
                model_population.parameters['ERK'].value = 0.273618423688386 * omega
                model_population.parameters['FAK'].value = 1.075820621939 * omega
                model_population.parameters['AKT'].value = 0.585067267459499 * omega
        
        #running the simulation with the experiment parameters (run_time)
        (cells, cell_count) = runAKTERK(init, run_time, model_concentration, model_population)
        
        if cells == 0:
            return 10000
        
        IMT_eoe_counts = retrieve_cellState(run_time, init_cell, cells)
        all_events_simulated = all_events(run_time, IMT_eoe_counts)
        
        #finding the parameters for the EMG and fractional proliferation model
        
        mitotic_lifespans = all_events_simulated['Lifespan'][all_events_simulated['event']=='mitosis']
        censored_lifespans = all_events_simulated['Lifespan'][all_events_simulated['event']=='End_of_Expt']
        
        number_of_deaths = len(all_events_simulated) - (all_events_simulated['event']=='apoptosis').value_counts()[0]
        cell_time_total = sum(all_events_simulated['Lifespan'])
        
        m = fp.mle(mitotic_lifespans, censored_lifespans)
        
        parameters_sim[condition] ={}

        parameters_sim[condition]['a'] = number_of_deaths/cell_time_total
        parameters_sim[condition]['mu'] = m['x'][0]
        parameters_sim[condition]['sigma'] = m['x'][1]
        parameters_sim[condition]['lambda'] = m['x'][2]
        parameters_sim[condition].update(fp.q_rates(m))
        
        #comparing those parameters:
        
        for name in parameter_names:
            print(name)
            cost += (parameters_sim[condition][name] - parameters_exp[condition][name])**2
        
    return cost


