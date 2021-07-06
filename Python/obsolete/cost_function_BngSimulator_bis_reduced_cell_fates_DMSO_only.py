# -*- coding: utf-8 -*-
"""
Created on Wed May 20 10:42:26 2020

@author: Grim.M
"""

path="/home/jeremie"

path=r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code"

import sys
sys.path.insert(0, path+"/cellular_cycle") #add a pathway to the model
sys.path.insert(0, path+"/fit") #add a pathway to the model

from Austin_proliferation_model_population import model
model_population=model
from Austin_proliferation_model_concentration_bis import model
model_concentration=model

from pysb.simulator import BngSimulator
from runAKTERK import runAKTERK
from retrieve_cellState_dw import retrieve_cellState
from plot_frac_div_q_ap1 import all_events
from rbf_mmd import nmmd2_rbf
import random
import pylab as pl
import pandas as pd
import math
import numpy as np

tmp_t=pl.linspace(0, 500, 2500)
sim = BngSimulator(model_population, tspan=tmp_t, verbose=True)
traj = sim.run()

###############################################################################
#IMTs for 55 kPa DMSO

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

#0.6 kPa DMSO
    
run_time=102.3
omega=200

df = pd.read_csv(path+'/cell_fates/ds8_06_DMSO_fates.csv')
lines = {}
for col in df :
    lines[col]=df[col].dropna().tolist()

IMT_eoe_counts = pd.DataFrame(columns=['End_of_Expt', 'Death' ,'Birthtime','Lifespan'])

for (death, mitotic_time, end_time, lifetime) in zip(lines['Death'], lines['Mitotic Time'], lines['End Time'], lines['Lifetime']):
    tmp_counts = {'End_of_Expt' : end_time==102.3, 'Death' : death==1, 'Birthtime' : mitotic_time, 'Lifespan' : lifetime}
    IMT_eoe_counts = IMT_eoe_counts.append(tmp_counts, ignore_index=True)
    
all_events_exp = all_events(run_time, IMT_eoe_counts)


proportion_exp_mitosis = {}
proportion_exp_quiescent_cells = {}
proportion_exp_apoptosis = {}
proportion_exp_end_of_experiment = {}

proportion_exp_mitosis['0.6 kPa DMSO'] = len(all_events_exp[all_events_exp['event']=='mitosis'])/len(all_events_exp)
proportion_exp_mitosis['55 kPa DMSO'] = 0

proportion_exp_quiescent_cells['0.6 kPa DMSO'] = len(all_events_exp[all_events_exp['event']=='quiescence'])/len(all_events_exp)
proportion_exp_quiescent_cells['55 kPa DMSO'] = 0

proportion_exp_apoptosis['0.6 kPa DMSO'] = len(all_events_exp[all_events_exp['event']=='apoptosis'])/len(all_events_exp)
proportion_exp_apoptosis['55 kPa DMSO'] = 0

proportion_exp_end_of_experiment['0.6 kPa DMSO'] = len(all_events_exp[all_events_exp['event']=='End_of_Expt'])/len(all_events_exp)
proportion_exp_end_of_experiment['55 kPa DMSO'] = 0


###############################################################################
#parameters

pysb_sampled_parameter_names = ["V_s2skp","V_Fp","k_M","k_EMY","k_CDF","k_CDS","k_ED","K_S","K_MY","K_DF","d_CD","k_P1", "k_P2", "K_CD", "k_E", "K_M", "K_CE", "V_2E2F", "K_2E2F", "k_RE", "k_b"]

init = 3
run_time = 102.3 
conditions=['0.2 kPa DMSO','55 kPa DMSO']

##############################################################################


def cost_function_reduced_quiescence_DMSO_only_21_parameters(parameter_vector):    
    
    print("les paramètres sont ")
    print(parameter_vector)
    
    global experimental_G1s, experimental_SMs, model_concentration, model_population, tmp_t, init, run_time
    
    param_dict = {pname: pvalue for pname, pvalue in zip(pysb_sampled_parameter_names, parameter_vector)}
    
    
    for pname, pvalue in param_dict.items():
        
        #Change model parameter values to current location in parameter space
        model_population.parameters[pname].value = (pvalue)
        
    #Simulate experimentally IMTs
    n_mmd2_tot=[]
    simulated_G1s=[]
    simulated_SMs=[]
    omega=200
    
    for condition in conditions:
        print(condition)
        
        if 'DMSO' in condition:            
            if '0.2 kPa' in condition :
                #model_concentration.parameters['ERK'].value = 1.06970165050038 
                #model_concentration.parameters['FAK'].value = 0.337143388171177
                #model_concentration.parameters['AKT'].value = 0.616658111482869
                
                model_population.parameters['ERK'].value = 1.06970165050038 * omega
                model_population.parameters['FAK'].value = 0.337143388171177 * omega
                model_population.parameters['AKT'].value = 0.616658111482869 * omega
            elif '4 kPa' in condition :
                #model_concentration.parameters['ERK'].value = 0.897485095798471
                #model_concentration.parameters['FAK'].value = 0.594150207445117
                #model_concentration.parameters['AKT'].value = 1.01975420186633
                
                model_population.parameters['ERK'].value = 0.897485095798471 * omega
                model_population.parameters['FAK'].value = 0.594150207445117 * omega
                model_population.parameters['AKT'].value = 1.01975420186633 * omega
            else : #55 kPa
                #model_concentration.parameters['ERK'].value = 1
                #model_concentration.parameters['FAK'].value = 1
                #model_concentration.parameters['AKT'].value = 1
                
                model_population.parameters['ERK'].value = 1 * omega
                model_population.parameters['FAK'].value = 1 * omega
                model_population.parameters['AKT'].value = 1 * omega
        else : #ERL treated
            if '0.2 kPa' in condition :
                #model_concentration.parameters['ERK'].value = 0.186745227418267
                #model_concentration.parameters['FAK'].value = 0.397060244114162
                #model_concentration.parameters['AKT'].value = 0.245702711121333
                
                model_population.parameters['ERK'].value = 0.186745227418267 * omega
                model_population.parameters['FAK'].value = 0.397060244114162 * omega
                model_population.parameters['AKT'].value = 0.245702711121333 * omega
            elif '4 kPa' in condition :
                #model_concentration.parameters['ERK'].value = 0.211453830597783
                #model_concentration.parameters['FAK'].value = 0.709172297576571
                #model_concentration.parameters['AKT'].value = 0.598518145949655
                
                model_population.parameters['ERK'].value = 0.211453830597783 * omega
                model_population.parameters['FAK'].value = 0.709172297576571 * omega
                model_population.parameters['AKT'].value = 0.598518145949655 * omega
            else : #55 kPa
                #model_concentration.parameters['ERK'].value = 0.273618423688386
                #model_concentration.parameters['FAK'].value = 1.075820621939
                #model_concentration.parameters['AKT'].value = 0.585067267459499
                
                model_population.parameters['ERK'].value = 0.273618423688386 * omega
                model_population.parameters['FAK'].value = 1.075820621939 * omega
                model_population.parameters['AKT'].value = 0.585067267459499 * omega       
        
        
        (cells, cell_count) = runAKTERK(init, run_time, model_concentration, model_population)
        IMT_eoe_counts = retrieve_cellState(run_time, init_cell, cells)
        all_events_sim = all_events(run_time, IMT_eoe_counts)
        

        proportion_sim_mitosis = len(all_events_sim[all_events_sim['event']=='mitosis'])/len(all_events_sim)
        
        proportion_sim_quiescent_cells =  len(all_events_sim[all_events_sim['event']=='quiescence'])/len(all_events_sim)
        
        proportion_sim_apoptosis = len(all_events_sim[all_events_sim['event']=='apoptosis'])/len(all_events_sim)
        
        proportion_sim_end_of_experiment = len(all_events_sim[all_events_sim['event']=='End_of_Expt'])/len(all_events_sim)
                
        
        if b==[] or a==[] or len(a)<25 or len(b)<25: #same nber of times
            return 10000
        
        n2=len(IMT_tmp)
        p2=sum(np.array(IMT_tmp)>102.3)/n2
        
        a=random.sample(a, 25)
        b=random.sample(b, 25)
        g1=random.sample(experimental_G1s[condition], 25)
        sm=random.sample(experimental_SMs[condition], 25)
        
        cost_quiescence=80/(1+math.exp(-25*abs(p1[condition]-p2)+1))
        
        n_mmd2_tot+=[nmmd2_rbf(np.array([b]).transpose(), np.array([g1]).transpose(), 1)]
        n_mmd2_tot+=[nmmd2_rbf(np.array([a]).transpose(), np.array([sm]).transpose(), 1)]
    print("ok")
    print("cost_quiescence")
    return sum(n_mmd2_tot)+cost_quiescence
