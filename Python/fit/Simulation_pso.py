# -*- coding: utf-8 -*-
"""
Created on Mon Dec 19 11:59:15 2022

@author: Grim.M
"""

import time

import re
import os
import sys
path= r'C:\Users\Grim.M\Documents\git_repos\ECM_cell_cycle_model'

path= os.getcwd()

sys.path.insert(0, path + "/Python/cellular_cycle") #add a pathway to the model #not sure about how to use folder dependency, to change
#sys.path.insert(0, path + "/Python/test_treatments") #add a pathway to the model #not sure about how to use folder dependency, to change
sys.path.insert(0, path + "/Python/fit") #add a pathway to the fit functions #not sure about how to use folder dependency, to change

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import random
import sobol_seq

##Importing the model, models differ on their initial conditions, one being a concentration and the other a protein count
from Austin_proliferation_model_concentration_5 import model
model_concentration=model
from Austin_proliferation_model_population_5 import model
model_population=model

from cost_function_BngSimulator_low_medium_stiffness_pso import cost_function_30_parameters_pso
#from cost_function_BngSimulator_low_medium_stiffness_pso_2 import cost_function_30_parameters_pso_2

##Import the simulator: for ODE and stochastic simulation :
from pysb.simulator import BngSimulator
from pysb.simulator import ScipyOdeSimulator

# Directory containing cell fate data
from pathlib import Path
pysb_sampled_parameter_names = ["V_s2skp","V_Fp","k_M","k_EMY","k_CDF","k_CDS","k_ED","K_S","K_MY","K_DF","d_CD","k_P1", "k_P2", "K_CD", "k_E", "K_M", "K_CE", "V_2E2F", "K_2E2F", "k_RE", "kf1", "kf2", "kb", "kasbaxmbcl2", "kdsbaxmbcl2", "k1sbh3", "k2sbh3", "kdbh3", "kasbh3bcl2", "kdsbh3bcl2"]



print(model_population.parameters["k_M"].value)
print(model_population.parameters["k_EMY"].value)
print(model_population.parameters["k_CDS"].value)
key = "k_M" #0.11 to 0.27
key = "k_EMY" #from 0.18 to 0.37
key = "k_CDS" #from 0.11 to 0.25
print(model_population.parameters[key].value)

#change the model parameters of the model if needed
omega = 200
new_param = pd.read_csv(path + "/umap_window/parameters.csv")
del new_param['Stress']
for key in new_param.keys():
    if model_population.parameters[key].value != new_param[key][0]:
        model_population.parameters[key].value = new_param[key]
    if model_concentration.parameters[key].value != new_param[key][0]:
        model_concentration.parameters[key].value = new_param[key]

model_concentration.parameters['omega'].value = 1
model_population.parameters['omega'].value

loc=[]
for parameter in pysb_sampled_parameter_names:
    #print((model.parameters[parameter].value))
    loc+=[np.log10(model.parameters[parameter].value)]


t = np.linspace(0, 500, 50000)
#simres = ScipyOdeSimulator(model_concentration, t, verbose=False).run();
#sim = BngSimulator(model_population, tspan=t, verbose=True).run()


original_parameters = []
for name in pysb_sampled_parameter_names:
    original_parameters += [np.log10(model_concentration.parameters[name].value)]



parameter_vector = 10**np.array(original_parameters)
idx = 'original_parameters'

one = time.time()
cost_function_30_parameters_pso_2(parameter_vector, 'original_parameters')
two = time.time()

three = time.time()
cost_function_30_parameters_pso(parameter_vector, 'original_parameters')
four = time.time()

print(four - three)
print(two - one)



#solved trouble shooting
"""

param_dict = {pname: pvalue for pname, pvalue in zip(pysb_sampled_parameter_names, parameter_vector)}   
for pname, pvalue in param_dict.items():
    #Change model parameter values to current location in parameter space
    model_population.parameters[pname].value = (pvalue)
    model_concentration.parameters[pname].value = (pvalue)
    
    print(model_population.parameters[pname].value )
    print(model_concentration.parameters[pname].value )
    print(pvalue)



conditions = ['0.6 kPa DMSO', '0.6 kPa ERL', '4 kPa DMSO', '4 kPa ERL']
inputs = {}
inputs['FAK'] = {'4 kPa DMSO': 0.6, '4 kPa ERL': 0.6, '0.6 kPa DMSO': 0.4, '0.6 kPa ERL': 0.4}
inputs['ERK'] = {'4 kPa DMSO': 1.0, '4 kPa ERL': 0.2, '0.6 kPa DMSO': 1.0, '0.6 kPa ERL': 0.2}
inputs['AKT'] = {'4 kPa DMSO': 1.0, '4 kPa ERL': 0.6, '0.6 kPa DMSO': 0.6, '0.6 kPa ERL': 0.3}

condition = '0.6 kPa ERL'
num_condition=0
parameters_sim={} 

t_sim = np.arange(0,120,1/3)

run_time = observed_cells_time = t_sim[-1]
omega=200
sigma=5
parameter_names = ['a','d','q']
init=1 #40
idx= 'test_p'

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


cost = {}
from runAKTERK_2 import runAKTERK
(cells, cell_count) = runAKTERK(init, run_time, observed_cells_time, model_concentration, model_population, idx, condition);


thr=0.26
sys_size=200
t = np.linspace(0, 1000, 1000)

from peakfinder import peakfinder


simres_2 = BngSimulator(model_population, t, verbose=False).run();

(a,max_tmp_cdc) = peakfinder(simres_2.observables['obsCdc20'], thr*sys_size, thr*sys_size, 1, includeEndpoints=1);
(b,max_tmp_cyce) = peakfinder(simres_2.observables['obsCycE'], thr*sys_size, thr*sys_size, 1, includeEndpoints=1);

max_tmp_cdc
"""

import pyswarms as ps


all_parameters_4kPa = pd.read_csv(path + '/Python/fit_data/sobol_seq/all_parameters_4kPa.csv', index_col='index')

all_parameters_4kPa = all_parameters_4kPa.sort_values(by=['cost 4kPa'])

#parameter_vector = np.log10(all_parameters_4kPa[all_parameters_4kPa.index == all_parameters_4kPa.index.iloc[9]][pysb_sampled_parameter_names].to_numpy())[0]


initial_positions = np.log10(all_parameters_4kPa.iloc[0:30][pysb_sampled_parameter_names].to_numpy())



def cost_function_f(parameter_vectors):
    cost = np.zeros(np.shape(parameter_vectors)[0])
    for i in range(len(cost)) :
        try :
            parameter_vector = parameter_vectors[i]
            print(parameter_vector)
            parameter_vector=[10**i for i in parameter_vector]
            cost[i] = + cost_function_30_parameters_pso(parameter_vector, i)
        except :
            cost[i] = + np.inf
    return cost


from multiprocessing.pool import ThreadPool as Pool


def cost_function_f_2(parameter_vectors):        
    print('letsgo')
    with Pool() as pool:
        #freeze_support()
        cost = pool.starmap(cost_function_30_parameters_pso, zip(parameter_vectors, range(len(parameter_vectors))))   
    return np.array(cost)


max_bound = np.array(loc)+2.5
min_bound = np.array(loc)-2.5
bounds = (min_bound, max_bound)




initial_positions[0:30]<max_bound
initial_positions[0:30]>min_bound

#options = {'c1': 0.5, 'c2': 0.3, 'w': 0.9, 'k': 3, 'p': 2} #hyperparameters 1
#options = {'c1': 0.05, 'c2': 0.1, 'w': 0.1, 'k': 6, 'p': 2} #hyperparameters 2
#options = {'c1': 0.1, 'c2': 0.3, 'w': 0.3, 'k': 5, 'p': 2} #hyperparameters 3
#options = {'c1': 0.2, 'c2': 0.3, 'w': 0.4, 'k': 5, 'p': 2} #hyperparameters 4
options = {'c1': 0.15, 'c2': 0.3, 'w': 0.4, 'k': 5, 'p': 2} #hyperparameters 5

#position = np.load("/home/jeremie/fit_data/pso/high_and_low_stiffness_cost/history/position_history.npy")


# Call instance of PSO
optimizer = ps.single.LocalBestPSO(n_particles=30, dimensions=30, options=options, bounds=bounds, init_pos = initial_positions[0:30])
#optimizer = ps.single.LocalBestPSO(n_particles=30, dimensions=30, options=options, bounds=bounds, init_pos = np.array(position[-1,:,:]))

# Provide your own init_pos
#optimizer.pos = np.array(initial_positions[0:10]).copy()

# Perform optimization
cost, pos = optimizer.optimize(cost_function_f_2, iters=120)



np.save(path + "/Python/fit_data/pso/history/position_history.npy", optimizer.pos_history)

np.save(path + "/fit_data/pso/history/cost_history.npy", optimizer.cost_history)


