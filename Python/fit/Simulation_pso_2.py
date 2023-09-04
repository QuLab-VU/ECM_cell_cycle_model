# -*- coding: utf-8 -*-
"""
Created on Sat May 13 16:35:22 2023

@author: Grim.M
"""

import os
path = os.getcwd()
import sys
sys.path.insert(0, path + "/../cellular_cycle") #add a pathway to the model
sys.path.insert(0, path + "/../fit") #add a pathway to the model

from cell_cycle_death_model import model
from cost_function_BngSimulator_low_medium_stiffness_pso import log_likelihood, log_prior, log_posterior
import numpy as np
from scipy.stats import multivariate_normal
import pandas as pd
from pysb.simulator import ScipyOdeSimulator, BngSimulator
import matplotlib.pyplot as plt

###############################################################################

# import pysb sampled parameters
pysb_sampled_parameter_names = np.loadtxt(path + "/../model_data/sampled_parameter_names.txt", dtype=str)

# Importing the parameters:
new_param = pd.read_csv(path + "/../model_data/parameters.csv")
del new_param['Stress']
for key in new_param.keys():
    try:
        if model.parameters[key].value != new_param[key][0]:
            model.parameters[key].value = new_param[key]
    except: print(key)
model.parameters['omega'].value = 200

# Setting initial concentrations '/../model_data/initial_concentrations.csv'
initial_concentrations = pd.read_csv(path + '/../model_data/initial_concentrations.csv')
for key in initial_concentrations.keys():
    try:
        model.parameters[key].value = initial_concentrations[key][0]/model.parameters['omega'].value
    except: print(key)


###############################################################################
# Setting original parameters

original_parameters = []
for name in pysb_sampled_parameter_names:
    original_parameters += [np.log10(model.parameters[name].value)]
parameter_vector = 10**np.array(original_parameters)
idx = 'original_parameters'
mean = 10**np.array(original_parameters)


# Test simulation with original parameters
time = np.linspace(0, 100, 10001)
sim_result_ODE = ScipyOdeSimulator(model, time, verbose=True).run()
sim_result_SSA = BngSimulator(model, time, verbose=True).run(seed=33)

observables=['obsCycE','obsCdc20', 'obsCycA', 'obsCycB', 'obsp21']
for i, obs in enumerate(observables):
    print(sim_result_ODE.observables[obs])
    print(len(sim_result_ODE.observables[obs]))
    print(time.shape)
    plt.plot(time, sim_result_ODE.observables['obsCdc20'])
    # plt.plot(time, sim_result_ODE.observables[i], label=obs)
plt.legend()
plt.show()

for i, obs in enumerate(observables):
    plt.plot(time, sim_result_SSA.observables[obs], label=obs)
plt.legend()
plt.show()



###############################################################################

conditions = ['0.6 kPa DMSO', '0.6 kPa ERL', '4 kPa DMSO', '4 kPa ERL', '55 kPa DMSO', '55 kPa ERL']
run_times = [120, 120, 120, 120, 60, 60]
init_cells = range(1, 11)


# for condition, run_time_tmp in zip(conditions, run_times):
#     for init_cell_tmp in init_cells:
#         name = condition.replace(' ', '_') + "_" + str(init_cell_tmp) + "_initial_cells"
#         print(name)
#
#         for iteration in range(10):
#             c = posterior(parameter_vector, 'test_initial_cell_number_' + name , pysb_sampled_parameter_names, [condition], [run_time_tmp], [init_cell_tmp],  [10], True, iteration)




###############################################################################

mean = np.log10(parameter_vector)
variance = 0.005 * np.ones(len(mean))
# Create a multivariate normal distribution with the given mean and variance
mvn = multivariate_normal(mean, np.diag(variance))

# Generate a random sample from the distribution
sample = mvn.rvs(size=10)

log_prior(parameter_vector, mean, variance)
[log_prior(sample[i], mean, variance) for i in range(sample.shape[0])]

parameter_vector
10**sample
###############################################################################

conditions = ['0.6 kPa DMSO', '0.6 kPa ERL', '4 kPa DMSO', '4 kPa ERL', '55 kPa DMSO', '55 kPa ERL']
run_times = [120, 120, 120, 120, 60, 60]
init_cells = (5*np.ones(6)).astype(int)
chain = 0
save=False


cost = log_posterior(np.log10(np.array(parameter_vector)), 'original_parameters', pysb_sampled_parameter_names, mean, variance, conditions, run_times, init_cells, save, 'original_parameters')

###############################################################################

def cost_function(parameter_vectors):
    global mean, variance, pysb_sampled_parameter_names
    conditions = ['0.6 kPa DMSO', '0.6 kPa ERL', '4 kPa DMSO', '4 kPa ERL', '55 kPa DMSO', '55 kPa ERL']
    run_times = [120, 120, 120, 120, 60, 60]
    init_cells = (5*np.ones(6)).astype(int)
    chain = 0
    save=False
    
    global iteration_swarm
    iteration_swarm += 1
    
    cost = np.zeros(np.shape(parameter_vectors)[0])
    
    for chain in range(len(cost)) :
        parameter_vector = parameter_vectors[chain]
        print(parameter_vector) 
        
        cost[chain] =+ log_posterior(parameter_vector, chain+1, pysb_sampled_parameter_names, mean, variance, conditions, run_times, init_cells, save, iteration_swarm)
        
    return cost

iteration_swarm = 0
c = cost_function(sample)

import time
#
# start = time.time()
# c_original = cost_function(np.log10(np.array([parameter_vector])))
# end = time.time()
#
# duration_one_cell = end - start
#
#
# ###############################################################################
#
# import pyswarms as ps
# options = {'c1': 0.1, 'c2': 0.1, 'w': 0.1, 'k': 4, 'p': 2} #hyperparameters 4
#
# n_particles = len(pysb_sampled_parameter_names)
# dimensions = len(pysb_sampled_parameter_names)
# initial_positions = mvn.rvs(size=n_particles)
#
#
# max_bound = np.array(mean)+2
# min_bound = np.array(mean)-2
# bounds = (min_bound, max_bound)
# iteration_swarm = 0
#
# optimizer = ps.single.LocalBestPSO(n_particles=n_particles, dimensions=dimensions, options=options, bounds=bounds, init_pos = initial_positions)
#
# cost, pos = optimizer.optimize(cost_function, iters=80)
#
#
