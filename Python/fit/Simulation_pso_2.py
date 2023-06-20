# -*- coding: utf-8 -*-
"""
Created on Sat May 13 16:35:22 2023

@author: Grim.M
"""


path="/home/jeremie"
#path=r"C:\Users\Grim.M\Documents\git_repos\ECM_cell_cycle_model"

import sys
sys.path.insert(0, path+"/Python/cellular_cycle") #add a pathway to the model
sys.path.insert(0, path+"/Python/fit") #add a pathway to the model

from Austin_proliferation_model_concentration_6 import model
model_concentration = model
from Austin_proliferation_model_population_6 import model
model_population = model

from cost_function_BngSimulator_low_medium_stiffness_pso import log_likelihood, log_prior, log_posterior

import numpy as np
from scipy.stats import multivariate_normal


###############################################################################

pysb_sampled_parameter_names = ["V_s2skp","V_Fp","k_M","k_EMY","k_CDF","k_CDS","k_ED","K_S","K_MY","K_DF","d_CD","k_P1", "k_P2", "K_CD", "k_E", "K_M", "K_CE", "V_2E2F", "K_2E2F", "k_RE", "kf1", "kf2", "kb", "kasbaxmbcl2", "kdsbaxmbcl2", "k1sbh3", "k2sbh3", "kdbh3", "kasbh3bcl2", "kdsbh3bcl2", 'v_s1p21', 'F_n', 'K_p21fak', 'a_stress', 'b_stress_AKT', 'b_stress_ERK']

###############################################################################
#Changing the parameters:

import pandas as pd
new_param = pd.read_csv(path + "/Matlab/v2/parameters.csv")
del new_param['Stress']
for key in new_param.keys():
    try:
        if model_population.parameters[key].value != new_param[key][0]:
            model_population.parameters[key].value = new_param[key]
        if model_concentration.parameters[key].value != new_param[key][0]:
            model_concentration.parameters[key].value = new_param[key]
    except: print(key)

model_population.parameters['omega'].value = 200
model_concentration.parameters['omega'].value = 1


model_concentration.parameters['v_s1p21'].value
model_concentration.parameters['v_s1p21'].value
# Setting initial concentrations

species_initials = ['Skp2_0', 'p21_0', 'p21p_0', 'cdh1act_0', 'cdh1inact_0', 'Myc_0', 'CycD_0', 'CycDp21_0', 'pRb_0', 'pRbp_0', 'E2F_0', 'pRbE2F_0', 'CycE_0', 'CycEp21_0', 'CycA_0', 'CycAp21_0', 'CycB_0', 'CycBp21_0', 'Cdc20_0', 'BAX_0', 'BAXmBCL2_0', 'BH3_0', 'BH3BCL2_0', 'BAXm_0', 'BCL2_0']
species_initials_values = [0.0199387776457489, 0.0212324326339528, 0.104907517295920, 0.970596980815461, 1-0.970596980815461, 0.321067821067829, 0.175953608509162, 0.0138939800293163, 0.301757475166483, 0.545620559818762, 0.0403080965235495, 0.554487916752975, 0.0640517424993558, 0.00487528808175496, 0.0182005140710040, 0.00234453838504852, 0.0234757871244827, 0.00435874152378834, 0.0204771947463057, 0.999000159839911, 0.000997752220497874, 5.81398691338436e-05, 0.0332751934642001, 1 - 0.999000159839911 - 0.000997752220497874, 0.95 - 0.0332751934642001 - 0.000997752220497874]

for species_initial, species_initials_value in zip(species_initials, species_initials_values):
    model_concentration.parameters[species_initial].value = species_initials_value
    model_population.parameters[species_initial].value = int(round(species_initials_value * model_population.parameters['omega'].value))

###############################################################################


original_parameters = []
for name in pysb_sampled_parameter_names:
    original_parameters += [np.log10(model_concentration.parameters[name].value)]

parameter_vector = 10**np.array(original_parameters)
idx = 'original_parameters'

mean = 10**np.array(original_parameters)

###############################################################################

conditions = ['0.6 kPa DMSO', '0.6 kPa ERL', '4 kPa DMSO', '4 kPa ERL', '55 kPa DMSO', '55 kPa ERL']
run_times = [120, 120, 120, 120, 60, 60]
init_cells = range(1, 11)


for condition, run_time_tmp in zip(conditions, run_times):
    for init_cell_tmp in init_cells:
        name = condition.replace(' ', '_') + "_" + str(init_cell_tmp) + "_initial_cells"
        print(name)
                
        for iteration in range(10):
            c = posterior(parameter_vector, 'test_initial_cell_number_' + name , pysb_sampled_parameter_names, [condition], [run_time_tmp], [init_cell_tmp],  [10], True, iteration)




###############################################################################

mean = np.log10(np.array([1.000e-02, 1.500e+00, 2.700e-01, 3.700e-01, 2.100e-01, 2.500e-01, 2.000e-01, 8.000e-01, 1.000e-01, 8.000e-01, 1.500e+00, 1.800e+01, 1.800e+01, 5.500e-01, 4.000e-01, 2.000e-01, 5.500e-01, 7.000e-01, 5.000e-01, 1.800e+02, 1.000e-02, 2.400e+00, 1.400e-01, 7.500e+01, 5.000e-02, 1.000e-03, 2.500e-02, 3.000e-02, 2.500e+01, 1.000e-02, 7.000e-01, 8.000e+00, 5.000e+00, 6.815e-01, 1.120e-01, 1.110e-02]))
variance = 0.005 * np.ones(len(mean))
# Create a multivariate normal distribution with the given mean and variance
mvn = multivariate_normal(mean, np.diag(variance))

# Generate a random sample from the distribution
sample = mvn.rvs(size=10)

log_prior(np.log10(parameter_vector))
[log_prior(sample[i]) for i in range(sample.shape[0])]

parameter_vector
10**sample
###############################################################################

def cost_function(parameter_vectors):
    pysb_sampled_parameter_names = ["V_s2skp","V_Fp","k_M","k_EMY","k_CDF","k_CDS","k_ED","K_S","K_MY","K_DF","d_CD","k_P1", "k_P2", "K_CD", "k_E", "K_M", "K_CE", "V_2E2F", "K_2E2F", "k_RE", "kf1", "kf2", "kb", "kasbaxmbcl2", "kdsbaxmbcl2", "k1sbh3", "k2sbh3", "kdbh3", "kasbh3bcl2", "kdsbh3bcl2", 'v_s1p21', 'F_n', 'K_p21fak', 'a_stress', 'b_stress_AKT', 'b_stress_ERK']
    conditions = ['0.6 kPa DMSO', '0.6 kPa ERL', '4 kPa DMSO', '4 kPa ERL', '55 kPa DMSO', '55 kPa ERL']
    run_times = [120, 120, 120, 120, 60, 60]
    init_cells = (5*np.ones(6)).astype(int)
    chain = 0
    save=True
    
    global iteration_swarm
    iteration_swarm += 1
    
    cost = np.zeros(np.shape(parameter_vectors)[0])
    
    for chain in range(len(cost)) :
        parameter_vector = parameter_vectors[chain]
        print(parameter_vector) 
        
        cost[chain] =+ log_posterior(parameter_vector, chain+1, pysb_sampled_parameter_names, conditions, run_times, init_cells, save, iteration_swarm)
        
    return cost

iteration_swarm = 0
c = cost_function(sample)

import time

start = time.time()
c_original = cost_function(np.log10(np.array([parameter_vector])))
end = time.time()

duration_one_cell = end - start


###############################################################################

import pyswarms as ps
options = {'c1': 0.1, 'c2': 0.1, 'w': 0.1, 'k': 4, 'p': 2} #hyperparameters 4

n_particles = len(pysb_sampled_parameter_names)
dimensions = len(pysb_sampled_parameter_names)
initial_positions = mvn.rvs(size=n_particles)


max_bound = np.array(mean)+2
min_bound = np.array(mean)-2
bounds = (min_bound, max_bound)
iteration_swarm = 0

optimizer = ps.single.LocalBestPSO(n_particles=n_particles, dimensions=dimensions, options=options, bounds=bounds, init_pos = initial_positions)

cost, pos = optimizer.optimize(cost_function, iters=80)


