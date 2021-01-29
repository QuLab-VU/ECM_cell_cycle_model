# -*- coding: utf-8 -*-
"""
Created on Fri Jan 29 11:43:11 2021

@author: Grim.M
"""



path="/home/jeremie"
#path=r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code"
import sys
sys.path.insert(0, path + "/cellular_cycle") #add a pathway to the model
sys.path.insert(0, path + "/fit") #add a pathway 

from pydream.core import run_dream
from pysb.integrate import Solver
import numpy as np
from pydream.parameters import SampledParam
from scipy.stats import norm, uniform
import os
import inspect
from pydream.convergence import Gelman_Rubin
import pandas as pd
import pylab as pl
import sobol_seq


from cost_function_BngSimulator_frac_prolif_quiescence_and_mmd_DMSO_only_30_parameters_fixed_initial_concentrations_mitosis import  cost_function_reduced_quiescence_DMSO_only_30_parameters


from Austin_proliferation_model_population_3 import model
model_population=model


###############################################################################


def cost_function_f(parameter_vector):
    parameter_vector=[10**i for i in parameter_vector]
    return -cost_function_reduced_quiescence_DMSO_only_30_parameters(parameter_vector)



pysb_sampled_parameter_names = ["V_s2skp","V_Fp","k_M","k_EMY","k_CDF","k_CDS","k_ED","K_S","K_MY","K_DF","d_CD","k_P1", "k_P2", "K_CD", "k_E", "K_M", "K_CE", "V_2E2F", "K_2E2F", "k_RE", "kf1", "kf2", "kb", "kasbaxmbcl2", "kdsbaxmbcl2", "k1sbh3", "k2sbh3", "kdbh3", "kasbh3bcl2", "kdsbh3bcl2"]

loc=[]
for parameter in pysb_sampled_parameter_names:
    print(np.log10(model.parameters[parameter].value))
    loc+=[np.log10(model.parameters[parameter].value)]


#######
    
l=sobol_seq.i4_sobol_generate(30, 4000, skip=0)
l2 = 4*l + (np.array(loc)-2)
cost=[]

for parameter_vector in l2 :
    c = cost_function_f(parameter_vector)
    print(cost)
    cost+=[c]

np.save(path + "/fit_data/cost_sobol_seq_0_4000.npy", cost)

#######
    
l=sobol_seq.i4_sobol_generate(30, 200, skip=4000)
l2 = 4*l + (np.array(loc)-2)
cost=[]

for parameter_vector in l2 :
    c = cost_function_f(parameter_vector)
    print(cost)
    cost+=[c]

np.save(path + "\fit_data\cost_sobol_seq_4000_4200.npy", cost)

#######
    
l=sobol_seq.i4_sobol_generate(30, 200, skip=4200)
l2 = 4*l + (np.array(loc)-2)
cost=[]

for parameter_vector in l2 :
    c = cost_function_f(parameter_vector)
    print(cost)
    cost+=[c]

np.save(path + "\fit_data\cost_sobol_seq_4200_4400.npy", cost)

#######
    
l=sobol_seq.i4_sobol_generate(30, 200, skip=4400)
l2 = 4*l + (np.array(loc)-2)
cost=[]

for parameter_vector in l2 :
    c = cost_function_f(parameter_vector)
    print(cost)
    cost+=[c]

np.save(path + "\fit_data\cost_sobol_seq_4400_4600.npy", cost)

#######
    
l=sobol_seq.i4_sobol_generate(30, 200, skip=4600)
l2 = 4*l + (np.array(loc)-2)
cost=[]

for parameter_vector in l2 :
    c = cost_function_f(parameter_vector)
    print(cost)
    cost+=[c]

np.save(path + "\fit_data\cost_sobol_seq_4600_4800.npy", cost)


#######
    
l=sobol_seq.i4_sobol_generate(30, 200, skip=4800)
l2 = 4*l + (np.array(loc)-2)
cost=[]

for parameter_vector in l2 :
    c = cost_function_f(parameter_vector)
    print(cost)
    cost+=[c]

np.save(path + "\fit_data\cost_sobol_seq_4800_5000.npy", cost)


#######
    
l=sobol_seq.i4_sobol_generate(30, 200, skip=5000)
l2 = 4*l + (np.array(loc)-2)
cost=[]

for parameter_vector in l2 :
    c = cost_function_f(parameter_vector)
    print(cost)
    cost+=[c]

np.save(path + "\fit_data\cost_sobol_seq_5000_5200.npy", cost)


#######
    
l=sobol_seq.i4_sobol_generate(30, 200, skip=5200)
l2 = 4*l + (np.array(loc)-2)
cost=[]

for parameter_vector in l2 :
    c = cost_function_f(parameter_vector)
    print(cost)
    cost+=[c]

np.save(path + "\fit_data\cost_sobol_seq_5200_5400.npy", cost)


#######
    
l=sobol_seq.i4_sobol_generate(30, 200, skip=5400)
l2 = 4*l + (np.array(loc)-2)
cost=[]

for parameter_vector in l2 :
    c = cost_function_f(parameter_vector)
    print(cost)
    cost+=[c]

np.save(path + "\fit_data\cost_sobol_seq_5400_5600.npy", cost)


#######
    
l=sobol_seq.i4_sobol_generate(30, 200, skip=5600)
l2 = 4*l + (np.array(loc)-2)
cost=[]

for parameter_vector in l2 :
    c = cost_function_f(parameter_vector)
    print(cost)
    cost+=[c]

np.save(path + "\fit_data\cost_sobol_seq_5600_5800.npy", cost)


#######
    
l=sobol_seq.i4_sobol_generate(30, 200, skip=5800)
l2 = 4*l + (np.array(loc)-2)
cost=[]

for parameter_vector in l2 :
    c = cost_function_f(parameter_vector)
    print(cost)
    cost+=[c]

np.save(path + "\fit_data\cost_sobol_seq_5800_6000.npy", cost)


#######
    
l=sobol_seq.i4_sobol_generate(30, 200, skip=6000)
l2 = 4*l + (np.array(loc)-2)
cost=[]

for parameter_vector in l2 :
    c = cost_function_f(parameter_vector)
    print(cost)
    cost+=[c]

np.save(path + "\fit_data\cost_sobol_seq_6000_6200.npy", cost)


#######
    
l=sobol_seq.i4_sobol_generate(30, 200, skip=6200)
l2 = 4*l + (np.array(loc)-2)
cost=[]

for parameter_vector in l2 :
    c = cost_function_f(parameter_vector)
    print(cost)
    cost+=[c]

np.save(path + "\fit_data\cost_sobol_seq_6200_6400.npy", cost)


#######
    
l=sobol_seq.i4_sobol_generate(30, 200, skip=6400)
l2 = 4*l + (np.array(loc)-2)
cost=[]

for parameter_vector in l2 :
    c = cost_function_f(parameter_vector)
    print(cost)
    cost+=[c]

np.save(path + "\fit_data\cost_sobol_seq_6400_6600.npy", cost)


#######
    
l=sobol_seq.i4_sobol_generate(30, 200, skip=6600)
l2 = 4*l + (np.array(loc)-2)
cost=[]

for parameter_vector in l2 :
    c = cost_function_f(parameter_vector)
    print(cost)
    cost+=[c]

np.save(path + "\fit_data\cost_sobol_seq_6600_6800.npy", cost)


#######
    
l=sobol_seq.i4_sobol_generate(30, 200, skip=6800)
l2 = 4*l + (np.array(loc)-2)
cost=[]

for parameter_vector in l2 :
    c = cost_function_f(parameter_vector)
    print(cost)
    cost+=[c]

np.save(path + "\fit_data\cost_sobol_seq_6800_7000.npy", cost)


#######
    
l=sobol_seq.i4_sobol_generate(30, 200, skip=7000)
l2 = 4*l + (np.array(loc)-2)
cost=[]

for parameter_vector in l2 :
    c = cost_function_f(parameter_vector)
    print(cost)
    cost+=[c]

np.save(path + "\fit_data\cost_sobol_seq_7000_7200.npy", cost)


#######
    
l=sobol_seq.i4_sobol_generate(30, 200, skip=7200)
l2 = 4*l + (np.array(loc)-2)
cost=[]

for parameter_vector in l2 :
    c = cost_function_f(parameter_vector)
    print(cost)
    cost+=[c]

np.save(path + "\fit_data\cost_sobol_seq_7200_7400.npy", cost)


