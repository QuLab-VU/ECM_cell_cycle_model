# -*- coding: utf-8 -*-
"""
Created on Tue Apr 21 11:07:42 2020

@author: Grim.M
"""

path="/home/jeremie"
#path=r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code"

from pysb.integrate import Solver
import numpy as np

from scipy.stats import norm, uniform
import os
import inspect

import pandas as pd
import pylab as pl

from simplepso.pso import PSO


import sys
sys.path.insert(0, path + "/cellular_cycle") #add a pathway to the model
sys.path.insert(0, path + "/fit") #add a pathway 


import math

from cost_function import cost_function

'''
#Initialize PySB solver object for running simulations.  Simulation timespan should match experimental data.
tmp_t = np.linspace(0, 1000, 10000)
sim = BngSimulator(model, tspan=tmp_t, verbose=True)
simres = sim.run(initials=Initial)
'''

#Load experimental data to which the model will be fit here.  The "experimental data" in this case is the IMTs distribution

experimental_SMs={}
experimental_G1s={}
experimental_IMTs={}



df = pd.read_csv(path+'/IMT_times/DS8_SM_0_6_4_55.csv')

for col in df.columns :
    experimental_SMs[col]=df[col].dropna().tolist()


df = pd.read_csv(path+'/IMT_times/DS8_G1_0_6_4_55.csv')

for col in df.columns :
    experimental_G1s[col]=df[col].dropna().tolist()

for k in experimental_G1s.keys() :
    minimum=min(len(experimental_G1s[k]),len(experimental_SMs[k]))
    experimental_IMTs[k]=[experimental_G1s[k][i]+experimental_SMs[k][i] for i in range(minimum)]



#Create lists of sampled pysb parameter names to use for subbing in parameter values in likelihood function.
pysb_sampled_parameter_names = ["V_s2skp","V_Fp","k_M","k_EMY","k_CDF","k_CDS","k_ED","K_S","K_MY","K_DF"]

# Add vector of PySB rate parameters to be sampled as unobserved random variables to DREAM with uniform priors.

original_params = [0.01 , 1.5 , 0.07 , 0.3 , 0.22 , 0.21 , 0.22 , 0.8 , 0.1 , 0.8 , 0.009]

#Define likelihood function to generate simulated data that corresponds to experimental time points.  
#This function should take as input a parameter vector (parameter values are in the order dictated by first argument to run_dream function below).
#The function returns a log probability value for the parameter vector given the experimental data.




#cost_function([1, 1, 1, 1, 1, 1, 1, 1, 1, 1])






pso = PSO(save_sampled=False, verbose=True)
pso.set_cost_function(cost_function)

pso.set_start_position(np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1]))

#pso.set_start_position(np.array([0.01 , 1.5 , 0.07 , 0.3 , 0.22 , 0.21 , 0.22 , 0.8 , 0.1 , 0.8 ]))
# allows particles to move +/- 2 orders of magnitude
#pso.set_bounds(np.array([0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001]))

pso.set_bounds(np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1]))

# sets maximum speed that a particle can travel

pso.set_speed(-.1, .1)

pso.run(num_particles=15, num_iterations=40, stop_threshold=1e-5)










