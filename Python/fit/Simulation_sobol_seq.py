# -*- coding: utf-8 -*-
"""
Created on Sat Nov  5 18:06:49 2022

@author: Grim.M
"""


import os
import sys
path= r'C:\Users\Grim.M\Documents\git_repos\ECM_cell_cycle_model\Python'

sys.path.insert(0, path + "/cellular_cycle") #add a pathway to the model #not sure about how to use folder dependency, to change
sys.path.insert(0, path + "/fit") #add a pathway to the fit functions #not sure about how to use folder dependency, to change

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

from cost_function_BngSimulator_frac_prolif_quiescence_and_mmd_ERL_only_30_parameters_fixed_initial_concentrations_mitosis_low_high_stiffness import cost_function_reduced_quiescence_ERL_only_30_parameters

##Import the simulator: for ODE and stochastic simulation :
from pysb.simulator import BngSimulator
from pysb.simulator import ScipyOdeSimulator

# Directory containing cell fate data
from pathlib import Path

pysb_sampled_parameter_names = ["V_s2skp","V_Fp","k_M","k_EMY","k_CDF","k_CDS","k_ED","K_S","K_MY","K_DF","d_CD","k_P1", "k_P2", "K_CD", "k_E", "K_M", "K_CE", "V_2E2F", "K_2E2F", "k_RE", "kf1", "kf2", "kb", "kasbaxmbcl2", "kdsbaxmbcl2", "k1sbh3", "k2sbh3", "kdbh3", "kasbh3bcl2", "kdsbh3bcl2"]
loc=[]
for parameter in pysb_sampled_parameter_names:
    #print((model.parameters[parameter].value))
    loc+=[np.log10(model.parameters[parameter].value)]

print(model_population.parameters["k_M"].value)
print(model.parameters["k_M"].value)
key = "k_M" #0.11 to 0.27
key = "k_EMY" #from 0.18 to 0.37
key = "k_CDS" #from 0.11 to 0.25

#change the model parameters of the model if needed
omega = 200
new_param = pd.read_csv(path + "/../umap_window/parameters.csv") 
del new_param['Stress']
for key in new_param.keys():
    if model_population.parameters[key].value != new_param[key][0]:
        model_population.parameters[key].value = new_param[key]
    if model_concentration.parameters[key].value != new_param[key][0]:
        model_concentration.parameters[key].value = new_param[key]
model_concentration.parameters['omega'].value = 1


l=sobol_seq.i4_sobol_generate(30, 110000, skip=0)
l2 = 4*l + (np.array(loc)-2)

t = np.linspace(0, 500, 50000)
simres = ScipyOdeSimulator(model_concentration, t, verbose=False).run(); #deterministic simulation


frac_prolif_parameters = pd.read_pickle(path + "/fit_data/sobol_seq/frac_prolif_parameters_all_parameters_006_55_kPa.pkl")



#indices in the sobol sequence for the 100 best parameters
indices = [int(i)-1 for i in frac_prolif_parameters.sort_values(by=['cost_total'])["index"].tolist()[0:100]] 
parameters = [l2[i] for i in indices]

#Change model parameter values to current location in parameter space
parameter_vector = 10**parameters[0]
idx = 'test'



cost_function_reduced_quiescence_ERL_only_30_parameters(parameter_vector, idx)











