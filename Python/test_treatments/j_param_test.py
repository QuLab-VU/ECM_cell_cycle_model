# -*- coding: utf-8 -*-
"""
Created on Tue Dec  6 01:36:21 2022

@author: olesk
"""

import os
import sys
sys.path.insert(0, os.getcwd() + "/cellular_cycle") #add a pathway to the model #not sure about how to use folder dependency, to change
sys.path.insert(0, os.getcwd() + "/fit") #add a pathway to the fit functions #not sure about how to use folder dependency, to change
sys.path.insert(0, os.getcwd() + "/fit_data/sobol_seq")

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

from pysb.simulator import BngSimulator
from pysb.simulator import ScipyOdeSimulator

params = pd.read_pickle("fit_data/sobol_seq/parameters_high_and_low_stiffness_ERL.pkl")

p_sort = params.sort_values(by=['cost_total_ERL'])
w_death = p_sort[p_sort["a_006_ERL"] > 0]
#p_sort = p_sort.iloc[:,10:]
parameter_names = ["V_s2skp","V_Fp","k_M","k_EMY","k_CDF","k_CDS","k_ED","K_S","K_MY","K_DF","d_CD","k_P1", "k_P2", "K_CD", "k_E", "K_M", "K_CE", "V_2E2F", "K_2E2F", "k_RE", "kf1", "kf2", "kb", "kasbaxmbcl2", "kdsbaxmbcl2", "k1sbh3", "k2sbh3", "kdbh3", "kasbh3bcl2", "kdsbh3bcl2"]

ii=0
for key in parameter_names:
    if model_population.parameters[key].value != p_sort[key].iloc[ii]:
        model_population.parameters[key].value = p_sort[key].iloc[ii]
    if model_concentration.parameters[key].value != p_sort[key].iloc[ii]:
        model_concentration.parameters[key].value = p_sort[key].iloc[ii]

"""
ii=4
for key in parameter_names:
    if model_population.parameters[key].value != w_death[key].iloc[ii]:
        model_population.parameters[key].value = w_death[key].iloc[ii]
    if model_concentration.parameters[key].value != w_death[key].iloc[ii]:
        model_concentration.parameters[key].value = w_death[key].iloc[ii]
"""

model_concentration.parameters['omega'].value = 1
omega=200

model_concentration.parameters['ERK'].value = 0.2 * 1
model_concentration.parameters['FAK'].value = 0.4 * 1
model_concentration.parameters['AKT'].value = 0.3 * 1

model_population.parameters['ERK'].value = 0.2 * omega
model_population.parameters['FAK'].value = 0.4 * omega
model_population.parameters['AKT'].value = 0.3 * omega


t = np.linspace(0, 500, 50000)
simres = ScipyOdeSimulator(model_concentration, t, verbose=False).run()
yout = simres.all;
observables=['obsCycD', 'obsCycE', 'obsCycA', 'obsCycB', 'obsCdc20']; #plot Cyclin E, Cdc20 and Cyclin A
for obs in range(5) :
    plt.plot(t[:],yout[observables[obs]][:], label=observables[obs], lw=2)
plt.legend(loc=1)
plt.ylim(0, 1)
plt.xlim(300, 500)
plt.show()

