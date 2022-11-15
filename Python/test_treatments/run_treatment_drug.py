# -*- coding: utf-8 -*-
"""
Created on Sat Oct 15 16:37:38 2022

@author: olesk
"""

import os
import sys
sys.path.insert(0, os.getcwd() + "/../cellular_cycle") #add a pathway to the model #not sure about how to use folder dependency, to change
sys.path.insert(0, os.getcwd() + "/../fit") #add a pathway to the fit functions #not sure about how to use folder dependency, to change

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import random
import sobol_seq

##Importing the model, models differ on their initial conditions, one being a concentration and the other a protein count
from Austin_proliferation_model_concentration_6 import model
model_concentration=model
from Austin_proliferation_model_population_6 import model
model_population=model

##Import the simulator: for ODE and stochastic simulation :
from pysb.simulator import BngSimulator
from pysb.simulator import ScipyOdeSimulator


params = {'axes.titley': 1,
         'axes.titlepad': -14,
         'axes.labelsize': 14,
         'axes.titlesize': 14,
         'xtick.labelsize': 12,
         'ytick.labelsize': 12}
plt.rcParams.update(params)

#change the model parameters of the model if needed
omega = 200
new_param = pd.read_csv("parameters_drugB.csv") 
del new_param['Stress']
for key in new_param.keys():
    if model_population.parameters[key].value != new_param[key][0]:
        model_population.parameters[key].value = new_param[key]
    if model_concentration.parameters[key].value != new_param[key][0]:
        model_concentration.parameters[key].value = new_param[key]
model_concentration.parameters['omega'].value = 1

t = np.linspace(0, 500, 50000)
simres = ScipyOdeSimulator(model_concentration, t, verbose=False).run()

model_population.parameters['ERK'].value = 0.2 * omega
model_population.parameters['FAK'].value = 1.0 * omega
model_population.parameters['AKT'].value = 0.6 * omega

model_concentration.parameters['ERK'].value = 0.2 * 1
model_concentration.parameters['FAK'].value = 1 * 1
model_concentration.parameters['AKT'].value = 0.6 * 1


omega = 200
new_param = pd.read_csv("parameters_drugD.csv") 
del new_param['Stress']
for key in new_param.keys():
    if model_population.parameters[key].value != new_param[key][0]:
        model_population.parameters[key].value = new_param[key]
    if model_concentration.parameters[key].value != new_param[key][0]:
        model_concentration.parameters[key].value = new_param[key]
model_concentration.parameters['omega'].value = 1
t = np.linspace(0, 500, 50000)
simres = ScipyOdeSimulator(model_concentration, t, verbose=False).run(); #deterministic simulation
yout = simres.all;
observables=['obsdrug', 'obsCycE', 'obsCdc20']; #plot Cyclin E, Cdc20 and Cyclin A
for obs in range(3) :
    plt.plot(t[:],yout[observables[obs]][:], label=observables[obs], lw=2)
plt.legend(loc=1)
plt.ylim(0, 1)
plt.show()


omega = 200
new_param = pd.read_csv("parameters_drugE.csv") 
del new_param['Stress']
for key in new_param.keys():
    if model_population.parameters[key].value != new_param[key][0]:
        model_population.parameters[key].value = new_param[key]
    if model_concentration.parameters[key].value != new_param[key][0]:
        model_concentration.parameters[key].value = new_param[key]
model_concentration.parameters['omega'].value = 1
t = np.linspace(0, 500, 50000)
simres = ScipyOdeSimulator(model_concentration, t, verbose=False).run(); #deterministic simulation
yout = simres.all;
observables=['obsdrug', 'obsCycE', 'obsCdc20', 'obsCycEdrug']; #plot Cyclin E, Cdc20 and Cyclin A
for obs in range(4) :
    plt.plot(t[:],yout[observables[obs]][:], label=observables[obs], lw=2)
plt.legend(loc=1)
plt.ylim(0, 1)
plt.show()


omega = 200
new_param = pd.read_csv("parameters_drugA.csv") 
del new_param['Stress']
for key in new_param.keys():
    if model_population.parameters[key].value != new_param[key][0]:
        model_population.parameters[key].value = new_param[key]
    if model_concentration.parameters[key].value != new_param[key][0]:
        model_concentration.parameters[key].value = new_param[key]
model_concentration.parameters['omega'].value = 1
t = np.linspace(0, 500, 50000)
simres = ScipyOdeSimulator(model_concentration, t, verbose=False).run(); #deterministic simulation
yout = simres.all;
observables=['obsdrug', 'obsCycA', 'obsCdc20', 'obsCycAdrug']; #plot Cyclin E, Cdc20 and Cyclin A
for obs in range(4) :
    plt.plot(t[:],yout[observables[obs]][:], label=observables[obs], lw=2)
plt.legend(loc=1)
plt.ylim(0, 1)
plt.show()


omega = 200
new_param = pd.read_csv("parameters_drugB.csv") 
del new_param['Stress']
for key in new_param.keys():
    if model_population.parameters[key].value != new_param[key][0]:
        model_population.parameters[key].value = new_param[key]
    if model_concentration.parameters[key].value != new_param[key][0]:
        model_concentration.parameters[key].value = new_param[key]
model_concentration.parameters['omega'].value = 1
t = np.linspace(0, 500, 50000)
simres = ScipyOdeSimulator(model_concentration, t, verbose=False).run(); #deterministic simulation
yout = simres.all;
observables=['obsdrug', 'obsCycB', 'obsCdc20', 'obsCycBdrug']; #plot Cyclin E, Cdc20 and Cyclin A
for obs in range(4) :
    plt.plot(t[:],yout[observables[obs]][:], label=observables[obs], lw=2)
plt.legend(loc=1)
plt.ylim(0, 1)
plt.show()


t = np.linspace(0, 500, 50000)
simres = ScipyOdeSimulator(model_concentration, t, verbose=False).run(); #deterministic simulation
yout = simres.all;
observables=['obsCycE','obsCycA', 'obsCycAdrug']; #plot Cyclin E, Cdc20 and Cyclin A
for obs in range(3) :
    plt.plot(t[:],yout[observables[obs]][:], label=observables[obs], lw=2)
plt.legend(loc=0)
plt.show()

t = np.linspace(0, 500, 50000)
simres = ScipyOdeSimulator(model_concentration, t, verbose=False).run();
yout = simres.all;
observables=['obsCycE','obsCycD']; #plot Cyclin E, Cdc20 and Cyclin A
for obs in range(2) :
    plt.plot(t[:],yout[observables[obs]][:], label=observables[obs], lw=2)
plt.legend(loc=0)
plt.show()

t = np.linspace(0, 500, 50000)
simres = ScipyOdeSimulator(model_concentration, t, verbose=False).run(); #deterministic simulation
yout = simres.all;
observables=['obsCycD','obsCdc20']; #plot Cyclin E, Cdc20 and Cyclin A

for obs in range(2) :
    plt.plot(t[:],yout[observables[obs]][:], label=observables[obs], lw=2)
plt.legend(loc=0)
plt.show()


t = np.linspace(0, 6000, 600000)
sim = BngSimulator(model_population, tspan=t, verbose=False); #stochastic simulator
simres = sim.run(method='ssa', verbose=False, seed=33); #running the simulation; remove seed for truly random; seed=33 is used for similarity to deterministic sim
yout = simres.observables;
observables=['obsCycE','obsCdc20', 'obsCycBdrug']; #plot Cyclin E, Cdc20 and Cyclin A
for obs in range(3) :    
    plt.plot(t[400000:],yout[observables[obs]][400000:], label=observables[obs], lw=2)
plt.legend(loc=1)
plt.show()

observables = [
    'obsSkp2', 'obsp21', 'obsp21p', 'obscdh1act', 'obsMyc', 'obsCycD', 
    'obsCycDp21', 'obspRb', 'obspRbp', 'obsE2F', 'obspRbE2F', 'obsCycE', 
    'obsCycEp21', 'obsCycA', 'obsCycAp21', 'obsCycB', 'obsCycBp21', 'obsCdc20'
    ]

all_prot = np.empty((500000,18))
i=0
for obs in observables:
    tmp = simres.observables[obs][100000:]
    tmp_prot = [x/200 for x in tmp]
    all_prot[:,i] = tmp_prot
    i+=1

drugb = simres.observables[100000:]
np.savetxt("C:\\Users\\olesk\\ECM_mod_output\\UMAP\\drugB_stoch.csv", all_prot, delimiter=",")




t = np.linspace(0, 6000, 600000)
sim = BngSimulator(model_population, tspan=t, verbose=False); #stochastic simulator
simres = sim.run(method='ssa', verbose=False, seed=33); #running the simulation; remove seed for truly random; seed=33 is used for similarity to deterministic sim
yout = simres.observables;
observables=['obsCycE','obsCdc20', 'obsCycEdrug']; #plot Cyclin E, Cdc20 and Cyclin A
for obs in range(3) :    
    plt.plot(t[550000:],yout[observables[obs]][550000:], label=observables[obs], lw=2)
plt.legend(loc=1)
plt.show()

observables = [
    'obsSkp2', 'obsp21', 'obsp21p', 'obscdh1act', 'obsMyc', 'obsCycD', 
    'obsCycDp21', 'obspRb', 'obspRbp', 'obsE2F', 'obspRbE2F', 'obsCycE', 
    'obsCycEp21', 'obsCycA', 'obsCycAp21', 'obsCycB', 'obsCycBp21', 'obsCdc20'
    ]

all_prot = np.empty((500000,18))
i=0
for obs in observables:
    tmp = simres.observables[obs][100000:]
    tmp_prot = [x/200 for x in tmp]
    all_prot[:,i] = tmp_prot
    i+=1

druge = simres.observables[100000:]
np.savetxt("C:\\Users\\olesk\\ECM_mod_output\\UMAP\\drugE_stoch.csv", all_prot, delimiter=",")