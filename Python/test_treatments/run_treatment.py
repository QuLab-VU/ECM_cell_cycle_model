# -*- coding: utf-8 -*-
"""
Created on Sat Oct 15 16:37:38 2022

@author: olesk
"""

import os
import sys
# sys.path.insert(0, os.getcwd() + "/Python/test_treatments/")#add a pathway to the model #not sure about how to use folder dependency, to change
# sys.path.insert(0, os.getcwd() + "/../fit") #add a pathway to the fit functions #not sure about how to use folder dependency, to change

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib; matplotlib.use('Qt5Agg')
import seaborn as sns
import random
import sobol_seq
from pathlib import Path

SAVE_DIR = Path(Path.home(),'ECM_mod_output')
if not SAVE_DIR.exists():
    SAVE_DIR.mkdir()

#Importing the model
from cell_cycle_death_model import model

##Import the simulator: for ODE and stochastic simulation :
from pysb.simulator import BngSimulator
from pysb.simulator import ScipyOdeSimulator


# params = {'axes.titley': 1,
#          'axes.titlepad': -14,
#          'axes.labelsize': 14,
#          'axes.titlesize': 14,
#          'xtick.labelsize': 12,
#          'ytick.labelsize': 12}
# plt.rcParams.update(params)

#change the model parameters of the model if needed

#new_param = pd.read_csv('parameters_cycei.csv')
#new_param = pd.read_csv('parameters_cycdi.csv')
new_param = pd.read_csv('parameters.csv')
for key in new_param.keys():
    try:
        if model.parameters[key].value != new_param[key][0]:
            model.parameters[key].value = new_param[key]
    except:
        print(key)
model.parameters['omega'].value = 200

model.parameters['AKT'].value = 0.3
model.parameters['ERK'].value = 0.2
model.parameters['FAK'].value = 0.4


t = np.linspace(0, 500, 50000)
simres = ScipyOdeSimulator(model, t, verbose=False).run()
yout = simres.all
observables=['obsCycE', 'obsCdc20', 'obsBAXm']; #plot Cyclin E, Cdc20 and Cyclin A
plt.figure()
for obs in range(3):
    plt.plot(t[:],yout[observables[obs]][:], label=observables[obs], lw=2)
plt.legend(loc=1)
#plt.ylim(0, 1)
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