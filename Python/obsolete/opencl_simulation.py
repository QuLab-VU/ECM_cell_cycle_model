# -*- coding: utf-8 -*-
"""
Created on Tue Apr  7 09:32:00 2020

@author: Grim.M
"""


path="/home/jeremie"
#path=r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code"
import sys
sys.path.insert(0, path+"/cellular_cycle") #add a pathway to the model
from Austin_proliferation_model_population import model

from pysb.simulator import OpenCLSSASimulator, ScipyOdeSimulator
from pysb.simulator import BngSimulator
import matplotlib.pyplot as plt
import numpy as np

###Bng ssa
t = np.linspace(0, 200, 1000)

sim = BngSimulator(model, tspan=t, verbose=True)
simres = sim.run(method='ssa', verbose=True)
yout = simres.observables

observables=['obsCdc20','obsCycD']

for obs in range(2) :    
    plt.plot(t,yout[observables[obs]], label=observables[obs], lw=2)
plt.legend(loc=0)
plt.show()

###Opencl

n_sim = 2
tspan = np.linspace(0, 200, 1000)
sim = OpenCLSSASimulator(model, tspan=tspan, verbose=True)
traj = sim.run(tspan=tspan, number_sim=n_sim)

result_Cdc20 = np.array(traj.observables)['obsCdc20']
result_CycE = np.array(traj.observables)['obsCycE']


np.save('/home/jeremie/result_obsCdc20.npy', result_Cdc20)
np.save('/home/jeremie/result_obsCycE.npy', result_CycE)


np.shape(result_obsCdc20)

x = np.array([tr[:] for tr in result_obsCdc20]).T

plt.plot(t,x, '0.5', lw=2, alpha=0.25)  # individual trajectories
plt.plot(t,x.mean(axis=1), 'b', lw=3, label='mean')
plt.plot(t,x.max(axis=1), 'k--', lw=2, label="min/max")
plt.plot(t,x.min(axis=1), 'k--', lw=2)
plt.show()
    

