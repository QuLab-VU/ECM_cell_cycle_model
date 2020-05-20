# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 11:52:01 2020

@author: Jérémie Marlhens
Title: Austin_proliferation_model_simulation
"""

path="/home/jeremie"

#path=r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code"
import sys
sys.path.insert(0, path+"/cellular_cycle") #add a pathway to the model
from Austin_proliferation_model_concentration import model

from Austin_proliferation_model_population import model



import sys
sys.path.insert(0,r"C:\Users\Grim.M\Documents\git_repos\ECM_cell_cycle_model\Python")



from pysb.simulator import ScipyOdeSimulator
from pysb.simulator import BngSimulator
from pysb.simulator import StochKitSimulator
from pysb.bng import generate_equations

from pysb.logging import EXTENDED_DEBUG

import pylab as pl
import numpy as np
from pysb.export import export 
import matplotlib.pyplot as plt


#Compare with matlab simulations

from scipy.io import loadmat
y = loadmat('C:/Users/Grim.M/Ecole/ens/Master/S2/Stage/Code/cellular_cycle/y.mat')
L=np.transpose(y['y'])
l=L.tolist()

tmp_x = loadmat('C:/Users/Grim.M/Ecole/ens/Master/S2/Stage/Code/cellular_cycle/tmp_x_stoch.mat')
L=np.transpose(tmp_x['tmp_x'])
tmp_x=L.tolist()

tmp_t = loadmat('C:/Users/Grim.M/Ecole/ens/Master/S2/Stage/Code/cellular_cycle/tmp_t_stoch.mat')
tmp_t=tmp_t['tmp_t']


'''
####Export the model to matlab
matlab_output = export(model, 'matlab')

with open('output_matlab.m', 'w') as f:
    f.write(matlab_output)

'''
'''
from pysb.kappa import contact_map
import networkx as nx
nx.draw_networkx(g, with_labels=True)

generate_equations(model)
print('\n'.join(f'{n}: {o}' for n,o in enumerate(model.odes)))
'''
####Deterministic simulation
t = pl.linspace(0, 1000, 3001)

m=ScipyOdeSimulator(model, t)

simres = ScipyOdeSimulator(model_concentration, t).run() 
yout = simres.all

observables=['obsCycE','obsCdc20']

for obs in range(2) :
    pl.plot(t,yout[observables[obs]], label=observables[obs], lw=2)
pl.legend(loc=0)
pl.show()

####Comparison with Matlab deterministic simulation

from pathlib import Path
fn = Path(r'C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code\cellular_cycle\comparison_pysb_matlab_0-10.pdf').expanduser()

obs = ['obsSkp2', 'obsp21', 'obsp21p', 'obscdh1act', 'obsMyc', 'obsCycD', 'obsCycDp21', 'obspRb', 'obspRbp', 'obsE2F', 'obspRbE2F', 'obsCycE', 'obsCycEp21', 'obsCycA', 'obsCycAp21', 'obsCycB', 'obsCycBp21', 'obsCdc20', 'obsBAX', 'obsBAXmBCL2', 'obsBH3','obsBH3BCL2','obsBAXm','obsBCL2']

fig, axs = plt.subplots(6, 4)


for i in range(4):
    for j in range(6):
        axs[j, i].plot(t[0:10],yout[obs[j+i*6]][0:10], label=obs[j+i*6] +' PySB', lw=0.01)
        axs[j, i].plot(t[0:10],l[j+i*6][0:10], label=obs[j+i*6] +' Matlab', lw=0.01)
        axs[j, i].legend(fontsize=2)

for ax in axs.flat:
    ax.label_outer()
    
plt.savefig(fn, bbox_inches = 'tight')


fn = Path(r'C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code\cellular_cycle\comparison_pysb_matlab.pdf').expanduser()

fig, axs = plt.subplots(6, 4)


for i in range(4):
    for j in range(6):
        axs[j, i].plot(t,yout[obs[j+i*6]], label=obs[j+i*6] +' PySB', lw=0.01)
        axs[j, i].plot(t,l[j+i*6], label=obs[j+i*6] +' Matlab', lw=0.01)
        axs[j, i].legend(fontsize=2)

for ax in axs.flat:
    ax.label_outer()
    
plt.savefig(fn, bbox_inches='tight')


####Stochastic simulation
t = pl.linspace(0, 100, 4000)

sim = BngSimulator(model_population, tspan=t, verbose=True)
simres = sim.run(method='ssa', verbose=True)
yout = simres.observables

observables=['obsCycE','obsCdc20']

for obs in range(2) :    
    plt.plot(t,yout[observables[obs]], label=observables[obs], lw=2)
plt.legend(loc=0)
plt.show()
    

simres2= sim.run(tspan = t[-1] + t, initials = simres.species[-1])
yout2 = simres2.observables

for obs in model.observables:    
    plt.plot(t[-1] + t,yout2[obs.name], label=obs.name, lw=2)
plt.legend(loc=0)
plt.show()

'''

####Comparing

from pathlib import Path
fn = Path(r'C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code\cellular_cycle\comparison_pysb_matlab_stoch.pdf').expanduser()

obs=['obsSkp2', 'obsp21', 'obsp21p', 'obscdh1act', 'obsMyc', 'obsCycD', 'obsCycDp21', 'obspRb', 'obspRbp', 'obsE2F', 'obspRbE2F', 'obsCycE', 'obsCycEp21', 'obsCycA', 'obsCycAp21', 'obsCycB', 'obsCycBp21', 'obsCdc20', 'obsBAX', 'obsBAXmBCL2', 'obsBH3','obsBH3BCL2','obsBAXm','obsBCL2']

fig, axs = plt.subplots(8, 6)


for i in range(6):
    for j in range(4):
        axs[2*j, i].plot(t,yout[obs[j+i*4]], label=obs[j+i*4] +' PySB', lw=0.01)
        axs[2*j, i].legend(fontsize=2)
        axs[2*j+1, i].plot(tmp_t,tmp_x[j+i*4], label=obs[j+i*4] +' Matlab', lw=0.01)
        axs[2*j+1, i].legend(fontsize=2)

for ax in axs.flat:
    ax.label_outer()
    
plt.savefig(fn, bbox_inches='tight')




