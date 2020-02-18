# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 11:52:01 2020

@author: Jérémie Marlhens
Title: Austin_proliferation_model_simulation
"""

import sys
sys.path.insert(0, r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code\cellular_cycle") #add a pathway to the model
from Austin_proliferation_model_population import model

from pysb.simulator import ScipyOdeSimulator
from pysb.simulator import BngSimulator
from pysb.bng import generate_equations

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

simres = ScipyOdeSimulator(model, t).run() 
yout = simres.all

observables=['obscdh1act','obsMyc']

for obs in range(2) :    
    pl.plot(t,yout[observables[obs]], label=observables[obs], lw=2)
pl.legend(loc=0)
pl.show()
    
####Comparison with Matlab simulation

from pathlib import Path
fn = Path(r'C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code\cellular_cycle\comparison_pysb_matlab_0-10.pdf').expanduser()

obs=['obsSkp2', 'obsp21', 'obsp21p', 'obscdh1act', 'obsMyc', 'obsCycD', 'obsCycDp21', 'obspRb', 'obspRbp', 'obsE2F', 'obspRbE2F', 'obsCycE', 'obsCycEp21', 'obsCycA', 'obsCycAp21', 'obsCycB', 'obsCycBp21', 'obsCdc20', 'obsBAX', 'obsBAXmBCL2', 'obsBH3','obsBH3BCL2','obsBAXm','obsBCL2']

fig, axs = plt.subplots(6, 4)


for i in range(4):
    for j in range(6):
        axs[j, i].plot(t[0:10],yout[obs[j+i*6]][0:10], label=obs[j+i*6] +' PySB', lw=0.01)
        axs[j, i].plot(t[0:10],l[j+i*6][0:10], label=obs[j+i*6] +' Matlab', lw=0.01)
        axs[j, i].legend(fontsize=2)

for ax in axs.flat:
    ax.label_outer()
    
plt.savefig(fn, bbox_inches='tight')


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
t = pl.linspace(0, 500, 450000)

sim = BngSimulator(model, tspan=t, verbose=True)
simres = sim.run(method='ssa', verbose=True)
yout = simres.observables

observables=['obsCycD','obsCdc20']

for obs in range(2) :    
    plt.plot(t,yout[observables[obs]], label=observables[obs], lw=2)
plt.legend(loc=0)
plt.show()
    
'''
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

# Hide x labels and tick labels for top plots and y ticks for right plots.
for ax in axs.flat:
    ax.label_outer()
    
plt.savefig(fn, bbox_inches='tight')

####


Skp2=l[0][0];
p21=l[1][0];
p21p=l[2][0];  
Cdh1=l[3][0]; 
Myc=l[4][0];
CycD=l[5][0];  
CycDp21=l[6][0];
Rb=l[7][0];
pRb=l[8][0];
E2F=l[9][0];
RbE2F=l[10][0];
CycE=l[11][0];
CycEp21=l[12][0];
CycA=l[13][0];
CycAp21=l[14][0];
CycB=l[15][0];
CycBp21=l[16][0];
Cdc20=l[17][0];
BAX = l[18][0];
BAXmBCL2 = l[19][0];
BH3=l[20][0];
BH3BCL2 = l[21][0];

Skp2=l[0][0];
p21=l[1][0];
p21p=l[2][0];  
Cdh1=l[3][0]; 
Myc=l[4][0];
CycD=l[5][0];  
CycDp21=l[6][0];
Rb=l[7][0];
pRb=l[8][0];
E2F=l[9][0];
RbE2F=l[10][0];
CycE=l[11][0];
CycEp21=l[12][0];
CycA=l[13][0];
CycAp21=l[14][0];
CycB=l[15][0];
CycBp21=l[16][0];
Cdc20=l[17][0];
BAX = l[18][0];
BAXmBCL2 = l[19][0];
BH3=l[20][0];
BH3BCL2 = l[21][0];




+ yout['obsE2F'][0]*model.parameters['v_s2p21'].value*(model.parameters['K_i13'].value / (model.parameters['K_i13'].value + yout['obspRb'][0])) * (model.parameters['K_i14'].value / (model.parameters['K_i14'].value + yout['obspRbp'][0]))+ model.parameters['v_s1p21'].value * model.parameters['omega'].value+ (yout['obsp21'][0]*yout['obsCycE'][0]*model.parameters['V_1p21'].value / (model.parameters['K_1p21'].value * model.parameters['omega'].value + yout['obsp21'][0]))*(-1) + (yout['obsp21'][0]*((model.parameters['V_Fp'].value * model.parameters['FAK'].value)/ ((model.parameters['K_Fp'].value * model.parameters['omega'].value) + yout['obsp21'][0])))*(-1) +  yout['obsp21p'][0]*((model.parameters['V_2p21'].value * model.parameters['omega'].value) / (model.parameters['K_2p21'].value * model.parameters['omega'].value + yout['obsp21p'][0]))+(yout['obsp21'][0]*model.parameters['k_ddp21'].value)*(-1)+(yout['obsp21'][0]*yout['obsCycD'][0]*(model.parameters['k_c1'].value / model.parameters['omega'].value))*(-1) + yout['obsCycDp21'][0]*model.parameters['k_c2'].value+ (yout['obsp21'][0]*yout['obsCycE'][0]*(model.parameters['k_c3'].value /model.parameters['omega'].value))*(-1) + yout['obsCycEp21'][0]*model.parameters['k_c4'].value+ (yout['obsp21'][0]*yout['obsCycA'][0]*(model.parameters['k_c5'].value / model.parameters['omega'].value))*(-1) + yout['obsCycAp21'][0]*model.parameters['k_c6'].value+ (yout['obsp21'][0]*yout['obsCycB'][0]*(model.parameters['k_c7'].value / model.parameters['omega'].value))*(-1)+ yout['obsCycBp21'][0]*model.parameters['k_c8'].value

 E2F*model.parameters['v_s2p21'].value*(model.parameters['K_i13'].value / (model.parameters['K_i13'].value + Rb)) * (model.parameters['K_i14'].value / (model.parameters['K_i14'].value + pRb))+ model.parameters['v_s1p21'].value * model.parameters['omega'].value+ (p21*CycE*model.parameters['V_1p21'].value / (model.parameters['K_1p21'].value * model.parameters['omega'].value + p21))*(-1) + (p21*((model.parameters['V_Fp'].value * model.parameters['FAK'].value)/ ((model.parameters['K_Fp'].value * model.parameters['omega'].value) + p21)))*(-1) +  p21p*((model.parameters['V_2p21'].value * model.parameters['omega'].value) / (model.parameters['K_2p21'].value * model.parameters['omega'].value + p21p))+(p21*model.parameters['k_ddp21'].value)*(-1)+(p21*CycD*(model.parameters['k_c1'].value / model.parameters['omega'].value))*(-1) + CycDp21*model.parameters['k_c2'].value+ (p21*CycE*(model.parameters['k_c3'].value /model.parameters['omega'].value))*(-1) + CycEp21*model.parameters['k_c4'].value+ (p21*CycA*(model.parameters['k_c5'].value / model.parameters['omega'].value))*(-1) + CycAp21*model.parameters['k_c6'].value+ (p21*CycB*(model.parameters['k_c7'].value / model.parameters['omega'].value))*(-1)+ CycBp21*model.parameters['k_c8'].value
