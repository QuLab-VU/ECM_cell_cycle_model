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

from Austin_proliferation_model_concentration_5 import model
model_concentration=model

from Austin_proliferation_model_population_5 import model
model_population=model


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
y = loadmat('C:/Users/Grim.M/Ecole/ens/Master/S2/Stage/Code/cellular_cycle/y_4.mat')
L=np.transpose(y['y'])
l=L.tolist()
#
#tmp_x = loadmat('C:/Users/Grim.M/Ecole/ens/Master/S2/Stage/Code/cellular_cycle/tmp_x_stoch.mat')
#L=np.transpose(tmp_x['tmp_x'])
#tmp_x=L.tolist()
#
#tmp_t = loadmat('C:/Users/Grim.M/Ecole/ens/Master/S2/Stage/Code/cellular_cycle/tmp_t_stoch.mat')
#tmp_t=tmp_t['tmp_t']


tmp_x = loadmat('C:/Users/Grim.M/Ecole/ens/Master/S2/Stage/Code/cellular_cycle/tmp_x_5.mat')
L=np.transpose(tmp_x['tmp_x'])
tmp_x=L.tolist()

tmp_t = loadmat('C:/Users/Grim.M/Ecole/ens/Master/S2/Stage/Code/cellular_cycle/tmp_t_5.mat')
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
#3rd version, changing the parameters:

import mat73

data_dict = mat73.loadmat(path+'/v1_test/params_8_20.mat')
data_dict_2 = mat73.loadmat(path+'/cellular_cycle/parameters__.mat')

data_dict_2 = mat73.loadmat(path+'/v1/params_7_11_21.mat')

new_param=data_dict['p']
del new_param['GF']
del new_param['Stress']

new_param_2=data_dict_2['p']
#del new_param_2['GF']
del new_param_2['Stress']



for key in new_param_2.keys():
    if new_param[key] != new_param_2[key]:
        print(key)



for key in new_param_2.keys():
    #print(key)
    if model_population.parameters[key].value == new_param_2[key]:
        print('')    
    else:
        print('pop')
        print(key)
        print(model_population.parameters[key].value)
        print(new_param_2[key])
    if model_concentration.parameters[key].value == new_param_2[key]:
        print('')    
    else:
        print('conc')
        print(key)
        print(model_concentration.parameters[key].value)
        print(new_param_2[key])
        
        


for key in new_param_2.keys():
    print(key)
    if model_population.parameters[key].value != new_param_2[key]:
        model_population.parameters[key].value = new_param_2[key]
    if model_concentration.parameters[key].value != new_param_2[key]:
        model_concentration.parameters[key].value = new_param_2[key]



for key in new_param.keys():
    print(key)
    if model_population.parameters[key].value != new_param[key]:
        model_population.parameters[key].value = new_param[key]
    if model_concentration.parameters[key].value != new_param[key]:
        model_concentration.parameters[key].value = new_param[key]
        
        
        
        

model_population.parameters['FAK'].value = (new_param['FAK']*200)
model_population.parameters['ERK'].value = (new_param['ERK']*200)
model_population.parameters['AKT'].value = (new_param['AKT']*200)
model_population.parameters['omega'].value 

#or 
model_population.parameters['k_EMY'].value

model_population.parameters['FAK'].value = (1*200)
model_population.parameters['ERK'].value = (1*200)
model_population.parameters['AKT'].value = (1*200)
model_population.parameters['omega'].value 

model_population.parameters['s'].value*model_population.parameters['omega'].value*model_population.parameters['V_sskp2'].value

model_population.parameters['s'].value*model_population.parameters['V_s2skp'].value*model_population.parameters['FAK'].value

model_population.parameters['s'].value * model_population.parameters['V_dskp2'].value * model_population.parameters['omega'].value * model_population.parameters['Skp2_0'].value / ( model_population.parameters['K_dskp2'].value * model_population.parameters['omega'].value + model_population.parameters['Skp2_0'].value ) * model_population.parameters['cdh1act_0'].value / ( model_population.parameters['K_cdh1'].value* model_population.parameters['omega'].value + model_population.parameters['cdh1act_0'].value)

model_concentration.parameters['m'].value*(-( model_concentration.parameters['kf1'].value + model_concentration.parameters['kf2'].value/model_concentration.parameters['omega'].value*model_concentration.parameters['BH3_0'].value)* model_concentration.parameters['BAX_0'].value + model_concentration.parameters['kb'].value * model_concentration.parameters['BAXm_0'].value + model_concentration.parameters['kb'].value * model_concentration.parameters['BAXmBCL2_0'].value)


model_population.parameters['s'].value * model_population.parameters['Skp2_0'].value * model_population.parameters['k_ddskp2'].value
    

print('\n'.join(f'{n}: {o}' for n,o in enumerate(model_concentration.odes)))
print('\n'.join(f'{n}: {o}' for n,o in enumerate(model.odes)))


model_concentration.species[19]




####Deterministic simulation


model_concentration.parameters['AKT'].value = 1
model_concentration.parameters['ERK'].value = 1
model_concentration.parameters['FAK'].value = 1
model_concentration.parameters['omega'].value = 1


t = pl.linspace(0, 1000, 3001)

simres = ScipyOdeSimulator(model_concentration, t).run() 
yout = simres.all

observables=['obsCycE','obsCdc20']

for obs in range(2) :
    pl.plot(t,yout[observables[obs]], label=observables[obs], lw=2)
pl.legend(loc=0)
pl.show()

pl.plot(t,yout['obsCdc20'], label='obsCdc20', lw=2)

pl.plot(t,yout['obsBAXmBCL2'] + yout['obsBAXm'] + yout['obsBAX'], label='obsBAX tot', lw=2)
pl.plot(t, yout['obsBAXmBCL2'] + yout['obsBH3BCL2'] + yout['obsBCL2'], label='obsBCL2 tot', lw=2)
pl.plot(t,yout['obsBAX'], label='obsBAX', lw=2)
pl.plot(t,yout['obsBCL2'], label='obsBCL2', lw=2)

pl.plot(l[23])

pl.plot(l[22])


####Comparison with Matlab deterministic simulation

from pathlib import Path
fn = Path(r'C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code\cellular_cycle\comparison_pysb_matlab_5_0-300.pdf').expanduser()

obs = ['obsSkp2', 'obsp21', 'obsp21p', 'obscdh1act', 'obsMyc', 'obsCycD', 'obsCycDp21', 'obspRb', 'obspRbp', 'obsE2F', 'obspRbE2F', 'obsCycE', 'obsCycEp21', 'obsCycA', 'obsCycAp21', 'obsCycB', 'obsCycBp21', 'obsCdc20', 'obsBAX', 'obsBAXmBCL2', 'obsBH3','obsBH3BCL2','obsBAXm','obsBCL2']

fig, axs = plt.subplots(6, 4)


for i in range(4):
    for j in range(6):
        axs[j, i].plot(t[0:300],yout[obs[j+i*6]][0:300], label=obs[j+i*6] +' PySB', lw=0.01)
        axs[j, i].plot(t[0:300],l[j+i*6][0:300], label=obs[j+i*6] +' Matlab', lw=0.01)
        axs[j, i].legend(fontsize=2)

for ax in axs.flat:
    ax.label_outer()
    
plt.savefig(fn, bbox_inches = 'tight')


fn = Path(r'C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code\cellular_cycle\comparison_pysb_matlab_4_.pdf').expanduser()

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

t = pl.linspace(0, 120, 170000)

sim = BngSimulator(model_population, tspan=t, verbose=True)
simres = sim.run(method='ssa', verbose=True)
yout = simres.observables

observables=['obsCycE','obsCdc20', 'obsCycA', 'obsCycB', 'obsp21']
#observables=['obsSkp2']


for obs in range(5) :    
    plt.plot(t,yout[observables[obs]], label=observables[obs], lw=2)
plt.legend(loc=0)
plt.show()

observables=['obsCycE','obsCdc20', 'obsCycD']
#observables=['obsSkp2']


for obs in range(3) :    
    plt.plot(t,yout[observables[obs]], label=observables[obs], lw=2)
plt.legend(loc=0)
plt.show()



###########
#save to matlab
from scipy.io import savemat
savemat(r'C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code\obsCycE.mat', {'obsCycE': yout['obsCycE']})

from peakfinder import peakfinder
(a,b)=peakfinder(yout['obsCycE'], [], 0.2)


simres2= sim.run(tspan = t[-1] + t, initials = simres.species[-1])
yout2 = simres2.observables

for obs in model.observables:    
    plt.plot(t[-1] + t,yout2[obs.name], label=obs.name, lw=2)
plt.legend(loc=0)
plt.show()


plt.plot(t[3100:],yout['obsCycE'][3100:], label='Cyclin E', lw=2)
plt.plot(t[3100:],yout['obsCdc20'][3100:], label='Cdc 20', lw=2)
plt.legend(loc=0)
plt.savefig(r'C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\report\figures\stoch_simulations\stoch_sim_param_0.pdf')


#############"


####Comparing

from pathlib import Path
fn = Path(r'C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code\cellular_cycle\comparison_pysb_matlab_stoch_5.pdf').expanduser()

obs=['obsSkp2', 'obsp21', 'obsp21p', 'obscdh1act', 'obsMyc', 'obsCycD', 'obsCycDp21', 'obspRb', 'obspRbp', 'obsE2F', 'obspRbE2F', 'obsCycE', 'obsCycEp21', 'obsCycA', 'obsCycAp21', 'obsCycB', 'obsCycBp21', 'obsCdc20', 'obsBAX', 'obsBAXmBCL2', 'obsBH3','obsBH3BCL2','obsBAXm','obsBCL2']

fig, axs = plt.subplots(8, 6)


for i in range(6):
    for j in range(4):
        axs[2*j, i].plot(t,yout[obs[j+i*4]], label=obs[j+i*4] +' PySB', lw=0.01, color='orange')
        axs[2*j, i].legend(fontsize=2)
        axs[2*j+1, i].plot(tmp_t,tmp_x[j+i*4], label=obs[j+i*4] +' Matlab', lw=0.01, color='blue')
        axs[2*j+1, i].legend(fontsize=2)

for ax in axs.flat:
    ax.label_outer()


plt.savefig(fn, bbox_inches='tight')

j=0
i=4
plt.plot(t,yout[obs[j+i*4]], label=obs[j+i*4] +' PySB', lw=1, color='orange')
plt.plot(tmp_t,tmp_x[j+i*4], label=obs[j+i*4] +' Matlab', lw=1, color='blue')
plt.show()

#############################################

import math

def RBS(x,x0,sigma) :
    observations , dimensions= np.shape(x0)
    s=0
    for i in range(observations):
        s+=math.exp(-np.linalg.norm(x-x0[i])**2/(2*sigma**2)) / (sigma*math.sqrt(2*math.pi))
    if observations != 0:
        return s/observations
    else:
        return s
    
    
from pysb.simulator import BngSimulator
from SMs_G1s_simulation_BngSimulator import IMTs
from runStochCell import runStochcell


tmp_t=pl.linspace(0, 4000, 8000)
sim = BngSimulator(model_population, tspan=tmp_t, verbose=True)
traj = sim.run()

variable_range = runStochcell(model_concentration)


ini={}
print(len(variable_range['obsBCL2']));
a=np.random.randint(len(variable_range['obsBCL2']))
for var in list(variable_range.keys()) :
    ini[var] = [round(variable_range[var][a]*200)]
    
Initial=[]
for j in ini:
    Initial=Initial+[ini[j][-1]]
initials=np.array(Initial)
        
(a,b) = IMTs(tmp_t, initials, model_population)

sigma=2
x=np.linspace(0,80,800)
y_sim=[RBS(i,np.transpose(np.array([b])),sigma) for i in x]    
# Plot the KDE as fitted using the default arguments
bins_list = np.linspace(0,90,25)
plt.plot(x, y_sim, lw=3, label = "G1 simulations")
plt.legend(loc='best')

plt.savefig(r'C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\report\figures\stoch_simulations\distrib_G1.pdf')


sigma=2
x=np.linspace(0,80,800)
y_sim=[RBS(i,np.transpose(np.array([a])),sigma) for i in x]    
# Plot the KDE as fitted using the default arguments
bins_list = np.linspace(0,90,25)
plt.plot(x, y_sim, lw=3, label = "S-M simulations")
plt.legend(loc='best')

plt.savefig(r'C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\report\figures\stoch_simulations\distrib_SM.pdf')





####################################
#find differences with matlab

from numpy import diff
observables=['obsCycE','obsCdc20']
dt = t[1] - t[0]

y = yout[observables[obs]]
dy = diff(y)/dt

derivate={}

for o in range(len(model_concentration.observables)):
    model_concentration.observables[o].name
    
    dt = t[1] - t[0]
    y = yout[model_concentration.observables[o].name]
    dy = diff(y)/dt
    derivate[model_concentration.observables[o].name]=dy[0]
    


V_1Mb=model_concentration.parameters['V_1Mb'].value
V_2Mb=model_concentration.parameters['V_2Mb'].value
K_ie=model_concentration.parameters['K_ie'].value
K_ib=model_concentration.parameters['K_ib'].value
b=model_concentration.parameters['b'].value
b2=model_concentration.parameters['b2'].value
K_1Mb=model_concentration.parameters['K_1Mb'].value
K_2Mb=model_concentration.parameters['K_2Mb'].value
V_1cdc20=model_concentration.parameters['V_1cdc20'].value
K_1cdc20=model_concentration.parameters['K_1cdc20'].value
V_2cdc20=model_concentration.parameters['V_2cdc20'].value
K_2cdc20=model_concentration.parameters['K_2cdc20'].value
CycB_tot=model_concentration.parameters['CycB_tot'].value
s=model_concentration.parameters['s'].value
k_c7=model_concentration.parameters['k_c7'].value
k_c8=model_concentration.parameters['k_c8'].value
s=model_concentration.parameters['s'].value

obsCycB=yout['obsCycB'][0]
obsCycE=yout['obsCycE'][0]
obsCycA=yout['obsCycA'][0]
obsCycBp21=yout['obsCycBp21'][0]
obsp21=yout['obsp21'][0]
obsCdc20=yout['obsCdc20'][0]



from scipy.io import loadmat
y = loadmat('C:/Users/Grim.M/Ecole/ens/Master/S2/Stage/Code/v1_test/cdc20.mat')
L=np.transpose(y['z'])
l=L.tolist()[0]

plt.plot(t, l)
plt.plot(t, yout['obsCdc20'], label='obsCdc20', lw=2)