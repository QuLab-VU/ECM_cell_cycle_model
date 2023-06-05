# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 11:52:01 2020

@author: Jérémie Marlhens
Title: Austin_proliferation_model_simulation
"""
###############################################################################

#path="/home/jeremie"
path=r"C:\Users\Grim.M\Documents\git_repos\ECM_cell_cycle_model"

import sys
sys.path.insert(0, path+"/Python/cellular_cycle") #add a pathway to the model
sys.path.insert(0, path+"/Python/fit") #add a pathway to the model

from Austin_proliferation_model_concentration_6 import model
model_concentration = model
from Austin_proliferation_model_population_6 import model
model_population = model

#model_population = model
#del model
"""
import copy

from Austin_proliferation_model_population_6 import model
model_concentration = copy.deepcopy(model)

del model
from Austin_proliferation_model_population_6 import model
model_population = copy.deepcopy(model)
"""

from pysb.simulator import ScipyOdeSimulator
from pysb.simulator import BngSimulator
from pysb.simulator import StochKitSimulator
from pysb.bng import generate_equations

from pysb.logging import EXTENDED_DEBUG

import numpy as np
from pysb.export import export
import matplotlib.pyplot as plt



###############################################################################
#Compare with matlab simulations

simulation_ode_matlab = {}
simulation_ode_matlab['55 kPa DMSO'] = np.genfromtxt(path + '/Matlab/v2/ode_100_100_100.csv', delimiter=',').T
simulation_ode_matlab['55 kPa ERL'] = np.genfromtxt(path + '/Matlab/v2/ode_060_020_100.csv', delimiter=',').T
simulation_ode_matlab['4 kPa DMSO'] = np.genfromtxt(path + '/Matlab/v2/ode_100_100_060.csv', delimiter=',').T
simulation_ode_matlab['0.6 kPa DMSO'] = np.genfromtxt(path + '/Matlab/v2/ode_060_100_040.csv', delimiter=',').T


stoch_run = {}
for file_tmp, condition in zip(['dmso_long_4_29_1_1', 'dmso_long_4_29_2_1', 'dmso_long_4_29_3_1', 'erl_long_4_29_1_1', 'erl_long_4_29_2_1', 'erl_long_4_29_3_1'], ['06 kPa DMSO', '4 kPa DMSO', '55 kPa DMSO', '06 kPa ERL', '4 kPa ERL', '55 kPa ERL']):
    stoch_run[condition] = np.genfromtxt(r"C:\Users\Grim.M\Ecole\Ens\Master\S2\Stage\Code\matlab_runs/" +  file_tmp + ".txt", delimiter=',', skip_footer=450000).T

plt.plot(stoch_run['55 kPa DMSO'][3, :10000])



###############################################################################
###############################################################################
###############################################################################
###############################################################################
#3rd version, changing the parameters:

import pandas as pd
#change the model parameters of the model if needed
omega = 200
#new_param = pd.read_csv(path + "/umap_window/parameters.csv")
new_param = pd.read_csv(path + "/Matlab/v2/parameters.csv")
del new_param['Stress']
for key in new_param.keys():
    try:
        if model_population.parameters[key].value != new_param[key][0]:
            model_population.parameters[key].value = new_param[key]
        if model_concentration.parameters[key].value != new_param[key][0]:
            model_concentration.parameters[key].value = new_param[key]

#        if model.parameters[key].value != new_param[key][0]:
#            model.parameters[key].value = new_param[key]

    except: print(key)



model_population.parameters['omega'].value = 200
model_concentration.parameters['omega'].value = 1
#model.parameters['omega'].value = 1


#model_concentration.parameters['v_s1p21'].value = 0.3
#model_population.parameters['v_s1p21'].value = 0.8
#model_population.parameters['F_n'].value = 8


###############################################################################
# Setting initial concentrations

species_initials = ['Skp2_0', 'p21_0', 'p21p_0', 'cdh1act_0', 'cdh1inact_0', 'Myc_0', 'CycD_0', 'CycDp21_0', 'pRb_0', 'pRbp_0', 'E2F_0', 'pRbE2F_0', 'CycE_0', 'CycEp21_0', 'CycA_0', 'CycAp21_0', 'CycB_0', 'CycBp21_0', 'Cdc20_0', 'BAX_0', 'BAXmBCL2_0', 'BH3_0', 'BH3BCL2_0', 'BAXm_0', 'BCL2_0']
species_initials_values = [0.0199387776457489, 0.0212324326339528, 0.104907517295920, 0.970596980815461, 1-0.970596980815461, 0.321067821067829, 0.175953608509162, 0.0138939800293163, 0.301757475166483, 0.545620559818762, 0.0403080965235495, 0.554487916752975, 0.0640517424993558, 0.00487528808175496, 0.0182005140710040, 0.00234453838504852, 0.0234757871244827, 0.00435874152378834, 0.0204771947463057, 0.999000159839911, 0.000997752220497874, 5.81398691338436e-05, 0.0332751934642001, 1 - 0.999000159839911 - 0.000997752220497874, 0.95 - 0.0332751934642001 - 0.000997752220497874]

for species_initial, species_initials_value in zip(species_initials, species_initials_values):
    
    # print(model_concentration.parameters[species_initial].value==model_concentration.parameters[species_initial].value)
    
    model_concentration.parameters[species_initial].value = species_initials_value
    model_population.parameters[species_initial].value = int(round(species_initials_value * model_population.parameters['omega'].value))
    

#
#model_concentration.initials[0].value.name
#
#
#for parameter in model.parameters:
#    print(parameter.name)
#    print(model.parameters[parameter.name].value==model_concentration.parameters[parameter.name].value)
#    if not model.parameters[parameter.name].value==model_concentration.parameters[parameter.name].value:
#        print(' AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA')


###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################

conditions = ['0.6 kPa DMSO', '0.6 kPa ERL', '4 kPa DMSO', '4 kPa ERL', '55 kPa DMSO', '55 kPa ERL']
inputs = {}
inputs['FAK'] = {'4 kPa DMSO': 0.6, '4 kPa ERL': 0.6, '0.6 kPa DMSO': 0.4, '0.6 kPa ERL': 0.4, '55 kPa DMSO': 1.0, '55 kPa ERL': 1.0}
inputs['ERK'] = {'4 kPa DMSO': 1.0, '4 kPa ERL': 0.2, '0.6 kPa DMSO': 1.0, '0.6 kPa ERL': 0.2, '55 kPa DMSO': 1.0, '55 kPa ERL': 0.2}
inputs['AKT'] = {'4 kPa DMSO': 1.0, '4 kPa ERL': 0.6, '0.6 kPa DMSO': 0.6, '0.6 kPa ERL': 0.3, '55 kPa DMSO': 1.0, '55 kPa ERL': 0.6}

###############################################################################

####Deterministic simulation

condition = '55 kPa DMSO'
condition = '0.6 kPa DMSO'
condition = '55 kPa ERL'

for input_signal in inputs.keys():
    model_concentration.parameters[input_signal].value = inputs[input_signal][condition]
    # model.parameters[input_signal].value = inputs[input_signal][condition]

t = np.linspace(0, 1000, 100001)
simres = ScipyOdeSimulator(model_concentration, t).run() 
yout = simres.all

observables=['obsCycE','obsCdc20', 'obsCycA', 'obsCycB', 'obsp21']
for obs in range(4) :
    plt.plot(t, yout[observables[obs]], label=observables[obs], lw=2)
    print(yout[observables[obs]][0])
plt.legend(loc=0)
plt.show()


####Comparison with Matlab deterministic simulation

nmax = 10000
from pathlib import Path
fn = Path(r'C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code\cellular_cycle\comparison_pysb_matlab_6_'+ condition.replace(' ', '_') +'_0-'+ str(t[nmax]) + '.pdf').expanduser()

obs = ['obsSkp2', 'obsp21', 'obsp21p', 'obscdh1act', 'obsMyc', 'obsCycD', 'obsCycDp21', 'obspRb', 'obspRbp', 'obsE2F', 'obspRbE2F', 'obsCycE', 'obsCycEp21', 'obsCycA', 'obsCycAp21', 'obsCycB', 'obsCycBp21', 'obsCdc20', 'obsBAX', 'obsBAXmBCL2', 'obsBH3','obsBH3BCL2','obsBAXm','obsBCL2']

fig, axs = plt.subplots(6, 4)

for i in range(4):
    for j in range(6):
        axs[j, i].plot(t[0:nmax], yout[obs[j+i*6]][0:nmax], label=obs[j+i*6] +' PySB', lw=0.01)
        axs[j, i].plot(t[0:nmax], simulation_ode_matlab[condition][j+i*6][0:nmax], label=obs[j+i*6] +' Matlab', lw=0.01)
        axs[j, i].legend(fontsize=2)

for ax in axs.flat:
    ax.label_outer()

plt.savefig(fn, bbox_inches = 'tight')

###############################################################################
###############################################################################

####Stochastic simulation

condition = '55 kPa DMSO'


for input_signal in inputs.keys():
    model_population.parameters[input_signal].value = inputs[input_signal][condition] * 200


t = np.linspace(0, 500, 50001)

sim = BngSimulator(model_population, tspan=t, verbose=True)
simres = sim.run(method='ssa', verbose=True)
yout = simres.observables

observables=['obsCycE','obsCdc20', 'obsCycD']
observables=['obsCycE','obsCdc20', 'obsCycA', 'obsCycB', 'obsp21']

for obs in range(len(observables)) :    
    plt.plot(t,yout[observables[obs]], label=observables[obs], lw=2)
plt.legend(loc=0)
plt.show()


####Comparison with Matlab Stochastic simulation

condition = '55 kPa DMSO'

from pathlib import Path
fn = Path(r'C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code\cellular_cycle\comparison_pysb_matlab_stoch_6_' + condition.replace(' ', '_') +'.pdf').expanduser()

obs=['obsSkp2', 'obsp21', 'obsp21p', 'obscdh1act', 'obsMyc', 'obsCycD', 'obsCycDp21', 'obspRb', 'obspRbp', 'obsE2F', 'obspRbE2F', 'obsCycE', 'obsCycEp21', 'obsCycA', 'obsCycAp21', 'obsCycB', 'obsCycBp21', 'obsCdc20', 'obsBAX', 'obsBAXmBCL2', 'obsBH3','obsBH3BCL2','obsBAXm','obsBCL2']

fig, axs = plt.subplots(8, 6)

for i in range(6):
    for j in range(4):
        axs[2*j, i].plot(t, yout[obs[j+i*4]], label=obs[j+i*4] +' PySB', lw=0.01, color='orange')
        axs[2*j, i].legend(fontsize=2)
        max_tmp = np.max(yout[obs[j+i*4]])
        
        if j+i*4<stoch_run[condition].shape[0]:
            axs[2*j+1, i].plot(t, stoch_run[condition][j+i*4]*200, label=obs[j+i*4] +' Matlab', lw=0.01, color='blue')
            axs[2*j+1, i].legend(fontsize=2)
            
            max_tmp = max(max_tmp, np.max(stoch_run[condition][j+i*4]))
            
        axs[2*j+1, i].set_ylim(0, 1.1*max_tmp)
        axs[2*j, i].set_ylim(0, 1.1*max_tmp)

for ax in axs.flat:
    ax.label_outer()

plt.savefig(fn, bbox_inches='tight')


###############################################################################
###############################################################################
###############################################################################
###############################################################################
"""

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


"""

"""
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
"""


###############################################################################
###############################################################################
###############################################################################
#GPU SSA test


PYOPENCL_CTX='0'

import numpy as np
from pysb.examples.earm_1_0 import model
from pysb.simulator import OpenCLSSASimulator, ScipyOdeSimulator

if __name__ == "__main__":
    n_sim = 100
    tspan = np.linspace(0, 5000, 101)
    name = 'tBid_total'
    sim = OpenCLSSASimulator(model, tspan=tspan, verbose=True)
    traj = sim.run(tspan=tspan, number_sim=n_sim)
    
    result = np.array(traj.observables)[name]
    
    x = np.array([tr[:] for tr in result]).T



from pysb.examples.schloegl import model



if __name__ == '__main__':
    
    n_sim=20
    x_0=100
    
    obs_name = 'X_total'
    tspan = np.linspace(0, 100, 101)
    x_0 = int(x_0)
    
    
    
    savename = 'schloegl_{}'.format(int(x_0))
    model.parameters['X_0'].value = x_0
    # create simulator and run simulations
    simulator = OpenCLSSASimulator(model, verbose=True,initials=[[10, 5] for i in range(10)]+ [[1, 2] for i in range(10)])
    
    traj = simulator.run(
        tspan,
        number_sim=n_sim
    )
    
    
    x = traj.dataframe[obs_name].unstack(0).values







# import the pysb module and all its methods and functions
from pysb import *
import numpy as np
from pysb.simulator import OpenCLSSASimulator, ScipyOdeSimulator


# instantiate a model
Model()

# declare monomers
Monomer('C8', ['b'])
Monomer('Bid', ['b', 'S'], {'S': ['u', 't']})

# input the parameter values
Parameter('kf', 1.0e-07)
Parameter('kr', 1.0e-03)
Parameter('kc', 1.0)

Expression('kg', kf/(kf+kc)) 

# now input the rules
Rule('C8_Bid_bind', C8(b=None) + Bid(b=None, S='u') | C8(b=1) % Bid(b=1, S='u'), kf, kr)
Rule('tBid_from_C8Bid', C8(b=1) % Bid(b=1, S='u') >> C8(b=None) + Bid(b=None, S='t'), kg)

model.monomers
model.parameters
model.rules

Observable('obsC8', C8(b=None))
Observable('obsBid', Bid(b=None, S='u'))
Observable('obstBid', Bid(b=None, S='t'))

Parameter('C8_0', 1000)
Parameter('Bid_0', 10000)
Initial(C8(b=None), C8_0)
Initial(Bid(b=None, S='u'), Bid_0)

from pysb.simulator import ScipyOdeSimulator

t = np.linspace(0, 20000)
simres = ScipyOdeSimulator(model, tspan=t).run()
yout = simres.all

n_sim=4
sim = OpenCLSSASimulator(model, tspan=t, verbose=True)
traj = sim.run(tspan=t, number_sim=n_sim)







# import the pysb module and all its methods and functions
from pysb import *
import numpy as np
from pysb.simulator import OpenCLSSASimulator, ScipyOdeSimulator
from pysb.simulator import BngSimulator


# instantiate a model
Model()

# declare monomers
Monomer('C8', ['b'])
Monomer('Bid', ['b', 'S'], {'S': ['u', 't']})

# input the parameter values
Parameter('kf', 1.0e-07)
Parameter('kr', 1.0e-03)
Parameter('kc', 1.0)



Observable('obsC8', C8(b=None))
Observable('obsBid', Bid(b=None, S='u'))
Observable('obstBid', Bid(b=None, S='t'))


#Expression('kg', (kf+kc)) #expression! OpenCLSSASimulator will work with an expression of this type
Expression('kg', kf/(kf+obsBid)) #expression! OpenCLSSASimulator won't work with an expression of this type

# now input the rules
Rule('C8_Bid_bind', C8(b=None) + Bid(b=None, S='u') | C8(b=1) % Bid(b=1, S='u'), kf, kr)
Rule('tBid_from_C8Bid', C8(b=1) % Bid(b=1, S='u') >> C8(b=None) + Bid(b=None, S='t'), kg)

model.monomers
model.parameters
model.rules

Parameter('C8_0', 1000)
Parameter('Bid_0', 10000)
Initial(C8(b=None), C8_0)
Initial(Bid(b=None, S='u'), Bid_0)

#Check if the simulation ODE works 
t = np.linspace(0, 20000)
simres = ScipyOdeSimulator(model, tspan=t).run()
yout = simres.all

#Check bng
sim = BngSimulator(model, tspan=t, verbose=True)
simres = sim.run(method='ssa', verbose=True)
yout = simres.observables

#And check openCLSSASimulator
n_sim=4
sim = OpenCLSSASimulator(model, tspan=t, verbose=True)
traj = sim.run(tspan=t, number_sim=n_sim)



###############################################################################
#test w cell cycle



from runStochCell import runStochcell
from peakfinder import peakfinder

from pysb.simulator import OpenCLSSASimulator #from gpu_ssa branch in Pysb github repos




pysb_sampled_parameter_names = ["V_s2skp","V_Fp","k_M","k_EMY","k_CDF","k_CDS","k_ED","K_S","K_MY","K_DF","d_CD","k_P1", "k_P2", "K_CD", "k_E", "K_M", "K_CE", "V_2E2F", "K_2E2F", "k_RE", "kf1", "kf2", "kb", "kasbaxmbcl2", "kdsbaxmbcl2", "k1sbh3", "k2sbh3", "kdbh3", "kasbh3bcl2", "kdsbh3bcl2"]
parameter_vector_original=[]
for pname in pysb_sampled_parameter_names:
    #Change model parameter values to current location in parameter space
    parameter_vector_original+=[model_population.parameters[pname].value]

v = 0.05;     
sys_size=200; #cell size
thr = 0.26
run_time = 120
observed_cells_time = run_time

idx = 1
condition =''

variable_range, check, tmp_score = runStochcell(model_population, model_concentration, idx, condition, thr);

check_ode = check==False or np.any(np.array([variable_range[key] for key in variable_range.keys()])<0) or tmp_score>0.2 #the model doesn't oscillate or the 

cells={}
cells['molecules']=[]; #Number of molecules taken randomly at one point during the cell cycle
cells['time']=[]; #
cells['division']=[]; #
cells['living_cell']=[]; #living cell?
cells['recording']=[];

#plt.plot(variable_range['obsCdc20'])

init = 2

for i in range(init) :
    cells['molecules']=cells['molecules']+[{}]
    a=np.random.randint(len(variable_range['obsSkp2']))
    for var in list(variable_range.keys()) :
        cells['molecules'][i][var] = [round(variable_range[var][a]*200)]
    cells['time']=cells['time'] + [[0]]
    cells['division']=cells['division'] + [0]
    cells['living_cell']=cells['living_cell'] + [True] #living cell = True
    cells['recording']=cells['recording'] + [True] #recording the cell = True (if it isn't dead and if no event didn't occur after observed_cells_time)

cell_count=np.array([[0., len(cells['molecules'])]]);   #start the cell counter at time zero and with number of initial cells
current_time=np.transpose(np.concatenate([np.zeros((1,len(cells['molecules']))), np.ones((1,len(cells['molecules'])))])); #keeps track of the current time of each cell and whether it's dead or not 


#index of cell we still have to simulate
index_cell_recording = [i for i, x in enumerate(cells['recording']) if x]

#index of cells that are below the 'observed time', in which we take into account new cell births
cell_below_observed_time = current_time[:,0]<observed_cells_time
index_cell_below_observed_time = [i for i, x in enumerate(cell_below_observed_time) if x]

#Initial molecule numbers for the different cells to be simulated
Initials = []
for ii in index_cell_recording:
    Initial=[]
    for j in cells['molecules'][ii] :
        Initial=Initial+[cells['molecules'][ii][j][-1]]
    Initials+=[Initial]

cell_count_tmp=np.zeros((1,2)); #temporary cell count with current cell time and whether a cell was added, taken away, or not


T = 1
tmp_t = np.linspace(0, T, int(T*10))


simulator = OpenCLSSASimulator(model_population, verbose=3, initials=Initials)
simres = simulator.run(tspan = tmp_t, number_sim=len(index_cell_recording));



simres.dataframe['obsBAX'].unstack(0).values


