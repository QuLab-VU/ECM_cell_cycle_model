"""
Created on Wed Feb  5 11:52:01 2020
@author: Jérémie Marlhens
"""
###############################################################################
import os
import sys
path = os.path.dirname(os.path.abspath(__file__)) #add a pathway to the model
sys.path.insert(0, path + '/../') #add a pathway to the model
sys.path.insert(0, path + '/../model_data') #add a pathway to the model
# import model
sys.path.insert(0, path + '/../cellular_cycle') #add a pathway to the model

from cell_cycle_death_model import model
import pandas as pd

from pysb.simulator import ScipyOdeSimulator
from pysb.simulator import BngSimulator
# from pysb.simulator import StochKitSimulator
# from pysb.bng import generate_equations
# from pysb.logging import EXTENDED_DEBUG
# from pysb.export import export

import numpy as np
import matplotlib.pyplot as plt

from Python.obsolete.Austin_proliferation_model_6 import model as model_previous

###############################################################################
#Compare with matlab simulations

simulation_ode_matlab = {}
simulation_ode_matlab['55 kPa DMSO'] = np.genfromtxt('/Users/jeremie/Downloads/ode_100_100_100.csv', delimiter=',').T

# simulation_ode_matlab['55 kPa ERL'] = np.genfromtxt(path + '/../../Matlab/v2/ode_060_020_100.csv', delimiter=',').T
# simulation_ode_matlab['4 kPa DMSO'] = np.genfromtxt(path + '/../../Matlab/v2/ode_100_100_060.csv', delimiter=',').T
# simulation_ode_matlab['0.6 kPa DMSO'] = np.genfromtxt(path + '/../../Matlab/v2/ode_060_100_040.csv', delimiter=',').T


# stoch_run = {}
# for file_tmp, condition in zip(['dmso_long_4_29_1_1', 'dmso_long_4_29_2_1', 'dmso_long_4_29_3_1', 'erl_long_4_29_1_1', 'erl_long_4_29_2_1', 'erl_long_4_29_3_1'], ['06 kPa DMSO', '4 kPa DMSO', '55 kPa DMSO', '06 kPa ERL', '4 kPa ERL', '55 kPa ERL']):
#     stoch_run[condition] = np.genfromtxt(r"C:\Users\Grim.M\Ecole\Ens\Master\S2\Stage\Code\matlab_runs/" +  file_tmp + ".txt", delimiter=',', skip_footer=450000).T
#
# plt.plot(stoch_run['55 kPa DMSO'][3, :10000])

###############################################################################

# import pysb sampled parameters
sampled_parameter_names = np.loadtxt(path + "/../model_data/sampled_parameter_names.txt", dtype=str)

# Importing the parameters:
new_param = pd.read_csv(path + "/../model_data/parameters.csv")
del new_param['Stress']
for key in new_param.keys():
    try:
        if model.parameters[key].value != new_param[key][0]:
            model.parameters[key].value = new_param[key]
            model_previous.parameters[key].value = new_param[key]
    except: print(key)
model.parameters['omega'].value = 200
model_previous.parameters['omega'].value = 200

# Setting initial concentrations '/../model_data/initial_concentrations.csv'
initial_concentrations = pd.read_csv(path + '/../model_data/initial_concentrations.csv')
for key in initial_concentrations.keys():
    try:
        model.parameters[key].value = initial_concentrations[key][0] / model.parameters['omega'].value
        model_previous.parameters[key].value = initial_concentrations[key][0]
    except: print(key)

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
condition = '55 kPa DMSO'
# condition = '0.6 kPa ERL'

for input_signal in inputs.keys():
    model.parameters[input_signal].value = inputs[input_signal][condition]
    model_previous.parameters[input_signal].value = inputs[input_signal][condition] * 200

t = np.linspace(0, 1000, 100001)
t = np.linspace(0, 100, 10001)
simres = ScipyOdeSimulator(model, t).run()
yout = simres.all

observables=['obsCycE','obsCdc20', 'obsCycA', 'obsCycB', 'obsp21']
for obs in range(4) :
    plt.plot(t, yout[observables[obs]], label=observables[obs], lw=2)
    print(yout[observables[obs]][0])
plt.legend(loc=0)
plt.show()

# Simulation with previous model
t = np.linspace(0, 100, 10001)
simres = ScipyOdeSimulator(model_previous, t).run()
yout_previous = simres.all


####Comparison with Matlab deterministic simulation

def comparaison_Matlab_Pysb_plot(yout, sim_matlab, omega):
    nmax = 10000
    obs = ['obsSkp2', 'obsp21', 'obsp21p', 'obscdh1act', 'obsMyc', 'obsCycD', 'obsCycDp21', 'obspRb', 'obspRbp', 'obsE2F', 'obspRbE2F', 'obsCycE', 'obsCycEp21', 'obsCycA', 'obsCycAp21', 'obsCycB', 'obsCycBp21', 'obsCdc20', 'obsBAX', 'obsBAXmBCL2', 'obsBH3','obsBH3BCL2','obsBAXm','obsBCL2']

    fig, axs = plt.subplots(6, 4)

    for i in range(4):
        for j in range(6):
            axs[j, i].plot(t[0:nmax], yout[obs[j+i*6]][0:nmax]/omega, label=obs[j+i*6] +' PySB', lw=0.1)
            axs[j, i].plot(t[0:nmax], sim_matlab[j+i*6][0:nmax], label=obs[j+i*6] +' Matlab', lw=0.1)
            axs[j, i].legend(fontsize=2)

    for ax in axs.flat:
        ax.label_outer()

    return fig, axs


def comparison_two_simulations(sim1, sim2):
    nmax = 10000


    obs = ['obsSkp2', 'obsp21', 'obsp21p', 'obscdh1act', 'obsMyc', 'obsCycD', 'obsCycDp21', 'obspRb', 'obspRbp', 'obsE2F', 'obspRbE2F', 'obsCycE', 'obsCycEp21', 'obsCycA', 'obsCycAp21', 'obsCycB', 'obsCycBp21', 'obsCdc20', 'obsBAX', 'obsBAXmBCL2', 'obsBH3','obsBH3BCL2','obsBAXm','obsBCL2']
    fig, axs = plt.subplots(6, 4)

    for i in range(4):
        for j in range(6):
            axs[j, i].plot(t[0:nmax], sim1[obs[j+i*6]][0:nmax], label=obs[j+i*6] +' sim new', lw=0.1)
            axs[j, i].plot(t[0:nmax], sim2[obs[j+i*6]][0:nmax], label=obs[j+i*6] +' sim old', lw=0.1)
            axs[j, i].legend(fontsize=2)

    for ax in axs.flat:
        ax.label_outer()

    return fig, axs

fig, axs = comparison_two_simulations(yout, yout_previous)
plt.savefig('/Users/jeremie/Downloads/comparison_pysb_pysb_55_kPa_DMSO.pdf', bbox_inches='tight')
plt.show()

sim_matlab = simulation_ode_matlab["55 kPa DMSO"]
fig, axs = comparaison_Matlab_Pysb_plot(yout, sim_matlab, model.parameters['omega'].value)
plt.savefig('/Users/jeremie/Downloads/comparison_matlab_pysb_55_kPa_DMSO.pdf', bbox_inches='tight')
plt.show()

###############################################################################
###############################################################################

####Stochastic simulation

condition = '55 kPa DMSO'


for input_signal in inputs.keys():
    model.parameters[input_signal].value = inputs[input_signal][condition]


t = np.linspace(0, 500, 50001)

sim = BngSimulator(model, tspan=t, verbose=True)
simres = sim.run(method='ssa', verbose=True)
yout = simres.observables

observables=['obsCycE','obsCdc20', 'obsCycD']
observables=['obsCycE','obsCdc20', 'obsCycA', 'obsCycB', 'obsp21']

for obs in range(len(observables)) :
    plt.plot(t,yout[observables[obs]], label=observables[obs], lw=2)
plt.legend(loc=0)
plt.show()


####Comparison with Matlab Stochastic simulation

def comparaison_Matlab_PySB_stoch():
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

#
# import math
#
# def RBS(x,x0,sigma) :
#     observations , dimensions= np.shape(x0)
#     s=0
#     for i in range(observations):
#         s+=math.exp(-np.linalg.norm(x-x0[i])**2/(2*sigma**2)) / (sigma*math.sqrt(2*math.pi))
#     if observations != 0:
#         return s/observations
#     else:
#         return s
#
#
# from pysb.simulator import BngSimulator
# from SMs_G1s_simulation_BngSimulator import IMTs
# from runStochCell import runStochcell
#
#
# tmp_t=pl.linspace(0, 4000, 8000)
# sim = BngSimulator(model, tspan=tmp_t, verbose=True)
# traj = sim.run()
#
# variable_range = runStochcell(model_concentration)
#
#
# ini={}
# print(len(variable_range['obsBCL2']));
# a=np.random.randint(len(variable_range['obsBCL2']))
# for var in list(variable_range.keys()) :
#     ini[var] = [round(variable_range[var][a]*200)]
#
# Initial=[]
# for j in ini:
#     Initial=Initial+[ini[j][-1]]
# initials=np.array(Initial)
#
# (a,b) = IMTs(tmp_t, initials, model_population)
#
# sigma=2
# x=np.linspace(0,80,800)
# y_sim=[RBS(i,np.transpose(np.array([b])),sigma) for i in x]
# # Plot the KDE as fitted using the default arguments
# bins_list = np.linspace(0,90,25)
# plt.plot(x, y_sim, lw=3, label = "G1 simulations")
# plt.legend(loc='best')
#
# plt.savefig(r'C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\report\figures\stoch_simulations\distrib_G1.pdf')
#
#
# sigma=2
# x=np.linspace(0,80,800)
# y_sim=[RBS(i,np.transpose(np.array([a])),sigma) for i in x]
# # Plot the KDE as fitted using the default arguments
# bins_list = np.linspace(0,90,25)
# plt.plot(x, y_sim, lw=3, label = "S-M simulations")
# plt.legend(loc='best')
#
# plt.savefig(r'C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\report\figures\stoch_simulations\distrib_SM.pdf')
#

####################################
#find differences with matlab
#
# from numpy import diff
# observables=['obsCycE','obsCdc20']
# dt = t[1] - t[0]
#
# y = yout[observables[obs]]
# dy = diff(y)/dt
#
# derivate={}
#
# for o in range(len(model_concentration.observables)):
#     model_concentration.observables[o].name
#
#     dt = t[1] - t[0]
#     y = yout[model_concentration.observables[o].name]
#     dy = diff(y)/dt
#     derivate[model_concentration.observables[o].name]=dy[0]
#
#
#
# V_1Mb=model_concentration.parameters['V_1Mb'].value
# V_2Mb=model_concentration.parameters['V_2Mb'].value
# K_ie=model_concentration.parameters['K_ie'].value
# K_ib=model_concentration.parameters['K_ib'].value
# b=model_concentration.parameters['b'].value
# b2=model_concentration.parameters['b2'].value
# K_1Mb=model_concentration.parameters['K_1Mb'].value
# K_2Mb=model_concentration.parameters['K_2Mb'].value
# V_1cdc20=model_concentration.parameters['V_1cdc20'].value
# K_1cdc20=model_concentration.parameters['K_1cdc20'].value
# V_2cdc20=model_concentration.parameters['V_2cdc20'].value
# K_2cdc20=model_concentration.parameters['K_2cdc20'].value
# CycB_tot=model_concentration.parameters['CycB_tot'].value
# s=model_concentration.parameters['s'].value
# k_c7=model_concentration.parameters['k_c7'].value
# k_c8=model_concentration.parameters['k_c8'].value
# s=model_concentration.parameters['s'].value
#
# obsCycB=yout['obsCycB'][0]
# obsCycE=yout['obsCycE'][0]
# obsCycA=yout['obsCycA'][0]
# obsCycBp21=yout['obsCycBp21'][0]
# obsp21=yout['obsp21'][0]
# obsCdc20=yout['obsCdc20'][0]
#
#
#
# from scipy.io import loadmat
# y = loadmat('C:/Users/Grim.M/Ecole/ens/Master/S2/Stage/Code/v1_test/cdc20.mat')
# L=np.transpose(y['z'])
# l=L.tolist()[0]
#
# plt.plot(t, l)
# plt.plot(t, yout['obsCdc20'], label='obsCdc20', lw=2)



###############################################################################
###############################################################################