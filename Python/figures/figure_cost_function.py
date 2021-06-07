# -*- coding: utf-8 -*-
"""
Created on Tue Apr 13 12:32:04 2021

@author: Grim.M
"""

path=r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code"

import sys
sys.path.insert(0, path + "/cellular_cycle") #add a pathway to the model
sys.path.insert(0, path + "/fit") #add a pathway 

import matplotlib.gridspec as gridspec
import numpy as np
import pylab as pl
import matplotlib.pyplot as plt
from pathlib import Path
import pandas as pd
from peakfinder import peakfinder


index=1
ode={}
ssa={}
try:
    ode["obsCycD"] = np.load(path + r"\fit_data\sobol_seq\cost_function_model_4\biochemical_data\personal_computer\ODE_timecourse_CycD_index_" + str(index) + ".npy")
    ode["obsCycE"] = np.load(path + r"\fit_data\sobol_seq\cost_function_model_4\biochemical_data\personal_computer\ODE_timecourse_CycE_index_" + str(index) + ".npy")
    ode["obsCdc20"] = np.load(path + r"\fit_data\sobol_seq\cost_function_model_4\biochemical_data\personal_computer\ODE_timecourse_Cdc20_index_" + str(index) + ".npy")
    ssa["obsCycD"] = np.load(path + r"\fit_data\sobol_seq\cost_function_model_4\biochemical_data\personal_computer\ssa_timecourse_CycD_index_" + str(index) + ".npy")
    ssa["obsCycE"] = np.load(path + r"\fit_data\sobol_seq\cost_function_model_4\biochemical_data\personal_computer\ssa_timecourse_CycE_index_" + str(index) + ".npy")
    ssa["obsCdc20"] = np.load(path + r"\fit_data\sobol_seq\cost_function_model_4\biochemical_data\personal_computer\ssa_timecourse_Cdc20_index_" + str(index) + ".npy")
except :
    ode["obsCycD"] = np.load(path + r"\fit_data\sobol_seq\cost_function_model_4\biochemical_data\Tempo\ODE_timecourse_CycD_index_" + str(index) + ".npy")
    ode["obsCycE"] = np.load(path + r"\fit_data\sobol_seq\cost_function_model_4\biochemical_data\Tempo\ODE_timecourse_CycE_index_" + str(index) + ".npy")
    ode["obsCdc20"] = np.load(path + r"\fit_data\sobol_seq\cost_function_model_4\biochemical_data\Tempo\ODE_timecourse_Cdc20_index_" + str(index) + ".npy")
    ssa["obsCycD"] = np.load(path + r"\fit_data\sobol_seq\cost_function_model_4\biochemical_data\Tempo\ssa_timecourse_CycD_index_" + str(index) + ".npy")
    ssa["obsCycE"] = np.load(path + r"\fit_data\sobol_seq\cost_function_model_4\biochemical_data\Tempo\ssa_timecourse_CycE_index_" + str(index) + ".npy")
    ssa["obsCdc20"] = np.load(path + r"\fit_data\sobol_seq\cost_function_model_4\biochemical_data\Tempo\ssa_timecourse_Cdc20_index_" + str(index) + ".npy")
        
omega=200
t = np.linspace(0, 800, 800)

T = t[1] - t[0]  # sampling interval 
N = ode['obsCycE'][600:].size
# 1/T = frequency
f = np.linspace(0, 1 / T, N)
###############################################################################

fn = Path(r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\report\figures\cost_function\cost_function_check.pdf").expanduser()

fig = plt.figure(figsize=(16, 20))
outer = gridspec.GridSpec(2, 3, wspace=0.1, hspace=0.15)
for i in range(3) :
    inner = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=outer[i], wspace=0.1, hspace=0.1)
    obs = ['obsCycD', 'obsCycE', 'obsCdc20'][i]
    ax = plt.Subplot(fig, inner[0])
    (a,max_tmp) = peakfinder(ssa[obs][400:520], 0.15*omega, 0.15*omega, 1, includeEndpoints=1);
    ax.plot(t[600:], ode[obs][600:], linewidth=1, color='orange')
    ax.set_xticks([])
    ax.set_yticks([])
    #ax.legend(fontsize=6)
#    ax.set_ylabel("Concentration")
#    ax.set_xlabel("time (hours)")
#    ax.set_ylabel("Concentration")
#    ax.set_xlabel("time (hours)")
    fig.add_subplot(ax)
    ax = plt.Subplot(fig, inner[1])
    fft = np.fft.fft(ode[obs][600:])
    ax.legend(fontsize=1)
#    ax.set_ylabel("Amplitude")
#    ax.set_xlabel("Frequency [Hz]")
    ax.bar(f[(f>0.005)*(f<0.1)], np.abs(fft)[(f>0.005)*(f<0.1)] * 1 / N, width=(1/750), color=(0.8, 0.4, 0.6, 0.6))  # 1 / N is a normalization factor
    ax.set_xticks([])
    ax.set_yticks([])
    fig.add_subplot(ax)

inner = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=outer[3], wspace=0.1, hspace=0.1)
obs = ['obsCycD', 'obsCycE', 'obsCdc20'][1]
ax = plt.Subplot(fig, inner[0])
(a,max_tmp_CycE) = peakfinder(ssa[obs][400:600], 0.15*omega, 0.15*omega, 1, includeEndpoints=1);
ax.plot(t[400:600], ssa[obs][400:600], linewidth=1, color='orange')
ax.scatter(np.array(max_tmp_CycE)+400, a, color='green', s=20)
#ax.set_ylabel("Pop")
#ax.set_xlabel("time (hours)")
#ax.set_yticklabels([])
#ax.set_xticklabels([])
ax.set_xticks([])
ax.set_yticks([])
fig.add_subplot(ax)

inner = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=outer[4], wspace=0.1, hspace=0.1)
obs = ['obsCycD', 'obsCycE', 'obsCdc20'][2]
ax = plt.Subplot(fig, inner[0])
(a,max_tmp_Cdc20) = peakfinder(ssa[obs][400:600], 0.15*omega, 0.15*omega, 1, includeEndpoints=1);
ax.plot(t[400:600], ssa[obs][400:600], linewidth=1, color='orange')
ax.scatter(np.array(max_tmp_Cdc20)+400, a, color='red', s=20)
#ax.set_ylabel("Pop")
#ax.set_xlabel("time (hours)")
ax.set_xticks([])
ax.set_yticks([])
fig.add_subplot(ax)

inner = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=outer[5], wspace=0.1, hspace=0.1)
ax = plt.Subplot(fig, inner[0])
ax.scatter(np.array(max_tmp_Cdc20)+400, [0 for i in range(len(max_tmp_Cdc20))], color='red', s=30)
ax.scatter(np.array(max_tmp_CycE)+400, [0 for i in range(len(max_tmp_CycE))], color='green', s=30)
ax.set_xticks([])
ax.set_yticks([])
fig.add_subplot(ax)


fig.show()

plt.savefig(fn, bbox_inches='tight')


###############################################################################
#cost calculation
from plot_frac_div_q_ap1 import all_events


import fracprolif as fp

######
#parameters for the simulation
init = 50
run_time = 102.3
observed_cells_time = 30
omega=200


#Experimental
path=r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code"
all_events_exp={}

IMT_eoe_counts = pd.DataFrame(columns=['End_of_Expt', 'Death' ,'Birthtime','Lifespan'])

df = pd.read_csv(path + '/cell_fates/ds8_06_DMSO_fates.csv')
lines = {}
for col in df :
    lines[col]=df[col].dropna().tolist()

for (death, mitotic_time, end_time, lifetime) in zip(lines['Death'], lines['Mitotic Time'], lines['End Time'], lines['Lifetime']):
    tmp_counts = {'End_of_Expt' : end_time==102.3, 'Death' : death==1, 'Birthtime' : mitotic_time, 'Lifespan' : lifetime}
    IMT_eoe_counts = IMT_eoe_counts.append(tmp_counts, ignore_index=True)

all_events_exp['0.6 kPa DMSO'] = all_events(run_time, IMT_eoe_counts)


#calculating rate of fractional proliferation model

mitotic_lifespans_exp = all_events_exp['0.6 kPa DMSO']['Lifespan'][all_events_exp['0.6 kPa DMSO']['event']=='mitosis']
censored_lifespans_exp = all_events_exp['0.6 kPa DMSO']['Lifespan'][all_events_exp['0.6 kPa DMSO']['event']=='End_of_Expt']

number_of_deaths = len(all_events_exp['0.6 kPa DMSO']) - (all_events_exp['0.6 kPa DMSO']['event']=='apoptosis').value_counts()[0]
cell_time_total = sum(all_events_exp['0.6 kPa DMSO']['Lifespan'])


m = fp.mle(mitotic_lifespans_exp, censored_lifespans_exp)

parameters_exp={}
parameters_exp['0.6 kPa DMSO'] ={}

parameters_exp['0.6 kPa DMSO']['a'] = number_of_deaths/cell_time_total
parameters_exp['0.6 kPa DMSO']['mu'] = m['x'][0]
parameters_exp['0.6 kPa DMSO']['sigma'] = m['x'][1]
parameters_exp['0.6 kPa DMSO']['lambda'] = m['x'][2]
parameters_exp['0.6 kPa DMSO'].update(fp.q_rates(m))

parameters_exp['0.6 kPa DMSO']['q']=1e-5
#######################################
#################################################################################
def RBS(x,x0,sigma) :
    observations , dimensions= np.shape(x0)
    s=0
    for i in range(observations):
        s+=math.exp(-np.linalg.norm(x-x0[i])**2/(2*sigma**2)) / (sigma*math.sqrt(2*math.pi))
    if observations != 0:
        return s/observations
    else:
        return s

################################################################################
mitotic_lifespans_exp.tolist()

sigma=5

x=np.linspace(0,120,400)
experimental_IMT_distribution_03kPa_DMSO=[RBS(i, np.transpose([mitotic_lifespans_exp.tolist()]), sigma) for i in x]


simulated_IMT_distribution_03kPa_DMSO = np.load(path + r'\fit_data\sobol_seq\cost_function_model_4\cost_and_IMT_distribution\personal_computer\IMT_distribution_'  + str(index) +"_50_initial_cells.npy")



experimental_SMs={}
experimental_G1s={}
experimental_IMTs={}


df = pd.read_csv(path+'/IMT_times/DS8_SM_0_2_55.csv')

for col in df.columns :
    experimental_SMs[col]=df[col].dropna().tolist()


df = pd.read_csv(path+'/IMT_times/DS8_G1_0_2_55.csv')

for col in df.columns :
    experimental_G1s[col]=df[col].dropna().tolist()

for k in experimental_G1s.keys() :
    minimum=min(len(experimental_G1s[k]),len(experimental_SMs[k]))
    experimental_IMTs[k]=[experimental_G1s[k][i]+experimental_SMs[k][i] for i in range(minimum)]


index=1

g1_distribution = np.load(path + r"\fit_data\sobol_seq\cost_function_model_4\cost_and_IMT_distribution\55_kPa\g1_distribution_" + str(index) + "_30_initial_cells_55kPa_DMSO.npy")
sm_distribution = np.load(path + r"\fit_data\sobol_seq\cost_function_model_4\cost_and_IMT_distribution\55_kPa\sm_distribution_" + str(index) + "_30_initial_cells_55kPa_DMSO.npy")

sigma=4
   
g1_distribution_sim = [RBS(i, np.transpose([g1_distribution]), sigma) for i in np.linspace(0,120,400)]
sm_distribution_sim = [RBS(i, np.transpose([sm_distribution]), sigma) for i in np.linspace(0,120,400)]
 
g1_distribution_exp = [RBS(i, np.transpose([experimental_G1s['55 kPa DMSO']]), sigma) for i in np.linspace(0,120,400)]
sm_distribution_exp = [RBS(i, np.transpose([experimental_SMs['55 kPa DMSO']]), sigma) for i in np.linspace(0,120,400)]
    

fn = Path(r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\report\figures\cost_function\cost_function_cost.pdf").expanduser()

fig, axs = plt.subplots(2, 2)

axs[0, 0].plot(np.linspace(0,120,400), g1_distribution_sim, linewidth=1, color='orange')
axs[0, 0].plot(np.linspace(0,120,400), g1_distribution_exp, linewidth=1, color='blue')
axs[0, 0].set_xticks([])
axs[0, 0].set_yticks([])


axs[0, 1].plot(np.linspace(0,120,400), sm_distribution_sim, linewidth=1, color='orange')
axs[0, 1].plot(np.linspace(0,120,400), sm_distribution_exp, linewidth=1, color='blue')
axs[0, 1].set_xticks([])
axs[0, 1].set_yticks([])


for ax in axs.flat:
    ax.label_outer()
    
    

axs[1, 1].spines["top"].set_visible(False)
axs[1, 1].spines["right"].set_visible(False)
axs[1, 1].spines["left"].set_visible(False)
axs[1, 1].spines["bottom"].set_visible(False)
axs[1, 1].set_xticks([])
axs[1, 1].set_yticks([])




axs[1, 0].plot(np.linspace(0,120,400), experimental_IMT_distribution_03kPa_DMSO, linewidth=1, color='blue')
axs[1, 0].plot(np.linspace(0,120,360), simulated_IMT_distribution_03kPa_DMSO, linewidth=1, color='orange')

axs[1, 0].set_xticks([])
axs[1, 0].set_yticks([])

fig.show()
plt.savefig(fn, bbox_inches='tight')

