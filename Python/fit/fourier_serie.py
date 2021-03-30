# -*- coding: utf-8 -*-
"""
Created on Thu Apr 30 10:14:32 2020

@author: Grim.M
"""

path="/home/jeremie"

#path=r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code"
import sys
sys.path.insert(0, path + "/cellular_cycle") #add a pathway to the model
sys.path.insert(0, path + "/fit") #add a pathway 

##

from numpy.fft import fft
import numpy as np
import math
from matplotlib.pyplot import *
import pylab as pl
from pysb.simulator import BngSimulator
from pysb.simulator import ScipyOdeSimulator
import matplotlib.pyplot as plt
from pathlib import Path
import pandas as pd


import sobol_seq

from Austin_proliferation_model_concentration_4 import model
model_concentration=model
from Austin_proliferation_model_population_4 import model
model_population=model

from runStochCell import runStochcell_check
from runStochCell import runStochcell

from peakfinder import peakfinder


t = np.linspace(0, 800, 800)

pysb_sampled_parameter_names = ["V_s2skp","V_Fp","k_M","k_EMY","k_CDF","k_CDS","k_ED","K_S","K_MY","K_DF","d_CD","k_P1", "k_P2", "K_CD", "k_E", "K_M", "K_CE", "V_2E2F", "K_2E2F", "k_RE", "kf1", "kf2", "kb", "kasbaxmbcl2", "kdsbaxmbcl2", "k1sbh3", "k2sbh3", "kdbh3", "kasbh3bcl2", "kdsbh3bcl2"]
loc=[]
for parameter in pysb_sampled_parameter_names:
    print((model.parameters[parameter].value))
    loc+=[np.log10(model.parameters[parameter].value)]

#
#cost = np.load(path+"/fit_data/cost_sobol_seq_0000_44000.npy")
#index = [idx for idx, element in enumerate(cost) if element>-6]
#cost_best = [element for idx, element in enumerate(cost) if element>-6]

l=sobol_seq.i4_sobol_generate(30, 120000, skip=0)
l2 = 4*l + (np.array(loc)-2)

index=1

parameter_vector=l2[index-1]
l2[0]
#parameter_vector = Cdc20_or_CycD_not_cycling_parameters_idx_16279
parameter_vector = [10**i for i in parameter_vector]

omega=200
model_concentration.parameters['ERK'].value = 1.06970165050038 
model_concentration.parameters['FAK'].value = 0.337143388171177
model_concentration.parameters['AKT'].value = 0.616658111482869

model_population.parameters['ERK'].value = 1.06970165050038 *200
model_population.parameters['FAK'].value = 0.337143388171177 *200
model_population.parameters['AKT'].value = 0.616658111482869 *200

param_dict = {pname: pvalue for pname, pvalue in zip(pysb_sampled_parameter_names, parameter_vector)}
for pname, pvalue in param_dict.items(): 
    model_concentration.parameters[pname].value = (pvalue)
    model_population.parameters[pname].value = (pvalue)
    
simres = ScipyOdeSimulator(model_concentration, t, verbose=True).run();
simres_2 = BngSimulator(model_population, t, verbose=True).run();

plt.plot(simres.observables['obsCdc20'])
plt.plot(simres.observables['obsCycE'])

yout=simres.all


############################

#saving  biochemical timecourse

"""
omega=200
model_concentration.parameters['ERK'].value = 1.06970165050038 
model_concentration.parameters['FAK'].value = 0.337143388171177
model_concentration.parameters['AKT'].value = 0.616658111482869

model_population.parameters['ERK'].value = 1.06970165050038 *200
model_population.parameters['FAK'].value = 0.337143388171177 *200
model_population.parameters['AKT'].value = 0.616658111482869 *200

t = np.linspace(0, 800, 800)

for index in frac_prolif_parameters['index'].tolist():
    index=int(index)
    parameter_vector=l2[index-1]
    parameter_vector = [10**i for i in parameter_vector]
    param_dict = {pname: pvalue for pname, pvalue in zip(pysb_sampled_parameter_names, parameter_vector)}
    for pname, pvalue in param_dict.items(): 
        model_concentration.parameters[pname].value = (pvalue)
        model_population.parameters[pname].value = (pvalue)
        
    simres = ScipyOdeSimulator(model_concentration, t, verbose=True).run();
    simres_2 = BngSimulator(model_population, t, verbose=True).run();

    np.save(r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code\fit_data\sobol_seq\cost_function_model_4\biochemical_data\ODE_timecourse_CycE_index_" + str(index), simres.observables['obsCycE'])
    np.save(r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code\fit_data\sobol_seq\cost_function_model_4\biochemical_data\ODE_timecourse_CycD_index_" + str(index), simres.observables['obsCycD'])
    np.save(r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code\fit_data\sobol_seq\cost_function_model_4\biochemical_data\ODE_timecourse_Cdc20_index_" + str(index) , simres.observables['obsCdc20'])
    
    np.save(r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code\fit_data\sobol_seq\cost_function_model_4\biochemical_data\ssa_timecourse_CycE_index_" + str(index), simres_2.observables['obsCycE'])
    np.save(r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code\fit_data\sobol_seq\cost_function_model_4\biochemical_data\ssa_timecourse_CycD_index_" + str(index) , simres_2.observables['obsCycD'])
    np.save(r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code\fit_data\sobol_seq\cost_function_model_4\biochemical_data\ssa_timecourse_Cdc20_index_" + str(index) , simres_2.observables['obsCdc20'])


"""
##########################################
#fft

T = t[1] - t[0]  # sampling interval 
N = yout['obsCycE'][600:].size
# 1/T = frequency
f = np.linspace(0, 1 / T, N)



fn = Path(r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\report\figures\fourier_transform\fft_CycD_Cdc20_sobol_seq_idx_" + str(index+1) + ".pdf").expanduser()

fig, axs = plt.subplots(3, 3)
for obs,i in zip(['obsCycD', 'obsCycE', 'obsCdc20'],range(3)):
    
    (a,max_tmp) = peakfinder(simres_2.observables[obs], 0.15*omega, 0.15*omega, 1, includeEndpoints=1);
    axs[i, 2].plot(t[0:], simres_2.observables[obs][0:], linewidth=0.5)
    axs[i, 2].scatter(max_tmp, a, color='red', s=0.1)
    #for peak in max_tmp:
    #    axs[i, 2].plot([peak for i in range(int(max(simres_2.observables[obs])))], range(int(max(simres_2.observables[obs]))), color='red')
    #axs[i, 2].legend(fontsize=6)
    axs[i, 2].set_ylabel("Pop")
    axs[i, 2].set_xlabel("time (hours)")
    
    
    axs[i, 1].plot(t[600:], yout[obs][600:], label = obs, linewidth=0.5)
    axs[i, 1].legend(fontsize=6)
    axs[i, 1].set_ylabel("Concentration")
    axs[i, 1].set_xlabel("time (hours)")
    
    fft = np.fft.fft(yout[obs][600:])
    axs[i, 0].legend(fontsize=1)
    axs[i, 0].set_ylabel("Amplitude")
    axs[i, 0].set_xlabel("Frequency [Hz]")
    axs[i, 0].bar(f[(f>0.005)*(f<0.1)], np.abs(fft)[(f>0.005)*(f<0.1)] * 1 / N, width=(1/750), color=(0.8, 0.4, 0.6, 0.6))  # 1 / N is a normalization factor
    
st = fig.suptitle("Sobol_seq_idx_" + str(index+1), fontsize="x-large")

for ax in axs.flat:
    ax.label_outer()
    
plt.savefig(fn, bbox_inches='tight')



##########################################################################
##########################################################################
#loading the whole dataset
path=r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code"

frac_prolif_parameters = pd.read_pickle(path + r"/fit_data/sobol_seq/cost_function_model_4/cost_and_IMT_distribution/frac_prolif_parameters_all_parameters_50_initial_cells.pkl")

IMT_distribution = pd.read_pickle(path + r"/fit_data/sobol_seq/cost_function_model_4/cost_and_IMT_distribution/IMT_distribution_all_parameters_50_initial_cells.pkl")

frac_prolif_parameters.sort_values(by=['cost'])[0:16]['index']

t = np.linspace(0, 800, 800)


import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec



fn = Path(r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\report\figures\fourier_transform\fft_CycD_Cdc20_sobol_seq_idx_1-100_best.pdf").expanduser()


fig = plt.figure(figsize=(50, 40))
outer = gridspec.GridSpec(10, 10, wspace=0.1, hspace=0.15)

for i in range(100):
    #setting parameters
    index = int(frac_prolif_parameters.sort_values(by=['cost'])[0:100]['index'].tolist()[i])
    #or
    #index = int(frac_prolif_parameters[300:400]['index'].tolist()[i])
    parameter_vector=l2[index-1]
    #or
    #index = int(frac_prolif_parameters.sort_values(by=['cost'])[-37:-1]['index'].tolist()[i])
    
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
        
        
#    #simulation
#    simres = ScipyOdeSimulator(model_concentration, t, verbose=True).run();
#    simres_2 = BngSimulator(model_population, t, verbose=True).run();

    #fft
    T = t[1] - t[0]  # sampling interval 
    N = ode['obsCycE'][600:].size
    # 1/T = frequency
    f = np.linspace(0, 1 / T, N)

    inner = gridspec.GridSpecFromSubplotSpec(3, 3,
                    subplot_spec=outer[i], wspace=0.1, hspace=0.1)
    
    (a,max_tmp_cdc) = peakfinder(ssa['obsCdc20'], 0.15*omega, 0.15*omega, 1, includeEndpoints=1);
    (a,max_tmp_cycE) = peakfinder(ssa['obsCycE'], 0.15*omega, 0.15*omega, 1, includeEndpoints=1);
    
    
    (check, tmp_score) = runStochcell_check(ssa, ode, index)
    
    for j in range(3):
        obs = ['obsCycD', 'obsCycE', 'obsCdc20'][j]
        ax = plt.Subplot(fig, inner[3*j+2])
        (a,max_tmp) = peakfinder(ssa[obs][400:520], 0.15*omega, 0.15*omega, 1, includeEndpoints=1);
        
        ax.plot(t[400:520], ssa[obs][400:520], linewidth=0.5)
        ax.scatter(np.array(max_tmp)+400, a, color='red', s=0.1)
#        ax.set_ylabel("Pop")
#        ax.set_xlabel("time (hours)")
        ax.set_yticklabels([])
        ax.set_xticklabels([])
        ax.set_xticks([])
        ax.set_yticks([])
        fig.add_subplot(ax)
        
        ax = plt.Subplot(fig, inner[3*j+1])
        if obs== 'obsCycD':
            ax.set_title("index_"+str(index) + " Check:" + str(check) + " Tmp score:"+ str(int(100*tmp_score)))
        ax.plot(t[600:], ode[obs][600:], label = obs, linewidth=0.5)
        ax.legend(fontsize=6)
        ax.set_yticklabels([])
        ax.set_xticklabels([])
#        ax.set_ylabel("Concentration")
#        ax.set_xlabel("time (hours)")
        fig.add_subplot(ax)

        ax = plt.Subplot(fig, inner[3*j+0])
        fft = np.fft.fft(ode[obs][600:])
        ax.legend(fontsize=1)
        ax.set_yticklabels([])
        ax.set_xticklabels([])
#        ax.set_ylabel("Amplitude")
#        ax.set_xlabel("Frequency [Hz]")
        ax.bar(f[(f>0.005)*(f<0.1)], np.abs(fft)[(f>0.005)*(f<0.1)] * 1 / N, width=(1/750), color=(0.8, 0.4, 0.6, 0.6))  # 1 / N is a normalization factor
        #st = ax.suptitle("Sobol_seq_idx_" + str(index+1), fontsize="x-large")
        fig.add_subplot(ax)

plt.savefig(fn, bbox_inches='tight')


fig.show()


#####################

#Remove bad tmp scores (old)
path=r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code"


bad=[]

for i in range(len(frac_prolif_parameters)):
    #setting parameters
    #index = int(frac_prolif_parameters.sort_values(by=['cost'])[0:100]['index'].tolist()[i])
    #or
    
    if not np.isnan(frac_prolif_parameters[:]['index'].tolist()[i]) :
        index = int(frac_prolif_parameters[:]['index'].tolist()[i])
        
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
            try :
                ode["obsCycD"] = np.load(path + r"\fit_data\sobol_seq\cost_function_model_4\biochemical_data\Tempo\ODE_timecourse_CycD_index_" + str(index) + ".npy")
                ode["obsCycE"] = np.load(path + r"\fit_data\sobol_seq\cost_function_model_4\biochemical_data\Tempo\ODE_timecourse_CycE_index_" + str(index) + ".npy")
                ode["obsCdc20"] = np.load(path + r"\fit_data\sobol_seq\cost_function_model_4\biochemical_data\Tempo\ODE_timecourse_Cdc20_index_" + str(index) + ".npy")
                ssa["obsCycD"] = np.load(path + r"\fit_data\sobol_seq\cost_function_model_4\biochemical_data\Tempo\ssa_timecourse_CycD_index_" + str(index) + ".npy")
                ssa["obsCycE"] = np.load(path + r"\fit_data\sobol_seq\cost_function_model_4\biochemical_data\Tempo\ssa_timecourse_CycE_index_" + str(index) + ".npy")
                ssa["obsCdc20"] = np.load(path + r"\fit_data\sobol_seq\cost_function_model_4\biochemical_data\Tempo\ssa_timecourse_Cdc20_index_" + str(index) + ".npy")
            except:
                print('n')
         
        if ssa!={}:
            (check, tmp_score) = runStochcell_check(ssa, ode, index)
        
            if tmp_score>0.1 :
                bad+=[index]

#remove them fro IMT and fracprolifparameters



for index in bad :
    frac_prolif_parameters = frac_prolif_parameters[frac_prolif_parameters['index'] != index]
    IMT_distribution.drop(index, axis=1)



###########################
#clustering



from tslearn.metrics import dtw
from tslearn.clustering import TimeSeriesKMeans
from tslearn.metrics import soft_dtw

dtw_score = dtw(ssa["obsCdc20"], ssa["obsCycD"])



model = TimeSeriesKMeans(n_clusters=5, metric="dtw", max_iter=10)
model.fit(Data)


Data = []

for i in range(len(frac_prolif_parameters)):

    index = int(frac_prolif_parameters['index'].tolist()[i])
        
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
    
    #(check, tmp_score) = runStochcell_check(ssa, ode, index)
    (a,max_tmp_cdc) = peakfinder(ssa["obsCdc20"][60:], 0.15*omega, 0.15*omega, 1, includeEndpoints=1);
    max_tmp_cdc = 60+np.array(max_tmp_cdc)
    
    if len(max_tmp_cdc)>2 and max_tmp_cdc[0]+100<800:
        Data += [np.array([ssa["obsCdc20"][max_tmp_cdc[0]:max_tmp_cdc[0]+80]]).T]


Data = np.array(Data)

model = TimeSeriesKMeans(n_clusters=40, metric="dtw", max_iter=10)

model = TimeSeriesKMeans(n_clusters=10, metric="softdtw", max_iter=10, metric_params={"gamma": .1})


model.fit(Data)


y_pred=model.predict(Data)



fn = Path(r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\report\figures\fourier_transform\clusters_softdtw_005_10_clusters.pdf").expanduser()

fig, axs = plt.subplots(2, 5)
for i in range(2):
    for j in range(5):
        
        yi=i*5+j
        for xx in Data[y_pred == yi][0:10]:
            axs[i, j].plot(xx.ravel(), "k-", alpha=.2,  linewidth=0.1)
        axs[i, j].plot(model.cluster_centers_[yi].ravel(), "r-",  linewidth=0.5)
        axs[i, j].set_ylabel("Pop")
        axs[i, j].set_xlabel("time (hours)")
    
for ax in axs.flat:
    ax.label_outer()
    
plt.savefig(fn, bbox_inches='tight')


"""

from tslearn.metrics import dtw
from tslearn.clustering import TimeSeriesKMeans
from tslearn.metrics import soft_dtw


dtw_score = dtw(ssa["obsCdc20"], ssa["obsCycD"])



model = TimeSeriesKMeans(n_clusters=4, metric="euclidean", max_iter=10)
model.fit(ffts_obsCdc20.to_numpy().T)


y_pred=model.predict(ffts_obsCdc20.T)





fn = Path(r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\report\figures\fourier_transform\clusters_softdtw_005_10_clusters.pdf").expanduser()

fig, axs = plt.subplots(2, 5)
for i in range(2):
    for j in range(5):
        
        yi=i*5+j
        for xx in ffts_obsCdc20.T[y_pred == yi][0:10]:
            axs[i, j].plot(xx.ravel(), "k-", alpha=.2,  linewidth=0.1)
        axs[i, j].plot(model.cluster_centers_[yi].ravel(), "r-",  linewidth=0.5)
        axs[i, j].set_ylabel("Pop")
        axs[i, j].set_xlabel("time (hours)")
    
for ax in axs.flat:
    ax.label_outer()
    
plt.savefig(fn, bbox_inches='tight')


##############################

#Fourier series on stochastic data -> clustering??




fn = Path(r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\report\figures\fourier_transform\fft_stoch_CycD_Cdc20_sobol_seq_idx_1-100_best.pdf").expanduser()


fig = plt.figure(figsize=(50, 40))
outer = gridspec.GridSpec(10, 10, wspace=0.1, hspace=0.15)

for i in range(100):
    #setting parameters
    index = int(frac_prolif_parameters.sort_values(by=['cost'])[0:100]['index'].tolist()[i])
    #or
    #index = int(frac_prolif_parameters[300:400]['index'].tolist()[i])
    parameter_vector=l2[index-1]
    #or
    #index = int(frac_prolif_parameters.sort_values(by=['cost'])[-37:-1]['index'].tolist()[i])
    
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
        
        
#    #simulation
#    simres = ScipyOdeSimulator(model_concentration, t, verbose=True).run();
#    simres_2 = BngSimulator(model_population, t, verbose=True).run();

    #fft
    T = t[1] - t[0]  # sampling interval 
    N = ssa['obsCycE'][200:].size
    # 1/T = frequency
    f = np.linspace(0, 1 / T, N)

    inner = gridspec.GridSpecFromSubplotSpec(3, 3,
                    subplot_spec=outer[i], wspace=0.1, hspace=0.1)
    
    (a,max_tmp_cdc) = peakfinder(ssa['obsCdc20'], 0.15*omega, 0.15*omega, 1, includeEndpoints=1);
    (a,max_tmp_cycE) = peakfinder(ssa['obsCycE'], 0.15*omega, 0.15*omega, 1, includeEndpoints=1);
    
    
    (check, tmp_score) = runStochcell_check(ssa, ode, index)
    
    for j in range(3):
        obs = ['obsCycD', 'obsCycE', 'obsCdc20'][j]
        ax = plt.Subplot(fig, inner[3*j+2])
        (a,max_tmp) = peakfinder(ssa[obs][400:520], 0.15*omega, 0.15*omega, 1, includeEndpoints=1);
        
        ax.plot(t[400:520], ssa[obs][400:520], linewidth=0.5)
        ax.scatter(np.array(max_tmp)+400, a, color='red', s=0.1)
#        ax.set_ylabel("Pop")
#        ax.set_xlabel("time (hours)")
        ax.set_yticklabels([])
        ax.set_xticklabels([])
        ax.set_xticks([])
        ax.set_yticks([])
        fig.add_subplot(ax)
        
        ax = plt.Subplot(fig, inner[3*j+1])
        if obs== 'obsCycD':
            ax.set_title("index_"+str(index) + " Check:" + str(check) + " Tmp score:"+ str(int(100*tmp_score)))
        ax.plot(t[600:], ode[obs][600:], label = obs, linewidth=0.5)
        ax.legend(fontsize=6)
        ax.set_yticklabels([])
        ax.set_xticklabels([])
#        ax.set_ylabel("Concentration")
#        ax.set_xlabel("time (hours)")
        fig.add_subplot(ax)

        ax = plt.Subplot(fig, inner[3*j+0])
        fft = np.fft.fft(ssa[obs][200:])
        ax.legend(fontsize=1)
        ax.set_yticklabels([])
        ax.set_xticklabels([])
#        ax.set_ylabel("Amplitude")
#        ax.set_xlabel("Frequency [Hz]")
        ax.bar(f[(f>0.005)*(f<0.1)], np.abs(fft)[(f>0.005)*(f<0.1)] * 1 / N, width=(1/750), color=(0.8, 0.4, 0.6, 0.6))  # 1 / N is a normalization factor
        #st = ax.suptitle("Sobol_seq_idx_" + str(index+1), fontsize="x-large")
        fig.add_subplot(ax)

plt.savefig(fn, bbox_inches='tight')



##############################

#save fft on stochastic data


for i in range(len(frac_prolif_parameters)):
    #setting parameters
    #or
    if not np.isnan(frac_prolif_parameters[:]['index'].tolist()[i]) :
    
        index = int(frac_prolif_parameters['index'].tolist()[i])
        
        parameter_vector=l2[index-1]
        #or
        #index = int(frac_prolif_parameters.sort_values(by=['cost'])[-37:-1]['index'].tolist()[i])
        
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
            
            
    #    #simulation
    #    simres = ScipyOdeSimulator(model_concentration, t, verbose=True).run();
    #    simres_2 = BngSimulator(model_population, t, verbose=True).run();
    
        #fft
        T = t[1] - t[0]  # sampling interval 
        N = ssa['obsCycE'][200:].size
        # 1/T = frequency
        f = np.linspace(0, 1 / T, N)
        obs = ['obsCycD', 'obsCycE', 'obsCdc20']
        
        for j in range(3):
            
            try:
                np.load(r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code\fit_data\sobol_seq\cost_function_model_4\fft\stochastic_simulation\fft_stochastic_simulation_" + obs[j] + "_200_to_800_hours_index_" + str(index) + ".npy")
                
            except:
                
                fft = np.fft.fft(ssa[obs[j]][200:])
                fft = np.abs(fft)
                
                np.save(r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code\fit_data\sobol_seq\cost_function_model_4\fft\stochastic_simulation\fft_stochastic_simulation_" + obs[j] + "_200_to_800_hours_index_" + str(index), fft)










