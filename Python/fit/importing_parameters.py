# -*- coding: utf-8 -*-
"""
Created on Tue Mar 30 15:17:10 2021

@author: Grim.M
"""

import numpy as np
import pandas as pd
import sobol_seq
import pylab as pl

##
from Austin_proliferation_model_concentration_5 import model
model_concentration=model
from Austin_proliferation_model_population_5 import model
model_population=model

from pysb.simulator import BngSimulator
from pysb.simulator import StochKitSimulator
from pysb.simulator import ScipyOdeSimulator

#import the model and its parameters

path=r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code"

#frac_prolif_parameters = pd.read_pickle(path + r"/fit_data/sobol_seq/cost_function_model_4/cost_and_IMT_distribution/frac_prolif_parameters_all_parameters_50_initial_cells.pkl")

frac_prolif_parameters = pd.read_pickle(path + r"/fit_data/sobol_seq/cost_function_model_5/cost_and_IMT_distribution/frac_prolif_parameters_all_parameters_006_55_kPa.pkl")


IMT_distribution = pd.read_pickle(path + r"/fit_data/sobol_seq/cost_function_model_5/cost_and_IMT_distribution/IMT_distribution_all_parameters_50_initial_cells.pkl")

#import parameters and their costs

frac_prolif_parameters.sort_values(by=['cost_total'])

#sort the parameters according to their glocal cost
#frac_prolif_parameters.sort_values(by=['cost_mmd'])

##########

#ploting cost distribution

#setting the global cost
sum_column = frac_prolif_parameters["cost_006"] +frac_prolif_parameters["cost_frac_prolif_model_006"] + 0 * frac_prolif_parameters["cost_mmd_55"]/2

frac_prolif_parameters["cost_total"] = sum_column

#different plots:
sns.distplot(frac_prolif_parameters['cost_mmd_006'])
#sns.distplot(frac_prolif_parameters['cost_frac_prolif_model_006'])
sns.distplot(frac_prolif_parameters['cost_006'])

plt.scatter(frac_prolif_parameters['cost_mmd_55']/2,frac_prolif_parameters['cost_006'])

sns.distplot(frac_prolif_parameters['cost_total'])
sns.distplot(frac_prolif_parameters['cost_mmd_55']*3/8)
sns.distplot(frac_prolif_parameters['cost_006'])

np.mean(frac_prolif_parameters['cost_mmd_55']*3/8)
np.mean(frac_prolif_parameters['cost_006'])

##########
#generate th sobol sequence

pysb_sampled_parameter_names = ["V_s2skp","V_Fp","k_M","k_EMY","k_CDF","k_CDS","k_ED","K_S","K_MY","K_DF","d_CD","k_P1", "k_P2", "K_CD", "k_E", "K_M", "K_CE", "V_2E2F", "K_2E2F", "k_RE", "kf1", "kf2", "kb", "kasbaxmbcl2", "kdsbaxmbcl2", "k1sbh3", "k2sbh3", "kdbh3", "kasbh3bcl2", "kdsbh3bcl2"]

loc=[]

for parameter in pysb_sampled_parameter_names:
    print((model.parameters[parameter].value))
    loc+=[np.log10(model.parameters[parameter].value)]



l = sobol_seq.i4_sobol_generate(30, 140000, skip = 0)
l2 = 4*l + (np.array(loc)-2)

######################
#selecting the best parameters
#example: the ten best

frac_prolif_parameters.sort_values(by=['cost_total'])[0:100]

indices = [int(i)-1 for i in frac_prolif_parameters.sort_values(by=['cost_total'])["index"].tolist()[0:100]] #indices in the sobol sequence

parameters = [l2[i] for i in indices]


#########################
#importing parameters from the pso space exploration

particle_positions_pso = pd.read_pickle(r"C:/Users/Grim.M/Ecole/ens/Master/S2/Stage/Code/fit_data/pso/high_and_low_stiffness_cost/compilation/particle_positions_pso.pkl")

frac_prolif_parameters_pso = pd.read_pickle(r"C:/Users/Grim.M/Ecole/ens/Master/S2/Stage/Code/fit_data/pso/high_and_low_stiffness_cost/compilation/frac_prolif_parameters_all_hyperparameters_pso.pkl")

position_history = np.load(r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code\fit_data\pso\high_and_low_stiffness_cost\compilation\position_history.npy", allow_pickle=True).tolist()

#Parameters sampled from the simulation

parameters_log_best=[]
cost_limit = frac_prolif_parameters.sort_values(by=['cost_total'])[0:1]['cost_total'].tolist()[-1] #set the cost limit to the best found in sobol sequence, can be changed

frac_prolif_parameters_pso['cost_total'].tolist()

for idx in frac_prolif_parameters_pso['index']:
    if frac_prolif_parameters_pso.loc[frac_prolif_parameters_pso['index']==idx].iloc[0]['cost_total']<cost_limit: #select the cost  lower than cost limit
        print(idx)
        frac_prolif_parameters_pso_temp = frac_prolif_parameters_pso.loc[frac_prolif_parameters_pso['index']==idx]
        parameter = particle_positions_pso[particle_positions_pso['identification']==idx]['parameter_vector'].iloc[0]
        pos = np.where(position_history[:,:,:] == parameter)
        position_history[pos[0][0],pos[1][0],:]
        parameters_log_best += [position_history[pos[0][0],pos[1][0],:].tolist()]


parameters_log_best = np.array(parameters_log_best) #here is the set of parameters, note it's a log scale
len(parameters_log_best)

parameter_vector = 10**parameters_log_best[1]

param_dict = {pname: pvalue for pname, pvalue in zip(pysb_sampled_parameter_names, parameter_vector)}   
for pname, pvalue in param_dict.items():
    #Change model parameter values to current location in parameter space
    model_population.parameters[pname].value = (pvalue)
    model_concentration.parameters[pname].value = (pvalue)

t = np.linspace(0, 800, 800)

simres_2 = BngSimulator(model_population, t, verbose=True).run();

pl.plot(simres_2.observables['obsCdc20'])


###########################################################################
#Usign clusters parameters

all_fft = np.load(r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code\fit_data\pso\high_and_low_stiffness_cost\hyperparameters_5\fft\stochastic_simulation\all_fft.npy")
indices = np.load(r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code\fit_data\pso\high_and_low_stiffness_cost\hyperparameters_5\fft\stochastic_simulation\indices.npy")

particle_positions_pso = pd.read_pickle(r"C:/Users/Grim.M/Ecole/ens/Master/S2/Stage/Code/fit_data/pso/high_and_low_stiffness_cost/compilation/particle_positions_pso.pkl")

frac_prolif_parameters_pso = pd.read_pickle(r"C:/Users/Grim.M/Ecole/ens/Master/S2/Stage/Code/fit_data/pso/high_and_low_stiffness_cost/compilation/frac_prolif_parameters_all_hyperparameters_pso.pkl")

position_history = np.load(r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code\fit_data\pso\high_and_low_stiffness_cost\compilation\position_history.npy", allow_pickle=True)

hyp=4
#hyperparameters 5


#Import the ffts calculated from stchastic data
from sklearn.cluster import KMeans

n_clusters=3

kmeans = KMeans(n_clusters, random_state=0).fit(all_fft)
kmeans.labels_
kmeans.cluster_centers_

parameters_log_clusters=[[] for c in range(n_clusters)]



for chain in range(np.shape(position_history[4])[1]) :
    for position_ind in range(np.shape(position_history[4])[0]) :    
        try :
            a = np.reshape(np.concatenate(particle_positions_pso[particle_positions_pso['hyperparameters']==hyp+1]['parameter_vector'].to_numpy()),(len(particle_positions_pso[particle_positions_pso['hyperparameters']==hyp+1]),30))
            parameter_vector = particle_positions_pso[particle_positions_pso['hyperparameters']==hyp+1][a==position_history[4][position_ind,chain,:].tolist()].iloc[0]['parameter_vector']
            idx = particle_positions_pso[particle_positions_pso['hyperparameters']==hyp+1][a==position_history[4][position_ind,chain,:].tolist()].iloc[0]['identification']
            

#           if frac_prolif_parameters_pso.loc[frac_prolif_parameters_pso['index']==idx].iloc[0]['cost_total']<cost_limit : #select the cost  lower than cost limit
            
            kmeans.labels_[indices==int(idx)][0]
            parameters_log_clusters[kmeans.labels_[np.array(indices)==int(idx)][0]] += [position_history[hyp][position_ind,chain,:].tolist()]

        except :
            print(idx)
            

            
#len(indices)         
#len(parameters_log_clusters[0])+len(parameters_log_clusters[1])+len(parameters_log_clusters[2])

parameters_log_clusters[0] = np.array(parameters_log_clusters[0])
parameters_log_clusters[1] = np.array(parameters_log_clusters[1])
parameters_log_clusters[2] = np.array(parameters_log_clusters[2])

#Cluster 1 :

parameter_vector = 10**parameters_log_clusters[0][0]

param_dict = {pname: pvalue for pname, pvalue in zip(pysb_sampled_parameter_names, parameter_vector)}   
for pname, pvalue in param_dict.items():
    #Change model parameter values to current location in parameter space
    model_population.parameters[pname].value = (pvalue)
    model_concentration.parameters[pname].value = (pvalue)

t = np.linspace(0, 100, 100)
simres_2 = BngSimulator(model_population, t, verbose=True).run();
pl.plot(simres_2.observables['obsCdc20'])


#Cluster 2 :

parameter_vector = 10**parameters_log_clusters[1][0]

param_dict = {pname: pvalue for pname, pvalue in zip(pysb_sampled_parameter_names, parameter_vector)}   
for pname, pvalue in param_dict.items():
    #Change model parameter values to current location in parameter space
    model_population.parameters[pname].value = (pvalue)
    model_concentration.parameters[pname].value = (pvalue)

t = np.linspace(0, 100, 100)
simres_2 = BngSimulator(model_population, t, verbose=True).run();
pl.plot(simres_2.observables['obsCdc20'])


#Cluster 3 :

parameter_vector = 10**parameters_log_clusters[2][0]

param_dict = {pname: pvalue for pname, pvalue in zip(pysb_sampled_parameter_names, parameter_vector)}   
for pname, pvalue in param_dict.items():
    #Change model parameter values to current location in parameter space
    model_population.parameters[pname].value = (pvalue)
    model_concentration.parameters[pname].value = (pvalue)

t = np.linspace(0, 100, 100)
simres_2 = BngSimulator(model_population, t, verbose=True).run();
pl.plot(simres_2.observables['obsCdc20'])

