# -*- coding: utf-8 -*-
"""
Created on Sat Nov  5 18:06:49 2022

@author: Grim.M
"""

import re

import os
import sys
path= r'C:\Users\Grim.M\Documents\git_repos\ECM_cell_cycle_model'

path= os.getcwd()

sys.path.insert(0, path + "/Python/cellular_cycle") #add a pathway to the model #not sure about how to use folder dependency, to change
#sys.path.insert(0, path + "/Python/test_treatments") #add a pathway to the model #not sure about how to use folder dependency, to change
sys.path.insert(0, path + "/Python/fit") #add a pathway to the fit functions #not sure about how to use folder dependency, to change

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import random
import sobol_seq

##Importing the model, models differ on their initial conditions, one being a concentration and the other a protein count
from Austin_proliferation_model_concentration_5 import model
model_concentration=model
from Austin_proliferation_model_population_5 import model
model_population=model

from cost_function_BngSimulator_frac_prolif_quiescence_and_mmd_ERL_only_30_parameters_fixed_initial_concentrations_mitosis_low_high_stiffness import cost_function_reduced_quiescence_ERL_only_30_parameters

from cost_function_BngSimulator_cell_count_30_parameters_middle_stiffness import cost_function_4kPa_30_parameters

##Import the simulator: for ODE and stochastic simulation :
from pysb.simulator import BngSimulator
from pysb.simulator import ScipyOdeSimulator

# Directory containing cell fate data
from pathlib import Path

pysb_sampled_parameter_names = ["V_s2skp","V_Fp","k_M","k_EMY","k_CDF","k_CDS","k_ED","K_S","K_MY","K_DF","d_CD","k_P1", "k_P2", "K_CD", "k_E", "K_M", "K_CE", "V_2E2F", "K_2E2F", "k_RE", "kf1", "kf2", "kb", "kasbaxmbcl2", "kdsbaxmbcl2", "k1sbh3", "k2sbh3", "kdbh3", "kasbh3bcl2", "kdsbh3bcl2"]
loc=[]
for parameter in pysb_sampled_parameter_names:
    #print((model.parameters[parameter].value))
    loc+=[np.log10(model.parameters[parameter].value)]


print(model_population.parameters["k_M"].value)
print(model.parameters["k_M"].value)
key = "k_M" #0.11 to 0.27
key = "k_EMY" #from 0.18 to 0.37
key = "k_CDS" #from 0.11 to 0.25
print(model_population.parameters[key].value)

#change the model parameters of the model if needed
omega = 200
new_param = pd.read_csv(path + "/umap_window/parameters.csv")
del new_param['Stress']
for key in new_param.keys():
    if model_population.parameters[key].value != new_param[key][0]:
        model_population.parameters[key].value = new_param[key]
    if model_concentration.parameters[key].value != new_param[key][0]:
        model_concentration.parameters[key].value = new_param[key]

model_concentration.parameters['omega'].value = 1


l=sobol_seq.i4_sobol_generate(30, 110000, skip=0)
l2 = 4*l + (np.array(loc)-2)

t = np.linspace(0, 500, 50000)
simres = ScipyOdeSimulator(model_concentration, t, verbose=False).run();
sim = BngSimulator(model_population, tspan=t, verbose=True)
traj = sim.run()



frac_prolif_parameters = pd.read_pickle(path  + "/Python/fit_data/sobol_seq/frac_prolif_parameters_all_parameters_006_55_kPa.pkl")

#indices in the sobol sequence for the 100 best parameters
indices = [int(i)-1 for i in frac_prolif_parameters.sort_values(by=['cost_total'])["index"].tolist()[0:]] 
parameters = [l2[i] for i in indices]


################






run_time=102.3
omega=200
observed_cells_time = 30
sigma=5
parameter_names = ['a','d','q']

print("The parameters are ")
print(parameter_vector)
#np.save(path + "/fit_data/sobol_seq/cost_function_model_4/particle_position/parameter_vector_" + str(idx) + ".npy", parameter_vector)

global model_concentration, model_population

param_dict = {pname: pvalue for pname, pvalue in zip(pysb_sampled_parameter_names, parameter_vector)}   
for pname, pvalue in param_dict.items():
    #Change model parameter values to current location in parameter space
    model_population.parameters[pname].value = (pvalue)
    model_concentration.parameters[pname].value = (pvalue)


"""
###############################################################################
###############################################################################
#Setting 1st condition
parameters_sim={}
thr=0.26
cost_mmd_006_ERL=0
idx = 'r'
cost_006_ERL=0    
condition = '0.6 kPa ERL'
init=1 #40
parameters_sim[condition] = {}

model_population.parameters['ERK'].value = 0.2 * omega
model_population.parameters['FAK'].value = 0.4 * omega
model_population.parameters['AKT'].value = 0.3 * omega

model_concentration.parameters['ERK'].value = 1.0 * 1 ###
model_concentration.parameters['FAK'].value = 0.4 * 1 ###
model_concentration.parameters['AKT'].value = 0.6 * 1 ###


model_population.parameters['ERK'].value = 0.2 * omega
model_population.parameters['FAK'].value = 1.0 * omega
model_population.parameters['AKT'].value = 0.6 * omega

model_concentration.parameters['ERK'].value = 1 * 1
model_concentration.parameters['FAK'].value = 1 * 1
model_concentration.parameters['AKT'].value = 1 * 1
###############################################################################
from runStochCell import runStochcell

variable_range, check, tmp_score = runStochcell(model_population, model_concentration, idx, condition, thr);
 
#try:
from runAKTERK_2 import runAKTERK

(cells, cell_count) = runAKTERK(init, 100, 20, model_concentration, model_population, idx, condition);


"""






################
#original parameters
original_parameters = []
for name in pysb_sampled_parameter_names:
    original_parameters += [np.log10(model_concentration.parameters[name].value)]

cost_function_reduced_quiescence_ERL_only_30_parameters(10**np.array(original_parameters), 'original_parameters')

cost_function_4kPa_30_parameters(10**np.array(original_parameters), 'original_parameters')

parameter_vector = 10**np.array(original_parameters)
idx='original_parameters'

from runAKTERK_2 import runAKTERK


parameters_original_ERL = pd.read_pickle(r"C:\Users\Grim.M\Ecole\Ens\Master\S2\Stage\Code\fit_data\sobol_seq\version 6\cost\parameters_original_parameters_high_and_low_stiffness_ERL.pkl")

parameters_original_4kpa = pd.read_pickle(r"C:\Users\Grim.M\Ecole\Ens\Master\S2\Stage\Code\fit_data\sobol_seq\version 6\cost\parameters_original_parameters_4kPa_DMSO_ERL.pkl")

################
#Change model parameter values to current location in parameter space
#parameter_vector = 10**parameters[0]
#idx = 'test'



###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
#Low stiff simulation

cost=[]
for parameter_vector, idx in zip(parameters[:1456], indices[:1456]):
    cost+=[cost_function_reduced_quiescence_ERL_only_30_parameters(10**parameter_vector, idx)]



i=1456
cost = cost_function_reduced_quiescence_ERL_only_30_parameters(10**parameters[i], indices[i])

i=1456
parameter_vector = 10**parameters[i]
idx = indices[i]

###############################################################################
from os import listdir
from os.path import isfile, join
#mypath = r"C:\Users\Grim.M\Ecole\Ens\Master\S2\Stage\Code\fit_data\sobol_seq\version 6\cost"
mypath = path + "/Python/fit_data/sobol_seq/cost"

onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f)) if 'high_and_low_stiffness_ERL' in f]



frac_prolif_parameters_2 = pd.DataFrame(columns = ["index", "a_006_ERL", "d_006_ERL", "q_006_ERL", "cost_mmd_006_ERL", "a_55_ERL", "d_55_ERL", "q_55_ERL", "cost_mmd_55_ERL", "cost_total_ERL"])


for file_temp in onlyfiles:
    frac_prolif_parameters_tmp = pd.read_pickle(mypath + "/"  + file_temp)
    frac_prolif_parameters_2 = pd.concat([frac_prolif_parameters_2, frac_prolif_parameters_tmp])

frac_prolif_parameters_2 = frac_prolif_parameters_2[frac_prolif_parameters_2['cost_total_ERL']!=np.inf]



frac_prolif_parameters_2.to_csv(r'C:\Users\Grim.M\Documents\git_repos\ECM_cell_cycle_model\Python\fit_data\sobol_seq\all_parameters_high_and_low_stiffness_ERL.csv', index=False)

frac_prolif_parameters_2['index']

onlyfiles = [int(f.split("_")[1]) for f in listdir(mypath) if isfile(join(mypath, f)) if '4kPa' in f]

for index_tmp in frac_prolif_parameters_2['index']:
    if index_tmp not in onlyfiles:
        try:
            parameter_vector = parameters[np.where(np.array(indices) == index_tmp)[0][0]]
            cost_function_4kPa_30_parameters(10**parameter_vector, index_tmp)
        
        except:
            print(index_tmp)


for file_temp in onlyfiles:
    try:
        frac_prolif_parameters_tmp = pd.read_pickle(mypath + "/"  + file_temp)
        index_tmp = int(file_temp.split("_")[1])
        parameter_vector = parameters[np.where(np.array(indices) == index_tmp)[0][0]]
        cost_function_4kPa_30_parameters(10**parameter_vector, index_tmp)
    
    except: print(file_temp)
    



###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
#data processing

DATADIR = path + "/Data/Exp_Cell_Counts"


filepaths = os.listdir(DATADIR)
filenames = [f for f in filepaths if re.match(r'^DIP_DS8', f)]
filepaths = [Path(DATADIR, f) for f in filenames]
filepaths.sort()
filepaths

all_response = []
for f in filepaths:
    all_response.append(np.loadtxt(f,delimiter=','))

means=[]
for vals in range(len(all_response)):
    means.append(np.mean(all_response[vals],axis=1))
    
stdevs=[]
for vals in range(len(all_response)):
    stdevs.append(np.std(all_response[vals],axis=1))
  
t = np.arange(0,420/3,1/3)

#all_response.keys()


slope_4kPa_DMSO, _, _, _ = np.linalg.lstsq(t[:,np.newaxis], means[2])

plt.plot(t, means[2], 'bo')
plt.plot(t, slope_4kPa_DMSO*t, 'r-')
plt.show()

slope_4kPa_ERL, _, _, _ = np.linalg.lstsq(t[:,np.newaxis], means[3])

plt.plot(t, means[3], 'bo')
plt.plot(t, slope_4kPa_ERL*t, 'r-')
plt.show()


##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################


###############################################################################

from os import listdir
from os.path import isfile, join
mypath_cc = r"C:\Users\Grim.M\Ecole\Ens\Master\S2\Stage\Code\fit_data\sobol_seq\version 6\cell_count"
mypath_p = r"C:\Users\Grim.M\Ecole\Ens\Master\S2\Stage\Code\fit_data\sobol_seq\version 6\cost"

onlyfiles = [f for f in listdir(mypath_cc) if isfile(join(mypath_cc, f))]

all_parameters = []


file_temp = 'cell_count_3955_4kPa_DMSO.npy'
file_temp = 'cell_count_41712_4kPa_DMSO.npy'
file_temp = 'cell_count_9606_4kPa_DMSO.npy'

file_temp = 'cell_count_637_4kPa_DMSO.npy'
file_temp = 'cell_count_33711_4kPa_DMSO.npy'

index = 637
index = 37007

index = 637
index = 33711

index=68336

for index in frac_prolif_parameters_2['index']:
    try:
        cost_ERL = cost_DMSO = cost_slope_ERL = cost_slope_DMSO = np.inf
        #file_name = file_temp.split('_')
        parameters_tmp = pd.read_pickle(mypath_p + '/parameters_' + str(index) + '_4kPa_DMSO_ERL.pkl')
        
        try :
            parameters_tmp = pd.read_pickle(mypath_p + '/parameters_' + str(index) + '_4kPa_ERL.pkl')
            
            parameters_tmp[['a_4_DMSO', 'q_4_DMSO', 'd_4_DMSO']] = pd.read_pickle(mypath_p + '/parameters_' + str(index) + '_4kPa_DMSO_ERL.pkl')[['a_4_DMSO', 'q_4_DMSO', 'd_4_DMSO']]
            
            #print(file_temp)
            #print('corr')
        except: pass
        
        
        if np.isnan(parameters_tmp['a_4_DMSO']).values[0]:
            print('no oscillations DMSO')
        else:
            cell_count_tmp = np.load(mypath_cc + '/' + '/cell_count_' + str(index) + '_4kPa_DMSO.npy')
            cell_count_tmp[:, 0]
            cell_count_t = np.zeros(len(t))
            time_tmp_m = 0
            
            for time_tmp_p, cell_count_time_tmp in cell_count_tmp:
                cell_count_t[(time_tmp_m <= t) & (t < time_tmp_p)] = cell_count_time_tmp
                time_tmp_m = time_tmp_p
            


    #        plt.plot(t, cell_count_t)
    #        plt.plot(cell_count_tmp[:, 0], cell_count_tmp[:, 1])
    #        plt.show()
            DIP_sim = np.log2(cell_count_t/cell_count_t[0])
            cost_DMSO = np.sqrt(np.sum((DIP_sim - means[2])[1:]**2 / stdevs[2][1:]**2))/len(t)
            parameters_tmp['slope_4kPa_DMSO'], _, _, _ = np.linalg.lstsq(t[:,np.newaxis], DIP_sim)
            cost_slope_DMSO = (parameters_tmp['slope_4kPa_DMSO'] - slope_4kPa_DMSO)[0]**2

            plt.scatter(t, DIP_sim, color = 'blue', s=1)
            plt.plot(t, parameters_tmp['slope_4kPa_DMSO'][0]*t, color = 'blue')
            plt.scatter(t, means[2], color = 'red', s=1)
            plt.plot(t, slope_4kPa_DMSO*t, color = 'red')
            plt.show()

        
        
        #file_name[-1] = 'ERL.npy'
        #cell_count_tmp = np.load(mypath_cc + '/' + '_'.join(file_name))
        
        if np.isnan(parameters_tmp['a_4_ERL']).values[0]:
            print('no oscillations ERL')

        else:
            cell_count_tmp = np.load(mypath_cc + '/' + '/cell_count_' + str(index) + '_4kPa_ERL.npy')
            t_ERL = t[t<cell_count_tmp[-1,0]]
            cell_count_tmp[:, 0]
            cell_count_t = np.zeros(len(t))
            time_tmp_m = 0
            
            for time_tmp_p, cell_count_time_tmp in cell_count_tmp:
                cell_count_t[(time_tmp_m <= t) & (t < time_tmp_p)] = cell_count_time_tmp
                time_tmp_m = time_tmp_p
            
    #        plt.plot(t_ERL, cell_count_t)
    #        plt.plot(cell_count_tmp[:, 0], cell_count_tmp[:, 1])
    #        plt.show()
            
            DIP_sim = np.log(cell_count_t/cell_count_t[0])
            cost_ERL = np.sqrt(np.sum((DIP_sim - means[3][:])[1:len(t)]**2 / stdevs[3][1:len(t)]**2))/len(t)
            parameters_tmp['slope_4kPa_ERL'], _, _, _ = np.linalg.lstsq(t[:,np.newaxis], DIP_sim)

            cost_slope_ERL = (parameters_tmp['slope_4kPa_ERL'] - slope_4kPa_ERL)[0]**2


            plt.scatter(t, DIP_sim, color = 'blue', s=1)
            plt.plot(t, parameters_tmp['slope_4kPa_ERL'][0]*t, color = 'blue')
            plt.scatter(t, means[3], color = 'red', s=1)
            plt.plot(t, slope_4kPa_ERL*t, color = 'red')
            plt.show()
    
        parameters_tmp['cost 4kPa DMSO'] = cost_DMSO
        parameters_tmp['cost 4kPa ERL'] = cost_ERL
        parameters_tmp['cost 4kPa'] = cost_ERL + cost_DMSO
        
        parameters_tmp['cost slope difference 4kPa'] = cost_slope_ERL + cost_slope_DMSO
        parameters_tmp['cost slope difference 4kPa ERL'] = cost_slope_ERL
        parameters_tmp['cost slope difference 4kPa DMSO'] = cost_slope_DMSO
        
        all_parameters += [parameters_tmp]
    except:
        print(index)
        print('error')

all_parameters_4kPa = pd.concat(all_parameters)


plt.scatter(all_parameters_4kPa['cost 4kPa ERL'], all_parameters_4kPa['cost 4kPa DMSO'])
plt.scatter(all_parameters_4kPa['cost 4kPa ERL'][all_parameters_4kPa.index=='original_parameters'], all_parameters_4kPa['cost 4kPa DMSO'][all_parameters_4kPa.index=='original_parameters'], c='red')
plt.show()

plt.scatter(all_parameters_4kPa['cost slope difference 4kPa ERL'], all_parameters_4kPa['cost slope difference 4kPa DMSO'])
plt.scatter(all_parameters_4kPa[all_parameters_4kPa.index=='original_parameters']['cost slope difference 4kPa ERL'], all_parameters_4kPa[all_parameters_4kPa.index=='original_parameters']['cost slope difference 4kPa DMSO'], c='red')
plt.show()

plt.scatter(all_parameters_4kPa['cost slope difference 4kPa'], all_parameters_4kPa['cost 4kPa'])
plt.scatter(all_parameters_4kPa[all_parameters_4kPa.index=='original_parameters']['cost slope difference 4kPa'], all_parameters_4kPa[all_parameters_4kPa.index=='original_parameters']['cost 4kPa'], c='red')
plt.show()

all_parameters_4kPa.to_csv(r'C:\Users\Grim.M\Documents\git_repos\ECM_cell_cycle_model\Python\fit_data\Sobol_seq\all_parameters_4kPa.csv', index=False)

###############################################################################





