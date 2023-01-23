# -*- coding: utf-8 -*-
"""
Created on Mon Dec 19 11:59:15 2022

@author: Grim.M
"""

import time
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

from matplotlib.axes._axes import _log as matplotlib_axes_logger
matplotlib_axes_logger.setLevel('ERROR')
 
import matplotlib.colors as mcolors
import matplotlib.cm as cm

##Importing the model, models differ on their initial conditions, one being a concentration and the other a protein count
from Austin_proliferation_model_concentration_5 import model
model_concentration=model
from Austin_proliferation_model_population_5 import model
model_population=model

#from cost_function_BngSimulator_low_medium_stiffness_pso import cost_function_30_parameters_pso
#from cost_function_BngSimulator_low_medium_stiffness_pso_2 import cost_function_30_parameters_pso_2
from cost_function_BngSimulator_low_medium_stiffness_pso_longer_simulation import cost_function_30_parameters_pso

##Import the simulator: for ODE and stochastic simulation :
from pysb.simulator import BngSimulator
from pysb.simulator import ScipyOdeSimulator

# Directory containing cell fate data
from pathlib import Path
pysb_sampled_parameter_names = ["V_s2skp","V_Fp","k_M","k_EMY","k_CDF","k_CDS","k_ED","K_S","K_MY","K_DF","d_CD","k_P1", "k_P2", "K_CD", "k_E", "K_M", "K_CE", "V_2E2F", "K_2E2F", "k_RE", "kf1", "kf2", "kb", "kasbaxmbcl2", "kdsbaxmbcl2", "k1sbh3", "k2sbh3", "kdbh3", "kasbh3bcl2", "kdsbh3bcl2"]



print(model_population.parameters["k_M"].value)
print(model_population.parameters["k_EMY"].value)
print(model_population.parameters["k_CDS"].value)
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
model_population.parameters['omega'].value

loc=[]
for parameter in pysb_sampled_parameter_names:
    #print((model.parameters[parameter].value))
    loc+=[np.log10(model.parameters[parameter].value)]


t = np.linspace(0, 500, 50000)
#simres = ScipyOdeSimulator(model_concentration, t, verbose=False).run();
#sim = BngSimulator(model_population, tspan=t, verbose=True).run()


original_parameters = []
for name in pysb_sampled_parameter_names:
    original_parameters += [np.log10(model_concentration.parameters[name].value)]



parameter_vector = 10**np.array(original_parameters)
idx = 'original_parameters'

#
#one = time.time()
#cost_function_30_parameters_pso_2(parameter_vector, 'original_parameters')
#two = time.time()
#
#three = time.time()
cost_function_30_parameters_pso(parameter_vector, 'original_parameters')
#four = time.time()
#
#print(four - three)
#print(two - one)



#solved trouble shooting
"""

param_dict = {pname: pvalue for pname, pvalue in zip(pysb_sampled_parameter_names, parameter_vector)}   
for pname, pvalue in param_dict.items():
    #Change model parameter values to current location in parameter space
    model_population.parameters[pname].value = (pvalue)
    model_concentration.parameters[pname].value = (pvalue)
    
    print(model_population.parameters[pname].value )
    print(model_concentration.parameters[pname].value )
    print(pvalue)



conditions = ['0.6 kPa DMSO', '0.6 kPa ERL', '4 kPa DMSO', '4 kPa ERL']
inputs = {}
inputs['FAK'] = {'4 kPa DMSO': 0.6, '4 kPa ERL': 0.6, '0.6 kPa DMSO': 0.4, '0.6 kPa ERL': 0.4}
inputs['ERK'] = {'4 kPa DMSO': 1.0, '4 kPa ERL': 0.2, '0.6 kPa DMSO': 1.0, '0.6 kPa ERL': 0.2}
inputs['AKT'] = {'4 kPa DMSO': 1.0, '4 kPa ERL': 0.6, '0.6 kPa DMSO': 0.6, '0.6 kPa ERL': 0.3}

condition = '0.6 kPa ERL'
num_condition=0
parameters_sim={} 

t_sim = np.arange(0,120,1/3)

run_time = observed_cells_time = t_sim[-1]
omega=200
sigma=5
parameter_names = ['a','d','q']
init=1 #40
idx= 'test_p'

cond_split = condition.split(' ')
cond_split_DMSO = condition.split(' ')
cond_split_DMSO[-1] = 'DMSO'
parameters_sim[condition] = {}


model_population.parameters['ERK'].value = inputs['ERK'][condition] * omega
model_population.parameters['FAK'].value = inputs['FAK'][condition] * omega
model_population.parameters['AKT'].value = inputs['AKT'][condition] * omega

model_concentration.parameters['ERK'].value = inputs['ERK'][' '.join(cond_split_DMSO)] * 1
model_concentration.parameters['FAK'].value = inputs['FAK'][' '.join(cond_split_DMSO)] * 1
model_concentration.parameters['AKT'].value = inputs['AKT'][' '.join(cond_split_DMSO)] * 1


cost = {}
from runAKTERK_2 import runAKTERK
(cells, cell_count) = runAKTERK(init, run_time, observed_cells_time, model_concentration, model_population, idx, condition);


thr=0.26
sys_size=200
t = np.linspace(0, 1000, 1000)

from peakfinder import peakfinder


simres_2 = BngSimulator(model_population, t, verbose=False).run();

(a,max_tmp_cdc) = peakfinder(simres_2.observables['obsCdc20'], thr*sys_size, thr*sys_size, 1, includeEndpoints=1);
(b,max_tmp_cyce) = peakfinder(simres_2.observables['obsCycE'], thr*sys_size, thr*sys_size, 1, includeEndpoints=1);

max_tmp_cdc
"""

import pyswarms as ps


all_parameters_4kPa = pd.read_csv(path + '/Python/fit_data/sobol_seq/all_parameters_4kPa.csv', index_col='index')

all_parameters_4kPa = all_parameters_4kPa.sort_values(by=['cost 4kPa'])

#parameter_vector = np.log10(all_parameters_4kPa[all_parameters_4kPa.index == all_parameters_4kPa.index.iloc[9]][pysb_sampled_parameter_names].to_numpy())[0]


initial_positions = np.log10(all_parameters_4kPa.iloc[0:30][pysb_sampled_parameter_names].to_numpy())


all_parameters_4kPa.iloc[8] #original

def cost_function_f(parameter_vectors):
    cost = np.zeros(np.shape(parameter_vectors)[0])
    for i in range(len(cost)) :
        parameter_vector = parameter_vectors[i]
        print(parameter_vector)
        parameter_vector=[10**i for i in parameter_vector]
        cost[i] = + cost_function_30_parameters_pso(parameter_vector, i)
    return cost

#
#from multiprocessing.pool import ThreadPool as Pool
#from multiprocessing import freeze_support
#
#
#def cost_function_f_2(parameter_vectors):        
#    print('letsgo')
#    with Pool() as pool:
#        freeze_support()
#        cost = pool.starmap(cost_function_30_parameters_pso, zip(parameter_vectors, range(len(parameter_vectors))))   
#    return np.array(cost)
#

max_bound = np.array(loc)+2.5
min_bound = np.array(loc)-2.5
bounds = (min_bound, max_bound)

#cost_function_f(initial_positions)

















#
#
#parameter_vector=[10**i for i in initial_positions[2]]
#
#param_dict = {pname: pvalue for pname, pvalue in zip(pysb_sampled_parameter_names, parameter_vector)}
#for pname, pvalue in param_dict.items():
#    #Change model parameter values to current location in parameter space
#    model_population.parameters[pname].value = (pvalue)
#    model_concentration.parameters[pname].value = (pvalue)
#    
#
#conditions = ['0.6 kPa DMSO', '0.6 kPa ERL', '4 kPa DMSO', '4 kPa ERL', '55 kPa DMSO', '55 kPa ERL']
#inputs = {}
#inputs['FAK'] = {'4 kPa DMSO': 0.6, '4 kPa ERL': 0.6, '0.6 kPa DMSO': 0.4, '0.6 kPa ERL': 0.4, '55 kPa DMSO': 1.0, '55 kPa ERL': 1.0}
#inputs['ERK'] = {'4 kPa DMSO': 1.0, '4 kPa ERL': 0.2, '0.6 kPa DMSO': 1.0, '0.6 kPa ERL': 0.2, '55 kPa DMSO': 1.0, '55 kPa ERL': 0.2}
#inputs['AKT'] = {'4 kPa DMSO': 1.0, '4 kPa ERL': 0.6, '0.6 kPa DMSO': 0.6, '0.6 kPa ERL': 0.3, '55 kPa DMSO': 1.0, '55 kPa ERL': 0.6}
#
#
#parameters_sim={} 
#cost = {}
#
#t_sim = np.arange(0,120,1/3)
#
#run_times = [120, 120, 120, 120, 60, 60] 
#omega = 200
#sigma = 5
#parameter_names = ['a','d','q']
#init = 1 #40
#
#condition = '55 kPa DMSO'
#condition = '0.6 kPa ERL'
#run_time = 120
#condition = '4 kPa ERL'
#
#run_time = 60
#condition = '55 kPa ERL'
################################################################################
################################################################################
#
#
#t_sim = np.arange(0,run_time,1/3)
#observed_cells_time = run_time
#cond_split = condition.split(' ')
#cond_split_DMSO = condition.split(' ')
#cond_split_DMSO[-1] = 'DMSO'
#parameters_sim[condition] = {}
#
#model_population.parameters['ERK'].value = inputs['ERK'][condition] * omega
#model_population.parameters['FAK'].value = inputs['FAK'][condition] * omega
#model_population.parameters['AKT'].value = inputs['AKT'][condition] * omega
#
#model_concentration.parameters['ERK'].value = inputs['ERK'][' '.join(cond_split_DMSO)] * 1
#model_concentration.parameters['FAK'].value = inputs['FAK'][' '.join(cond_split_DMSO)] * 1
#model_concentration.parameters['AKT'].value = inputs['AKT'][' '.join(cond_split_DMSO)] * 1
#
################################################################################
#    
##try:
#
#(cells, cell_count) = runAKTERK(init, run_time, run_time, model_concentration, model_population, idx, condition);
#
#print('Cell count')
#print(cell_count)
#








parameter_vector=[10**i for i in initial_positions[2]]
cost_function_30_parameters_pso(parameter_vector, 2)

initial_positions[0:30]<max_bound
initial_positions[0:30]>min_bound

#options = {'c1': 0.5, 'c2': 0.3, 'w': 0.9, 'k': 3, 'p': 2} #hyperparameters 1
#options = {'c1': 0.05, 'c2': 0.1, 'w': 0.1, 'k': 6, 'p': 2} #hyperparameters 2
#options = {'c1': 0.1, 'c2': 0.3, 'w': 0.3, 'k': 5, 'p': 2} #hyperparameters 3
#options = {'c1': 0.2, 'c2': 0.3, 'w': 0.4, 'k': 5, 'p': 2} #hyperparameters 4
options = {'c1': 0.15, 'c2': 0.3, 'w': 0.4, 'k': 5, 'p': 2} #hyperparameters 5
#options = {'c1': 0.16, 'c2': 0.32, 'w': 0.45, 'k': 6, 'p': 2} #hyperparameters 6

#position = np.load("/home/jeremie/fit_data/pso/high_and_low_stiffness_cost/history/position_history.npy")


# Call instance of PSO
optimizer = ps.single.LocalBestPSO(n_particles=30, dimensions=30, options=options, bounds=bounds, init_pos = initial_positions[0:30])
#optimizer = ps.single.LocalBestPSO(n_particles=30, dimensions=30, options=options, bounds=bounds, init_pos = np.array(position[-1,:,:]))

# Provide your own init_pos
#optimizer.pos = np.array(initial_positions[0:10]).copy()

# Perform optimization
cost, pos = optimizer.optimize(cost_function_f, iters=60)



np.save(path + "/Python/fit_data/pso/history/position_history.npy", optimizer.pos_history)

np.save(path + "/Python/fit_data/pso/history/cost_history.npy", optimizer.cost_history)



###############################################################################
###############################################################################
###############################################################################
#Longer simulations

frac_prolif_parameters_pso = pd.read_pickle(path + "/Python/fit_data/pso/all_parameters_pso.pkl")
frac_prolif_parameters_original_parameters = pd.read_pickle(path + "/Python/fit_data/pso/parameters_original_parameters_pso.pkl")



parameters_longer_sim = frac_prolif_parameters_pso[frac_prolif_parameters_pso['Cost']<   frac_prolif_parameters_original_parameters['Cost'][0]]

parameters_longer_sim = parameters_longer_sim.sort_values(by=['Cost'])


parameter_vector = parameters_longer_sim[pysb_sampled_parameter_names].to_numpy()[0]
idx = parameters_longer_sim.index[0]
cost_function_30_parameters_pso(parameter_vector, '0.' + idx)

for parameter_vector, idx in zip(parameters_longer_sim[pysb_sampled_parameter_names].to_numpy()[1:], parameters_longer_sim.index[1:]):
    print(parameter_vector)
    print('0.' + idx)
    cost_function_30_parameters_pso(parameter_vector, '0.' + idx)




frac_prolif_parameters_pso = pd.read_pickle(path + "/Python/fit_data/pso/all_parameters_pso.pkl")
frac_prolif_parameters_original_parameters = pd.read_pickle(path + "/Python/fit_data/pso/parameters_original_parameters_pso.pkl")

parameters_longer_sim = frac_prolif_parameters_pso[frac_prolif_parameters_pso['Cost']<   np.inf]
parameters_longer_sim = parameters_longer_sim.sort_values(by=['Cost'], ascending=False)


for parameter_vector, idx in zip(parameters_longer_sim[pysb_sampled_parameter_names].to_numpy()[0:60], parameters_longer_sim.index[0:60]):
    print(parameter_vector)
    print('0.' + idx)
    cost_function_30_parameters_pso(parameter_vector, '0.' + idx)



###############################################################################
###############################################################################
#Data procesing


path_cost = r"C:\Users\Grim.M\Ecole\Ens\Master\S2\Stage\Code\fit_data\pso/cost_3_bis_bis_bis/longer_sim/cost"
frames=[]
import os
os.chdir(path_cost)
for filename in os.listdir(os.getcwd()):
    if 'parameters_' in filename:
        print(filename)
        frac_prolif_parameters_temp = pd.read_pickle(path_cost + "\\" + filename)
        frac_prolif_parameters_temp['index'] = re.findall('\d+', filename)[-1] 
#        frac_prolif_parameters_temp['particle'] = int(re.search('parameters_(.*)_0.', filename).group(1))
        
#        frac_prolif_parameters_temp['step'] = np.where(pos_history[:, int(re.search('parameters_(.*)_0.', filename).group(1)), 0] == np.log10(frac_prolif_parameters_temp[pysb_sampled_parameter_names].to_numpy()[0][0]))[0][0]
        frames += [frac_prolif_parameters_temp]

frac_prolif_parameters_pso_longer_sim = pd.concat(frames, sort=False)
frac_prolif_parameters_pso_longer_sim.set_index('index', inplace=True)
frac_prolif_parameters_pso_longer_sim['index'] = frac_prolif_parameters_pso_longer_sim.index

frac_prolif_parameters_pso_longer_sim['Cost'] = frac_prolif_parameters_pso_longer_sim['Cost 0.6 kPa DMSO'] + frac_prolif_parameters_pso_longer_sim['Cost 0.6 kPa ERL'] + frac_prolif_parameters_pso_longer_sim['Cost 4 kPa DMSO'] + frac_prolif_parameters_pso_longer_sim['Cost 4 kPa ERL'] + frac_prolif_parameters_pso_longer_sim['Cost 55 kPa DMSO'] + frac_prolif_parameters_pso_longer_sim['Cost 55 kPa ERL']

###############################################################################
###############################################################################
###############################################################################
#Data processing

pos_history = np.load(r"C:\Users\Grim.M\Ecole\Ens\Master\S2\Stage\Code\fit_data\pso\cost_3_bis_bis_bis\history\position_history.npy")
cost_history = np.load(r"C:\Users\Grim.M\Ecole\Ens\Master\S2\Stage\Code\fit_data\pso\cost_3_bis_bis_bis\history\cost_history.npy")



path_cost = r"C:\Users\Grim.M\Ecole\Ens\Master\S2\Stage\Code\fit_data\pso\cost_3_bis_bis_bis\cost"
frames=[]
import os
os.chdir(path_cost)
for filename in os.listdir(os.getcwd()):
    if 'parameters_' in filename:
        #print(filename)
        frac_prolif_parameters_temp = pd.read_pickle(path_cost + "\\" + filename)
        frac_prolif_parameters_temp['index'] = re.findall('\d+', filename)[-1] 
        frac_prolif_parameters_temp['particle'] = int(re.search('parameters_(.*)_0.', filename).group(1))
        
        frac_prolif_parameters_temp['step'] = np.where(pos_history[:, int(re.search('parameters_(.*)_0.', filename).group(1)), 0] == np.log10(frac_prolif_parameters_temp[pysb_sampled_parameter_names].to_numpy()[0][0]))[0][0]
        frames += [frac_prolif_parameters_temp]


frac_prolif_parameters_pso = pd.concat(frames, sort=False)
frac_prolif_parameters_pso.set_index('index', inplace=True)
frac_prolif_parameters_pso['index'] = frac_prolif_parameters_pso.index

frac_prolif_parameters_pso.to_pickle(r"")

frac_prolif_parameters_pso['Cost'] = frac_prolif_parameters_pso['Cost 0.6 kPa DMSO'] + frac_prolif_parameters_pso['Cost 0.6 kPa ERL'] + frac_prolif_parameters_pso['Cost 4 kPa DMSO'] + frac_prolif_parameters_pso['Cost 4 kPa ERL'] + frac_prolif_parameters_pso['Cost 55 kPa DMSO'] + frac_prolif_parameters_pso['Cost 55 kPa ERL']



#original parameter
frac_prolif_parameters_original_parameters = pd.read_pickle(r'C:\Users\Grim.M\Ecole\Ens\Master\S2\Stage\Code\fit_data\pso\cost_3_bis_bis_bis\original_parameters\parameters_original_parameters_0.43973113833887045_pso.pkl')
frac_prolif_parameters_original_parameters['Cost'] = frac_prolif_parameters_original_parameters['Cost 0.6 kPa DMSO'] + frac_prolif_parameters_original_parameters['Cost 0.6 kPa ERL'] + frac_prolif_parameters_original_parameters['Cost 4 kPa DMSO'] + frac_prolif_parameters_original_parameters['Cost 4 kPa ERL'] + frac_prolif_parameters_original_parameters['Cost 55 kPa DMSO'] + frac_prolif_parameters_original_parameters['Cost 55 kPa ERL']


#dim reduction
from sklearn import decomposition
import umap as umap

pca_1 = decomposition.PCA(n_components=2)
pca_1.fit(np.log10(frac_prolif_parameters_pso[pysb_sampled_parameter_names].to_numpy()))
pca_tmp = pca_1.transform(np.log10(frac_prolif_parameters_pso[pysb_sampled_parameter_names].to_numpy()))
pca_original_param = pca_1.transform(np.log10(frac_prolif_parameters_original_parameters[pysb_sampled_parameter_names].to_numpy()))

umap_1 = umap.UMAP().fit(np.log10(frac_prolif_parameters_pso[pysb_sampled_parameter_names].to_numpy()))

umap_tmp = umap_1.transform(np.log10(frac_prolif_parameters_pso[pysb_sampled_parameter_names].to_numpy()))
umap_original_param = umap_1.transform(np.log10(frac_prolif_parameters_original_parameters[pysb_sampled_parameter_names].to_numpy()))

frac_prolif_parameters_pso['pca 1'] = pca_tmp[:,0]
frac_prolif_parameters_pso['pca 2'] = pca_tmp[:,1]
frac_prolif_parameters_pso['umap 1'] = umap_tmp[:,0]
frac_prolif_parameters_pso['umap 2'] = umap_tmp[:,1]


frac_prolif_parameters_original_parameters['pca 1'] = pca_original_param[:,0]
frac_prolif_parameters_original_parameters['pca 2'] = pca_original_param[:,1]
frac_prolif_parameters_original_parameters['umap 1'] = umap_original_param[:,0]
frac_prolif_parameters_original_parameters['umap 2'] = umap_original_param[:,1]



frac_prolif_parameters_pso = frac_prolif_parameters_pso.sort_values(by=['step'])


###############################################################################
###############################################################################
###############################################################################
#Saving

frac_prolif_parameters_pso.to_pickle(r"C:\Users\Grim.M\Documents\git_repos\ECM_cell_cycle_model\Python\fit_data\pso\all_parameters_pso.pkl") 

frac_prolif_parameters_original_parameters.to_pickle(r"C:\Users\Grim.M\Documents\git_repos\ECM_cell_cycle_model\Python\fit_data\pso\parameters_original_parameters_pso.pkl")

frac_prolif_parameters_pso_longer_sim.to_pickle(r"C:\Users\Grim.M\Documents\git_repos\ECM_cell_cycle_model\Python\fit_data\pso\all_parameters_pso_long_sim.pkl") 


###############################################################################
###############################################################################
###############################################################################
#Loading


frac_prolif_parameters_pso = pd.read_pickle(path + "/Python/fit_data/pso/all_parameters_pso.pkl")
frac_prolif_parameters_original_parameters = pd.read_pickle(path + "/Python/fit_data/pso/parameters_original_parameters_pso.pkl")
pos_history = np.load(r"C:\Users\Grim.M\Ecole\Ens\Master\S2\Stage\Code\fit_data\pso\cost_3_bis_bis_bis\history\position_history.npy")


frac_prolif_parameters_pso_longer_sim = pd.read_pickle(r"C:\Users\Grim.M\Documents\git_repos\ECM_cell_cycle_model\Python\fit_data\pso\all_parameters_pso_long_sim.pkl") 

###############################################################################
###############################################################################
###############################################################################
#fig 1
dim_reduction = 'umap'
dim_reduction = 'pca'

colors = sns.color_palette(n_colors=np.shape(pos_history)[1])

fig = plt.figure(figsize=(3,3),dpi=200)
ax = fig.add_subplot(111)
for chain in range(np.shape(pos_history)[1]) :
    x_tmp = frac_prolif_parameters_pso[frac_prolif_parameters_pso['particle']==chain][[dim_reduction + ' 1', dim_reduction + ' 2']].to_numpy()
    ax.scatter(x_tmp[:, 0], x_tmp[:, 1], c = colors[chain], alpha = 0.25,s = 0.1);
    ax.plot(x_tmp[:,0], x_tmp[:,1], c = colors[chain], alpha = 0.5, lw = 0.25, ls = '--');
    ax.plot(x_tmp[0,0], x_tmp[0,1],c=colors[chain], marker='+')
    ax.plot(x_tmp[-1,0], x_tmp[-1,1],c=colors[chain], marker='o')

x_tmp = frac_prolif_parameters_original_parameters[[dim_reduction + ' 1', dim_reduction + ' 2']].to_numpy()
ax.scatter(x_tmp[0,0], x_tmp[0,1],c='black', marker="*", s=400)

plt.savefig(path + "/../report/figures/pso/all_particles_running_pca_version_3_hyperparameters_1.pdf")


############################################################
#fig 2

colormap = cm.rainbow
normalize = mcolors.Normalize(vmin = min(frac_prolif_parameters_pso["Cost"]), vmax = max(frac_prolif_parameters_pso["Cost"][frac_prolif_parameters_pso["Cost"]<np.inf]))

frac_prolif_parameters_pso['color'] = list(map(tuple, colormap(normalize(frac_prolif_parameters_pso['Cost']))))
frac_prolif_parameters_pso['color'][frac_prolif_parameters_pso['Cost'] == np.inf] = 'grey'


fig = plt.figure(figsize=(3,3),dpi=200)
ax = fig.add_subplot(111)

for chain in range(np.shape(pos_history)[1]) :
    print(chain)
    x_tmp = frac_prolif_parameters_pso[frac_prolif_parameters_pso['particle']==chain][[dim_reduction + ' 1', dim_reduction + ' 2']].to_numpy()
    color_tmp = frac_prolif_parameters_pso[frac_prolif_parameters_pso['particle']==chain]['color']

    ax.plot(x_tmp[:,0], x_tmp[:,1], c =  'grey', alpha = 0.5, lw = 0.25, ls = '--');
    ax.scatter(x_tmp[:, 0], x_tmp[:, 1], c = color_tmp, alpha = 0.9,s = 0.3);
    ax.scatter(x_tmp[0,0], x_tmp[0,1], c=color_tmp[0], marker='+', s=30)
    ax.scatter(x_tmp[-1,0], x_tmp[-1,1], c=color_tmp[-1], marker='o', s=30)

scalarmappaple = cm.ScalarMappable(norm=normalize, cmap=colormap)
scalarmappaple.set_array([])

cb = fig.colorbar(scalarmappaple, shrink=0.8)
cb.set_label(label='Cost', fontsize=10)
cb.ax.tick_params(labelsize=6)

x_tmp = frac_prolif_parameters_original_parameters[[dim_reduction + ' 1', dim_reduction + ' 2']].to_numpy()
ax.scatter(x_tmp[0,0], x_tmp[0,1],c='black', marker="*", s=400)


pl.savefig(path + "/../report/figures/pso/all_particles_running_pca_hyperparameters_cost.pdf", bbox_inches='tight')



############################################################
#fig 3

df = frac_prolif_parameters_pso[frac_prolif_parameters_pso['Cost']<np.inf][['index', 'Cost', 'particle', 'step']]
df = df.sort_values(by=['Cost'])

fig = plt.figure(figsize=(3,3),dpi=200)
ax = fig.add_subplot(111)
ax.bar(x='index', height='Cost', data=df, color=colormap(normalize(df.Cost.values)))
ax.set_xticklabels([]) 
ax.set_xticks([]) 
ax.plot([frac_prolif_parameters_original_parameters['Cost'][0] for n in range(len(frac_prolif_parameters_pso))], '--', color = 'black', linewidth=1) #original parameters

ax.set_xlabel('Sets of parameters', fontsize = 16)
ax.set_ylabel('Cost', fontsize = 16)


############################################################
#fig 4

conditions = ['0.6 kPa DMSO', '0.6 kPa ERL', '4 kPa DMSO', '4 kPa ERL', '55 kPa DMSO', '55 kPa ERL']

frac_prolif_parameters_pso_longer_sim = frac_prolif_parameters_pso_longer_sim[frac_prolif_parameters_pso_longer_sim['Cost']<np.inf]


frac_prolif_parameters_pso_longer_sim['Cost longer sim'] = frac_prolif_parameters_pso_longer_sim['Cost']
for condition in conditions:
    frac_prolif_parameters_pso_longer_sim['Cost ' + condition + ' longer sim'] = frac_prolif_parameters_pso_longer_sim['Cost ' + condition]


frac_prolif_parameters_pso[frac_prolif_parameters_pso['index'].isin(frac_prolif_parameters_pso_longer_sim.index)]['Cost']


frac_prolif_parameters_pso_full = pd.merge(frac_prolif_parameters_pso_longer_sim[['Cost longer sim']], frac_prolif_parameters_pso[frac_prolif_parameters_pso['index'].isin(frac_prolif_parameters_pso_longer_sim.index)], right_index=True, left_index=True)



plt.scatter(frac_prolif_parameters_pso_full['Cost'],frac_prolif_parameters_pso_full['Cost longer sim'])
plt.scatter(frac_prolif_parameters_pso_full['Cost longer sim'],frac_prolif_parameters_pso_full['Cost'])



fig, axs = plt.subplots(nrows=3, ncols=2, figsize=(10,14))
for ax, condition in zip(axs.flatten(), conditions):
    ax




###############################################################################
###############################################################################
###############################################################################

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
run_times = [120, 120, 120, 120, 60, 60]
conditions = ['0.6 kPa DMSO', '0.6 kPa ERL', '4 kPa DMSO', '4 kPa ERL', '55 kPa DMSO', '55 kPa ERL']
#all_response.keys()

slope_exp = {}
for run_time, condition, num_condition in zip(run_times, conditions, range(6)):
    t = np.arange(0,len(means[num_condition])/3,1/3)
    slope, _, _, _ = np.linalg.lstsq(t[:,np.newaxis], means[num_condition])
    slope_exp[condition]= slope

###############################################################################
###############################################################################
###############################################################################



mypath_cc = r"C:\Users\Grim.M\Ecole\Ens\Master\S2\Stage\Code\fit_data\pso\cost_3_bis_bis_bis\cell_count"

slope_sim = {}

i=0
chain_tmp = str(frac_prolif_parameters_pso[(frac_prolif_parameters_pso.step == 0) & (frac_prolif_parameters_pso.Cost == np.inf)].particle[i])
index_tmp = frac_prolif_parameters_pso[(frac_prolif_parameters_pso.step == 0) & (frac_prolif_parameters_pso.Cost == np.inf)].index[i]


chain_tmp = str(df[df.Cost == df.Cost.sort_values()[0]].particle[0])
index_tmp = df[df.Cost == df.Cost.sort_values()[0]].index[0]


chain_tmp = str(frac_prolif_parameters_pso[(frac_prolif_parameters_pso['particle']==8) & (frac_prolif_parameters_pso['step']==0)].particle[0])
index_tmp = frac_prolif_parameters_pso[(frac_prolif_parameters_pso['particle']==8) & (frac_prolif_parameters_pso['step']==0)].index[0]

cost = {}

fig, axs = plt.subplots(nrows=3, ncols=2, figsize=(10,14))
for ax, run_time, condition, num_condition in zip(axs.flatten(), run_times, conditions, range(6)):
    cond_split = condition.split(' ')
    cell_count_tmp = np.load(mypath_cc + '/' + '/cell_count_' + chain_tmp +'_0.' + index_tmp + '_' + '_'.join(cond_split) + '.npy')
    t = np.arange(0,run_time,1/3)
    cell_count_tmp[:, 0]
    cell_count_t = np.zeros(len(t))
    time_tmp_m = 0
    for time_tmp_p, cell_count_time_tmp in cell_count_tmp:
        cell_count_t[(time_tmp_m <= t) & (t < time_tmp_p)] = cell_count_time_tmp
        time_tmp_m = time_tmp_p
    
    DIP_sim = np.log2(cell_count_t/cell_count_t[0])
    slope_sim[condition], _, _, _ = np.linalg.lstsq(t[:,np.newaxis], DIP_sim)
    #cost_slope_DMSO = (parameters_tmp['slope_4kPa_DMSO'] - slope_4kPa_DMSO)[0]**2
    
    ax.set_title(condition)
    ax.scatter(t, DIP_sim, color = 'blue', s=1)
    ax.plot(t, slope_sim[condition][0]*t, color = 'blue')
    ax.scatter(t, means[num_condition][:len(t)], color = 'red', s=1)
    ax.plot(t, slope_exp[condition]*t, color = 'red')
    
###############################################################################         
    #cost computation: rmse dip rates
    cost['Cost ' + condition] = np.sqrt(np.sum((DIP_sim[1:len(t)] - means[num_condition][1:len(t)])**2 / stdevs[num_condition][1:len(t)]**2))/len(t)  
###############################################################################         

frac_prolif_parameters_pso[frac_prolif_parameters_pso.index==index_tmp][['Cost']]
tmp = frac_prolif_parameters_pso[frac_prolif_parameters_pso.index==index_tmp][['Cost 0.6 kPa DMSO', 'Cost 0.6 kPa ERL', 'Cost 4 kPa DMSO', 'Cost 4 kPa ERL', 'Cost 55 kPa DMSO', 'Cost 55 kPa ERL']]

plt.tight_layout()




###############################################################################
###############################################################################






df = frac_prolif_parameters_pso[frac_prolif_parameters_pso['Cost']<np.inf][['index', 'Cost', 'particle', 'step']]
df = df.sort_values(by=['Cost'])

fig = plt.figure(figsize=(3,3),dpi=200)
ax = fig.add_subplot(111)
ax.bar(x='index', height='Cost', data=df, color=colormap(normalize(df.Cost.values)))
ax.set_xticklabels([]) 
ax.set_xticks([]) 
ax.plot([frac_prolif_parameters_original_parameters['Cost'][0] for n in range(len(frac_prolif_parameters_pso))], '--', color = 'black', linewidth=1) #original parameters

ax.set_xlabel('Sets of parameters', fontsize = 16)
ax.set_ylabel('Cost', fontsize = 16)

