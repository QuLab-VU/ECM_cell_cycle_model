# -*- coding: utf-8 -*-
"""
Created on Mon Jul 27 11:58:46 2020

@author: Grim.M
"""

path=r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code"
import sys
sys.path.insert(0, path + "/cellular_cycle") #add a pathway to the model
sys.path.insert(0, path + "/fit") #add a pathway 
from pydream.core import run_dream
from pysb.integrate import Solver
import numpy as np
from pydream.parameters import SampledParam
from scipy.stats import norm, uniform
import os
import inspect
from pydream.convergence import Gelman_Rubin
import pandas as pd
import pylab as pl

import seaborn as sns
from matplotlib import pyplot as plt

##################################################################################
#Merging the different iterations:


logps={}
sampled_params_dict={}

for chain in range(10):   
    logps_temp = np.array([])
    params_temp = np.array([])
    for it in range(int(20/10)):
        if logps_temp.size==0:
            print(str(10*it))
            logps_temp = vars()['cellular_cycle_dreamzs_5chain_logps_chain_' + str(chain) + '_' + str(10*it)]
            params_temp = vars()['cellular_cycle_dreamzs_5chain_sampled_params_chain_' + str(chain) + '_' + str(10*it)]
        else:
            print(str(10*it))
            logps_temp = np.concatenate([logps_temp, vars()['cellular_cycle_dreamzs_5chain_logps_chain_' + str(chain) + '_' + str(10*it)]])
            params_temp = np.concatenate([ params_temp, vars()['cellular_cycle_dreamzs_5chain_sampled_params_chain_' + str(chain) + '_' + str(10*it)]])
    logps['cellular_cycle_dreamzs_5chain_logps_chain_' + str(chain)] = logps_temp
    sampled_params_dict['cellular_cycle_dreamzs_5chain_sampled_params_chain_' + str(chain)] = params_temp


for name_logps,name in zip(logps.keys(), sampled_params_dict.keys()): 
    np.save(path + "/fit_data/pydream/frac_prolif_quiescence_DMSO_only_30_parameters/gamma10/" + name + "_0_2010", sampled_params_dict[name])
    np.save(path + "/fit_data/pydream/frac_prolif_quiescence_DMSO_only_30_parameters/gamma10/" + name_logps + "_0_2010", logps[name_logps])


#################
    

logps={}
sampled_params_dict={}

for chain in range(10):     
    logps_temp = np.array([])
    params_temp = np.array([])
    for it in range(int(2000/50)):
        if logps_temp.size==0:
            print(str((61+it)*50))
            logps_temp = vars()['cellular_cycle_dreamzs_5chain_logps_chain_' + str(chain) + '_' + str((61+it)*50)]
            params_temp = vars()['cellular_cycle_dreamzs_5chain_sampled_params_chain_' + str(chain) + '_' + str((61+it)*50)]
        else:
            print(str((61+it)*50))
            logps_temp = np.concatenate([logps_temp, vars()['cellular_cycle_dreamzs_5chain_logps_chain_' + str(chain) + '_' + str((61+it)*50)]])
            params_temp = np.concatenate([ params_temp, vars()['cellular_cycle_dreamzs_5chain_sampled_params_chain_' + str(chain) + '_' + str((61+it)*50)]])            
    logps['cellular_cycle_dreamzs_5chain_logps_chain_' + str(chain)] = logps_temp
    sampled_params_dict['cellular_cycle_dreamzs_5chain_sampled_params_chain_' + str(chain)] = params_temp

for name_logps,name in zip(logps.keys(), sampled_params_dict.keys()): 
    np.save(path + "/fit_data/pydream/bis_log_DMSO_only_quiescence_21_parameters/" + name + "_3000-5000", sampled_params_dict[name])
    np.save(path + "/fit_data/pydream/bis_log_DMSO_only_quiescence_21_parameters/" + name_logps + "_3000-5000", logps[name_logps])
    

####################


logps={}
sampled_params_dict={}
for chain in range(10):   
    logps_temp = np.concatenate([vars()['cellular_cycle_dreamzs_5chain_logps_chain_' + str(chain) + '_014000'], vars()['cellular_cycle_dreamzs_5chain_logps_chain_' + str(chain) + '_1450036500']])
    params_temp = np.concatenate([vars()['cellular_cycle_dreamzs_5chain_sampled_params_chain_' + str(chain) + '_014000'],  vars()['cellular_cycle_dreamzs_5chain_sampled_params_chain_' + str(chain) + '_1450036500']])   
    logps['cellular_cycle_dreamzs_5chain_logps_chain_' + str(chain)] = logps_temp
    sampled_params_dict['cellular_cycle_dreamzs_5chain_sampled_params_chain_' + str(chain)] = params_temp

for name_logps,name in zip(logps.keys(), sampled_params_dict.keys()): 
    np.save(path + "/fit_data/pydream/bis_log_DMSO_only_quiescence_21_parameters/" + name + "_0-36500", sampled_params_dict[name])
    np.save(path + "/fit_data/pydream/bis_log_DMSO_only_quiescence_21_parameters/" + name_logps + "_0-36500", logps[name_logps])
    
    
    
###################################################################################
    
#loading these parameters

number_of_chains = 40
number_of_chains = 10
number_of_chains = 14


for chain in range(10): 
    np.load(path + r"/fit_data/pydream/frac_prolif_quiescence_DMSO_only_30_parameters/cellular_cycle_dreamzs_5chain_sampled_params_chain_" + str(chain) + "_0_1000.npy")
    np.load(path + "/fit_data/pydream/frac_prolif_quiescence_DMSO_only_30_parameters/cellular_cycle_dreamzs_5chain_logps_chain_" + str(chain) + "_0_1000.npy")




logps={}
sampled_params_dict={}

for chain in range(number_of_chains):   
    logps_temp = np.array([])
    params_temp = np.array([])
    
    logps_temp = vars()['cellular_cycle_dreamzs_5chain_logps_chain_' + str(chain) + "_0"]
    params_temp = vars()['cellular_cycle_dreamzs_5chain_sampled_params_chain_' + str(chain) + "_0"]          
    logps['cellular_cycle_dreamzs_5chain_logps_chain_' + str(chain)] = logps_temp
    sampled_params_dict['cellular_cycle_dreamzs_5chain_sampled_params_chain_' + str(chain)] = params_temp



sample_parameters=[]
logps_tot=[]

for (logps_chain, chain) in zip(logps.keys(), sampled_params_dict.keys()):
    sample_parameters+=[sampled_params_dict[chain]]
    logps_tot+=[logps[logps_chain].T[0]]


GR =[]

for chain in range(number_of_chains):   
    sampled_parameters_temp = vars()['cellular_cycle_dreamzs_5chain_sampled_params_chain_' + str(chain) + '_0']
    
    GR += [[Gelman_Rubin(sample_parameters)]]

    pl.plot(logps_temp.transpose()[0], label = 'chain_' +str(chain))
    pl.legend() 
    



old_samples = sample_parameters

#np.save(path + "/fit_data/pydream/bis_log/sample_parameters.npy", sample_parameters)


#####################################################################
#plot log evolution


for chain in range(number_of_chains):   
    logps_temp = vars()['cellular_cycle_dreamzs_5chain_logps_chain_' + str(chain) + '_0']
    
    pl.plot(logps_temp.transpose()[0], label = 'chain_' +str(chain))
    pl.legend() 
    
pl.savefig(path + "/../report/figures/pydream/frac_prolif/logps_10chains_gamma10.pdf")


for chain in range(number_of_chains):   
    logps_temp = vars()['cellular_cycle_dreamzs_5chain_logps_chain_' + str(chain) + '_0']
    
    pl.plot(logps_temp.transpose()[0][1500:], label = 'chain_' +str(chain))
    pl.legend() 
pl.savefig(path + "/../report/figures/pydream/frac_prolif/logps_10chains_gamma10_.pdf")


#############################################################################"""
#
#pysb_sampled_parameter_names = ["V_s2skp","V_Fp","k_M","k_EMY","k_CDF","k_CDS","k_ED","K_S","K_MY","K_DF"]
#
#total_iterations = len(old_samples[0])
#burnin = int(total_iterations / 2)
#
#samples = np.concatenate([old_samples[i][burnin:, :] for i in range(10)])
#logps_tot_p = np.concatenate([logps_tot[i][burnin:] for i in range(10)])



samples_after_burning = []
logps_after_burning = []

for chain in range(number_of_chains):   
    logps_temp = vars()['cellular_cycle_dreamzs_5chain_logps_chain_' + str(chain) + '_0']
    sampled_parameters_temp = vars()['cellular_cycle_dreamzs_5chain_sampled_params_chain_' + str(chain) + '_0']
    
    samples_after_burning += [sampled_parameters_temp[np.where(logps_temp>-100)[0][0]:]]
    logps_after_burning += [logps_temp[np.where(logps_temp>-100)[0][0]:]]
    

#before burning
samples = np.concatenate([sample_parameters[i] for i in range(10)])
all_logps = np.concatenate([logps_tot[i] for i in range(10)])


######################################################

#separate death and mitosis parameters

samples_death = np.array([sample[20:30].tolist() for sample in  samples])
samples_mitosis = np.array([sample[0:20].tolist() for sample in  samples])

sample_parameters_death = []
sample_parameters_mitosis = []

for chain in range(number_of_chains):
    samples_death_temp = np.array([sample[20:30].tolist()for sample in  sample_parameters[chain]])
    samples_mitosis_temp = np.array([sample[0:20].tolist()for sample in  sample_parameters[chain]])

    sample_parameters_death+=[samples_death_temp]
    sample_parameters_mitosis+=[samples_mitosis_temp]



#################################################################################
#Plot walk in the space parameters

new_array = [tuple(row) for row in samples]
unique_vectors=np.unique(new_array, axis=0)

new_array = [tuple(row) for row in samples_death]
unique_vectors_death=np.unique(new_array, axis=0)

new_array = [tuple(row) for row in samples_mitosis]
unique_vectors_mitosis=np.unique(new_array, axis=0)


from sklearn import decomposition
import umap

colors = sns.color_palette(n_colors=10)

#####################################
#before burning
#all sampled parameters

#pca
pca = decomposition.PCA(n_components=2)
pca.fit(unique_vectors)
#pca.fit(samples)

X_pca = []
for chain in range(10):   
    X_pca += [pca.transform(sample_parameters[chain])]

fig = pl.figure(figsize=(3,3),dpi=200)
ax = fig.add_subplot(111)
for chain in range(10) :
    ax.scatter(X_pca[chain][:,0], X_pca[chain][:,1],c = colors[chain], alpha = 0.25,s = 0.05);
    ax.plot(X_pca[chain][:,0], X_pca[chain][:,1], c = colors[chain], alpha = 0.5, lw = 0.25, ls = '--');
    ax.plot(X_pca[chain][0,0], X_pca[chain][0,1],c=colors[chain], marker='+')
    ax.plot(X_pca[chain][-1,0], X_pca[chain][-1,1],c=colors[chain], marker='o')

pl.savefig(path + "/../report/figures/pydream/frac_prolif/chains_running_pca_10chains_gamma10_before_burning_unique.pdf")

#umap   
trans = umap.UMAP().fit(unique_vectors)

X_umap = []
for chain in range(10):   
    X_umap += [trans.transform(sample_parameters[chain])]


fig = pl.figure(figsize=(3,3),dpi=200)
ax = fig.add_subplot(111)
for chain in range(10) :
    ax.scatter(X_umap[chain][:,0], X_umap[chain][:,1],c = colors[chain], alpha = 0.25,s = 0.05);
    ax.plot(X_umap[chain][:,0], X_umap[chain][:,1], c = colors[chain], alpha = 0.5, lw = 0.25, ls = '--');
    ax.plot(X_umap[chain][0,0], X_umap[chain][0,1],c=colors[chain], marker='+')
    ax.plot(X_umap[chain][-1,0], X_umap[chain][-1,1],c=colors[chain], marker='o')
pl.savefig(path + "/../report/figures/pydream/frac_prolif/chains_running_umap_10chains_gamma10_before_burning.pdf")


#tsne

from sklearn.manifold import TSNE
X_embedded = TSNE(n_components=2)
X_embedded.fit_transform(samples)

#####
#only death parameters

#pca
pca = decomposition.PCA(n_components=2)
pca.fit(unique_vectors_death)
#pca.fit(samples)


X_pca = []
for chain in range(10):   
    X_pca += [pca.transform(sample_parameters_death[chain])]


fig = pl.figure(figsize=(4,4),dpi=200)
ax = fig.add_subplot(111)

for chain in range(10) :
    ax.scatter(X_pca[chain][:,0], X_pca[chain][:,1],c = colors[chain], alpha = 0.25,s = 0.05);
    ax.plot(X_pca[chain][:,0], X_pca[chain][:,1], c = colors[chain], alpha = 0.5, lw = 0.25, ls = '--');
    ax.plot(X_pca[chain][0,0], X_pca[chain][0,1],c=colors[chain], marker='+')
    ax.plot(X_pca[chain][-1,0], X_pca[chain][-1,1],c=colors[chain], marker='o')

pl.savefig(path + "/../report/figures/pydream/frac_prolif/chains_running_pca_10chains_gamma10_before_burning_unique_death_parameters_only.pdf")

#####
#only mitosis parameters

#pca
pca = decomposition.PCA(n_components=2)
pca.fit(unique_vectors_mitosis)
#pca.fit(samples)


X_pca = []
for chain in range(10):   
    X_pca += [pca.transform(sample_parameters_mitosis[chain])]


fig = pl.figure(figsize=(4,4),dpi=200)
ax = fig.add_subplot(111)

for chain in range(10) :
    ax.scatter(X_pca[chain][:,0], X_pca[chain][:,1],c = colors[chain], alpha = 0.25,s = 0.05);
    ax.plot(X_pca[chain][:,0], X_pca[chain][:,1], c = colors[chain], alpha = 0.5, lw = 0.25, ls = '--');
    ax.plot(X_pca[chain][0,0], X_pca[chain][0,1],c=colors[chain], marker='+')
    ax.plot(X_pca[chain][-1,0], X_pca[chain][-1,1],c=colors[chain], marker='o')

pl.savefig(path + "/../report/figures/pydream/frac_prolif/chains_running_pca_10chains_gamma10_before_burning_unique_mitosis_parameters_only.pdf")



#################################################
#after burning

#pca
pca = decomposition.PCA(n_components=2)
pca.fit(unique_vectors)

X_pca = []
for chain in range(10):   
    X_pca += [pca.transform(samples_after_burning[chain])]

fig = pl.figure(figsize=(3,3),dpi=200)
ax = fig.add_subplot(111)
for chain in range(10) :
    ax.scatter(X_pca[chain][:,0], X_pca[chain][:,1],c = colors[chain], alpha = 0.25,s = 0.05);
    ax.plot(X_pca[chain][:,0], X_pca[chain][:,1], c = colors[chain], alpha = 0.5, lw = 0.25, ls = '--');
    ax.plot(X_pca[chain][0,0], X_pca[chain][0,1],c=colors[chain], marker='+')
    ax.plot(X_pca[chain][-1,0], X_pca[chain][-1,1],c=colors[chain], marker='o')
pl.savefig(path + "/../report/figures/pydream/frac_prolif/chains_running_pca_10chains_gamma10_after_burning.pdf")


#umap   
trans = umap.UMAP().fit(samples)

X_umap = []
for chain in range(10):   
    X_umap += [trans.transform(samples_after_burning[chain])]


fig = pl.figure(figsize=(3,3),dpi=200)
ax = fig.add_subplot(111)
for chain in range(10) :
    ax.scatter(X_umap[chain][:,0], X_umap[chain][:,1],c = colors[chain], alpha = 0.25,s = 0.05);
    ax.plot(X_umap[chain][:,0], X_umap[chain][:,1], c = colors[chain], alpha = 0.5, lw = 0.25, ls = '--');
    ax.plot(X_umap[chain][0,0], X_umap[chain][0,1],c=colors[chain], marker='+')
    ax.plot(X_umap[chain][-1,0], X_umap[chain][-1,1],c=colors[chain], marker='o')
pl.savefig(path + "/../report/figures/pydream/frac_prolif/chains_running_umap_10chains_gamma10_after_burning.pdf")

#tsne

from sklearn.manifold import TSNE

X_embedded = TSNE(n_components=2)
X_embedded.fit_transform(samples)



############################################################

#Test
samples 
all_logps 

max(all_logps)

10**samples[np.argmax(all_logps)]

max(cellular_cycle_dreamzs_5chain_logps_chain_8_0)
10**cellular_cycle_dreamzs_5chain_sampled_params_chain_0_0[np.argmax(cellular_cycle_dreamzs_5chain_logps_chain_0_0)]

min(all_logps)
np.argmin(all_logps)


from cost_function_BngSimulator_bis_reduced import cost_function_reduced
#or
from cost_function_BngSimulator_bis_reduced_DMSO_only import cost_function_reduced_DMSO_only
#or
from cost_function_BngSimulator_bis_reduced_quiescence_DMSO_only import cost_function_reduced_quiescence_DMSO_only

from cost_function_BngSimulator_bis_reduced_quiescence_DMSO_only_21_parameters import cost_function_reduced_quiescence_DMSO_only_21_parameters

from cost_function_BngSimulator_frac_prolif_quiescence_DMSO_only_21_parameters import cost_function_reduced_quiescence_DMSO_only_21_parameters

from cost_function_BngSimulator_frac_death_rate_DMSO_only_switch_death_parameters import  cost_function_reduced_death_rate_DMSO_switch_death_parameters

from cost_function_BngSimulator_frac_prolif_quiescence_DMSO_only_30_parameters import  cost_function_reduced_quiescence_DMSO_only_30_parameters


parameter_vector = 10**samples[15]

parameter_vector = 10**samples[np.argmax(all_logps)]

parameter_vector = 10**samples[np.argmin(all_logps)]

##bis_log
#parameter_vector = np.array([0.76062416, 0.02256117, 0.30398523, 0.94528072, 1.34990567, 0.30460912, 0.11683466, 0.10229111, 0.82909919, 0.74736298])
#parameter_vector = np.array([0.55205378, 1.71817805, 0.86381869, 1.01441367, 1.76416909, 1.3477477 , 0.01771722, 0.67373353, 0.01335147, 1.55034297])
#parameter_vector = np.array([5.20103250e-02, 4.42837548e-02, 2.39078059e-01, 3.91191419e-02, 1.71177745e-04, 1.23105202e-04, 4.61787386e-01, 1.12304753e+00, 1.59868274e-02, 4.27401052e-01])


##bis_log_DMSO_only:
#best param (before 24600)
#parameter_vector = np.array([4.11575939e-04, 1.15914699e-01, 1.45504306e-01, 1.21348177e-02, 1.43575822e-03, 8.87765894e-03, 5.55602576e-01, 1.03522823e-03,1.69709926e-01, 2.92699744e-03])

#bis_log_quiescence_DMSO_only
#best param before 6000:
#parameter_vector= np.array([1.80358494e+00, 2.14465000e-03, 4.40430589e-01, 3.49731647e-03, 1.09110080e-03, 4.00638043e-01, 5.31736197e-03, 3.80273671e-03,       2.28068478e+00, 1.89408293e+00])

#bis_log_quiescence_DMSO_only
#best param before 13000:
#parameter_vector= np.array([6.65561930e-02, 3.78074446e-03, 2.04759230e-01, 2.53518251e-03, 1.63049830e-02, 4.10932785e-01, 5.43363731e-02, 1.97787556e-03,       3.71571282e+00, 2.01144088e+00]) #test with paramters 2

#bis_log_quiescence_DMSO_only wiht 21 parameters
#parameter_vector = np.array([2.88095027e-02, 1.08653936e-04, 1.06190007e-03, 5.91535448e-02, 9.63052769e-04, 1.12765757e-04, 3.35453630e-02, 1.18449277e-04, 1.24550538e-01, 2.00076132e-03, 1.28792094e-01, 6.23484593e-03, 8.39769019e-02, 2.02874490e-02, 2.16102188e+00, 6.79877479e-04, 1.18330213e-03, 8.82600428e+00, 3.06767094e-04, 1.46608884e-04, 3.84589838e-02])

#frac prolif
#parameter_vector = np.array([5.08434023e-03, 4.14619045e-01, 1.58080261e+00, 1.87548356e-02, 2.43834294e-04, 2.67242878e-03, 9.17533143e-02, 2.39374299e-04, 3.14708506e-03, 1.37142120e-01, 7.69247749e-01, 1.41305571e-04, 1.90911690e-04, 1.03374254e-03, 3.77969796e+00, 1.22245523e-03, 8.07754038e-04, 6.38585894e+00, 4.09374334e-02, 1.73584598e-03, 4.19911872e-01])

#frac prolif 30
#parameter_vector = np.array([1.35981351e-03, 3.27621894e-04, 3.33605499e-01, 1.67110971e+01, 1.33333794e+01, 7.26881136e-04, 7.24728198e+01, 1.56434290e+01, 8.28685189e-02, 3.77250858e+01, 1.46225845e+00, 4.07452470e-03, 2.57708238e+00, 1.69128821e-02, 4.33971772e-01, 3.08585275e+01, 1.58804859e-02, 3.60589231e-03, 8.64121317e-02, 5.59717739e-01, 3.23506196e-01, 2.41717684e-04, 2.24990105e+01, 3.05306899e-01, 3.12691267e+00, 4.82800772e-02, 1.72436648e-03, 1.11840705e-02, 2.60040648e-02, 3.36539927e-02])

#parameter_vector = np.array([9.73436639e+01, 1.20482457e+00, 9.44180958e-03, 8.09900119e+00, 2.78845243e-01, 3.87094008e+01, 1.65148700e-03, 1.46029110e+00, 1.81932279e-03, 8.32673504e+01, 1.31658430e+00, 9.29026540e-04, 4.17686561e+00, 4.23690936e-01, 1.62391120e-01, 1.90052578e-03, 1.63257207e-03, 1.12756428e-04, 2.55672808e-02, 1.47315792e-04, 2.15972787e-04, 8.14881900e+01, 5.62618511e+00, 4.78496152e-02, 2.84263099e+01, 2.24898042e-03, 1.87618954e-03, 8.22937520e-01, 4.10033319e+00, 4.31574492e-04])

#parameter_vector = np.array([1.14637063e-01, 2.13614373e-04, 2.08256410e-02, 4.84990660e-02, 6.03643430e+00, 1.74932077e-04, 2.65252075e-04, 4.81838284e+00, 7.25610453e-03, 3.67227612e-02, 9.11806494e-02, 4.50512863e+01, 5.21316141e-03, 5.09559069e-03, 2.85905519e+00, 1.36298838e-02, 2.10806043e-03, 7.26355984e+00, 1.16149358e-03, 1.14447665e-03, 4.05081200e-03, 7.21596812e+01, 4.07609420e+00, 4.40626117e-04, 1.54503854e+00, 1.22937820e-03, 8.39269962e-03, 5.11573206e+00, 1.34492066e-02, 4.18501642e-01])

    
logps_tot_p[350]

test0 = cost_function_reduced(parameter_vector)
test = cost_function_reduced_DMSO_only(parameter_vector)
test_q = cost_function_reduced_quiescence_DMSO_only(parameter_vector)
test_21 = cost_function_reduced_quiescence_DMSO_only_21_parameters(parameter_vector)

test_death = cost_function_reduced_death_rate_DMSO_switch_death_parameters(parameter_vector)

test = cost_function_reduced_quiescence_DMSO_only_30_parameters(parameter_vector)
##############################################################################
#plot distributions

###plot


ndims = len(old_samples[0][0])
colors = sns.color_palette(n_colors=ndims)

for dim in range(ndims):
    fig = plt.figure()
    sns.distplot(samples[:, dim], color=colors[dim])
    fig.savefig(r'C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\report\figures\pydream\bis_log_quiescence_DMSO_only\PyDREAM_example_cellular_cycle_dimension_' + str(dim) + "=" + pysb_sampled_parameter_names[dim] + ".pdf")
    
    
dict_param={}
#pysb_sampled_parameter_names = ["V_s2skp","V_Fp","k_M","k_EMY","k_CDF","k_CDS","k_ED","K_S","K_MY","K_DF"]
#or
#pysb_sampled_parameter_names = ["V_s2skp","V_Fp","k_M","k_EMY","k_CDF","k_CDS","k_ED","K_S","K_MY","K_DF","d_CD","k_P1", "k_P2", "K_CD", "k_E", "K_M", "K_CE", "V_2E2F", "K_2E2F", "k_RE", "k_b"]

#pysb_sampled_parameter_names = ["kf1","kf2","kb","kasbaxmbcl2","kdsbaxmbcl2","k1sbh3","k2sbh3","kdbh3","kasbh3bcl2","kdsbh3bcl2"]


for p in range(10):
    dict_param[pysb_sampled_parameter_names[p]]=samples[:, p]

df=df = pd.DataFrame(dict_param)

fig = pl.figure()
sns.violinplot(data=df)
fig.savefig(r'C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\report\figures\pydream\bis_log_quiescence_DMSO_only_21_parameters\PyDREAM_parameter_distribution_violin.pdf')



from pathlib import Path
fn = Path(r'C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\report\figures\pydream\bis_log_quiescence_DMSO_only_21_parameters\PyDREAM_parameter_distribution.pdf').expanduser()

a=5
b=5
fig, axs = pl.subplots(ncols=a, nrows=b, sharex=True, sharey=True, figsize=(15,6))
for dim, ax in zip(range(10), axs.flat):
        sns.distplot(samples[:, dim], color=colors[dim], ax=ax, kde_kws={"lw": 0.5, "label": "variable: " + pysb_sampled_parameter_names[dim]}, bins=np.linspace(-4,1,50))
        ax.legend(loc='upper left', fontsize = 6)
pl.savefig(fn, bbox_inches='tight')





##############################################################################
#Dimensionnality reduction and clustering

#PCA

from sklearn import decomposition

pca = decomposition.PCA(n_components=2)
pca.fit(samples)
X1 = pca.transform(samples)


dimension_reduction="PCA"

#umap
import umap

trans = umap.UMAP().fit(samples)
X1=trans.transform(samples)
dimension_reduction="umap"


#tsne
from sklearn.manifold import TSNE

X_embedded = TSNE(n_components=2).fit_transform(X)


dimension_reduction="tsne"

############plot dimensionality reduction

max(logps_tot_p)
np.argmax(logps_tot_p)

samples[np.argmin(logps_tot_p)]

best_params_5p = np.array([])

threshold=0.75*(max(logps_tot_p)-min(logps_tot_p))+min(logps_tot_p)


np.where(logps_tot_p >threshold)



best_param_thresh=samples[np.where(logps_tot_p >threshold)[0],:]
best_logps_thresh=logps_tot_p[np.where(logps_tot_p >threshold)[0]]
best_X1_thresh=X1[np.where(logps_tot_p >threshold)[0]]


plt.scatter(X1[:, 0], X1[:, 1], c=logps_tot_p, cmap='rainbow')
plt.colorbar()

df = pd.DataFrame(X1, columns=["x", "y"])
best_X1_thresh = pd.DataFrame(best_X1_thresh, columns=["x", "y"])



sns.kdeplot(df['x'], df['y'], cmap='magma_r')
#plt.scatter(best_X1_thresh['x'], best_X1_thresh['y'], marker='x', c=best_logps_thresh, alpha=0.5)
#plt.colorbar()
plt.xlabel(dimension_reduction + ' 1')
plt.ylabel(dimension_reduction + ' 2')
plt.title("Samples")
plt.savefig(path + "/fit_data/pydream/bis_log/density_plot_"+ dimension_reduction+'_best_samples_thresh_75p.pdf')



####Clustering

#KMEANS
from sklearn.cluster import KMeans

y_pred = KMeans(n_clusters=6).fit_predict(X1)

clustering_method="KMeans"

############
#Process clusters
sample_plot={}
sample_clusters={}
logps_clusters={}

for g in np.unique(y_pred):
    ix = np.where(y_pred == g)
    df=pd.DataFrame(X1[ix], columns=["x", "y"])
    sample_plot[str(g)]=pd.DataFrame(X1[ix], columns=["x", "y"])
    sample_clusters[str(g)]=samples[ix]
    logps_clusters[str(g)]=logps_tot_p[ix]
    
########################
#plot clusters
    
plt.scatter(X1[:, 0], X1[:, 1], c=y_pred, s=2)


color=["Blues_r", "Blues", "Greens", "Greys", "Oranges", "Purples", "Reds", "Wistia", "autumn", "bone", "Greens_r", "Greys_r", "Oranges_r", "Purples", "Reds", "Wistia", "autumn", "bone", "bwr", "cividis","cool", "coolwarm", "coolwarm_r", "copper", "copper_r", "cubehelix", "cubehelix_r", "flag", "flag_r", "gist_earth", "gist_earth_r", "gist_gray", "gist_gray_r", "gist_heat", "gist_heat_r", "gist_ncar", "gist_ncar_r", "gist_rainbow", "gist_rainbow_r", "gist_stern", "gist_stern_r", "gist_yarg", "gist_yarg_r", "gnuplot", "gnuplot2", "gnuplot2_r", "gnuplot_r", "gray", "gray_r", "hot", "hot_r", "hsv", "hsv_r", "icefire", "icefire_r", "inferno", "inferno_r", "jet", "jet_r", "magma", "magma_r", "mako", "mako_r", "nipy_spectral", "nipy_spectral_r", "ocean", "ocean_r", "pink", "pink_r", "plasma", "plasma_r", "prism", "prism_r", "rainbow", "rainbow_r", "rocket", "rocket_r", "seismic", "seismic_r", "spring", "spring_r", "summer", "summer_r", "tab10", "tab10_r", "tab20", "tab20_r", "tab20b", "tab20b_r", "tab20c", "tab20c_r", "terrain", "terrain_r", "twilight", "twilight_r", "twilight_shifted", "twilight_shifted_r", "viridis", "viridis_r", "vlag", "vlag_r", "winter", "winter_r"]


plt.ioff()
f, ax = plt.subplots(figsize=(15, 15))
plt.ioff()
f, ax = plt.subplots(figsize=(15, 15))
for (cluster,number) in  zip(sample_plot.keys(),range(len(sample_plot.keys()))):
    ax.tick_params(axis='both', which='major', labelsize=22)
    ax = sns.kdeplot(sample_plot[cluster].x, sample_plot[cluster].y, cmap=color[number],   shade=True, shade_lowest=False, alpha=0.66, legend=False, cbar=False)
plt.axis('equal')
plt.xlabel(dimension_reduction + " 1", fontsize=22)
plt.ylabel(dimension_reduction + " 2", fontsize=22)
plt.savefig(path + "/fit_data/pydream/bis_log/clustering_parameters_dimension_reduction_" + dimension_reduction + "_clustering_" + clustering_method + ".pdf")




###############################################################################
#extracting parameters

###############
#best point among the clusters
logps_clusters
sample_clusters.keys()

best_parameters=[]
best_parameters_logps=[]
best_parameters_plot=[]


for cluster in logps_clusters.keys():
    print(cluster)
    print(max(logps_clusters[cluster]))
    best_parameters_logps+=[max(logps_clusters[cluster])]
    best_parameters+=[sample_clusters[cluster][np.argmax(logps_clusters[cluster])]]


#X2 = trans.transform(best_parameters)
X2 = pca.transform(best_parameters)

df2 = pd.DataFrame(np.concatenate([X2.T, -np.array([best_parameters_logps])]).T, columns=["x", "y","size"])

df2["size"]=df2["size"]*10



plt.ioff()
f, ax = plt.subplots(figsize=(15, 15))
for (cluster,number) in  zip(sample_plot.keys(),range(len(sample_plot.keys()))):
    ax.tick_params(axis='both', which='major', labelsize=22)
    ax = sns.kdeplot(sample_plot[cluster].x, sample_plot[cluster].y, cmap=color[number],   shade=True, shade_lowest=False, alpha=0.66, legend=False, cbar=False)
sns.scatterplot(x="x", y="y", size="size", data=df2, legend=False)
plt.axis('equal')
plt.xlabel("PCA 1", fontsize=22)
plt.ylabel("PCA 2", fontsize=22)
plt.savefig(path + "/fit_data/pydream/bis_log/bestparam_bestlogps_" + dimension_reduction + "_clustering_" + clustering_method + ".pdf")


np.save(path + "/fit_data/pydream/bis/best_param_bestlog.npy", best_parameters)

############
#mean of different points

best=[]
for cluster in sample_clusters.keys():
    best+=[np.mean(sample_clusters[cluster],axis=0)]

X3 = pca.transform(best)
#X3 = trans.transform(best)


df3 = pd.DataFrame(X3, columns=["x", "y"])

plt.ioff()
f, ax = plt.subplots(figsize=(15, 15))

for (cluster,number) in  zip(sample_plot.keys(),range(len(sample_plot.keys()))):
    ax.tick_params(axis='both', which='major', labelsize=22)
    ax = sns.kdeplot(sample_plot[cluster].x, sample_plot[cluster].y, cmap=color[number],   shade=True, shade_lowest=False, alpha=0.66, legend=False, cbar=False)
ax = sns.scatterplot(x="x", y="y", data=df3)
plt.axis('equal')
plt.xlabel("PCA 1", fontsize=22)
plt.ylabel("PCA 2", fontsize=22)
plt.savefig(path + "/fit_data/pydream/bis/clustering_parameters_best_mean_" + dimension_reduction + "_clustering_" + clustering_method + ".pdf")


np.save(path + "/fit_data/pydream/bis/params_from_mean.npy", best)

##############
#weighted average

weighted={}
best_parameters_weighted=[]

for cluster in logps_clusters.keys():
    print(cluster)
    print(max(logps_clusters[cluster]))
    
    weighted[cluster]= (logps_clusters[cluster] - min(logps_clusters[cluster])) / (max(logps_clusters[cluster]) - min(logps_clusters[cluster]))
    
    best_parameters_weighted+=[sum(sample_clusters[cluster]*weighted[cluster][:,np.newaxis])/(len(weighted[cluster]))]


X4 = pca.transform(best_parameters_weighted)
#X4 = trans.transform(best_parameters_weighted)

df4 = pd.DataFrame(X4, columns=["x", "y"])

plt.ioff()
f, ax = plt.subplots(figsize=(15, 15))

for (cluster,number) in  zip(sample_plot.keys(),range(len(sample_plot.keys()))):
    ax.tick_params(axis='both', which='major', labelsize=22)
    ax = sns.kdeplot(sample_plot[cluster].x, sample_plot[cluster].y, cmap=color[number],   shade=True, shade_lowest=False, alpha=0.66, legend=False, cbar=False)
ax = sns.scatterplot(x="x", y="y", data=df4)
plt.axis('equal')
plt.xlabel("PCA 1", fontsize=22)
plt.ylabel("PCA 2", fontsize=22)
plt.savefig(path + "/fit_data/pydream/bis/clustering_parameters_best_weighted_mean_" + dimension_reduction + "_clustering_" + clustering_method + ".pdf")


np.save(path + "/fit_data/pydream/bis/params_from_mean.npy", best)


#####################
#max density



def RBS(x,x0,sigma) :
    observations , dimensions= np.shape(x0)
    s=0
    for i in range(observations):
        s+=math.exp(-np.linalg.norm(x-x0[i])**2/(2*sigma**2)) / (sigma*math.sqrt(2*math.pi))
    if observations != 0:
        return s/observations
    else:
        return s


#########""
best_params=[]
density_rbs={}

sigma=50
p = np.linspace(-4,1,2000)
pysb_sampled_parameter_names = ["V_s2skp","V_Fp","k_M","k_EMY","k_CDF","k_CDS","k_ED","K_S","K_MY","K_DF"]

for (dim, param) in zip(range(10), pysb_sampled_parameter_names):
    density_rbs[param]=[]
    density_rbs[param]=[RBS(i,np.array([samples[:, dim]]),sigma) for i in p]
    best_params+=[p[np.argmax(density_rbs[param])]]

plt.plot(p,density_rbs[param])

a=5
b=2
fig, axs = pl.subplots(ncols=a, nrows=b, sharex=True, sharey=True, figsize=(15,6))
for param, dim, ax in zip(pysb_sampled_parameter_names, range(a*b), axs.flat):
        ax.plot(p,density_rbs[param], color=colors[dim], label=param)
        ax.scatter(p[np.argmax(density_rbs[param])], max(density_rbs[param]))
        ax.legend(loc='upper left', fontsize = 6)



sigma=2


x=np.linspace(0,2,200)
sample_clusters["0"][:,0]


y={}
for cluster in sample_clusters.keys():
    cluster_temp=pd.DataFrame()
    for param in range(10):
        rbs_temp=[RBS(i,np.transpose(np.array([sample_clusters[cluster][:,param]])),sigma) for i in x]
        cluster_temp[pysb_sampled_parameter_names[param]] = pd.Series(np.array(rbs_temp))
    y[cluster]=cluster_temp



best_param={}

a=2
b=5
for cluster in sample_clusters.keys():
    best_param[cluster]={}
    fig, axs = pl.subplots(ncols=a, nrows=b)
    for i in range(a):
        for j in range(b):
            param=pysb_sampled_parameter_names[b*i+j]
            density=y[cluster][param]
            bb = np.argmax(density)
            print(x[bb])
            best_param[cluster][param]=[x[bb]]
            aa = y[cluster][param][bb]
            axs[j,i].plot(x,density, label=param)
            axs[j,i].scatter(x[bb], aa)
            axs[j,i].legend(fontsize=5)
            
"""  
    best_param[cluster] = pd.DataFrame(data=best_param[cluster])
    np.save(path + "/fit_data/pydream/bis/best_parameters_max_density_"+cluster+".npy",best_param[cluster].values[0])
    for ax in axs.flat:
        ax.label_outer()
"""


pd.DataFrame(data=best_param)


max_density_cluster=[]

for cluster in sample_clusters.keys():
    max_density_cluster+=[np.array(list(best_param[cluster].values())).T[0]]



#X5 = pca.transform(max_density_cluster)
X5 = trans.transform(max_density_cluster)

df5 = pd.DataFrame(X5, columns=["x", "y"])

plt.ioff()
f, ax = plt.subplots(figsize=(15, 15))

for (cluster,number) in  zip(sample_plot.keys(),range(len(sample_plot.keys()))):
    ax.tick_params(axis='both', which='major', labelsize=22)
    ax = sns.kdeplot(sample_plot[cluster].x, sample_plot[cluster].y, cmap=color[number],   shade=True, shade_lowest=False, alpha=0.66, legend=False, cbar=False)
ax = sns.scatterplot(x="x", y="y", data=df4)
plt.axis('equal')
plt.xlabel("PCA 1", fontsize=22)
plt.ylabel("PCA 2", fontsize=22)
plt.savefig(path + "/fit_data/pydream/bis/clustering_parameters_best_max_density_" + dimension_reduction + "_clustering_" + clustering_method + ".pdf")





###################
#test

from cost_function_BngSimulator_bis import cost_function


lower_limits = [0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0]
parameters_to_sample = SampledParam(uniform, loc=lower_limits, scale=2) #loc is the lower bdry and scale is the size of the distribution
sampled_parameter_names = [parameters_to_sample]
niterations = 5
converged = False
total_iterations = niterations
nchains = 3

sampled_params, log_ps = run_dream(sampled_parameter_names, cost_function_f, niterations=niterations, nchains=nchains, multitry=False, gamma_levels=4, adapt_gamma=True, history_thin=1, model_name='teeeest', verbose=True)


parameter_vector=sampled_params[-1][-1]
log_ps[-1][-1]


cost_function(parameter_vector)


