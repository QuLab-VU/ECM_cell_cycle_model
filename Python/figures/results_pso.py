# -*- coding: utf-8 -*-
"""
Created on Mon Jun  7 20:57:16 2021

@author: Grim.M
"""


import numpy as np
import pandas as pd
import pylab as pl
import random 
import sobol_seq as sbs
import seaborn as sns
import re
from sklearn import decomposition
import matplotlib.cm as cm
import matplotlib.colors as mcolors


################################################

particle_positions_pso = pd.read_pickle(r"C:/Users/Grim.M/Ecole/ens/Master/S2/Stage/Code/fit_data/pso/high_and_low_stiffness_cost/compilation/particle_positions_pso.pkl")

frac_prolif_parameters_pso = pd.read_pickle(r"C:/Users/Grim.M/Ecole/ens/Master/S2/Stage/Code/fit_data/pso/high_and_low_stiffness_cost/compilation/frac_prolif_parameters_all_hyperparameters_pso.pkl")

position_history = np.load(r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code\fit_data\pso\high_and_low_stiffness_cost\compilation\position_history.npy", allow_pickle=True)

frac_prolif_parameters_006_55 = pd.read_pickle(r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code\fit_data\sobol_seq\cost_function_model_4\cost_and_IMT_distribution\frac_prolif_parameters_all_parameters_006_55_kPa.pkl")
frac_prolif_parameters = frac_prolif_parameters_006_55


#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################

title=["Hyperparameters 1: \n {'c1': 0.5, 'c2': 0.3, 'w': 0.9, 'k': 3, 'p': 2} \n 8 particles \n 100 steps", "Hyperparameters 2: \n {'c1': 0.05, 'c2': 0.1, 'w': 0.1, 'k': 6, 'p': 2} \n 16 particles \n 40 steps","Hyperparameters 3: \n {'c1': 0.1, 'c2': 0.3, 'w': 0.3, 'k': 5, 'p': 2} \n 16 particles \n 20 steps", "Hyperparameters 4: \n {'c1': 0.2, 'c2': 0.3, 'w': 0.4, 'k': 5, 'p': 2} \n 16 particles \n 30 steps", "Hyperparameters 5: \n {'c1': 0.15, 'c2': 0.3, 'w': 0.4, 'k': 5, 'p': 2} \n 30 particles \n 39 steps"]

unique_vectors=[]
unique_vectors += [np.concatenate(position_history[0])]
unique_vectors += [np.concatenate(position_history[1])]
unique_vectors += [np.concatenate(position_history[2])]
unique_vectors += [np.concatenate(position_history[3])]
unique_vectors += [np.concatenate(position_history[4])]



#pca&
pca_1 = decomposition.PCA(n_components=2)
pca_1.fit(unique_vectors[0])
pca_2 = decomposition.PCA(n_components=2)
pca_2.fit(unique_vectors[1])
pca_3 = decomposition.PCA(n_components=2)
pca_3.fit(unique_vectors[2])
pca_4 = decomposition.PCA(n_components=2)
pca_4.fit(unique_vectors[3])
pca_5 = decomposition.PCA(n_components=2)
pca_5.fit(unique_vectors[4])

X_pca = []
X_pca+=[[]]
for particle in range(np.shape(position_history[0])[1]):   
    X_pca[0] += [pca_1.transform(position_history[0][:,particle,:])]
X_pca+=[[]]
for particle in range(np.shape(position_history[1])[1]):   
    X_pca[1] += [pca_2.transform(position_history[1][:,particle,:])]
X_pca+=[[]]
for particle in range(np.shape(position_history[2])[1]):   
    X_pca[2] += [pca_3.transform(position_history[2][:,particle,:])]
X_pca+=[[]]
for particle in range(np.shape(position_history[3])[1]):   
    X_pca[3] += [pca_4.transform(position_history[3][:,particle,:])]
X_pca+=[[]]
for particle in range(np.shape(position_history[4])[1]):   
    X_pca[4] += [pca_5.transform(position_history[4][:,particle,:])]



fig, axs = pl.subplots(nrows=2, ncols=5, figsize=(30,12))

for hyp in range(5) :
    colors = sns.color_palette(n_colors=np.shape(position_history[hyp])[1])
    for chain in range(np.shape(position_history[hyp])[1]) :
        axs[0, hyp].scatter(X_pca[hyp][chain][:,0], X_pca[hyp][chain][:,1],c = colors[chain], alpha = 0.25,s = 0.05);
        axs[0, hyp].plot(X_pca[hyp][chain][:,0], X_pca[hyp][chain][:,1], c = colors[chain], alpha = 0.5, lw = 0.25, ls = '--');
        axs[0, hyp].plot(X_pca[hyp][chain][0,0], X_pca[hyp][chain][0,1],c=colors[chain], marker='+')
        axs[0, hyp].plot(X_pca[hyp][chain][-1,0], X_pca[hyp][chain][-1,1],c=colors[chain], marker='o')
        #axs[0, hyp].title.set_text(title[hyp])
        axs[0, hyp].set_title(title[hyp],fontsize= 16)
        
    axs[0, hyp].set_xlabel('PCA 1', fontsize = 16) 

    a = np.sort(frac_prolif_parameters_pso[frac_prolif_parameters_pso['hyperparameters']==hyp+1]['cost_total'].tolist())
    axs[1, hyp].bar(range(len(a)), a)
    axs[1, hyp].plot([frac_prolif_parameters['cost_total'].loc[frac_prolif_parameters['index'] == 1] for n in range(len(a))], 'red')
    axs[1, hyp].set_xlabel('Sets of parameters', fontsize = 16) 

axs[0, 0].set_ylabel('PCA 2', fontsize = 16) 
axs[1, 0].set_ylabel('Cost', fontsize = 16) 


pl.savefig(path + "/../report/figures/pso/all_particles_running_pca_all_hyperparameters.pdf", bbox_inches='tight')

#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################

title=["Hyperparameters 1: \n {'c1': 0.5, 'c2': 0.3, 'w': 0.9, 'k': 3, 'p': 2} \n 8 particles \n 100 steps", "Hyperparameters 2: \n {'c1': 0.05, 'c2': 0.1, 'w': 0.1, 'k': 6, 'p': 2} \n 16 particles \n 40 steps","Hyperparameters 3: \n {'c1': 0.1, 'c2': 0.3, 'w': 0.3, 'k': 5, 'p': 2} \n 16 particles \n 20 steps", "Hyperparameters 4: \n {'c1': 0.2, 'c2': 0.3, 'w': 0.4, 'k': 5, 'p': 2} \n 16 particles \n 30 steps", "Hyperparameters 5: \n {'c1': 0.15, 'c2': 0.3, 'w': 0.4, 'k': 5, 'p': 2} \n 30 particles \n 39 steps"]


colormap = cm.rainbow
normalize = mcolors.Normalize(vmin = min(frac_prolif_parameters_pso["cost_total"]), vmax = max(frac_prolif_parameters_pso["cost_total"]))

fig, axs = pl.subplots(nrows=1, ncols=5, figsize=(34,5))

for hyp in range(5) :
    axs[hyp].title.set_text(title[hyp])
    for chain in range(np.shape(position_history[hyp])[1]) :
        axs[hyp].plot(X_pca[hyp][chain][:,0], X_pca[hyp][chain][:,1], c = 'grey', alpha = 0.5, lw = 0.25, ls = '--');
        axs[hyp].plot(X_pca[hyp][chain][0,0], X_pca[hyp][chain][0,1], c = 'grey', marker='+', alpha = 0.6)
        axs[hyp].plot(X_pca[hyp][chain][-1,0], X_pca[hyp][chain][-1,1], c = 'grey', marker='o', alpha = 0.6)
        for position_ind in range(np.shape(position_history[hyp])[0]) :    
            try :
                identification_temp = particle_positions_pso[particle_positions_pso['parameter_vector']==position_history[hyp][position_ind,chain,0]].iloc[0]['identification']
                
                frac_prolif_parameters_temp = frac_prolif_parameters_pso[frac_prolif_parameters_pso['index']==identification_temp]
                
                cost_temp = frac_prolif_parameters_temp.iloc[0]['cost_total']
                
                axs[hyp].scatter(X_pca[hyp][chain][position_ind,0], X_pca[hyp][chain][position_ind,1],c = colormap(normalize(cost_temp)) , alpha = 1, s = 5);
            except :
                print(identification_temp)
                axs[hyp].scatter(X_pca[hyp][chain][position_ind,0], X_pca[hyp][chain][position_ind,1],c = 'white', alpha = 0.25, s = 5, edgecolor='black', linewidth=0.3);

scalarmappaple = cm.ScalarMappable(norm=normalize, cmap=colormap)
scalarmappaple.set_array([])

#fig.subplots_adjust(right=0.8)
#cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
cb = fig.colorbar(scalarmappaple, ax=axs.ravel().tolist(), shrink=0.8)
cb.set_label(label='cost_total', fontsize=24)
cb.ax.tick_params(labelsize=12)

pl.savefig(path + "/../report/figures/pso/all_particles_running_pca_hyperparameters_cost.pdf", bbox_inches='tight')


#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################
###############################################################################################
###############################################################################################
#plotting best costs


from plot_frac_div_q_ap1 import all_events
import fracprolif as fp
import math

###############################################################################################
###############################################################################################
######
#parameters for the simulation
init = 50
run_time = 102.3
observed_cells_time = 30
omega=200

#parameters of the frac prolif model to compare
parameter_names = ['a','d','q']

#######################################################################################
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

x=np.linspace(0,120,360)
experimental_IMT_distribution=[RBS(i, np.transpose([mitotic_lifespans_exp.tolist()]), sigma) for i in x]



sigma=4
 
g1_distribution_exp = [RBS(i, np.transpose([experimental_G1s['55 kPa DMSO']]), sigma) for i in np.linspace(0,120,400)]
sm_distribution_exp = [RBS(i, np.transpose([experimental_SMs['55 kPa DMSO']]), sigma) for i in np.linspace(0,120,400)]


##################################################################################
hyp = 4
hyperparameters_folder = "/hyperparameters_5"


g1_distribution_55 = pd.read_pickle(r"C:/Users/Grim.M/Ecole/ens/Master/S2/Stage/Code/fit_data/pso/high_and_low_stiffness_cost" + hyperparameters_folder + "/cost_and_IMT_distribution/g1_distribution_all_parameters_50_initial_cells_55.pkl")

sm_distribution_55 = pd.read_pickle(r"C:/Users/Grim.M/Ecole/ens/Master/S2/Stage/Code/fit_data/pso/high_and_low_stiffness_cost" + hyperparameters_folder + "/cost_and_IMT_distribution/sm_distribution_all_parameters_50_initial_cells_55.pkl")


IMT_distribution = pd.read_pickle(r"C:/Users/Grim.M/Ecole/ens/Master/S2/Stage/Code/fit_data/pso/high_and_low_stiffness_cost" + hyperparameters_folder + "/cost_and_IMT_distribution/IMT_distribution_all_parameters_50_initial_cells.pkl")

####################################################################################

from pathlib import Path

fn = Path(path + "/../report/figures/pso/results_pca_hyperparameters_5.pdf")

fig, axs = pl.subplots(ncols=2, nrows=2, sharex=False, sharey=False, figsize=(40,20), constrained_layout=False)

############### 1st subplot
minimum = frac_prolif_parameters_pso.sort_values(by=['cost_total'])[0:50]['cost_total'].tolist()[2]
#cost_limit = frac_prolif_parameters_pso.sort_values(by=['cost_total'])[0:60]['cost_total'].tolist()[-1]
cost_limit = frac_prolif_parameters.sort_values(by=['cost_total'])[0:1]['cost_total'].tolist()[-1]


for ax, col in zip(axs[0], ['55 kPa DMSO','0.6 kPa DMSO']):
    ax.annotate(col, xy=(0.5, 1), xytext=(0, 5),
                xycoords='axes fraction', textcoords='offset points',
                fontsize=36, ha='center', va='baseline')

x=np.linspace(0,120,360)
colormap = cm.RdYlBu
normalize = mcolors.Normalize(vmin = minimum, vmax = cost_limit)
#or 

for idx in frac_prolif_parameters_pso[frac_prolif_parameters_pso['hyperparameters']==hyp+1]['index']:
    print(idx)
    #idx = int(idx)
    if frac_prolif_parameters_pso.loc[frac_prolif_parameters_pso['index']==idx].iloc[0]['cost_total']<cost_limit :
    #if frac_prolif_parameters_pso.loc[frac_prolif_parameters_pso['index']==idx].iloc[0]['cost_mmd']< inf :
        axs[0, 1].plot(x, IMT_distribution[idx].tolist(), color=colormap(normalize(frac_prolif_parameters_pso.loc[frac_prolif_parameters_pso['index']==idx].iloc[0]['cost_total'])), lw=0.4)
#        axs[0, 1].grid(True)
        axs[0, 1].tick_params(axis="y", labelsize=24)
        axs[0, 1].tick_params(axis="x", labelsize=24)
#plotting original parameters and expt
axs[0, 1].plot(x, experimental_IMT_distribution, color='black', lw=6)
#axs[0, 1].plot(x, IMT_distribution[1].tolist(), lw=5, color=colormap(normalize(frac_prolif_parameters_pso.loc[frac_prolif_parameters_pso['index']==1.0].iloc[0]['cost_total'])), path_effects=[mpe.withStroke(linewidth=6, foreground='black')])


axs[0, 1].set_ylabel('Density', fontsize = 24) 
axs[0, 1].set_xlabel('Intermitotic times', fontsize = 24) 

############### 2nd subplot


for idx in frac_prolif_parameters_pso[frac_prolif_parameters_pso['hyperparameters']==hyp+1]['index']:
    if frac_prolif_parameters_pso.loc[frac_prolif_parameters_pso['index']==idx].iloc[0]['cost_total']<cost_limit : 
        frac_prolif_parameters_pso_temp = frac_prolif_parameters_pso.loc[frac_prolif_parameters_pso['index']==idx]
        parameters_temp = frac_prolif_parameters_pso_temp[['a_006','d_006','q_006']].iloc[0].as_matrix()
        
        axs[1, 1].plot(parameter_names, np.log10(parameters_temp), color=colormap(normalize(frac_prolif_parameters_pso.loc[frac_prolif_parameters_pso['index']==idx].iloc[0]['cost_total'])), lw=0.4)
#        axs[1, 1].grid(True)
        axs[1, 1].tick_params(axis="y", labelsize=24)
        axs[1, 1].tick_params(axis="x", labelsize=24)

#plotting original parameters and expt
axs[1, 1].plot(parameter_names,  np.log10([parameters_exp['0.6 kPa DMSO']['a'], parameters_exp['0.6 kPa DMSO']['d'], parameters_exp['0.6 kPa DMSO']['q']]), color='black', lw=6)
#axs[1, 1].plot(parameter_names, np.log10(frac_prolif_parameters_pso.loc[frac_prolif_parameters_pso['index']==1.0][['a_006','d_006','q_006']].iloc[0].as_matrix()), lw=5, color=colormap(normalize(frac_prolif_parameters_pso.loc[frac_prolif_parameters_pso['index']==1.0].iloc[0]['cost_total'])), path_effects=[mpe.withStroke(linewidth=6, foreground='black')])


axs[1, 1].set_ylabel('log(parameters)', fontsize = 24) 
axs[1, 1].set_xlabel('parameter from fractional proliferation model', fontsize = 24) 


############### 3rd subplot

x=np.linspace(0,120,400)

for idx in frac_prolif_parameters_pso[frac_prolif_parameters_pso['hyperparameters']==hyp+1]['index']:
    if frac_prolif_parameters_pso.loc[frac_prolif_parameters_pso['index']==idx].iloc[0]['cost_total']<cost_limit:
    #if frac_prolif_parameters_pso.loc[frac_prolif_parameters_pso['index']==idx].iloc[0]['cost_mmd']< inf :
        axs[0, 0].plot(x, g1_distribution_55[idx].tolist()
,color=colormap(normalize(frac_prolif_parameters_pso.loc[frac_prolif_parameters_pso['index']==idx].iloc[0]['cost_total'])), lw=0.4)
#        axs[0, 1].grid(True)
        axs[0, 0].tick_params(axis="y", labelsize=24)
        axs[0, 0].tick_params(axis="x", labelsize=24)
#plotting original parameters and expt
axs[0, 0].plot(x, g1_distribution_exp, color='black', lw=6)
#axs[0, 0].plot(x, g1_distribution_55[1].tolist(), lw=4,color=colormap(normalize(frac_prolif_parameters_pso.loc[frac_prolif_parameters_pso['index']==1.0].iloc[0]['cost_total'])), path_effects = [mpe.Stroke(linewidth=8, foreground='black')])
#axs[0,0].plot(x, g1_distribution_55[1].tolist(), lw=5, color='red', path_effects=[mpe.withStroke(linewidth=8, foreground='black')], label="stroke no alpha")
#axs[0,0].plot(x, g1_distribution_55[1].tolist(), lw=5, color=colormap(normalize(frac_prolif_parameters_pso.loc[frac_prolif_parameters_pso['index']==1.0].iloc[0]['cost_total'])), path_effects=[mpe.withStroke(linewidth=6, foreground='black')])



axs[0, 0].set_ylabel('Density', fontsize = 24) 
axs[0, 0].set_xlabel('G1 times', fontsize = 24) 

############### 4th subplot


x=np.linspace(0,120,400)

for idx in frac_prolif_parameters_pso[frac_prolif_parameters_pso['hyperparameters']==hyp+1]['index']:
    if frac_prolif_parameters_pso.loc[frac_prolif_parameters_pso['index']==idx].iloc[0]['cost_total']<cost_limit:
    #if frac_prolif_parameters_pso.loc[frac_prolif_parameters_pso['index']==idx].iloc[0]['cost_mmd']< inf :
        axs[1, 0].plot(x, sm_distribution_55[idx].tolist()
,color=colormap(normalize(frac_prolif_parameters_pso.loc[frac_prolif_parameters_pso['index']==idx].iloc[0]['cost_total'])), lw=0.4)
#        axs[0, 1].grid(True)
        axs[1, 0].tick_params(axis="y", labelsize=24)
        axs[1, 0].tick_params(axis="x", labelsize=24)
#plotting original parameters and expt
axs[1, 0].plot(x, sm_distribution_exp, color='black', lw=6)
#axs[1, 0].plot(x, sm_distribution_55[1].tolist(), lw=4, color=colormap(normalize(frac_prolif_parameters_pso.loc[frac_prolif_parameters_pso['index']==1.0].iloc[0]['cost_total'])), path_effects=[mpe.withStroke(linewidth=6, foreground='black')])


axs[1, 0].set_ylabel('Density', fontsize = 24) 
axs[1, 0].set_xlabel('S to M times', fontsize = 24) 



scalarmappaple = cm.ScalarMappable(norm=normalize, cmap=colormap)
scalarmappaple.set_array([])

#fig.subplots_adjust(right=0.8)
#cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
cb = fig.colorbar(scalarmappaple, ax=axs.ravel().tolist(), shrink=0.8)
cb.set_label(label='cost_total', fontsize=24)
cb.ax.tick_params(labelsize=12)

#pl.tight_layout()

#pl.tight_layout()

pl.savefig(fn, bbox_inches='tight')


#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################
parameters_log_best=[]
cost_limit = frac_prolif_parameters.sort_values(by=['cost_total'])[0:1]['cost_total'].tolist()[-1]


for idx in frac_prolif_parameters_pso[frac_prolif_parameters_pso['hyperparameters']==hyp+1]['index']:
    if frac_prolif_parameters_pso.loc[frac_prolif_parameters_pso['index']==idx].iloc[0]['cost_total']<cost_limit : 
        frac_prolif_parameters_pso_temp = frac_prolif_parameters_pso.loc[frac_prolif_parameters_pso['index']==idx]
        parameter = particle_positions_pso[particle_positions_pso[particle_positions_pso['hyperparameters']==hyp+1]['identification']==idx]['parameter_vector'].iloc[0]
        pos = np.where(position_history[hyp][:,:,0] == parameter)
        position_history[hyp][pos[0],pos[1],:]
        parameters_log_best += position_history[hyp][pos[0],pos[1],:].tolist()

parameters_log_best = np.array(parameters_log_best)
parameters_log_original = np.array(l2[0])

colors = sns.color_palette(n_colors=np.shape(parameters_log_best)[1])

a=3
b=10
fn = Path(r"C:\Users\Grim.M\Ecole\Ens\Master\S2\Stage\report\figures\pso\results_distribution_best_parameters_hyperparameters_5.pdf")
fig, axs = pl.subplots(ncols=a, nrows=b, sharex=False, sharey=True, figsize=(15,30))
for dim, ax in zip(range(30), axs.flat):
    sns.distplot(parameters_log_best[:, dim], color=colors[dim], ax=ax, kde_kws={"lw": 0.5, "label": "variable: " + pysb_sampled_parameter_names[dim]}, bins=np.linspace(loc[dim]-2,loc[dim]+2,50))
    ax.plot([parameters_log_original[dim] for i in range (0,8)],range (0,8), color = 'red', lw=2)
    ax.plot([parameters_log_original[dim]-2 for i in range (0,8)],range (0,8), '--', color = 'black', lw=1)
    ax.plot([parameters_log_original[dim]+2 for i in range (0,8)],range (0,8), '--', color = 'black', lw=1)
    ax.set_ylim([0, 7])
    ax.legend(loc='upper left', fontsize = 12)
    ax.set_xlabel('log(parameter)', fontsize = 12) 
    ax.set_ylabel('Density', fontsize = 12) 

for ax in axs.flat:
    ax.label_outer()
    
pl.savefig(fn, bbox_inches='tight')


pl.scatter(parameters_log_best[:, 0], parameters_log_best[:, 1])
