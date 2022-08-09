# -*- coding: utf-8 -*-
"""
Created on Thu Sep  2 17:25:15 2021

@author: Grim.M
"""


path=r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code"
import sys
sys.path.insert(0, path + "/cellular_cycle") #add a pathway to the model
from Austin_proliferation_model_concentration_5 import model
model_concentration=model

import math
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
import umap as umap
import sobol_seq
from peakfinder import peakfinder
from pathlib import Path
import matplotlib.gridspec as gridspec
from pathlib import Path
import matplotlib.pyplot as plt
from scipy import stats


########################################################################################
pysb_sampled_parameter_names = ["V_s2skp","V_Fp","k_M","k_EMY","k_CDF","k_CDS","k_ED","K_S","K_MY","K_DF","d_CD","k_P1", "k_P2", "K_CD", "k_E", "K_M", "K_CE", "V_2E2F", "K_2E2F", "k_RE", "kf1", "kf2", "kb", "kasbaxmbcl2", "kdsbaxmbcl2", "k1sbh3", "k2sbh3", "kdbh3", "kasbh3bcl2", "kdsbh3bcl2"]
loc=[]
for parameter in pysb_sampled_parameter_names:
    print((model.parameters[parameter].value))
    loc+=[np.log10(model.parameters[parameter].value)]

l=sobol_seq.i4_sobol_generate(30, 100000, skip=0)
l2 = 4*l + (np.array(loc)-2)
#######################################################################################

position_history = np.load(r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code\fit_data\pso\high_and_low_stiffness_cost_2\history\position_history.npy", allow_pickle=True)

#position_history = np.load("/home/jeremie/fit_data/pso/high_and_low_stiffness_cost/history/position_history.npy", allow_pickle=True)

particle_positions_pso = pd.read_pickle(r"C:/Users/Grim.M\Ecole\Ens\Master\S2\Stage\Code\fit_data\pso\high_and_low_stiffness_cost_2\particle_position/particle_position_all_parameters_pso.pkl")

frac_prolif_parameters_pso = pd.read_pickle(r"C:/Users/Grim.M/Ecole/ens/Master/S2/Stage/Code/fit_data/pso/high_and_low_stiffness_cost_2/cost_and_IMT_distribution/frac_prolif_parameters_pso.pkl")

#particle_positions_pso = pd.read_pickle("/home/jeremie/fit_data/pso/high_and_low_stiffness_cost/particle_position/particle_position_all_parameters_pso.pkl")



particle_positions_pso = pd.read_pickle(r"C:/Users/Grim.M/Ecole/ens/Master/S2/Stage/Code/fit_data/pso/high_and_low_stiffness_cost_2/particle_position/particle_position_all_parameters_pso.pkl")


frac_prolif_parameters = pd.read_pickle(path + r"/fit_data/sobol_seq/cost_function_model_5/cost_and_IMT_distribution/frac_prolif_parameters_all_parameters_006_55_kPa.pkl")


IMT_distribution = pd.read_pickle(r"C:/Users/Grim.M/Ecole/ens/Master/S2/Stage/Code/fit_data/pso/high_and_low_stiffness_cost_2/cost_and_IMT_distribution/IMT_distribution_all_parameters_50_initial_cells.pkl")

g1_distribution_55 = pd.read_pickle(r"C:/Users/Grim.M/Ecole/ens/Master/S2/Stage/Code/fit_data/pso/high_and_low_stiffness_cost_2/cost_and_IMT_distribution/g1_distribution_all_parameters_50_initial_cells_55.pkl")

sm_distribution_55 = pd.read_pickle(r"C:/Users/Grim.M/Ecole/ens/Master/S2/Stage/Code/fit_data/pso/high_and_low_stiffness_cost_2/cost_and_IMT_distribution/sm_distribution_all_parameters_50_initial_cells_55.pkl")


#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
### Figure 0 ###

parameters_log_best=[]
cost_limit = frac_prolif_parameters.sort_values(by=['cost_total'])[0:1]['cost_total'].tolist()[-1]
cost_limit = frac_prolif_parameters_pso.sort_values(by=['cost_total'])[0:100]['cost_total'].tolist()[-1]

n=0
for idx in frac_prolif_parameters_pso['index']:
    if frac_prolif_parameters_pso.loc[frac_prolif_parameters_pso['index']==idx].iloc[0]['cost_total'] < cost_limit : 
        n+=1
        print(n)
        frac_prolif_parameters_pso_temp = frac_prolif_parameters_pso.loc[frac_prolif_parameters_pso['index']==idx]
        parameter = particle_positions_pso[particle_positions_pso['identification']==idx]['parameter_vector'].iloc[0]
        pos = np.where(position_history[:,:,:] == parameter)
        position_history[pos[0][0],pos[1][0],:]
        print(len(position_history[pos[0][0],pos[1][0],:]))
        parameters_log_best += [position_history[pos[0][0],pos[1][0],:].tolist()]

parameters_log_best = np.array(parameters_log_best)
parameters_log_original = np.array(l2[0])

colors = sns.color_palette(n_colors = np.shape(parameters_log_best)[1])


a=3
b=10
fn = Path(r"C:\Users\Grim.M\Ecole\Ens\Master\S2\Stage\report\figures\pso\results_distribution_best_sobolseq_2.pdf")
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



#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#figure 1

#umap

unique_vectors=[]
unique_vectors += [np.concatenate(position_history)]

umap_1 = umap.UMAP().fit(unique_vectors[0])
X_umap = []
for particle in range(np.shape(position_history)[1]):   
    X_umap += [umap_1.transform(position_history[:,particle,:])]

##########################################

fig, axs = pl.subplots(nrows=1, ncols=5, figsize=(24,6), gridspec_kw={"width_ratios":[2, 2, 0.03, 0.03, 2]})

#### subplot 1 ####

colors = sns.color_palette(n_colors=np.shape(position_history)[1])
for chain in range(np.shape(position_history)[1]) :
    axs[0].scatter(X_umap[chain][:,0], X_umap[chain][:,1],c = colors[chain], alpha = 0.25,s = 0.05);
    axs[0].plot(X_umap[chain][:,0], X_umap[chain][:,1], c = colors[chain], alpha = 0.5, lw = 0.25, ls = '--');
    axs[0].plot(X_umap[chain][0,0], X_umap[chain][0,1],c=colors[chain], marker="^")
    axs[0].plot(X_umap[chain][-1,0], X_umap[chain][-1,1],c=colors[chain], marker='o')
    #axs[0, hyp].title.set_text(title)
    axs[0].set_title('A.',fontsize= 26, loc='left')
axs[0].set_xlabel('umap 1', fontsize = 16)
axs[0].set_ylabel('umap 2', fontsize = 16)


#### subplot 2 ####
vmin = min(frac_prolif_parameters_pso["cost_total"])
#vmax = max(frac_prolif_parameters_pso["cost_total"])
#vmax = frac_prolif_parameters.sort_values(by=['cost_total'])['cost_total'].tolist()[100]
vmax = 20
vmax = frac_prolif_parameters['cost_total'].loc[frac_prolif_parameters['index']==1][0]

colormap = cm.viridis
normalize = mcolors.Normalize(vmin = vmin, vmax = vmax)

for chain in range(np.shape(position_history)[1]) :
    axs[1].plot(X_umap[chain][:,0], X_umap[chain][:,1], c = 'grey', alpha = 0.5, lw = 0.25, ls = '--');
    axs[1].plot(X_umap[chain][0,0], X_umap[chain][0,1], c = 'grey', marker="^", alpha = 0.6)
    axs[1].plot(X_umap[chain][-1,0], X_umap[chain][-1,1], c = 'grey', marker='o', alpha = 0.6)
    a = np.reshape(np.concatenate(particle_positions_pso['parameter_vector'].to_numpy()),(len(particle_positions_pso),30))
    for position_ind in range(np.shape(position_history)[0]) :    
        try :
            identification_temp = particle_positions_pso[a==position_history[position_ind,chain,:].tolist()].iloc[0]['identification']
            frac_prolif_parameters_temp = frac_prolif_parameters_pso[frac_prolif_parameters_pso['index']==identification_temp]
            cost_temp = frac_prolif_parameters_temp.iloc[0]['cost_total']
            axs[1].scatter(X_umap[chain][position_ind,0], X_umap[chain][position_ind,1],c = colormap(normalize(cost_temp)) , alpha = 1, s = 5);
            print(position_ind)
        except :
            print(identification_temp)
            axs[1].scatter(X_umap[chain][position_ind,0], X_umap[chain][position_ind,1],c = 'white', alpha = 0.25, s = 5, edgecolor='black', linewidth=0.3);
axs[1].set_title('B.',fontsize= 26, loc='left')
scalarmappaple = cm.ScalarMappable(norm=normalize, cmap=colormap)
scalarmappaple.set_array([])
axs[1].set_xlabel('umap 1', fontsize = 16)
axs[1].set_ylabel('umap 2', fontsize = 16)

#CBar
cbar = fig.colorbar(scalarmappaple, cax=axs[2])
cbar.set_clim(vmin, vmax)
cbar.set_label('Predicted activity')
axs[3].axis('off')

#### subplot 3 ####
#Cost distribution
import pandas as pd

a = np.sort(frac_prolif_parameters_pso['cost_total'].tolist())         
#axs[3].bar(range(len(a)), a, color=colormap(normalize(cost_temp)))
frac_prolif_parameters_pso = frac_prolif_parameters_pso.sort_values(['cost_total'])
df = frac_prolif_parameters_pso[['index','cost_total']]
axs[4].bar(x='index', height='cost_total', data=df, color=colormap(normalize(df.cost_total.values)))
axs[4].set_xticks([])
axs[4].plot([frac_prolif_parameters['cost_total'].loc[frac_prolif_parameters['index']==1] for n in range(len(a))], '--', color = 'black', linewidth=2)
axs[4].plot([frac_prolif_parameters.sort_values(by=['cost_total'])['cost_total'].tolist()[0] for n in range(len(a))], color = 'black', linewidth=2)
axs[4].plot([a[99] for n in range(len(a))],':', color = 'black', linewidth=2)
axs[4].set_xlabel('Sets of parameters', fontsize = 16) 
axs[4].set_ylabel('Cost', fontsize = 16) 

pl.tight_layout()
pl.savefig(path + "/../report/figures/pso/all_particles_running_2_umap.pdf", bbox_inches='tight')



#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#figure 1'

a=3
b=10
fn = Path(r"C:\Users\Grim.M\Ecole\Ens\Master\S2\Stage\report\figures\pso\results_distribution_parameters_umap_2.pdf")
fig, axs = pl.subplots(ncols=a, nrows=b, sharex=True, sharey=True, figsize=(15,30))
for dim, ax in zip(range(30), axs.flat):
    ax.title.set_text(str(pysb_sampled_parameter_names[dim]))
    for chain in range(np.shape(position_history)[1]) :
        scatt_temp = ax.scatter(X_umap[chain][:,0], X_umap[chain][:,1], c = position_history[:,chain,dim], alpha = 1,s = 0.5, vmin=loc[dim]-2, vmax=loc[dim]+2)
        
#    cbar = plt.colorbar(scatt_temp, ax=ax, ticks=loc[dim] + range(-2,3))
#    cbar.ax.set_yticklabels(['log('+str(pysb_sampled_parameter_names[dim])+ ') - 2', 'log('+str(pysb_sampled_parameter_names[dim])+ ') - 1', 'log('+str(pysb_sampled_parameter_names[dim])+ ')', 'log('+str(pysb_sampled_parameter_names[dim])+ ') + 1', 'log('+str(pysb_sampled_parameter_names[dim])+ ') + 2'])
#    ax.set_xlabel('umap 1', fontsize = 16)
#    ax.set_ylabel('umap 2', fontsize = 16)
#pl.tight_layout()
fig.text(0.5, 0.04, 'umap 1', ha='center')
fig.text(0.04, 0.5, 'umap 2', va='center', rotation='vertical')
for ax in axs.flat:
    ax.label_outer()


#pl.tight_layout()

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
cbar = fig.colorbar(scatt_temp, cax=cbar_ax, ticks=loc[dim] + range(-2,3))
cbar.ax.set_yticklabels(['log(original parameter) - 2', 'log(original parameter) - 1', 'log(original parameter)', 'log(original parameter) + 1', 'log(original parameter) + 2'], rotation=-45)

#pl.tight_layout()


pl.savefig(fn, bbox_inches='tight')

#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
####################################################################################
#figure 2
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
#################################################################################

from pathlib import Path

fn = Path(path + "/../report/figures/pso/results_model_5.pdf")

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
colormap = cm.viridis
normalize = mcolors.Normalize(vmin = minimum, vmax = cost_limit)
#or 

for idx in frac_prolif_parameters_pso['index']:
    print(idx)
    if frac_prolif_parameters_pso.loc[frac_prolif_parameters_pso['index']==idx].iloc[0]['cost_total']<cost_limit :
        axs[0, 1].plot(x, IMT_distribution[idx].tolist(), color=colormap(normalize(frac_prolif_parameters_pso.loc[frac_prolif_parameters_pso['index']==idx].iloc[0]['cost_total'])), lw=0.4)
        axs[0, 1].tick_params(axis="y", labelsize=24)
        axs[0, 1].tick_params(axis="x", labelsize=24)
#plotting original parameters and expt
axs[0, 1].plot(x, experimental_IMT_distribution, color='black', lw=6)

axs[0, 1].set_ylabel('Density', fontsize = 24) 
axs[0, 1].set_xlabel('Intermitotic times', fontsize = 24) 

############### 2nd subplot

for idx in frac_prolif_parameters_pso['index']:
    if frac_prolif_parameters_pso.loc[frac_prolif_parameters_pso['index']==idx].iloc[0]['cost_total']<cost_limit : 
        frac_prolif_parameters_pso_temp = frac_prolif_parameters_pso.loc[frac_prolif_parameters_pso['index']==idx]
        parameters_temp = frac_prolif_parameters_pso_temp[['a_006','d_006','q_006']].iloc[0].to_numpy() 
        axs[1, 1].plot(parameter_names, np.log10(parameters_temp), color=colormap(normalize(frac_prolif_parameters_pso.loc[frac_prolif_parameters_pso['index']==idx].iloc[0]['cost_total'])), lw=0.4)
        axs[1, 1].tick_params(axis="y", labelsize=24)
        axs[1, 1].tick_params(axis="x", labelsize=24)
#plotting original parameters and expt
axs[1, 1].plot(parameter_names,  np.log10([parameters_exp['0.6 kPa DMSO']['a'], parameters_exp['0.6 kPa DMSO']['d'], parameters_exp['0.6 kPa DMSO']['q']]), color='black', lw=6)

axs[1, 1].set_ylabel('log(parameters)', fontsize = 24) 
axs[1, 1].set_xlabel('parameter from fractional proliferation model', fontsize = 24) 


############### 3rd subplot

x=np.linspace(0,120,400)
for idx in frac_prolif_parameters_pso['index']:
    if frac_prolif_parameters_pso.loc[frac_prolif_parameters_pso['index']==idx].iloc[0]['cost_total']<cost_limit:
        axs[0, 0].plot(x, g1_distribution_55[idx].tolist()
,color=colormap(normalize(frac_prolif_parameters_pso.loc[frac_prolif_parameters_pso['index']==idx].iloc[0]['cost_total'])), lw=0.4)
        axs[0, 0].tick_params(axis="y", labelsize=24)
        axs[0, 0].tick_params(axis="x", labelsize=24)
#plotting original parameters and expt
axs[0, 0].plot(x, g1_distribution_exp, color='black', lw=6)

axs[0, 0].set_ylabel('Density', fontsize = 24) 
axs[0, 0].set_xlabel('G1 times', fontsize = 24) 

############### 4th subplot

x=np.linspace(0,120,400)
for idx in frac_prolif_parameters_pso['index']:
    if frac_prolif_parameters_pso.loc[frac_prolif_parameters_pso['index']==idx].iloc[0]['cost_total']<cost_limit:
        axs[1, 0].plot(x, sm_distribution_55[idx].tolist()
,color=colormap(normalize(frac_prolif_parameters_pso.loc[frac_prolif_parameters_pso['index']==idx].iloc[0]['cost_total'])), lw=0.4)
        axs[1, 0].tick_params(axis="y", labelsize=24)
        axs[1, 0].tick_params(axis="x", labelsize=24)
#plotting original parameters and expt
axs[1, 0].plot(x, sm_distribution_exp, color='black', lw=6)
axs[1, 0].set_ylabel('Density', fontsize = 24) 
axs[1, 0].set_xlabel('S to M times', fontsize = 24) 


scalarmappaple = cm.ScalarMappable(norm=normalize, cmap=colormap)
scalarmappaple.set_array([])
cb = fig.colorbar(scalarmappaple, ax=axs.ravel().tolist(), shrink=0.8)
cb.set_label(label='cost_total', fontsize=24)
cb.ax.tick_params(labelsize=12)

pl.savefig(fn, bbox_inches='tight')

#############################################################################
#############################################################################
#############################################################################
#############################################################################
#figure 3
"""

all_fft=[]
indices=[]

for index in frac_prolif_parameters_pso['index']:
        
    fft=[]
    obs = ['obsCycE', 'obsCdc20']
    try :
        for j in range(2):
            
            fft+= np.load(r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code\fit_data\pso\high_and_low_stiffness_cost_2\fft\stochastic_simulation\fft_stochastic_simulation_" + obs[j] + "_50_to_800_hours_index_" + str(index)+"_55_kPa.npy").tolist()
            
            fft+= np.load(r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code\fit_data\pso\high_and_low_stiffness_cost_2\fft\stochastic_simulation\fft_stochastic_simulation_" + obs[j] + "_50_to_800_hours_index_" + str(index)+"_006_kPa.npy").tolist()
        
        all_fft+=[fft]
        indices+=[index]
        print('ok')
    except:
        print(index)
        
all_fft = np.array(all_fft)    
indices = np.array(indices)    
np.shape(all_fft)

indices = [int(numeric_string) for numeric_string in indices]

np.save(r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code\fit_data\pso\high_and_low_stiffness_cost_2\fft\stochastic_simulation\all_fft.npy", all_fft)
np.save(r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code\fit_data\pso\high_and_low_stiffness_cost_2\fft\stochastic_simulation\indices.npy", indices)
"""

all_fft = np.load(r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code\fit_data\pso\high_and_low_stiffness_cost_2\fft\stochastic_simulation\all_fft.npy")
indices = np.load(r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code\fit_data\pso\high_and_low_stiffness_cost_2\fft\stochastic_simulation\indices.npy")

###
umap_fft = umap.UMAP().fit(all_fft)
#trans = umap.UMAP().fit(l2)
X_umap_fft = umap_fft.transform(all_fft)
pl.scatter(X_umap_fft[:,0], X_umap_fft[:,1], s=0.8)
from sklearn.cluster import KMeans
import numpy as np

n_clusters=10
kmeans = KMeans(n_clusters, random_state=0).fit(all_fft)
kmeans.labels_
kmeans.cluster_centers_

#umap parameter space#

unique_vectors=[]
unique_vectors += [np.concatenate(position_history)]
umap_1 = umap.UMAP().fit(unique_vectors[0])
X_umap = []
for particle in range(np.shape(position_history)[1]):   
    X_umap += [umap_1.transform(position_history[:,particle,:])]


#############################################################################
#plot#

fig, axs = pl.subplots(nrows=1, ncols=2, figsize=(12,5))
colors = [plt.cm.Spectral(each) for each in np.linspace(0, 1, n_clusters)]


axs[0].title.set_text('Umap on spectral data')
for chain in range(np.shape(position_history)[1]) :
    axs[0].plot(X_umap[chain][:,0], X_umap[chain][:,1], c = 'grey', alpha = 0.5, lw = 0.25, ls = '--');
    axs[0].plot(X_umap[chain][0,0], X_umap[chain][0,1], c = 'grey', marker="^", alpha = 0.6)
    axs[0].plot(X_umap[chain][-1,0], X_umap[chain][-1,1], c = 'grey', marker='o', alpha = 0.6)
    for position_ind in range(np.shape(position_history)[0]) :    
        try :
            a = np.reshape(np.concatenate(particle_positions_pso['parameter_vector'].to_numpy()),(len(particle_positions_pso),30))
            identification_temp = particle_positions_pso[a==position_history[position_ind,chain,:].tolist()].iloc[0]['identification']
            frac_prolif_parameters_temp = frac_prolif_parameters_pso[frac_prolif_parameters_pso['index']==identification_temp]
            cost_temp = frac_prolif_parameters_temp.iloc[0]['cost_total']
            axs[0].scatter(X_umap[chain][position_ind,0], X_umap[chain][position_ind,1], c = colors[kmeans.labels_[indices==int(identification_temp)][0]] , alpha = 1, s = 5);
            print(kmeans.labels_[np.array(indices)==int(identification_temp)])

        except :
            print(identification_temp)
            axs[0].scatter(X_umap[chain][position_ind,0], X_umap[chain][position_ind,1],c = 'white', alpha = 0.25, s = 5, edgecolor='black', linewidth=0.3);


axs[1].title.set_text('Umap on spectral data')

for k in range(len(X_umap_fft)):
    # Black used for noise.
    axs[1].scatter(X_umap_fft[k,0], X_umap_fft[k,1], s=0.8, c = colors[kmeans.labels_[k]])

pl.savefig(path + "/../report/figures/pso/all_particles_running_umap_all_hyperparameters_clustering_fft_kmeans_" + str(n_clusters) + "_clusters_model_5.pdf", bbox_inches='tight')


######################################################################################
######################################################################################

colors = [plt.cm.Spectral(each) for each in np.linspace(0, 1, n_clusters)]


cost_clusters = [[] for c in range(n_clusters)]
idx_clusters = [[] for c in range(n_clusters)]

for chain in range(np.shape(position_history)[1]) :
    for position_ind in range(np.shape(position_history)[0]) :    
        try :
            a = np.reshape(np.concatenate(particle_positions_pso['parameter_vector'].to_numpy()),(len(particle_positions_pso),30))
            idx = particle_positions_pso[a==position_history[position_ind,chain,:].tolist()].iloc[0]['identification']
            
#           if frac_prolif_parameters_pso.loc[frac_prolif_parameters_pso['index']==idx].iloc[0]['cost_total']<cost_limit : #select the cost  lower than cost limit
            
            kmeans.labels_[indices==int(idx)][0]
            cost_clusters[kmeans.labels_[np.array(indices)==int(idx)][0]] += [frac_prolif_parameters_pso.loc[frac_prolif_parameters_pso['index']==idx].iloc[0]['cost_total']]
            idx_clusters[kmeans.labels_[np.array(indices)==int(idx)][0]] += [idx]
        except :
            print(idx)

mean_clusters=[]
stdv_clusters=[]
for clust_i in range(n_clusters) :
    mean_clusters += [np.mean(cost_clusters[clust_i])]
    stdv_clusters += [np.std(cost_clusters[clust_i])]
fig, ax = plt.subplots(1)
x=range(n_clusters)
for pos, y, stdv, color in zip(x, mean_clusters, 
                               stdv_clusters, colors):
    ax.errorbar(pos, y, stdv, lw = 2,
                marker="o", linestyle="none", 
                color = color)
plt.sca(ax)
plt.xticks(x, range(n_clusters), fontsize=16)
plt.xlabel('Cluster')
#plt.xlim(right=3) 

plt.ylabel('Cost total', fontsize=16)

plt.savefig(path + "/../report/figures/pso/Cost_comparison_between_clusters_" + str(n_clusters) + "_clusters_model_5.pdf", bbox_inches='tight')




#Welch's t-test

stats.ttest_ind(cost_clusters[0], cost_clusters[1], equal_var=False)
stats.ttest_ind(cost_clusters[0], cost_clusters[2], equal_var=False)
stats.ttest_ind(cost_clusters[1], cost_clusters[2], equal_var=False)

(stdv_clusters[0]/stdv_clusters[1])**2
(stdv_clusters[0]/stdv_clusters[2])**2
(stdv_clusters[2]/stdv_clusters[1])**2


#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#figure 4


fn = Path(path + "/../report/figures/pso/results_model_5_" + str(n_clusters) + "_clusters_nolim.pdf")

fig, axs = pl.subplots(ncols=2, nrows=2, sharex=False, sharey=False, figsize=(40,20), constrained_layout = False)

############### 1st subplot
minimum = frac_prolif_parameters_pso.sort_values(by=['cost_total'])[0:50]['cost_total'].tolist()[2]
#cost_limit = frac_prolif_parameters_pso.sort_values(by=['cost_total'])[0:60]['cost_total'].tolist()[-1]
cost_limit = frac_prolif_parameters.sort_values(by=['cost_total'])[0:1]['cost_total'].tolist()[-1]
cost_limit = math.inf


for ax, col in zip(axs[0], ['55 kPa DMSO','0.6 kPa DMSO']):
    ax.annotate(col, xy=(0.5, 1), xytext=(0, 5),
                xycoords='axes fraction', textcoords='offset points',
                fontsize=36, ha='center', va='baseline')

x=np.linspace(0,120,360)
colormap = cm.RdYlBu
normalize = mcolors.Normalize(vmin = minimum, vmax = cost_limit)
#or 
for cluster_i in range(n_clusters):
    for idx in idx_clusters[cluster_i]:
        print(idx)
        #idx = int(idx)
        if frac_prolif_parameters_pso.loc[frac_prolif_parameters_pso['index']==idx].iloc[0]['cost_total']<cost_limit  :

            axs[0, 1].plot(x, IMT_distribution[idx].tolist(), color=colors[cluster_i], lw=0.4)
    #        axs[0, 1].grid(True)
            axs[0, 1].tick_params(axis="y", labelsize=24)
            axs[0, 1].tick_params(axis="x", labelsize=24)

axs[0, 1].plot(x, experimental_IMT_distribution, color='black', lw=6)

axs[0, 1].set_ylabel('Density', fontsize = 24) 
axs[0, 1].set_xlabel('Intermitotic times', fontsize = 24) 

############### 2nd subplot

for cluster_i in range(n_clusters):
    for idx in idx_clusters[cluster_i]:
        if frac_prolif_parameters_pso.loc[frac_prolif_parameters_pso['index']==idx].iloc[0]['cost_total']<cost_limit : 
            frac_prolif_parameters_pso_temp = frac_prolif_parameters_pso.loc[frac_prolif_parameters_pso['index']==idx]
            parameters_temp = frac_prolif_parameters_pso_temp[['a_006','d_006','q_006']].iloc[0].to_numpy() 
            
            axs[1, 1].plot(parameter_names, np.log10(parameters_temp), color=colors[cluster_i], lw=0.4)
            axs[1, 1].tick_params(axis="y", labelsize=24)
            axs[1, 1].tick_params(axis="x", labelsize=24)

axs[1, 1].plot(parameter_names,  np.log10([parameters_exp['0.6 kPa DMSO']['a'], parameters_exp['0.6 kPa DMSO']['d'], parameters_exp['0.6 kPa DMSO']['q']]), color='black', lw=6)


axs[1, 1].set_ylabel('log(parameters)', fontsize = 24) 
axs[1, 1].set_xlabel('parameter from fractional proliferation model', fontsize = 24) 


############### 3rd subplot

x=np.linspace(0,120,400)

for cluster_i in range(n_clusters):
    for idx in idx_clusters[cluster_i]:
        if frac_prolif_parameters_pso.loc[frac_prolif_parameters_pso['index']==idx].iloc[0]['cost_total']<cost_limit:
            axs[0, 0].plot(x, g1_distribution_55[idx].tolist()
    , color=colors[cluster_i], lw=0.4)
            axs[0, 0].tick_params(axis="y", labelsize=24)
            axs[0, 0].tick_params(axis="x", labelsize=24)

axs[0, 0].set_ylabel('Density', fontsize = 24) 
axs[0, 0].set_xlabel('G1 times', fontsize = 24) 

axs[0, 0].plot(x, g1_distribution_exp, color='black', lw=6)


############### 4th subplot

x=np.linspace(0,120,400)
for cluster_i in range(n_clusters):
    for idx in idx_clusters[cluster_i]:
        if frac_prolif_parameters_pso.loc[frac_prolif_parameters_pso['index']==idx].iloc[0]['cost_total']<cost_limit:
            axs[1, 0].plot(x, sm_distribution_55[idx].tolist(), color=colors[cluster_i], lw=0.4)
            axs[1, 0].tick_params(axis="y", labelsize=24)
            axs[1, 0].tick_params(axis="x", labelsize=24)

axs[1, 0].plot(x, sm_distribution_exp, color='black', lw=6)

axs[1, 0].set_ylabel('Density', fontsize = 24) 
axs[1, 0].set_xlabel('S to M times', fontsize = 24) 


#pl.tight_layout()

pl.savefig(fn, bbox_inches='tight')

#################################################################################
#################################################################################
#################################################################################
#################################################################################
#figure 4

for cluster_i in range(n_clusters):
    fn = Path(path + "/../report/figures/pso/results_model_5_cluster_n_"+ str(cluster_i) + "_cluster_number_" + str(n_clusters) + "_clusters_nolim.pdf")
    
    fig, axs = pl.subplots(ncols=2, nrows=2, sharex=False, sharey=False, figsize=(40,20), constrained_layout = False)
    
    ############### 1st subplot
    minimum = frac_prolif_parameters_pso.sort_values(by=['cost_total'])[0:50]['cost_total'].tolist()[2]
    #cost_limit = frac_prolif_parameters_pso.sort_values(by=['cost_total'])[0:60]['cost_total'].tolist()[-1]
    cost_limit = frac_prolif_parameters.sort_values(by=['cost_total'])[0:1]['cost_total'].tolist()[-1]
    cost_limit = math.inf
    
    
    for ax, col in zip(axs[0], ['55 kPa DMSO','0.6 kPa DMSO']):
        ax.annotate(col, xy=(0.5, 1), xytext=(0, 5),
                    xycoords='axes fraction', textcoords='offset points',
                    fontsize=36, ha='center', va='baseline')
    
    x=np.linspace(0,120,360)
    colormap = cm.RdYlBu
    normalize = mcolors.Normalize(vmin = minimum, vmax = cost_limit)
    #or 
    #for cluster_i in range(n_clusters):
    for idx in idx_clusters[cluster_i]:
        print(idx)
        #idx = int(idx)
        if frac_prolif_parameters_pso.loc[frac_prolif_parameters_pso['index']==idx].iloc[0]['cost_total']<cost_limit  :

            axs[0, 1].plot(x, IMT_distribution[idx].tolist(), color=colors[cluster_i], lw=0.4)
    #        axs[0, 1].grid(True)
            axs[0, 1].tick_params(axis="y", labelsize=24)
            axs[0, 1].tick_params(axis="x", labelsize=24)
    
    axs[0, 1].plot(x, experimental_IMT_distribution, color='black', lw=6)
    
    axs[0, 1].set_ylabel('Density', fontsize = 24) 
    axs[0, 1].set_xlabel('Intermitotic times', fontsize = 24) 
    
    ############### 2nd subplot
    
    #for cluster_i in range(n_clusters):
    for idx in idx_clusters[cluster_i]:
        if frac_prolif_parameters_pso.loc[frac_prolif_parameters_pso['index']==idx].iloc[0]['cost_total']<cost_limit : 
            frac_prolif_parameters_pso_temp = frac_prolif_parameters_pso.loc[frac_prolif_parameters_pso['index']==idx]
            parameters_temp = frac_prolif_parameters_pso_temp[['a_006','d_006','q_006']].iloc[0].to_numpy() 
            
            axs[1, 1].plot(parameter_names, np.log10(parameters_temp), color=colors[cluster_i], lw=0.4)
            axs[1, 1].tick_params(axis="y", labelsize=24)
            axs[1, 1].tick_params(axis="x", labelsize=24)

    axs[1, 1].plot(parameter_names,  np.log10([parameters_exp['0.6 kPa DMSO']['a'], parameters_exp['0.6 kPa DMSO']['d'], parameters_exp['0.6 kPa DMSO']['q']]), color='black', lw=6)
    
    
    axs[1, 1].set_ylabel('log(parameters)', fontsize = 24) 
    axs[1, 1].set_xlabel('parameter from fractional proliferation model', fontsize = 24) 
    
    
    ############### 3rd subplot
    
    x=np.linspace(0,120,400)
    
    #for cluster_i in range(n_clusters):
    for idx in idx_clusters[cluster_i]:
        if frac_prolif_parameters_pso.loc[frac_prolif_parameters_pso['index']==idx].iloc[0]['cost_total']<cost_limit:
            axs[0, 0].plot(x, g1_distribution_55[idx].tolist()
    , color=colors[cluster_i], lw=0.4)
            axs[0, 0].tick_params(axis="y", labelsize=24)
            axs[0, 0].tick_params(axis="x", labelsize=24)
    
    axs[0, 0].set_ylabel('Density', fontsize = 24) 
    axs[0, 0].set_xlabel('G1 times', fontsize = 24) 
    
    axs[0, 0].plot(x, g1_distribution_exp, color='black', lw=6)
    
    
    ############### 4th subplot
    
    x=np.linspace(0,120,400)
    #for cluster_i in range(n_clusters):
    for idx in idx_clusters[cluster_i]:
        if frac_prolif_parameters_pso.loc[frac_prolif_parameters_pso['index']==idx].iloc[0]['cost_total']<cost_limit:
            axs[1, 0].plot(x, sm_distribution_55[idx].tolist(), color=colors[cluster_i], lw=0.4)
            axs[1, 0].tick_params(axis="y", labelsize=24)
            axs[1, 0].tick_params(axis="x", labelsize=24)
    
    axs[1, 0].plot(x, sm_distribution_exp, color='black', lw=6)
    
    axs[1, 0].set_ylabel('Density', fontsize = 24) 
    axs[1, 0].set_xlabel('S to M times', fontsize = 24) 
    
    
    #pl.tight_layout()
    
    pl.savefig(fn, bbox_inches='tight')

#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################



fn = Path(r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\report\figures\pso\fft_kmeans_model_5_" + str(n_clusters) + "_clusters.pdf").expanduser()
t = np.linspace(0, 800, 800)

obs = ['obsCycE', 'obsCdc20']
stiffness = ['0.6 kPa', '55 kPa']

fig = plt.figure(figsize=(30, 16))
outer = gridspec.GridSpec(n_clusters, 1, wspace=1.8, hspace=0.3)

for i in range(n_clusters):
    np.array(indices)[kmeans.labels_==i]
    index = random.choice(np.array(indices)[kmeans.labels_==i])
    index = np.array(indices)[kmeans.labels_==i][50]
    ssa={}
    ssa['0.6 kPa']={}
    ssa['55 kPa']={}
    fft={}
    fft['0.6 kPa']={}
    fft['55 kPa']={}
    ssa['0.6 kPa']["obsCycE"] = np.load(path + r"\fit_data\pso\high_and_low_stiffness_cost_2\biochemical_data\ssa_timecourse_CycE_index_" + str(index).replace('0.','') + "_006_kPa_DMSO.npy")
    ssa['0.6 kPa']["obsCdc20"] =  np.load(path + r"\fit_data\pso\high_and_low_stiffness_cost_2\biochemical_data\ssa_timecourse_Cdc20_index_" + str(index).replace('0.','') + "_006_kPa_DMSO.npy")
    ssa['55 kPa']["obsCycE"] = np.load(path + r"\fit_data\pso\high_and_low_stiffness_cost_2\biochemical_data\ssa_timecourse_CycE_index_" + str(index).replace('0.','') + "_55_kPa_DMSO.npy")
    ssa['55 kPa']["obsCdc20"] =  np.load(path + r"\fit_data\pso\high_and_low_stiffness_cost_2\biochemical_data\ssa_timecourse_Cdc20_index_" + str(index).replace('0.','') + "_55_kPa_DMSO.npy")
    fft['55 kPa']["obsCdc20"] = np.load(r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code\fit_data\pso\high_and_low_stiffness_cost_2\fft\stochastic_simulation\fft_stochastic_simulation_obsCdc20_50_to_800_hours_index_" + str(index).replace('0.','')+"_55_kPa.npy").tolist()   
    fft['55 kPa']["obsCycE"] = np.load(r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code\fit_data\pso\high_and_low_stiffness_cost_2\fft\stochastic_simulation\fft_stochastic_simulation_obsCycE_50_to_800_hours_index_" + str(index).replace('0.','')+"_55_kPa.npy").tolist()
    fft['0.6 kPa']["obsCdc20"] = np.load(r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code\fit_data\pso\high_and_low_stiffness_cost_2\fft\stochastic_simulation\fft_stochastic_simulation_obsCdc20_50_to_800_hours_index_" + str(index).replace('0.','')+"_006_kPa.npy").tolist()
    fft['0.6 kPa']["obsCycE"] = np.load(r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code\fit_data\pso\high_and_low_stiffness_cost_2\fft\stochastic_simulation\fft_stochastic_simulation_obsCycE_50_to_800_hours_index_" + str(index).replace('0.','')+"_006_kPa.npy").tolist()
    inner = gridspec.GridSpecFromSubplotSpec(2, 4, subplot_spec=outer[i], wspace=0.1, hspace=0.1)
    T = t[1] - t[0]  # sampling interval 
    N = ssa['55 kPa']['obsCycE'][50:].size
    f = np.linspace(0, 1 / T, N) # 1/T = frequency
    for j in range(2):
        for s in range(2):
            ax = plt.Subplot(fig, inner[j + 2*s])
            ax.plot(t[400:520], ssa[stiffness[s]][obs[j]][400:520], linewidth=0.5, color=colors[i])
            ax.set_title(obs[j] + ' ' + stiffness[s])
            fig.add_subplot(ax)
            ax = plt.Subplot(fig, inner[4+j + 2*s])
            ax.plot(f, fft[stiffness[s]][obs[j]], linewidth=0.5, color=colors[i])
            fig.add_subplot(ax)

plt.savefig(fn, bbox_inches='tight')



######################################################################################
######################################################################################


#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################

idx_temp = frac_prolif_parameters_pso.sort_values(by=['cost_total'])[0:1]['index'].tolist()[0]
particle_temp = frac_prolif_parameters_pso.sort_values(by=['cost_total'])[0:1]['particle'].tolist()[0]
IMT_best_pso = pd.read_pickle(r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code\fit_data\pso\high_and_low_stiffness_cost_2\all_events_simulated\all_events_simulated_particle_" + particle_temp + "_0." + idx + "_50_initial_cells.pkl")
IMT_best_pso = IMT_best_pso[IMT_best_pso["event"]=='mitosis']["Lifespan"].tolist()

G1_best_pso = np.load(r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code\fit_data\pso\high_and_low_stiffness_cost_2\cost_and_IMT_distribution\g1_distribution_particle_" + particle_temp + "_0." + idx + "_30_initial_cells_55kPa_DMSO.npy")
SM_best_pso = np.load(r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code\fit_data\pso\high_and_low_stiffness_cost_2\cost_and_IMT_distribution\sm_distribution_particle_" + particle_temp + "_0." + idx + "_30_initial_cells_55kPa_DMSO.npy")


idx_temp = int(frac_prolif_parameters.sort_values(by=['cost_total'])[0:1]['index'].tolist()[0])
IMT_best_sobol = pd.read_pickle(r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code\fit_data/sobol_seq/cost_function_model_5/all_events_simulated/all_events_simulated_" + str(idx_temp) + "_50_initial_cells.pkl")
IMT_best_sobol = IMT_best_sobol[IMT_best_sobol["event"]=='mitosis']["Lifespan"].tolist()
np.mean(IMT_best_sobol)
np.var(IMT_best_sobol)


G1_best_sobol = np.load(r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code\fit_data/sobol_seq/cost_function_model_5/cost_and_IMT_distribution\g1_distribution_" + str(idx_temp) + "_30_initial_cells_55kPa_DMSO.npy")
SM_best_sobol = np.load(r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code\fit_data/sobol_seq/cost_function_model_5/cost_and_IMT_distribution\sm_distribution_" + str(idx_temp) + "_30_initial_cells_55kPa_DMSO.npy")


experimental_G1s['55 kPa DMSO']
experimental_SMs['55 kPa DMSO']




fn = Path(path + "/../report/figures/pso/results_mean_var_params_times_model_5_" + str(n_clusters) + "_clusters_nolim.pdf")

fig, axs = pl.subplots(ncols=2, nrows=2, sharex=False, sharey=False, figsize=(20,10), constrained_layout = False)

############### 1st subplot

means=[np.mean(IMT_best_pso), np.mean(IMT_best_sobol), np.mean(mitotic_lifespans_exp)]
stdvs=[np.var(IMT_best_pso), np.var(IMT_best_sobol), np.var(mitotic_lifespans_exp)]


for pos, y, stdv in zip(x, means, stdvs):
    axs[0, 1].errorbar(pos, y, stdv, lw = 2, marker="o", linestyle="none")
plt.sca(axs[0, 1])
plt.xticks(range(3),['best pso', 'best sobol', 'exp'], fontsize=16)
axs[0, 1].set_ylabel('Intermitotic times', fontsize = 16) 

############### 2nd subplot


means=[np.mean(G1_best_pso), np.mean(G1_best_sobol), np.mean(experimental_G1s['55 kPa DMSO'])]
stdvs=[np.var(G1_best_pso), np.var(G1_best_sobol), np.var(experimental_G1s['55 kPa DMSO'])]


for pos, y, stdv in zip(x, means, stdvs):
    axs[0, 0].errorbar(pos, y, stdv, lw = 2, marker="o", linestyle="none")
plt.sca(axs[0, 0])
plt.xticks(range(3),['best pso', 'best sobol', 'exp'], fontsize=16)
axs[0, 0].set_ylabel('G1 times', fontsize = 16) 
axs[0, 0].plot(x, g1_distribution_exp, color='black', lw=6)


############### 3rd subplot


means=[np.mean(SM_best_pso), np.mean(SM_best_sobol), np.mean(experimental_SMs['55 kPa DMSO'])]
stdvs=[np.var(SM_best_pso), np.var(SM_best_sobol), np.var(experimental_SMs['55 kPa DMSO'])]


for pos, y, stdv in zip(x, means, stdvs):
    axs[1, 0].errorbar(pos, y, stdv, lw = 2, marker="o", linestyle="none")
plt.sca(axs[1, 0])
plt.xticks(range(3),['best pso', 'best sobol', 'exp'], fontsize=16)
axs[1, 0].set_ylabel('SM times', fontsize = 16) 


#pl.tight_layout()

pl.savefig(fn, bbox_inches='tight')







fn = Path(path + "/../report/figures/pso/results_mean_params_times_model_5_" + str(n_clusters) + "_clusters_nolim.pdf")

fig, axs = pl.subplots(ncols=2, nrows=2, sharex=False, sharey=False, figsize=(20,10), constrained_layout = False)

############### 1st subplot

means=[np.mean(IMT_best_pso), np.mean(IMT_best_sobol), np.mean(mitotic_lifespans_exp)]
stdvs=[np.var(IMT_best_pso), np.var(IMT_best_sobol), np.var(mitotic_lifespans_exp)]


axs[0, 1].scatter(range(3), means)
plt.sca(axs[0, 1])
plt.xticks(range(3),['best pso', 'best sobol', 'exp'], fontsize=16)
axs[0, 1].set_ylabel('Intermitotic times', fontsize = 16) 

############### 2nd subplot


means=[np.mean(G1_best_pso), np.mean(G1_best_sobol), np.mean(experimental_G1s['55 kPa DMSO'])]
stdvs=[np.var(G1_best_pso), np.var(G1_best_sobol), np.var(experimental_G1s['55 kPa DMSO'])]

axs[0, 0].scatter(range(3), means)
plt.sca(axs[0, 0])
plt.xticks(range(3),['best pso', 'best sobol', 'exp'], fontsize=16)
axs[0, 0].set_ylabel('G1 times', fontsize = 16) 


############### 3rd subplot


means=[np.mean(SM_best_pso), np.mean(SM_best_sobol), np.mean(experimental_SMs['55 kPa DMSO'])]
stdvs=[np.var(SM_best_pso), np.var(SM_best_sobol), np.var(experimental_SMs['55 kPa DMSO'])]

axs[1, 0].scatter(range(3), means)
plt.sca(axs[1, 0])
plt.xticks(range(3),['best pso', 'best sobol', 'exp'], fontsize=16)
axs[1, 0].set_ylabel('SM times', fontsize = 16) 


#pl.tight_layout()

pl.savefig(fn, bbox_inches='tight')




fn = Path(path + "/../report/figures/pso/results_var_params_times_model_5_" + str(n_clusters) + "_clusters_nolim.pdf")

fig, axs = pl.subplots(ncols=2, nrows=2, sharex=False, sharey=False, figsize=(20,10), constrained_layout = False)

############### 1st subplot

means=[np.mean(IMT_best_pso), np.mean(IMT_best_sobol), np.mean(mitotic_lifespans_exp)]
stdvs=[np.var(IMT_best_pso), np.var(IMT_best_sobol), np.var(mitotic_lifespans_exp)]


axs[0, 1].scatter(range(3), stdvs)
plt.sca(axs[0, 1])
plt.xticks(range(3),['best pso', 'best sobol', 'exp'], fontsize=16)
axs[0, 1].set_ylabel('Intermitotic times', fontsize = 16) 

############### 2nd subplot


means=[np.mean(G1_best_pso), np.mean(G1_best_sobol), np.mean(experimental_G1s['55 kPa DMSO'])]
stdvs=[np.var(G1_best_pso), np.var(G1_best_sobol), np.var(experimental_G1s['55 kPa DMSO'])]

axs[0, 0].scatter(range(3), stdvs)
plt.sca(axs[0, 0])
plt.xticks(range(3),['best pso', 'best sobol', 'exp'], fontsize=16)
axs[0, 0].set_ylabel('G1 times', fontsize = 16) 


############### 3rd subplot


means=[np.mean(SM_best_pso), np.mean(SM_best_sobol), np.mean(experimental_SMs['55 kPa DMSO'])]
stdvs=[np.var(SM_best_pso), np.var(SM_best_sobol), np.var(experimental_SMs['55 kPa DMSO'])]

axs[1, 0].scatter(range(3), stdvs)
plt.sca(axs[1, 0])
plt.xticks(range(3),['best pso', 'best sobol', 'exp'], fontsize=16)
axs[1, 0].set_ylabel('SM times', fontsize = 16) 


#pl.tight_layout()

pl.savefig(fn, bbox_inches='tight')



