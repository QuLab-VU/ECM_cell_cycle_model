#!/usr/bin/env python
# coding: utf-8

# Simulation of the cell cycle and apoptosis model biochemically and in a cell population
# =========
# 
# ## 1 Simulation with the original parameters
# 
# The parameters refer to the biochemical constants. The first step is importing the model and dependencies. The model and the simulations are done within the package PySB [1].
# 
# 

# In[1]:


import os
import sys
sys.path.insert(0, os.getcwd() + "/cellular_cycle") #add a pathway to the model #not sure about how to use folder dependency, to change
sys.path.insert(0, os.getcwd() + "/fit") #add a pathway to the fit functions #not sure about how to use folder dependency, to change

import numpy as np
import pandas as pd
import pylab as pl
import seaborn as sns
import random

##Importing the model, models differ on their initial conditions, one being a concentration and the other a protein count
from Austin_proliferation_model_concentration_5 import model
model_concentration=model
from Austin_proliferation_model_population_5 import model
model_population=model

##Import the simulator: for ODE and stochastic simulation :
from pysb.simulator import BngSimulator
from pysb.simulator import ScipyOdeSimulator


#change the model parameters of the model if needed
omega = 200
new_param = pd.read_csv(os.getcwd() + "/../umap_window/parameters.csv") 
del new_param['Stress']
#for key in new_param.keys():
#    if model_population.parameters[key].value != new_param[key][0]:
#        model_population.parameters[key].value = new_param[key]
#    if model_concentration.parameters[key].value != new_param[key][0]:
#        model_concentration.parameters[key].value = new_param[key]

# Directory containing cell fate data
from pathlib import Path

DATADIR = Path('../Data/Exp_Intermitotic_Times/')


# ## 2 Simulation of biochemical modem
# 
# Now we execute the model with the original parameters and the initial conditions. Let's start with simulating the biochemical model with a set of ordinary differential equations and then with the Gillespie's Stochastic Simulation Algorithm both available in the PySB package.
# 
# ### Deteriministic simulation
# 

# In[2]:


#ODE simulation

t = pl.linspace(0, 120, 120000)
simres = ScipyOdeSimulator(model_concentration, t, verbose=False).run(); #deterministic simulation
yout = simres.all;
observables=['obsCycE','obsCdc20', 'obsCycA']; #plot Cyclin E, Cdc20 and Cyclin A

for obs in range(3) :
    pl.plot(t,yout[observables[obs]], label=observables[obs], lw=2)
pl.legend(loc=0)
pl.show()


# ### Stochastic simulation
# Set random seed for repeatability.

# In[3]:


get_ipython().run_cell_magic('time', '', "\n### Stochastic simulation\n\nt = pl.linspace(0, 120, 170000);\nsim = BngSimulator(model_population, tspan=t, verbose=False); #stochastic simulator\nsimres = sim.run(method='ssa', verbose=False, seed=33); #running the simulation; remove seed for truly random; seed=33 is used for similarity to deterministic sim\nyout = simres.observables;\nobservables=['obsCycE','obsCdc20', 'obsCycA']; #plot Cyclin E, Cdc20 and Cyclin A\n\n\nfor obs in range(3) :    \n    pl.plot(t,yout[observables[obs]], label=observables[obs], lw=2)\npl.legend(loc=0)\npl.show()")


# ## 3 Simulation of a population of cells
# 
# In order to reproduce faithfully the experiments we run a population of cells. Each events is registered, if death occurs the simulation stops, if a mitosis occurs, a new simulation occurs from the division time.
# 
# ### Simulating the population
# 

# In[4]:


get_ipython().run_cell_magic('time', '', '#from runAKTERK_2_starting_at_mitosis import runAKTERK #Cell population simulation\n#from runStochCell import runStochcell #Simulate the biochemical model with ODEs\nfrom runAKTERK_2 import runAKTERK #Simulate a population of cells from independants biochemical models\nfrom runStochCell import runStochcell\n\nfrom retrieve_cellState_dw import retrieve_cellState #\nfrom call_cell_events import call_events #\n\nimport fracprolif as fp\nfrom rbf_mmd import nmmd2_rbf\n\nrun_time=120 #how long is the experiment\nobserved_cells_time=120 #how long are we looking at new cells\nnber_of_starting_cells=6 #size of initial population\nindex="_" #to be removed...\ncondition="_" #to be removed...\n\n\n#(cells, cell_count) = runAKTERK(nber_of_starting_cells, run_time, observed_cells_time, model_concentration, model_population, index, condition);\n(cells, cell_count) = runAKTERK(nber_of_starting_cells, run_time, observed_cells_time, model_concentration, model_population, index, condition);\n#index and conditions to save the identity of the simulation (to remove?)')


# In[5]:


from call_cell_events import call_events #

IMT_eoe_counts = retrieve_cellState(run_time, nber_of_starting_cells, cells); #Gives the cells life per cycle
all_events_simulated = call_events(run_time, IMT_eoe_counts); #Gives all events (apoptosis, mitosis, quiescence) and when they happen

#plot the results
grid = sns.JointGrid(x='Birthtime', y='Lifespan', data=all_events_simulated) #x axis = birthtime of a cell, y = endtime of this cell
g = grid.plot_joint(sns.scatterplot, hue='event', data=all_events_simulated) 
g.ax_joint.plot([run_time, 0], [0, run_time], 'b-', linewidth = 1, color='black') #Recorded events behind this line
sns.kdeplot(all_events_simulated.loc[all_events_simulated['event']=='mitosis', 'Lifespan'].to_list(), ax=g.ax_marg_y, vertical=True, legend=False, bw=0.3) #output distribution of intermitotic times
pl.ylim(0, 110)


# In[6]:


print(cells.keys()) #different components of 'cells', molecules are all the time courses
len(cells['molecules']) #number of cells life (starts at a mitosis event or time 0 of the experiment and ends at a mitosis, cell death or end of the experiment)

#example
#[cells['molecules'][i]['obsSkp2'] for i in range(len(cells['molecules']))] #all time courses of  Skp2 concatenated

#we iterates over all the cells for all the observable proteins and make it a Dataframe with column=proteins
prot_all = pd.DataFrame(np.array([np.concatenate([cells['molecules'][i][obs.name] for i in range(len(cells['molecules']))]) for obs in model_concentration.observables]).T, columns = [obs.name.replace('obs','') for obs in model_concentration.observables])

#total observed time
print(str(len(prot_all) * (cells['time'][0][1] - cells['time'][0][0])) + ' hours')

#np.array([np.concatenate([cells['molecules'][i][obs.name] for i in range(len(cells['molecules']))]) for obs in model_concentration.observables]).shape
prot_all.head()


# ### Calculating the parameters assuming the fractional proliferation model
# We can calculated the parameters associated with the fractional proliferation model
# 
# *Recall the formula for fractional proliferation model + ref
# 

# In[7]:



mitotic_lifespans = all_events_simulated['Lifespan'][all_events_simulated['event']=='mitosis']
censored_lifespans = all_events_simulated['Lifespan'][all_events_simulated['event']=='End_of_Expt']

number_of_deaths = len(all_events_simulated) - (all_events_simulated['event']=='apoptosis').value_counts()[0]
cell_time_total = sum(all_events_simulated['Lifespan'])
print("death rate = " + str(number_of_deaths/cell_time_total)) #death rate (death event/time)


m = fp.mle(mitotic_lifespans, censored_lifespans);
print("μ = " + str(m['x'][0])) #mu
print("σ = " + str(m['x'][1])) #sigma
print("λ = " + str(m['x'][2])) #lambda


# ### Calculating a cost for this set of parameters
# 
# Once the simulation is done we can calculate a cost associated with the parameters. Let's try with our initial parameters.
# We have two conditions: high and low stiffness (55 and 0.6 kPa). We change the signaling input FAK, ERK and AKT as entries in the model to simulate these different conditions. Then the parameters of the fractional proliferation model are calculated and the difference between the intermitotic distribution for the low stiffness and the difference between the G1 and S to M time distribution for the high stiffness.
# 
# ##### Importing experimental data

# In[8]:


#0.6 kPa low stiffness

all_events_exp = {}
IMT_eoe_counts = pd.DataFrame(columns=['End_of_Expt', 'Death' ,'Birthtime','Lifespan'])
dataframe = pd.read_csv(Path(DATADIR,'ds8_06_DMSO_fates.csv'))
lines = {}
for col in dataframe :
    lines[col]=dataframe[col].dropna().tolist()

for (death, mitotic_time, end_time, lifetime) in zip(lines['Death'], lines['Mitotic Time'], lines['End Time'], lines['Lifetime']):
    tmp_counts = {'End_of_Expt' : end_time==102.3, 'Death' : death==1, 'Birthtime' : mitotic_time, 'Lifespan' : lifetime}
    IMT_eoe_counts = IMT_eoe_counts.append(tmp_counts, ignore_index=True)

all_events_exp['0.6 kPa DMSO'] = call_events(run_time, IMT_eoe_counts)
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


#55 kPa high stiffness #importing G1 and S to M time (0.2 kPa to be removed)
experimental_SMs={}
experimental_G1s={}

dataframe = pd.read_csv(Path(DATADIR,'DS8_SM_0_2_55.csv')) 
for col in dataframe.columns :
    experimental_SMs[col]=dataframe[col].dropna().tolist()
dataframe = pd.read_csv(Path(DATADIR,'DS8_G1_0_2_55.csv'))
for col in dataframe.columns :
    experimental_G1s[col]=dataframe[col].dropna().tolist()


# ##### Low stiffness

# In[9]:


get_ipython().run_cell_magic('time', '', '\nfrom runAKTERK_2_starting_at_mitosis import runAKTERK #Cell population simulation\nimport random\n\nrun_time=102.3 #how long is the experiment\nobserved_cells_time=30 #how long are we looking at new cells\nnber_of_starting_cells=50 #size of initial population\nindex="_" #to be removed...\n\n#Input of the model: Activation of different signalling pathway measured with Western blot\ncondition = \'0.6 kPa DMSO\'\nmodel_concentration.parameters[\'ERK\'].value = 1.06970165050038 \nmodel_concentration.parameters[\'FAK\'].value = 0.337143388171177\nmodel_concentration.parameters[\'AKT\'].value = 0.616658111482869\nmodel_population.parameters[\'ERK\'].value = 1.06970165050038 * 200\nmodel_population.parameters[\'FAK\'].value = 0.337143388171177 * 200\nmodel_population.parameters[\'AKT\'].value = 0.616658111482869 * 200\n\n#Setting the cost to zero\nparameters_sim={}\ncost_frac_prolif_model_006=0\ncost_mmd_006=0\ncost_006=0\n\n#Simulation of the cell population\n(cells, cell_count) = runAKTERK(nber_of_starting_cells, run_time, observed_cells_time, model_concentration, model_population, index, condition);\nIMT_eoe_counts = retrieve_cellState(run_time, nber_of_starting_cells, cells)\nall_events_simulated = call_events(run_time, IMT_eoe_counts)\nall_events_simulated = all_events_simulated[np.array(all_events_simulated[\'Birthtime\']<30)] #We kept the events whose birthtime is inferior to 30 hours to compare to the experiments.\nmitotic_lifespans = all_events_simulated[\'Lifespan\'][all_events_simulated[\'event\']==\'mitosis\']\ncensored_lifespans = all_events_simulated[\'Lifespan\'][all_events_simulated[\'event\']==\'End_of_Expt\']\n\n#calculating fractional proliferation model parameters\nnumber_of_deaths = len(all_events_simulated) - (all_events_simulated[\'event\']==\'apoptosis\').value_counts()[0]\ncell_time_total = sum(all_events_simulated[\'Lifespan\'])\nm = fp.mle(mitotic_lifespans, censored_lifespans)\nparameters_sim[condition] = {}\nparameters_sim[condition][\'a\'] = number_of_deaths/cell_time_total\nparameters_sim[condition][\'mu\'] = m[\'x\'][0]\nparameters_sim[condition][\'sigma\'] = m[\'x\'][1]\nparameters_sim[condition][\'lambda\'] = m[\'x\'][2]\nparameters_sim[condition].update(fp.q_rates(m))\nfor name in [\'a\',\'d\',\'q\']: #apoptosis, division and quiescence rates\n    if parameters_exp[condition][name]<10**-5:\n        parameters_exp[condition][name] = 1e-5\n    if parameters_sim[condition][name]<10**-5:\n        parameters_sim[condition][name] = 1e-5\n    cost_frac_prolif_model_006 += abs(np.log(parameters_sim[condition][name]) - np.log(parameters_exp[condition][name]))\n            \n\n#comparing distributions\nnumber_of_IMT_to_sample=72\nsimulated_IMT = random.sample(mitotic_lifespans.to_list(), number_of_IMT_to_sample)\nexperimental_IMT = random.sample(mitotic_lifespans_exp.to_list(), number_of_IMT_to_sample)\ncost_mmd_006 += nmmd2_rbf(np.array([simulated_IMT]).transpose(), np.array([experimental_IMT]).transpose(), 1) #difference between IMT distribution\n\n\n#Adding the two costs\ncost_006 += cost_frac_prolif_model_006 + cost_mmd_006\n\nprint(cost_006) ')


# ##### High stiffness
# In this case we compare the distribution of G1 and S to M time to the experimental data.

# In[ ]:


get_ipython().run_cell_magic('time', '', '\ninit = 30\nrun_time = 102.3 #how long is the experiment\nobserved_cells_time = 30 #how long are we looking at new cells\nnber_of_starting_cells = 50 #size of initial population\nindex = "_" #to be removed...\ncost_mmd_55 = 0 #set the cost to zero\n\n#Input of the model: Activation of different signalling pathway measured with Western blot\ncondition = \'55 kPa DMSO\'\nmodel_concentration.parameters[\'ERK\'].value = 1\nmodel_concentration.parameters[\'FAK\'].value = 1\nmodel_concentration.parameters[\'AKT\'].value = 1\nmodel_population.parameters[\'ERK\'].value = 1 * 200\nmodel_population.parameters[\'FAK\'].value = 1 * 200\nmodel_population.parameters[\'AKT\'].value = 1 * 200\n\n#Simulation of the cell population\n(cells, cell_count) = runAKTERK(nber_of_starting_cells, run_time, observed_cells_time, model_concentration, model_population, index, condition);\nIMT_eoe_counts = retrieve_cellState(run_time, nber_of_starting_cells, cells)\nall_events_simulated = call_events(run_time, IMT_eoe_counts)\nall_events_simulated = all_events_simulated[np.array(all_events_simulated[\'Birthtime\']<30)] #We kept the events whose birthtime is inferior to 30 hours to compare to the experiments.\nmitotic_lifespans = all_events_simulated[\'Lifespan\'][all_events_simulated[\'event\']==\'mitosis\']\ncensored_lifespans = all_events_simulated[\'Lifespan\'][all_events_simulated[\'event\']==\'End_of_Expt\']\n\n#Comparing G1 and S to M distribution\nnumber_of_SM_to_sample=160\nnumber_of_G1_to_sample=160\nsimulated_G1s = random.sample(simulated_g1, number_of_G1_to_sample)\nexperimental_G1s_picked = random.sample(experimental_G1s[\'55 kPa DMSO\'], number_of_G1_to_sample)\nsimulated_SMs = random.sample(simulated_sm, number_of_SM_to_sample)\nexperimental_SMs_picked = random.sample(experimental_SMs[\'55 kPa DMSO\'], number_of_SM_to_sample)\ncost_mmd_55 += nmmd2_rbf(np.array([simulated_G1s]).transpose(), np.array([experimental_G1s_picked]).transpose(), 1) #difference between IMT distribution\ncost_mmd_55 += nmmd2_rbf(np.array([simulated_SMs]).transpose(), np.array([experimental_SMs_picked]).transpose(), 1) #difference between IMT distribution\n\nprint(cost_mmd_55)')


# ## 4 Simulation with the fitted parameters
# 
# #### Sobol sequence
# From the Sobol sequence

# In[8]:


frac_prolif_parameters = pd.read_pickle(os.getcwd() + "/fit_data/sobol_seq/frac_prolif_parameters_all_parameters_006_55_kPa.pkl")

#Original parameters
pysb_sampled_parameter_names = ["V_s2skp","V_Fp","k_M","k_EMY","k_CDF","k_CDS","k_ED","K_S","K_MY","K_DF","d_CD","k_P1", "k_P2", "K_CD", "k_E", "K_M", "K_CE", "V_2E2F", "K_2E2F", "k_RE", "kf1", "kf2", "kb", "kasbaxmbcl2", "kdsbaxmbcl2", "k1sbh3", "k2sbh3", "kdbh3", "kasbh3bcl2", "kdsbh3bcl2"]
loc=[]
for parameter in pysb_sampled_parameter_names:
    loc+=[np.log10(model.parameters[parameter].value)] 


#Generating the sobol sequence
import sobol_seq
l = sobol_seq.i4_sobol_generate(30, 140000, skip = 0)
l2 = 4*l + (np.array(loc)-2)


# In[9]:


#indices in the sobol sequence for the 100 best parameters
indices = [int(i)-1 for i in frac_prolif_parameters.sort_values(by=['cost_total'])["index"].tolist()[0:100]] 
parameters = [l2[i] for i in indices]

#Change model parameter values to current location in parameter space
parameter_vector = 10**parameters[0] #original parameters
param_dict = {pname: pvalue for pname, pvalue in zip(pysb_sampled_parameter_names, parameter_vector)}   
for pname, pvalue in param_dict.items():
    model_population.parameters[pname].value = (pvalue)
    model_concentration.parameters[pname].value = (pvalue)

sim = BngSimulator(model_population, tspan=t, verbose=False);
simres2 = sim.run(method='ssa', verbose=False, seed=30);
yout = simres2.observables;
observables=['obsCycE','obsCdc20', 'obsCycA'];


for obs in range(3) :    
    pl.plot(t,yout[observables[obs]], label=observables[obs], lw=2)
pl.legend(loc=0)
pl.show()


# #### Particle Swarm Optimization
# 

# In[6]:


frac_prolif_parameters = pd.read_pickle(os.getcwd() + "/fit_data/sobol_seq/frac_prolif_parameters_all_parameters_006_55_kPa.pkl")


# 

# In[7]:


#importing parameters from the pso space exploration

frac_prolif_parameters = pd.read_pickle(os.getcwd() + "/fit_data/sobol_seq/frac_prolif_parameters_all_parameters_006_55_kPa.pkl")
position_history = np.load(os.getcwd() + "/fit_data/pso/position_history.npy", allow_pickle=True)
particle_positions_pso = pd.read_pickle(os.getcwd() + "/fit_data/pso/particle_position_all_parameters_pso.pkl")
frac_prolif_parameters_pso = pd.read_pickle(os.getcwd() + "/fit_data/pso/frac_prolif_parameters_pso.pkl")

#Parameters sampled from the simulation

parameters_log_best=[]
cost_limit = frac_prolif_parameters.sort_values(by=['cost_total'])[0:1]['cost_total'].tolist()[-1] #set the cost limit to the best found in sobol sequence, can be changed
frac_prolif_parameters_pso['cost_total'].tolist()
for idx in frac_prolif_parameters_pso['index']:
    if frac_prolif_parameters_pso.loc[frac_prolif_parameters_pso['index']==idx].iloc[0]['cost_total']<cost_limit: #select the cost  lower than cost limit
        #print(idx)
        frac_prolif_parameters_pso_temp = frac_prolif_parameters_pso.loc[frac_prolif_parameters_pso['index']==idx]
        parameter = particle_positions_pso[particle_positions_pso['identification']==idx]['parameter_vector'].iloc[0]
        pos = np.where(position_history[:,:,:] == parameter)
        position_history[pos[0][0],pos[1][0],:]
        parameters_log_best += [position_history[pos[0][0],pos[1][0],:].tolist()]


parameters_log_best = np.array(parameters_log_best) #here is the set of parameters, note it's a log scale
len(parameters_log_best)


# In[10]:


#parameters_log_best is the list of parameters

parameter_vector = 10**parameters_log_best[1]

param_dict = {pname: pvalue for pname, pvalue in zip(pysb_sampled_parameter_names, parameter_vector)}   
for pname, pvalue in param_dict.items():
    #Change model parameter values to current location in parameter space
    model_population.parameters[pname].value = (pvalue)
    model_concentration.parameters[pname].value = (pvalue)


sim = BngSimulator(model_population, tspan=t, verbose=False);
simres2 = sim.run(method='ssa', verbose=False);
yout = simres2.observables;
observables=['obsCycE','obsCdc20', 'obsCycA'];


for obs in range(3) :    
    pl.plot(t,yout[observables[obs]], label=observables[obs], lw=2)
pl.legend(loc=0)
pl.show()


# ### Bibliography 
# 
# [1]
