# -*- coding: utf-8 -*-
"""
Created on Mon Feb 15 20:30:27 2021

@author: Grim.M
"""

#time stopper:

import signal

class TimeoutException(RuntimeError):
    """ Time out occurred! """
    pass

def handler(signum, frame):
    print('forever is over!')
    raise TimeoutException()


##############################################

path="/home/jeremie"
#path=r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code"
import sys
sys.path.insert(0, path + "/cellular_cycle") #add a pathway to the model
sys.path.insert(0, path + "/fit") #add a pathway 

##
from Austin_proliferation_model_concentration_4 import model
model_concentration=model
from Austin_proliferation_model_population_4 import model
model_population=model

##
import os
os.environ["PYOPENCL_CTX"] ="0"

##
from runAKTERK_2_starting_at_mitosis import runAKTERK
from runStochCell import runStochcell_check

#or
from retrieve_cellState_dw import retrieve_cellState
from plot_frac_div_q_ap1 import all_events

##
from pysb.simulator import BngSimulator
from pysb.simulator import ScipyOdeSimulator

##
import matplotlib.pyplot as plt
import fracprolif as fp
from rbf_mmd import nmmd2_rbf
import numpy as np
import pandas as pd
import pylab as pl
from math import inf
import math
import random
import seaborn as sns
from pathlib import Path
import sobol_seq


import matplotlib.colors as mcolors
import matplotlib.cm as cm
import matplotlib as mpl


#or

from cost_function_BngSimulator_frac_prolif_quiescence_and_mmd_DMSO_only_30_parameters_fixed_initial_concentrations_mitosis import  cost_function_reduced_quiescence_DMSO_only_30_parameters


from cost_function_BngSimulator_frac_prolif_quiescence_and_mmd_DMSO_only_30_parameters_fixed_initial_concentrations_mitosis_time_stopper import  cost_function_reduced_quiescence_DMSO_only_30_parameters





tmp_t=pl.linspace(0, 10, 100)
sim = BngSimulator(model_population, tspan=tmp_t, verbose=True)
traj = sim.run()


sim = ScipyOdeSimulator(model_concentration, tspan=tmp_t, verbose=True)
traj = sim.run()


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

"""
#######################################
#fracprolif birthtime X IMT Length

grid = sns.JointGrid(x='Birthtime', y='Lifespan', data=all_events_exp['0.6 kPa DMSO'])
g = grid.plot_joint(sns.scatterplot, hue='event', data=all_events_exp['0.6 kPa DMSO'])
#sns.kdeplot(all_events.loc[all_events['event']=='mitosis', 'Birthtime'], ax=g.ax_marg_x, legend=False)
sns.kdeplot(all_events_exp['0.6 kPa DMSO'].loc[all_events_exp['0.6 kPa DMSO']['event']=='mitosis', 'Lifespan'], ax=g.ax_marg_y, vertical=True, legend=False, bw=0.3)
plt.ylim(0, 110)



#best param for mmdcost
path=r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code"


idx = int(frac_prolif_parameters[frac_prolif_parameters['cost_mmd'] == frac_prolif_parameters['cost_mmd'].min()]['index'].iloc[0])


try:
    all_events_temp = pd.read_pickle(path + "/fit_data/sobol_seq/cost_function_model_4/all_events_simulated/Tempo/all_events_simulated_" + str(idx) + "_50_initial_cells.pkl")
except:
    all_events_temp = pd.read_pickle(path + "/fit_data/sobol_seq/cost_function_model_4/all_events_simulated/personnal_computer/all_events_simulated_" + str(idx) + "_50_initial_cells.pkl")
all_events_temp_2 = all_events_temp[np.array(all_events_temp['Birthtime']<30)]


grid = sns.JointGrid(x='Birthtime', y='Lifespan', data=all_events_temp_2)
g = grid.plot_joint(sns.scatterplot, hue='event', data=all_events_temp_2)
#sns.kdeplot(all_events.loc[all_events['event']=='mitosis', 'Birthtime'], ax=g.ax_marg_x, legend=False)
sns.kdeplot(all_events_temp_2.loc[all_events_temp_2['event']=='mitosis', 'Lifespan'], ax=g.ax_marg_y, vertical=True, legend=False, bw=0.3)
plt.ylim(0, 110)


#best param for frac prolif cost


idx = int(frac_prolif_parameters[frac_prolif_parameters['cost_frac_prolif_model'] == frac_prolif_parameters['cost_frac_prolif_model'].min()]['index'].iloc[0])

try:
    all_events_temp = pd.read_pickle(path + "/fit_data/sobol_seq/cost_function_model_4/all_events_simulated/Tempo/all_events_simulated_" + str(idx) + "_50_initial_cells.pkl")
except:
    all_events_temp = pd.read_pickle(path + "/fit_data/sobol_seq/cost_function_model_4/all_events_simulated/personal_computer/all_events_simulated_" + str(idx) + "_50_initial_cells.pkl")
    
all_events_temp_2 = all_events_temp[np.array(all_events_temp['Birthtime']<30)]


grid = sns.JointGrid(x='Birthtime', y='Lifespan', data=all_events_temp_2)
g = grid.plot_joint(sns.scatterplot, hue='event', data=all_events_temp_2)
#sns.kdeplot(all_events.loc[all_events['event']=='mitosis', 'Birthtime'], ax=g.ax_marg_x, legend=False)
sns.kdeplot(all_events_temp_2.loc[all_events_temp_2['event']=='mitosis', 'Lifespan'], ax=g.ax_marg_y, vertical=True, legend=False, bw=0.3)
plt.ylim(0, 110)


#best param for global cost


idx = int(frac_prolif_parameters[frac_prolif_parameters['cost'] == frac_prolif_parameters['cost'].min()]['index'].iloc[0])

try:
    all_events_temp = pd.read_pickle(path + "/fit_data/sobol_seq/cost_function_model_4/all_events_simulated/Tempo/all_events_simulated_" + str(idx) + "_50_initial_cells.pkl")
except:
    all_events_temp = pd.read_pickle(path + "/fit_data/sobol_seq/cost_function_model_4/all_events_simulated/personnal_computer/all_events_simulated_" + str(idx) + "_50_initial_cells.pkl")
all_events_temp_2 = all_events_temp[np.array(all_events_temp['Birthtime']<30)]
   


grid = sns.JointGrid(x='Birthtime', y='Lifespan', data=all_events_temp_2)
g = grid.plot_joint(sns.scatterplot, hue='event', data=all_events_temp_2)
#sns.kdeplot(all_events.loc[all_events['event']=='mitosis', 'Birthtime'], ax=g.ax_marg_x, legend=False)
sns.kdeplot(all_events_temp_2.loc[all_events_temp_2['event']=='mitosis', 'Lifespan'], ax=g.ax_marg_y, vertical=True, legend=False, bw=0.3)
plt.ylim(0, 110)



#first fittedparameters

all_events_temp = pd.read_pickle(path + "/fit_data/all_events_simulated_first_fitted_parameters.pkl")


all_events_temp_2 = all_events_temp[np.array(all_events_temp['Birthtime']<30)]


grid = sns.JointGrid(x='Birthtime', y='Lifespan', data=all_events_temp_2)
g = grid.plot_joint(sns.scatterplot, hue='event', data=all_events_temp_2)
#sns.kdeplot(all_events.loc[all_events['event']=='mitosis', 'Birthtime'], ax=g.ax_marg_x, legend=False)
sns.kdeplot(all_events_temp_2.loc[all_events_temp_2['event']=='mitosis', 'Lifespan'], ax=g.ax_marg_y, vertical=True, legend=False, bw=0.3)
plt.ylim(0, 110)


#Combined parameters

all_events_combined_parameters = pd.read_pickle(path + "/fit_data/sobol_seq/cost_function/events_for_best_parameters/all_events_simulated_combined_parameters_50_initial_cells.pkl")
all_events_temp_2 = all_events_combined_parameters[np.array(all_events_combined_parameters['Birthtime']<30)]


grid = sns.JointGrid(x='Birthtime', y='Lifespan', data=all_events_temp_2)
g = grid.plot_joint(sns.scatterplot, hue='event', data=all_events_temp_2)
#sns.kdeplot(all_events.loc[all_events['event']=='mitosis', 'Birthtime'], ax=g.ax_marg_x, legend=False)
sns.kdeplot(all_events_temp_2.loc[all_events_temp_2['event']=='mitosis', 'Lifespan'], ax=g.ax_marg_y, vertical=True, legend=False, bw=0.3)
plt.ylim(0, 110)

"""

##########################################################################
#simulation

####################################################################################
#Generating the parameters according to the sobol sequence

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


def cost_function_f(parameter_vector, idx):
    parameter_vector=[10**i for i in parameter_vector]
    return -cost_function_reduced_quiescence_DMSO_only_30_parameters(parameter_vector, idx)


cost_function_f(l2[0], 1)

###############################
"""
starts=np.linspace(70101,80101,101)
ends=np.linspace(70200,82200,101)

    l=sobol_seq.i4_sobol_generate(30, 22, skip = 10188)


for n in range(len(starts)):
    l=sobol_seq.i4_sobol_generate(30, 22, skip = int(starts[n]-1))
    l2 = 4*l + (np.array(loc)-2)
    cost=[]
    for i in range(len(l2)) :
        parameter_vector=l2[i]
        print(parameter_vector)
        c = cost_function_f(parameter_vector, int(starts[n] + i))
        print(cost)
        cost+=[c]
        
    np.save(path + r"\fit_data\sobol_seq\cost_function_model_4\sobol_sequence\cost_sobol_seq_test_" + str(int(starts[n])) + "_" + str(int(ends[n])) +".npy", cost)



starts=np.linspace(1401, 9801,43)
ends=np.linspace(1600, 10000,43)


for n in range(len(starts)):
    l=sobol_seq.i4_sobol_generate(30, 200, skip = int(starts[n]-1))
    l2 = 4*l + (np.array(loc)-2)
    cost=[]
    for i in range(len(l2)) :
        parameter_vector=l2[i]
        print(parameter_vector)
        c = cost_function_f(parameter_vector, int(starts[n] + i))
        print(cost)
        cost+=[c]
    np.save(path + r"/fit_data/sobol_seq/cost_function_model_4/sobol_sequence/cost_sobol_seq_test_" + str(int(starts[n])) + "_" + str(int(ends[n])) +".npy", cost)



starts=np.linspace(30101, 50101,101)
ends=np.linspace(30300, 50300,101)


starts=np.linspace(50301, 60101,50)
ends=np.linspace(50300, 60300,50)


starts=np.linspace(82001,99801,90)
ends=np.linspace(82200,100000,90)


starts=np.linspace(100001,119801,100)
ends=np.linspace(100200,120000,100)



l=sobol_seq.i4_sobol_generate(30, 200000, skip=0)
l2 = 4*l + (np.array(loc)-2)


"""
########################$

path=r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code"

frac_prolif_parameters = pd.read_pickle(path + r"/fit_data/sobol_seq/cost_function_model_4/cost_and_IMT_distribution/frac_prolif_parameters_all_parameters_50_initial_cells.pkl")


indexes_calculated = np.array([int(i) for i in frac_prolif_parameters['index'].tolist() if not np.isnan(i)])

starts=np.linspace(1,9801,50)
ends=np.linspace(200,10000,50)


starts = np.linspace(10701,19901,93)
ends = np.linspace(10800,20000,93)

starts=np.linspace(20001,29801,50)
ends=np.linspace(20200,30000,50)




starts=np.linspace(30001,39901,100)
ends=np.linspace(30100,40000,100)



starts=np.linspace(40001,42801,15)
ends=np.linspace(40200,43000,15)


starts=np.linspace(43001,49801,35)
ends=np.linspace(43200,50000,35)



starts=np.linspace(50001,59801,50)
ends=np.linspace(50200,60000,50)



l=sobol_seq.i4_sobol_generate(30, 863, skip = 69137)
l2 = 4*l + (np.array(loc)-2)
cost=[]
for i in range(len(l2)) :
    if np.any(indexes_calculated==int(69138 + i)):
        parameter_vector=l2[i]
        print(parameter_vector)
    else:
        try :
            parameter_vector=l2[i]
            print(parameter_vector)
            c = cost_function_f(parameter_vector, int(69138 + i))
            print(cost)
            cost+=[c]
        except :
            print(69138 + i)
   


starts=np.linspace(60001,69901,100)
ends=np.linspace(60100,70000,100)




for n in range(len(starts)):
    l=sobol_seq.i4_sobol_generate(30, 100, skip = int(starts[n]-1))
    l2 = 4*l + (np.array(loc)-2)
    cost=[]
    for i in range(len(l2)) :
        if np.any(indexes_calculated==int(starts[n] + i)):
            parameter_vector=l2[i]
            print(parameter_vector)
        else:
            try :
                parameter_vector=l2[i]
                print(parameter_vector)
                c = cost_function_f(parameter_vector, int(starts[n] + i))
                print(cost)
                cost+=[c]
            except :
                print('time')



"""
l=sobol_seq.i4_sobol_generate(30, 111, skip = 489)
l2 = 4*l + (np.array(loc)-2)
cost=[]
for i in range(len(l2)) :
    if np.any(indexes_calculated==int(490 + i)):
        parameter_vector=l2[i]
        print(parameter_vector)
    else:
        signal.alarm(1000)
        try :
            parameter_vector=l2[i]
            print(parameter_vector)
            c = cost_function_f(parameter_vector, int(490 + i))
            signal.alarm(0)
            print(cost)
            cost+=[c]
        except TimeoutException as exc:
            print(490 + i)
            signal.alarm(0)
"""


starts=np.linspace(70001,79801,50)
ends=np.linspace(70200,80000,50)


for n in range(len(starts)):
    l=sobol_seq.i4_sobol_generate(30, 200, skip = int(starts[n]-1))
    l2 = 4*l + (np.array(loc)-2)
    cost=[]
    for i in range(len(l2)) :
        if np.any(indexes_calculated==int(starts[n] + i)):
            parameter_vector=l2[i]
            print(parameter_vector)
        else:
            signal.alarm(1500)
            try :
                parameter_vector=l2[i]
                print(parameter_vector)
                c = cost_function_f(parameter_vector, int(starts[n] + i))
                signal.alarm(0)
                print(cost)
                cost+=[c]
            except TimeoutException as exc:
                print(starts[n] + i)




"""
######################
#best parameters

l=sobol_seq.i4_sobol_generate(30, 120000, skip=0)
l2 = 4*l + (np.array(loc)-2)


path=r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code"

frac_prolif_parameters = pd.read_pickle(path + "/fit_data/sobol_seq/cost_function_model_4/cost_and_IMT_distribution/frac_prolif_parameters_all_parameters_50_initial_cells.pkl")



for idx in frac_prolif_parameters.sort_values(by=['cost'])[0:10]['index'] :
    parameter_vector=l2[int(idx)-1]
    print(idx)
    print(parameter_vector)
    
    cost=[]
    for n in range(16,20):
        c = cost_function_f(parameter_vector, str(int(idx)) + '_iteration_n_'+str(n))
        cost+=[c]    
    
    
    np.save(path + r"/fit_data/sobol_seq/cost_function_model_4/cost_and_IMT_distribution/cost_variability_sobol_index_" + str(int(idx+1)) +".npy", np.array(cost))

#111131
#array([-0.36190796, -0.8420667 , -2.4096266 , -0.5112375 , -1.37202509,
#       -1.78584731, -2.54818265, -1.13673545, -1.06246948,  1.1843705 ,
#       -0.40078252,  3.23156035,  0.8024832 ,  0.31247573,  0.37192693,
#       -0.71071927, -1.86654039,  1.39682533, -0.79452975,  4.25371611,
#       -2.88876343,  0.1516651 , -1.39735951,  0.69625877, -1.44900973,
#       -3.16769409, -2.44248358, -1.92085848,  2.43538508, -3.36050415])

####################""
#Comparing costs
    


sigma=1.5

cost_range = np.linspace(-35,0,400)
cost_distribution = pd.DataFrame(cost_range)  


path=r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code\fit_data\sobol_seq\cost_function_model_4\cost_and_IMT_distribution\personal_computer"
os.chdir(path)
for filename in os.listdir(os.getcwd()) :
    if 'cost_variability_sobol_index' in filename:
        print(filename)
        cost_temp = np.load(path + '/' + filename)
        
        cost_distribution_temp = [RBS(i, np.transpose([cost_temp]), sigma) for i in np.linspace(-35,0,400)]
        cost_distribution[int(filename.replace("cost_variability_sobol_index_", "").replace(".npy", ""))] = pd.Series(cost_distribution_temp, index=cost_distribution.index)
        
path=r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code\fit_data\sobol_seq\cost_function_model_4\cost_and_IMT_distribution\Tempo"
os.chdir(path)
import os
for filename in os.listdir(os.getcwd()) :
    if 'cost_variability_sobol_index' in filename:
        print(filename)
        cost_temp = np.load(path + '/' + filename)
        cost_distribution_temp = [RBS(i, np.transpose([cost_temp]), sigma) for i in np.linspace(-35,0,400)]
        cost_distribution[int(filename.replace("cost_variability_sobol_index_", "").replace(".npy", ""))] = pd.Series(cost_distribution_temp, index=cost_distribution.index)
        
    

cost_range=np.linspace(-35,0,400)

for col in cost_distribution.columns:
    plt.plot(cost_range, cost_distribution[col])


"""



"""


#setting the parameters for the model simulations
condition= '0.6 kPa DMSO'

if 'DMSO' in condition:
    if '0.6 kPa' in condition :
        model_concentration.parameters['ERK'].value = 1.06970165050038 
        model_concentration.parameters['FAK'].value = 0.337143388171177
        model_concentration.parameters['AKT'].value = 0.616658111482869
        
        model_population.parameters['ERK'].value = 1.06970165050038 * omega
        model_population.parameters['FAK'].value = 0.337143388171177 * omega
        model_population.parameters['AKT'].value = 0.616658111482869 * omega
    elif '4 kPa' in condition :
        model_concentration.parameters['ERK'].value = 0.897485095798471
        model_concentration.parameters['FAK'].value = 0.594150207445117
        model_concentration.parameters['AKT'].value = 1.01975420186633
    
        model_population.parameters['ERK'].value = 0.897485095798471 * omega
        model_population.parameters['FAK'].value = 0.594150207445117 * omega
        model_population.parameters['AKT'].value = 1.01975420186633 * omega
    else : #55 kPa
        model_concentration.parameters['ERK'].value = 1
        model_concentration.parameters['FAK'].value = 1
        model_concentration.parameters['AKT'].value = 1
        
        model_population.parameters['ERK'].value = 1 * omega
        model_population.parameters['FAK'].value = 1 * omega
        model_population.parameters['AKT'].value = 1 * omega
else : #ERL treated
    if '0.6 kPa' in condition :
        model_concentration.parameters['ERK'].value = 0.186745227418267
        model_concentration.parameters['FAK'].value = 0.397060244114162
        model_concentration.parameters['AKT'].value = 0.245702711121333
        
        model_population.parameters['ERK'].value = 0.186745227418267 * omega
        model_population.parameters['FAK'].value = 0.397060244114162 * omega
        model_population.parameters['AKT'].value = 0.245702711121333 * omega
    elif '4 kPa' in condition :
        model_concentration.parameters['ERK'].value = 0.211453830597783
        model_concentration.parameters['FAK'].value = 0.709172297576571
        model_concentration.parameters['AKT'].value = 0.598518145949655
    
        model_population.parameters['ERK'].value = 0.211453830597783 * omega
        model_population.parameters['FAK'].value = 0.709172297576571 * omega
        model_population.parameters['AKT'].value = 0.598518145949655 * omega
    else : #55 kPa
        model_concentration.parameters['ERK'].value = 0.273618423688386
        model_concentration.parameters['FAK'].value = 1.075820621939
        model_concentration.parameters['AKT'].value = 0.585067267459499
        
        model_population.parameters['ERK'].value = 0.273618423688386 * omega
        model_population.parameters['FAK'].value = 1.075820621939 * omega
        model_population.parameters['AKT'].value = 0.585067267459499 * omega




init = 50

cost = np.load(path+"/fit_data/cost_sobol_seq_0000_44000.npy")
index = [idx for idx, element in enumerate(cost) if element<=-6 and element>-inf]


for idx in index[416:616]:
    print(idx)
    parameter_vector = [10**i for i in l2[idx]]
    print(parameter_vector)
    pysb_sampled_parameter_names = ["V_s2skp","V_Fp","k_M","k_EMY","k_CDF","k_CDS","k_ED","K_S","K_MY","K_DF","d_CD","k_P1", "k_P2", "K_CD", "k_E", "K_M", "K_CE", "V_2E2F", "K_2E2F", "k_RE", "kf1", "kf2", "kb", "kasbaxmbcl2", "kdsbaxmbcl2", "k1sbh3", "k2sbh3", "kdbh3", "kasbh3bcl2", "kdsbh3bcl2"]
    param_dict = {pname: pvalue for pname, pvalue in zip(pysb_sampled_parameter_names, parameter_vector)}
    for pname, pvalue in param_dict.items():
        model_population.parameters[pname].value = (pvalue)
        model_concentration.parameters[pname].value = (pvalue)
        
    signal.alarm(600)
    try:
        (cells, cell_count) = runAKTERK(init, run_time, observed_cells_time, model_concentration, model_population)
        IMT_eoe_counts = retrieve_cellState(run_time, init, cells)
        all_events_simulated = all_events(run_time, IMT_eoe_counts)
        all_events_simulated.to_pickle(path + "/fit_data/all_events_simulated_" + str(idx) + "_50_initial_cells.pkl")
        signal.alarm(0)
    except TimeoutException as exc:
        print('f')




for idx in index[861:900] :
    print(idx)
    parameter_vector = [10**i for i in l2[idx]]
    print(parameter_vector)
    pysb_sampled_parameter_names = ["V_s2skp","V_Fp","k_M","k_EMY","k_CDF","k_CDS","k_ED","K_S","K_MY","K_DF","d_CD","k_P1", "k_P2", "K_CD", "k_E", "K_M", "K_CE", "V_2E2F", "K_2E2F", "k_RE", "kf1", "kf2", "kb", "kasbaxmbcl2", "kdsbaxmbcl2", "k1sbh3", "k2sbh3", "kdbh3", "kasbh3bcl2", "kdsbh3bcl2"]
    param_dict = {pname: pvalue for pname, pvalue in zip(pysb_sampled_parameter_names, parameter_vector)}
    for pname, pvalue in param_dict.items():
        model_population.parameters[pname].value = (pvalue)
        model_concentration.parameters[pname].value = (pvalue)
        
    (cells, cell_count) = runAKTERK(init, run_time, observed_cells_time, model_concentration, model_population)
    IMT_eoe_counts = retrieve_cellState(run_time, init, cells)
    all_events_simulated = all_events(run_time, IMT_eoe_counts)
    all_events_simulated.to_pickle(path + "/fit_data/all_events_simulated_" + str(idx) + "_50_initial_cells.pkl")
    
    

################################################
    #simulation from the parameters obtained in the first fit

pysb_sampled_parameter_names = ["V_s2skp","V_Fp","k_M","k_EMY","k_CDF","k_CDS","k_ED","K_S","K_MY","K_DF"]

parameter_vector = np.array([5.20103250e-02, 4.42837548e-02, 2.39078059e-01, 3.91191419e-02, 1.71177745e-04, 1.23105202e-04, 4.61787386e-01, 1.12304753e+00, 1.59868274e-02, 4.27401052e-01]) #god
param_dict = {pname: pvalue for pname, pvalue in zip(pysb_sampled_parameter_names, parameter_vector)}

for pname, pvalue in param_dict.items():  #Change model parameter values to current location in parameter space
    model_population.parameters[pname].value = (pvalue)
    model_concentration.parameters[pname].value = (pvalue)


(cells, cell_count) = runAKTERK(init, run_time, observed_cells_time, model_concentration, model_population)
IMT_eoe_counts = retrieve_cellState(run_time, init, cells)
all_events_simulated = all_events(run_time, IMT_eoe_counts)
all_events_simulated.to_pickle(path + "/fit_data/all_events_simulated_first_fitted_parameters.pkl")


######################################################################

all_events_simulated_2 = all_events_simulated[np.array(all_events_simulated['Birthtime']<30)]


grid = sns.JointGrid(x = 'Birthtime', y = 'Lifespan', data = all_events_simulated_2)
g = grid.plot_joint(sns.scatterplot, hue = 'event', data = all_events_simulated_2)
#sns.kdeplot(all_events.loc[all_events['event']=='mitosis', 'Birthtime'], ax=g.ax_marg_x, legend=False)
sns.kdeplot(all_events_simulated_2.loc[all_events_simulated_2['event']=='mitosis', 'Lifespan'], ax=g.ax_marg_y, vertical=True, legend=False, bw=4)
plt.ylim(0, 110)


#######################################################################

index = [271, 645, 2552, 2734, 3072, 3189, 3438, 4531, 5728, 5844, 6462, 7031, 9326, 9714, 10165, 10248, 10673, 12970, 14492, 14638, 17291, 17567, 18206, 18358, 19054, 20294, 24059, 24090, 24823, 26001, 26755, 29368, 31644, 31979, 32114, 32427, 32726,32972, 33013, 33965, 34305, 35499, 35902, 36068, 36284, 36307, 38622, 40987, 42037, 43083, 43586] #todo

index_r = [271, 3072, 5844, 32114]


index_ = [2552, 2734, 3189,  4531, 5728, 6462, 7031, 9326, 9714, 10165, 10248, 10673, 12970, 14492, 14638, 17291, 17567, 18358, 19054, 20294, 24059, 24090, 24823, 26001, 26755, 31644, 31979, 32427, 32726, 32972, 33965, 34305, 35499, 35902, 36068, 36284, 36307, 38622, 40987, 43083, 43586]


##################################


cost = np.load(path+"/fit_data/cost_sobol_seq_0000_44000.npy")
index = [idx for idx, element in enumerate(cost) if element<=-6 and element>-inf]


index_ = [0, 38, 98, 136, 226, 252, 261, 341, 376, 399, 405, 434, 449, 481, 488, 522, 532, 595, 627, 637, 660, 692, 704, 742, 745, 764, 782, 813, 892, 900, 917, 977, 989, 1003, 1047, 1063, 1070, 1096, 1157, 1205, 1280, 1323, 1338, 1362, 1369, 1409, 1418, 1471, 1496, 1537, 1539, 1552, 1621, 1658, 1670, 1738, 1746, 1764, 1820, 1876, 1890, 1891, 1970, 2000, 2003, 2008, 2014, 2037, 2055, 2072, 2144, 2153, 2202, 2280, 2313, 2345, 2353, 2399, 2418, 2498, 2522, 2530, 2585, 2604, 2667, 2684, 2691, 2757, 2800, 2823, 2996, 3091, 3096, 3129, 3165, 3207, 3209, 3281, 3309, 3326, 3350, 3351, 3371, 3454, 3471, 3495, 3512, 3554, 3560, 3568, 3587, 3595, 3619, 3670, 3711, 3716, 3726, 3743, 3789]

##################################"

#load

not_found=[]
index__=[]
index_ = []
for idx in index:
    try :
        all_events_best_temp = pd.read_pickle(path + "/fit_data/sobol_seq/cost_function/events_for_best_parameters/Tempo/all_events_simulated_" + str(idx) + "_50_initial_cells.pkl")
        print("Tempo")
        index_+=["Tempo/all_events_simulated_" + str(idx) + "_50_initial_cells.pkl"]
        index__+=[idx]
    except:
        print('personnal?')
        try :
            all_events_best_temp = pd.read_pickle(path + "/fit_data/sobol_seq/cost_function/events_for_best_parameters/personnal_computer_j/all_events_simulated_" + str(idx) + "_50_initial_cells.pkl")
            print("personnal")
            index_+=["personnal_computer_j/all_events_simulated_" + str(idx) + "_50_initial_cells.pkl"]
            index__+=[idx]
        except :
            print("not found")
            not_found+=[idx]

#remove those that have been calculated:


for (idx_, idx__) in zip(index_, index__):
    print(idx_)
    try :
        all_events_best_temp = pd.read_pickle(path + "/fit_data/sobol_seq/cost_function/cost_and_IMT_distribution/parameters_" + str(idx__) + "_50_initial_cells.pkl")
        print("Tempo")
        index_.remove(idx_)
        index__.remove(idx__)
    except:
        print('personnal?')




for idx_file, idx in zip(index_,index__):
    
    all_events_best_temp = pd.read_pickle(path + "/fit_data/sobol_seq/cost_function/events_for_best_parameters/" + idx_file)
    all_events_best_temp_2 = all_events_best_temp[np.array(all_events_best_temp['Birthtime']<30)]
    
    
    #Calculating cost for fractional proliferation model
    mitotic_lifespans = all_events_best_temp_2['Lifespan'][all_events_best_temp_2['event']=='mitosis']
    censored_lifespans = all_events_best_temp_2['Lifespan'][all_events_best_temp_2['event']=='End_of_Expt']
    parameters_sim={}
    number_of_deaths = len(all_events_best_temp_2) - (all_events_best_temp_2['event']=='apoptosis').value_counts()[0]
    
    
    cell_time_total = sum(all_events_best_temp_2['Lifespan'])
    m = fp.mle(mitotic_lifespans, censored_lifespans)
        
    parameters_sim['0.6 kPa DMSO'] = {}
    parameters_sim['0.6 kPa DMSO']['a'] = number_of_deaths/cell_time_total
    parameters_sim['0.6 kPa DMSO']['mu'] = m['x'][0]
    parameters_sim['0.6 kPa DMSO']['sigma'] = m['x'][1]
    parameters_sim['0.6 kPa DMSO']['lambda'] = m['x'][2]
    parameters_sim['0.6 kPa DMSO'].update(fp.q_rates(m))
    
    cost_frac_prolif_model = 0
    for name in parameter_names:
        print(name)
        if parameters_exp['0.6 kPa DMSO'][name]<10**-5:
            parameters_exp['0.6 kPa DMSO'][name] = 1e-5
        if parameters_sim['0.6 kPa DMSO'][name]<10**-5:
            parameters_sim['0.6 kPa DMSO'][name] = 1e-5
        
        cost_frac_prolif_model += abs(np.log(parameters_sim['0.6 kPa DMSO'][name]) - np.log(parameters_exp['0.6 kPa DMSO'][name]))
        
    if len(mitotic_lifespans)>=50:
        #comparing those parameters and mmd:
        #number_of_IMT_to_sample = min (len(mitotic_lifespans), len(mitotic_lifespans_exp))
        number_of_IMT_to_sample = 50
        simulated_IMT = random.sample(mitotic_lifespans.to_list(), number_of_IMT_to_sample)
        experimental_IMT = random.sample(mitotic_lifespans_exp.to_list(), number_of_IMT_to_sample)
        
        cost_mmd = nmmd2_rbf(np.array([simulated_IMT]).transpose(), np.array([experimental_IMT]).transpose(), 1) #difference between IMT distribution
    
    else :
        cost_mmd=50
        
    #adding the costs
    cost = cost_mmd + cost_frac_prolif_model
    
    x=np.linspace(0,120,360)
    y=[RBS(i, np.transpose([mitotic_lifespans.tolist()]), sigma) for i in x]
    
    #updating
    
    
    frac_prolif_parameters_temp = {"index":[idx], "a":[parameters_sim['0.6 kPa DMSO']['a']], "d":[parameters_sim['0.6 kPa DMSO']['d']], "q":[parameters_sim['0.6 kPa DMSO']['q']], "cost":[cost], "cost_mmd":[cost_mmd], "cost_frac_prolif_model":[cost_frac_prolif_model]}
    
    frac_prolif_parameters_temp = pd.DataFrame(data=frac_prolif_parameters_temp)
    
    
    np.save(path + "/fit_data/sobol_seq/cost_function/cost_and_IMT_distribution/IMT_distribution_" + str(idx) + "_50_initial_cells.npy",np.array(y))
    
    frac_prolif_parameters_temp.to_pickle(path + "/fit_data/sobol_seq/cost_function/cost_and_IMT_distribution/parameters_" + str(idx) + "_50_initial_cells.pkl")
    
    
    
"""
####################################################
#testing

"""
#parameters to test

t = np.linspace(0, 600, 400)
parameters_to_be_tested = [l2[0]]

path=r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code\fit_data"
os.chdir(path)

for filename in os.listdir(os.getcwd())[0:99] :
    a = np.load(path+ "/" +  filename)
    parameters_to_be_tested = np.concatenate([parameters_to_be_tested, [a]], axis=0)

model_concentration.parameters['ERK'].value = 1.06970165050038 
model_concentration.parameters['FAK'].value = 0.337143388171177
model_concentration.parameters['AKT'].value = 0.616658111482869


parameter_vector = [10**i for i in parameters_to_be_tested[0]]
param_dict = {pname: pvalue for pname, pvalue in zip(pysb_sampled_parameter_names, parameter_vector)}  
for pname, pvalue in param_dict.items(): 
    model_concentration.parameters[pname].value = (pvalue) #Change model parameter values to current location in parameter space
simres = ScipyOdeSimulator(model_concentration, t, verbose=True).run();

plt.plot(simres.observables['obsCdc20'])
plt.plot(simres.observables['obsCycD'])


from pathlib import Path
fn = Path(r'C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code\fit_data\Cdc20_CycD.pdf').expanduser()


fig, axs = plt.subplots(10, 10)

axs[0, 0].plot(t,simres.observables['obsCdc20'], label="idx_1_Cdc20", color='orange', lw=0.05)
axs[0, 0].plot(t,simres.observables['obsCycD'], label="idx_1_CycD", color='blue', lw=0.05)
axs[0, 0].legend(fontsize=1)

for i in range(10):
    for j in range(10):
        if (i!=0 or j!=0) :
            parameter_vector = [10**i for i in parameters_to_be_tested[i+j*10]]
            param_dict = {pname: pvalue for pname, pvalue in zip(pysb_sampled_parameter_names, parameter_vector)}  
            for pname, pvalue in param_dict.items(): 
                model_concentration.parameters[pname].value = (pvalue) #Change model parameter values to current location in parameter space
            simres = ScipyOdeSimulator(model_concentration, t, verbose=True).run();
            
            filename = os.listdir(os.getcwd())[0:99][i-1+j*10]
            filename = filename.replace("Cdc20_or_CycD_not_cycling_parameters_idx_", "").replace(".npy", "")
            
            axs[i, j].plot(t,simres.observables['obsCdc20'], label= "idx_" + filename + "_Cdc20", color='orange', lw=0.05)
            axs[i, j].plot(t,simres.observables['obsCycD'], label= "idx_" + filename + "_CycD", color='blue', lw=0.05)
            axs[i, j].legend(fontsize=1)


for ax in axs.flat:
    ax.label_outer()
    
plt.savefig(fn, bbox_inches='tight')
"""

#####################################################################

#####################################################################

"""
path=r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code\fit_data\sobol_seq\cost_function_model_4\sobol_sequence\personal_computer"
os.chdir(path)

cost_sobol_seq_test_10001_20000=np.array([])

import os
for filename in os.listdir(os.getcwd())[0:200] :
    
    a = np.load(path+ "/" +  filename)
    
    cost_sobol_seq_test_10001_26300 = np.concatenate([cost_sobol_seq_test_10001_26300, a])
 

cost_sobol_seq_test_20001_30000=np.array([])

import os
for filename in os.listdir(os.getcwd())[201:272] :
    
    a = np.load(path+ "/" +  filename)
    
    cost_sobol_seq_test_20001_26300 = np.concatenate([cost_sobol_seq_test_20001_26300, a])




path=r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code\fit_data\sobol_seq\cost_function_model_4\sobol_sequence\Tempo"
os.chdir(path)

cost_sobol_seq_test_00001_10000=np.array([])

starts = np.linspace(1, 551, 12)

cost_sobol_seq_test_00001_30000=np.array([])


for s in starts:
    a = np.load(path+ "/" +  "cost_sobol_seq_test_" + str(int(s)) + "_" + str(int(s+49)) + ".npy")
    cost_sobol_seq_test_00001_10000 = np.concatenate([cost_sobol_seq_test_20001_26300, a])


for s in starts:
    a = np.load(path+ "/" +  "cost_sobol_seq_test_" + str(int(s)) + "_" + str(int(s+49)) + ".npy")
    cost_sobol_seq_test_00001_10000 = np.concatenate([cost_sobol_seq_test_20001_26300, a])
"""
#####################################################################


#####################################################################


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
        
        
################
        


x=np.linspace(0,120,360)
IMT_distribution = pd.DataFrame(x)  


frac_prolif_parameters = {"index":[], "a":[], "d":[], "q":[], "cost":[], "cost_mmd":[], "cost_frac_prolif_model":[]}
frac_prolif_parameters = pd.DataFrame(data=frac_prolif_parameters)

#Personal comp

path=r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code\fit_data\sobol_seq\cost_function_model_4\cost_and_IMT_distribution\personal_computer"
os.chdir(path)
import os
for filename in os.listdir(os.getcwd()) :
    if 'parameters' in filename and '_50_initial_cells.pkl' in filename:
        frac_prolif_parameters_temp = pd.read_pickle(path + "/" + filename)
        frames = [frac_prolif_parameters, frac_prolif_parameters_temp]
        frac_prolif_parameters = pd.concat(frames)
    
    elif 'IMT_distribution_' in filename:
        y = np.load(path + '/' + filename)
        IMT_distribution[int(filename.replace("IMT_distribution_", "").replace("_50_initial_cells.npy", ""))] = pd.Series(y, index=IMT_distribution.index)

#Tempo

path=r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code\fit_data\sobol_seq\cost_function_model_4\cost_and_IMT_distribution\Tempo"
os.chdir(path)
import os
for filename in os.listdir(os.getcwd()) :
    if 'parameters' in filename:
        frac_prolif_parameters_temp = pd.read_pickle(path + "/" + filename)
        frames = [frac_prolif_parameters, frac_prolif_parameters_temp]
        frac_prolif_parameters = pd.concat(frames)
    
    elif 'IMT_distribution_' in filename:
        y = np.load(path + '/' + filename)
        IMT_distribution[int(filename.replace("IMT_distribution_", "").replace("_50_initial_cells.npy", ""))] = pd.Series(y, index=IMT_distribution.index)




##########

t = np.linspace(0, 800, 800)

T = t[1] - t[0]  # sampling interval 
N = ssa['obsCycE'][200:].size
N=600
# 1/T = frequency
f = np.linspace(0, 1 / T, N)


ffts_obsCycE = pd.DataFrame(f)  


#ffts_obsCycE={}
ffts_obsCycD= pd.DataFrame(f)  
ffts_obsCdc20= pd.DataFrame(f)  

path=r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code\fit_data\sobol_seq\cost_function_model_4\fft\stochastic_simulation\\"
os.chdir(path)
import os
for filename in os.listdir(os.getcwd()) :
    if "obsCycE" in filename :
        fft_temp = np.load(path + filename)
        ffts_obsCycE[int(filename.replace("fft_stochastic_simulation_obsCycE_200_to_800_hours_index_", "").replace(".npy", ""))] = pd.Series(fft_temp, index=ffts_obsCycE.index)
    elif "obsCycD" in filename :
        fft_temp = np.load(path + filename)
        ffts_obsCycD[int(filename.replace("fft_stochastic_simulation_obsCycD_200_to_800_hours_index_", "").replace(".npy", ""))] = pd.Series(fft_temp, index=ffts_obsCycD.index)
    elif "obsCdc20" in filename :
        fft_temp = np.load(path + filename)
        ffts_obsCdc20[int(filename.replace("fft_stochastic_simulation_obsCdc20_200_to_800_hours_index_", "").replace(".npy", ""))] = pd.Series(fft_temp, index=ffts_obsCdc20.index)
        



###########################################################################

#Remove bad tmp scores (old)

path=r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code"

bad=[]

for i in range(len(frac_prolif_parameters)):
    #setting parameters
    #index = int(frac_prolif_parameters.sort_values(by=['cost'])[0:100]['index'].tolist()[i])
    #or
    index = int(frac_prolif_parameters[:]['index'].tolist()[i])
    parameter_vector=l2[index-1]
    
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
    ffts_obsCdc20.drop(index, axis=1)
    ffts_obsCycE.drop(index, axis=1)
    ffts_obsCycD.drop(index, axis=1)
    
    
#########################################################################




frac_prolif_parameters.to_pickle('C:\\Users\\Grim.M\\Ecole\\ens\\Master\\S2\\Stage\\Code\\fit_data\\sobol_seq\\cost_function_model_4\\cost_and_IMT_distribution\\frac_prolif_parameters_all_parameters_50_initial_cells.pkl')

IMT_distribution.to_pickle('C:\\Users\\Grim.M\\Ecole\\ens\\Master\\S2\\Stage\\Code\\fit_data\\sobol_seq\\cost_function_model_4\\cost_and_IMT_distribution\\IMT_distribution_all_parameters_50_initial_cells.pkl')

ffts_obsCdc20.to_pickle('C:\\Users\\Grim.M\\Ecole\\ens\\Master\\S2\\Stage\\Code\\fit_data\\sobol_seq\\cost_function_model_4\\fft\\stochastic_simulation\\all_ffts_obsCdc20.pkl')

ffts_obsCycE.to_pickle('C:\\Users\\Grim.M\\Ecole\\ens\\Master\\S2\\Stage\\Code\\fit_data\\sobol_seq\\cost_function_model_4\\fft\\stochastic_simulation\\all_ffts_obsCycE.pkl')

ffts_obsCycD.to_pickle('C:\\Users\\Grim.M\\Ecole\\ens\\Master\\S2\\Stage\\Code\\fit_data\\sobol_seq\\cost_function_model_4\\fft\\stochastic_simulation\\all_ffts_obsCycD.pkl')



"""

##################################"
#loading and concatenate
    
#load all the cost calculated:
x=np.linspace(0,120,360)


not_found=[]
index__=[]
index_ = []
for idx in index:
    print(idx)
    try :
        pd.read_pickle(path + "/fit_data/sobol_seq/cost_function/cost_and_IMT_distribution/parameters_" + str(idx) + "_50_initial_cells.pkl")
        print("found")
        index_+=["Tempo/all_events_simulated_" + str(idx) + "_50_initial_cells.pkl"]
        index__+=[idx]
    except:
        print('not found')


IMT_distribution = pd.DataFrame(x)  
   

frac_prolif_parameters = {"index":[], "a":[], "d":[], "q":[], "cost":[], "cost_mmd":[], "cost_frac_prolif_model":[]}
frac_prolif_parameters = pd.DataFrame(data=frac_prolif_parameters)
 
for idx in index__:
    
    y = np.load(path + r"/fit_data/sobol_seq/cost_function/cost_and_IMT_distribution/IMT_distribution_" + str(idx) + "_50_initial_cells.npy")
    
    IMT_distribution[idx] = pd.Series(y, index=IMT_distribution.index)
    
    frac_prolif_parameters_temp = pd.read_pickle(path + "/fit_data/sobol_seq/cost_function/cost_and_IMT_distribution/parameters_" + str(idx) + "_50_initial_cells.pkl")
    
    frames = [frac_prolif_parameters, frac_prolif_parameters_temp]
    frac_prolif_parameters = pd.concat(frames)
    
frac_prolif_parameters.to_pickle(path + "/fit_data/sobol_seq/cost_function/cost_and_IMT_distribution/frac_prolif_parameters_all_parameters_50_initial_cells.pkl")

IMT_distribution.to_pickle(path + "/fit_data/sobol_seq/cost_function/cost_and_IMT_distribution/IMT_distribution_all_parameters_50_initial_cells.pkl")

"""
#loading the whole dataset
path=r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code"

frac_prolif_parameters = pd.read_pickle(path + r"/fit_data/sobol_seq/cost_function_model_4/cost_and_IMT_distribution/frac_prolif_parameters_all_parameters_50_initial_cells.pkl")

IMT_distribution = pd.read_pickle(path + r"/fit_data/sobol_seq/cost_function_model_4/cost_and_IMT_distribution/IMT_distribution_all_parameters_50_initial_cells.pkl")





#
sum(frac_prolif_parameters["cost_mmd"]<inf)

sum(frac_prolif_parameters["cost_mmd"]==50)



#comapring original parameters and others

a = np.sort(frac_prolif_parameters['cost'].tolist())
plt.bar(range(len(a)), a)
plt.plot([frac_prolif_parameters['cost'].loc[frac_prolif_parameters['index']==1] for n in range(len(a))], 'red')


sum(frac_prolif_parameters['cost'].tolist()<frac_prolif_parameters['cost'].iloc[0])

"""


    #adding the new columns
    all_events_best_temp_2['index'] = pd.Series([idx for i in range(len(all_events_best_temp_2))], index=all_events_best_temp_2.index)
    all_events_best_temp_2['cost'] = pd.Series([cost for i in range(len(all_events_best_temp_2))], index=all_events_best_temp_2.index)
    all_events_best_temp_2['cost_mmd'] = pd.Series([cost_mmd for i in range(len(all_events_best_temp_2))], index=all_events_best_temp_2.index)
    all_events_best_temp_2['cost_frac_prolif_model'] = pd.Series([cost_frac_prolif_model for i in range(len(all_events_best_temp_2))], index=all_events_best_temp_2.index)
    


    #merging the dataframes
    frames = [all_events_best, all_events_best_temp_2]
    all_events_best = pd.concat(frames)
"""



##########################################################################
#original parameters 

"""
all_events_best_temp = pd.read_pickle(path + "/fit_data/sobol_seq/cost_function/events_for_best_parameters/Tempo/all_events_simulated_0_50_initial_cells.pkl")
all_events_best_original_parameters = all_events_best_temp[np.array(all_events_best_temp['Birthtime']<30)]


#Calculating cost for fractional proliferation model
mitotic_lifespans = all_events_best_original_parameters['Lifespan'][all_events_best_original_parameters['event']=='mitosis']
censored_lifespans = all_events_best_original_parameters['Lifespan'][all_events_best_original_parameters['event']=='End_of_Expt']
parameters_sim_original_parameters={}
number_of_deaths = len(all_events_best_temp_2) - (all_events_best_temp_2['event']=='apoptosis').value_counts()[0]

cell_time_total = sum(all_events_best_temp_2['Lifespan'])
m = fp.mle(mitotic_lifespans, censored_lifespans)

parameters_sim_original_parameters['0.6 kPa DMSO'] = {}
parameters_sim_original_parameters['0.6 kPa DMSO']['a'] = number_of_deaths/cell_time_total
parameters_sim_original_parameters['0.6 kPa DMSO']['mu'] = m['x'][0]
parameters_sim_original_parameters['0.6 kPa DMSO']['sigma'] = m['x'][1]
parameters_sim_original_parameters['0.6 kPa DMSO']['lambda'] = m['x'][2]
parameters_sim_original_parameters['0.6 kPa DMSO'].update(fp.q_rates(m))

cost_frac_prolif_model_original_parameters = 0
for name in parameter_names:
    print(name)
    if parameters_exp['0.6 kPa DMSO'][name]<10**-5:
        parameters_exp['0.6 kPa DMSO'][name] = 1e-5
    if parameters_sim_original_parameters['0.6 kPa DMSO'][name]<10**-5:
        parameters_sim_original_parameters['0.6 kPa DMSO'][name] = 1e-5
    
    cost_frac_prolif_model_original_parameters += abs(np.log(parameters_sim_original_parameters['0.6 kPa DMSO'][name]) - np.log(parameters_exp['0.6 kPa DMSO'][name]))
    

    
#comparing those parameters and mmd:
number_of_IMT_to_sample = min (len(mitotic_lifespans), len(mitotic_lifespans_exp))
number_of_IMT_to_sample = 20
simulated_IMT = random.sample(mitotic_lifespans.to_list(), number_of_IMT_to_sample)
experimental_IMT = random.sample(mitotic_lifespans_exp.to_list(), number_of_IMT_to_sample)

cost_mmd_original_parameters = nmmd2_rbf(np.array([simulated_IMT]).transpose(), np.array([experimental_IMT]).transpose(), 1) #difference between IMT distribution

#calculating IMT densities

IMT_distribution_original_parameters = [RBS(i, np.transpose([mitotic_lifespans.tolist()]), sigma) for i in x]


#adding the costs
cost_original_parameters = cost_mmd_original_parameters + cost_frac_prolif_model_original_parameters



#############################################
"""

x=np.linspace(0,120,360)



colormap = cm.RdYlBu
normalize = mcolors.Normalize(vmin = min(frac_prolif_parameters["cost_mmd"]), vmax = max(frac_prolif_parameters["cost_mmd"]))
#or 
normalize = mcolors.Normalize(vmin = min(frac_prolif_parameters["cost_mmd"]), vmax = frac_prolif_parameters.loc[frac_prolif_parameters['index']==1].iloc[0]['cost_mmd'])


fig = pl.figure(figsize=(6,4),dpi=200)
ax = fig.add_subplot(111)
for idx in frac_prolif_parameters['index']:
    idx = int(idx)
    if frac_prolif_parameters.loc[frac_prolif_parameters['index']==idx].iloc[0]['cost_mmd']< (frac_prolif_parameters.loc[frac_prolif_parameters['index']==1].iloc[0]['cost_mmd']-6) : #or inf
    #if frac_prolif_parameters.loc[frac_prolif_parameters['index']==idx].iloc[0]['cost_mmd']< inf : #or inf
        ax.plot(x, IMT_distribution[idx].tolist()
,color=colormap(normalize(frac_prolif_parameters.loc[frac_prolif_parameters['index']==idx].iloc[0]['cost_mmd'])), lw=0.4)
        ax.grid(True)
        ax.tick_params(axis="y", labelsize=8)
        ax.tick_params(axis="x", labelsize=8)

#plotting original parameters and expt
ax.plot(x, experimental_IMT_distribution, color='black', lw=2)
ax.plot(x, IMT_distribution[1].tolist(), color=colormap(normalize(frac_prolif_parameters.loc[frac_prolif_parameters['index']==1.0].iloc[0]['cost_mmd'])), lw=2)

plt.subplots_adjust(bottom=0.1, right=0.075, left=0.07, top=0.9, wspace=0.2, hspace=0.4)
scalarmappaple = cm.ScalarMappable(norm=normalize, cmap=colormap)
scalarmappaple.set_array([])
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
cb = fig.colorbar(scalarmappaple, cax=cbar_ax)
cb.set_label(label='cost_mmd', fontsize=12)
cb.ax.tick_params(labelsize=12)


#fig.savefig(r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\report\figures\sobol_seq\sobol_seq_cost_function\cost_mmd.pdf")

fig.savefig(r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\report\figures\sobol_seq\sobol_seq_cost_function\cost_mmd_better_than_original_4.pdf")



########################################################################


colormap = cm.rainbow

normalize = mcolors.Normalize(vmin=min(frac_prolif_parameters["cost_frac_prolif_model"])
, vmax=max(frac_prolif_parameters["cost_frac_prolif_model"]))

fig = pl.figure(figsize=(6.5,4),dpi=200)
ax = fig.add_subplot(111)
for idx in frac_prolif_parameters['index']:
    idx = int(idx)
    frac_prolif_parameters_temp = frac_prolif_parameters.loc[frac_prolif_parameters['index']==idx]
    parameters_temp = frac_prolif_parameters_temp[['a','d','q']].iloc[0].as_matrix()
    
    ax.plot(parameter_names, np.log10(parameters_temp), color=colormap(normalize(frac_prolif_parameters_temp.iloc[0]['cost_frac_prolif_model'])), lw=0.4)
    ax.grid(True)
    ax.tick_params(axis="y", labelsize=8)
    ax.tick_params(axis="x", labelsize=8)

#plotting original parameters and expt
ax.plot(parameter_names,  np.log10([parameters_exp['0.6 kPa DMSO']['a'], parameters_exp['0.6 kPa DMSO']['d'], parameters_exp['0.6 kPa DMSO']['q']]), color='black', lw=2)
ax.plot(parameter_names, np.log10(frac_prolif_parameters.loc[frac_prolif_parameters['index']==1.0][['a','d','q']].iloc[0].as_matrix()), color=colormap(normalize(frac_prolif_parameters.loc[frac_prolif_parameters['index']==1.0].iloc[0]['cost_frac_prolif_model'])), lw=2)

plt.subplots_adjust(bottom=0.1, right=0.075, left=0.07, top=0.9, wspace=0.2, hspace=0.4)
scalarmappaple = cm.ScalarMappable(norm=normalize, cmap=colormap)
scalarmappaple.set_array([])
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
cb = fig.colorbar(scalarmappaple, cax=cbar_ax)
cb.set_label(label='cost_frac_prolif_model', fontsize=12)
cb.ax.tick_params(labelsize=10)


fig.savefig(r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\report\figures\sobol_seq\sobol_seq_cost_function\cost_frac_prolif_model.pdf")

#fig.savefig(r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\report\figures\sobol_seq\sobol_seq_cost_function\cost_frac_prolif_model_better_than_original_4.pdf")





##################################################################################
#Combining best apoptosis rate and best mmd cost

l = sobol_seq.i4_sobol_generate(30, 100000, skip = 0)
l2 = 4*l + (np.array(loc)-2)

frac_prolif_parameters[frac_prolif_parameters['cost_mmd'] < frac_prolif_parameters['cost_mmd'].min()+1]
    
#best mmd cost
idx = int(frac_prolif_parameters[frac_prolif_parameters['cost_mmd'] == frac_prolif_parameters['cost_mmd'].min()]['index'].iloc[0] -1)


frac_prolif_parameters.sort_values(by=['cost'])[0:15]['index'] 
idx = 44010 - 1
idx = 80080 - 1
idx= 74204 -1

parameter_vector_best_mmd = [10**i for i in l2[idx]]

print(parameter_vector_best_mmd)


#best apoptosis rate
best_idx=np.array([[1,1],[1,1]])

for idx in frac_prolif_parameters['index']:
    idx=int(idx)
    frac_prolif_parameters_temp =  frac_prolif_parameters.loc[frac_prolif_parameters['index']==idx]

    best_idx = np.concatenate([np.array([[abs(parameters_exp['0.6 kPa DMSO']['a']-frac_prolif_parameters_temp.iloc[0]['a']), idx]]), best_idx], axis=0)
    
parameter_vector_best_apoptosis_rate = [10**i for i in l2[int(best_idx[np.argmin(best_idx[:,0]),1])]]

idx=int(best_idx[np.argmin(best_idx[:,0]),1])


#combination and simulation

parameter_vector = parameter_vector_best_mmd[0:20] + parameter_vector_best_apoptosis_rate[20:]

pysb_sampled_parameter_names = ["V_s2skp","V_Fp","k_M","k_EMY","k_CDF","k_CDS","k_ED","K_S","K_MY","K_DF","d_CD","k_P1", "k_P2", "K_CD", "k_E", "K_M", "K_CE", "V_2E2F", "K_2E2F", "k_RE", "kf1", "kf2", "kb", "kasbaxmbcl2", "kdsbaxmbcl2", "k1sbh3", "k2sbh3", "kdbh3", "kasbh3bcl2", "kdsbh3bcl2"]
param_dict = {pname: pvalue for pname, pvalue in zip(pysb_sampled_parameter_names, parameter_vector)}
for pname, pvalue in param_dict.items():
    model_population.parameters[pname].value = (pvalue)
    model_concentration.parameters[pname].value = (pvalue)
    


(cells, cell_count) = runAKTERK(init, run_time, observed_cells_time, model_concentration, model_population)
IMT_eoe_counts = retrieve_cellState(run_time, init, cells)
all_events_simulated = all_events(run_time, IMT_eoe_counts)
all_events_simulated.to_pickle(path + "/fit_data/all_events_simulated_combined_parameters_50_initial_cells.pkl")



############################################################################


all_events_combined_parameters = pd.read_pickle(path + "/fit_data/sobol_seq/cost_function/events_for_best_parameters/all_events_simulated_combined_parameters_50_initial_cells.pkl")
all_events_combined_parameters = all_events_combined_parameters[np.array(all_events_combined_parameters['Birthtime']<30)]


#Calculating cost for fractional proliferation model
mitotic_lifespans = all_events_combined_parameters['Lifespan'][all_events_combined_parameters['event']=='mitosis']
censored_lifespans = all_events_combined_parameters['Lifespan'][all_events_combined_parameters['event']=='End_of_Expt']
parameters_combined_parameters = {}
number_of_deaths = len(all_events_combined_parameters) - (all_events_combined_parameters['event']=='apoptosis').value_counts()[0]

cell_time_total = sum(all_events_combined_parameters['Lifespan'])
m = fp.mle(mitotic_lifespans, censored_lifespans)

parameters_combined_parameters['0.6 kPa DMSO'] = {}
parameters_combined_parameters['0.6 kPa DMSO']['a'] = number_of_deaths/cell_time_total
parameters_combined_parameters['0.6 kPa DMSO']['mu'] = m['x'][0]
parameters_combined_parameters['0.6 kPa DMSO']['sigma'] = m['x'][1]
parameters_combined_parameters['0.6 kPa DMSO']['lambda'] = m['x'][2]
parameters_combined_parameters['0.6 kPa DMSO'].update(fp.q_rates(m))

cost_frac_prolif_model_combined_parameters = 0
for name in parameter_names:
    print(name)
    if parameters_exp['0.6 kPa DMSO'][name]<10**-5:
        parameters_exp['0.6 kPa DMSO'][name] = 1e-5
    if parameters_combined_parameters['0.6 kPa DMSO'][name]<10**-5:
        parameters_combined_parameters['0.6 kPa DMSO'][name] = 1e-5
    
    cost_frac_prolif_model_combined_parameters += abs(np.log(parameters_combined_parameters['0.6 kPa DMSO'][name]) - np.log(parameters_exp['0.6 kPa DMSO'][name]))
    

    
#comparing those parameters and mmd:
number_of_IMT_to_sample = 20
simulated_IMT = random.sample(mitotic_lifespans.to_list(), number_of_IMT_to_sample)
experimental_IMT = random.sample(mitotic_lifespans_exp.to_list(), number_of_IMT_to_sample)

cost_mmd_combined_parameters = nmmd2_rbf(np.array([simulated_IMT]).transpose(), np.array([experimental_IMT]).transpose(), 1) #difference between IMT distribution

#calculating IMT densities

IMT_distribution_combined_parameters = [RBS(i, np.transpose([mitotic_lifespans.tolist()]), sigma) for i in x]


#adding the costs
cost_combined_parameters = cost_mmd_combined_parameters + cost_frac_prolif_model_combined_parameters


#compare with others

print(min(frac_prolif_parameters['cost_mmd']))
print(cost_mmd_combined_parameters)

print(min(frac_prolif_parameters['cost_frac_prolif_model']))
print(cost_frac_prolif_model_combined_parameters)

print(min(frac_prolif_parameters['cost']))
print(cost_combined_parameters)


print(parameter_vector)

##################################################################################
#multivariate normal probability from original parameters:

from scipy.stats import multivariate_normal



mean = l2[0]
cov = 4 * np.identity(30)

rv = multivariate_normal(mean, cov)
rv.pdf(l2[1])


sample = rv.rvs()


#or euclidian distance?

np.linalg.norm(l2[0]-l2[1])

#################################################################################

from sklearn import neighbors


##############################################################################


import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import ListedColormap
from sklearn import neighbors, datasets

n_neighbors = 5



l2
target = np.array([0 for i in range(120000)])
target[indices]=1


clf = neighbors.KNeighborsClassifier(n_neighbors, weights=weights)
clf.fit(l2, target)
    

predict=[]
for n in range(1000):
    sample = rv.rvs()

    predict+=[clf.predict([sample])]

    
    
h = .02  # step size in the mesh

# Create color maps
cmap_light = ListedColormap(['orange', 'cyan', 'cornflowerblue'])
cmap_bold = ['darkorange', 'c', 'darkblue']



###################################################################################
#Umap

frac_prolif_parameters[frac_prolif_parameters["cost_mmd"]!=50]


indices = [int(i) for i in frac_prolif_parameters["index"].tolist()]
#or

indices = [int(i) for i in frac_prolif_parameters[frac_prolif_parameters["cost_mmd"]!=50]["index"].tolist()]

indices = [int(i)-1 for i in frac_prolif_parameters["index"].tolist()]

#normal umap


l=sobol_seq.i4_sobol_generate(30, 120000, skip=0)
l2 = 4*l + (np.array(loc)-2)

parameters_log = np.array([l2[i] for i in indices])
parameters_log_not_working = np.array([l2[i] for i in range(len(l2)) if i not in indices])

parameters = 10**parameters_log
all_parameters = 10**l2

parameters_not_working = 10**parameters_log_not_working



from sklearn import decomposition
import umap as umap


#################################
clusters= np.array([y_pred[ffts_obsCdc20.columns==ind+1] for ind in indices]).T[0]

#number of clusters
N=4
# define the colormap
cmap = plt.cm.jet
# extract all colors from the .jet map
cmaplist = [cmap(i) for i in range(cmap.N)]
# create the new map
cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)

# define the bins and normalize
bounds = np.linspace(0,N,N+1)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)


####################################

#umap   
trans = umap.UMAP().fit(parameters)
#trans = umap.UMAP().fit(all_parameters)
X_umap = trans.transform(parameters)
X_umap_2 = trans.transform(parameters_not_working)

###
trans = umap.UMAP().fit(parameters_log)
#trans = umap.UMAP().fit(l2)
X_umap = trans.transform(parameters_log)
X_umap_2 = trans.transform(parameters_log_not_working)


###########################################
fig = pl.figure(figsize=(4,3),dpi=200)
ax = fig.add_subplot(111)
sc = ax.scatter(X_umap[:,0], X_umap[:,1], c=clusters, alpha = 0.8, s = 1, cmap=cmap, norm=norm);



# create the colorbar
cb = pl.colorbar(sc, spacing='proportional',ticks=bounds)
cb.set_label('Custom cbar')
ax.set_title('Discrete color mappings')
pl.show()

pl.savefig(r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\report\figures\sobol_seq\sobol_seq_cost_function\Umap_all_parameters_clusters.pdf")

###########################################"
#figure
fig = pl.figure(figsize=(4,3),dpi=200)
ax = fig.add_subplot(111)

sc = ax.scatter(X_umap[:,0], X_umap[:,1], c=frac_prolif_parameters[frac_prolif_parameters["cost_mmd"]!=50]["cost"], cmap="RdYlBu", alpha = 0.8, s = 1);


#or

#sc = ax.scatter(X_umap[:,0], X_umap[:,1], c= np.log(rv.pdf(parameters_log)/max(rv.pdf(parameters_log))), cmap="brg", alpha = 0.8, s = 1);

#sc = ax.scatter(X_umap[:,0], X_umap[:,1], c = [np.linalg.norm(l2[0]-i) for i in parameters_log] , cmap="seismic", alpha = 0.8, s = 1);

ax.scatter(X_umap_2[:,0], X_umap_2[:,1], c="black", alpha = 0.1,s = 0.1);
 
#ax.scatter(test_pts[:,0], test_pts[:,1], c=np.array(cost)<inf, alpha = 1, cmap="gnuplot",s = 1.5);
pl.colorbar(sc)

pl.savefig(r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\report\figures\sobol_seq\sobol_seq_cost_function\Umap_all_parameters.pdf")
#pl.savefig(r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\report\figures\sobol_seq\sobol_seq_cost_function\Umap_working_parameters_test_points.pdf")
pl.savefig(r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\report\figures\sobol_seq\sobol_seq_cost_function\Umap_20_first_log_parameters_plus_cost_mmd_grid_check.pdf")




###########################################"
#adding the general cost as a dimension
parameters_plus_cost = np.concatenate([parameters, np.array([frac_prolif_parameters[frac_prolif_parameters["cost_mmd"]!=50]["cost"]]).T], axis=1)

parameters_log_plus_cost = np.concatenate([parameters_log, np.array([frac_prolif_parameters[frac_prolif_parameters["cost_mmd"]!=50]["cost"]]).T], axis=1)
#parameters_log_not_working_plus_cost = 

#umap   
trans = umap.UMAP().fit(parameters_plus_cost)
#trans = umap.UMAP().fit(all_parameters)
X_umap = trans.transform(parameters_plus_cost)
#X_umap_2 = trans.transform(parameters_not_working)

###
trans = umap.UMAP().fit(parameters_log_plus_cost)
#trans = umap.UMAP().fit(l2)
X_umap = trans.transform(parameters_log_plus_cost)
#X_umap_2 = trans.transform(parameters_log_not_working)



######################################
#20 first parameters
parameters_log_20 = parameters_log[:,0:20]
parameters_log_not_working_20 = parameters_log_not_working[:,0:20]


parameters_20 = parameters[:,0:20]
parameters_not_working_20 = parameters_not_working[:,0:20]


trans = umap.UMAP().fit(parameters_log_20)
X_umap = trans.transform(parameters_log_20)
X_umap_2 = trans.transform(parameters_not_working_20)


trans = umap.UMAP().fit(parameters_20)
X_umap = trans.transform(parameters_20)
X_umap_2 = trans.transform(parameters_not_working_20)



######################################
#20 first parameters and cost

parameters_log_20_plus_cost = np.concatenate([parameters_log_20, np.array([frac_prolif_parameters[frac_prolif_parameters["cost_mmd"]!=50]["cost"]]).T], axis=1)

parameters_20_plus_cost = np.concatenate([parameters_20, np.array([frac_prolif_parameters[frac_prolif_parameters["cost_mmd"]!=50]["cost"]]).T], axis=1)



trans = umap.UMAP().fit(parameters_log_20_plus_cost)
X_umap = trans.transform(parameters_log_20_plus_cost)

trans = umap.UMAP().fit(parameters_20_plus_cost)
X_umap = trans.transform(parameters_20_plus_cost)


######################################
#20 first parameters and mmd cost

parameters_log_20_plus_cost = np.concatenate([parameters_log_20, np.array([frac_prolif_parameters[frac_prolif_parameters["cost_mmd"]!=50]["cost_mmd"]]).T], axis=1)

parameters_20_plus_cost = np.concatenate([parameters_20, np.array([frac_prolif_parameters[frac_prolif_parameters["cost_mmd"]!=50]["cost_mmd"]]).T], axis=1)



trans = umap.UMAP().fit(parameters_log_20_plus_cost)
X_umap = trans.transform(parameters_log_20_plus_cost)

trans = umap.UMAP().fit(parameters_20_plus_cost)
X_umap = trans.transform(parameters_20_plus_cost)



########################################

corners = np.array([
    [9, 11],  # 
    [9, 13],  # 
    [11, 13],  # 
    [11, 11],  # 
])
   

    
test_pts = np.array([
    (corners[0]*(1-x) + corners[1]*x)*(1-y) +
    (corners[3]*(1-x) + corners[2]*x)*y
    for y in np.linspace(0, 1, 10)
    for x in np.linspace(0, 1, 10)
])

   
#eliminate negative rates
    
inv_transformed_points = trans.inverse_transform(test_pts)


a = [vector for (_, vector) in enumerate(inv_transformed_points) if np.all(vector>0)]
test_pts = np.array([test_pts[position].tolist() for (position, vector) in enumerate(inv_transformed_points) if np.all(vector>0)])


np.array([parameter_vector_best_apoptosis_rate[20:]])
inv_transformed_points = 10**inv_transformed_points

a = [np.concatenate([inv_transformed_points[i,0:20], np.array(parameter_vector_best_apoptosis_rate[20:])],axis=0) for i in range(len(inv_transformed_points))]


np.save(path + "/fit_data/inv_transformed_points_log___.npy", a)

np.save(path + "/fit_data/inv_transformed_points_log.npy", inv_transformed_points)
np.save(path + "/fit_data/inv_transformed_points_2.npy", inv_transformed_points)


inv_transformed_points_log = np.load(path+"/fit_data/inv_transformed_points_log.npy")
inv_transformed_points = np.load(path+"/fit_data/inv_transformed_points_2.npy")


cost=[]

for parameter_vector in a :
    print(parameter_vector)
    c = cost_function_reduced_quiescence_DMSO_only_30_parameters(parameter_vector)
    print(cost)
    cost+=[c]

sum(np.array(cost)<inf)

#########################################################################



path=r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code"

frac_prolif_parameters = pd.read_pickle(path + "/fit_data/sobol_seq/cost_function_model_4/cost_and_IMT_distribution/frac_prolif_parameters_all_parameters_50_initial_cells.pkl")


indices = [int(i)-1 for i in frac_prolif_parameters["index"].tolist()]


indices = [int(i)-1 for i in frac_prolif_parameters["index"].tolist()]

parameters_log

parameters

colors = sns.color_palette(n_colors=30)


from pathlib import Path
#fn = Path(r'C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\report\figures\pydream\bis_log_quiescence_DMSO_only_21_parameters\PyDREAM_parameter_distribution.pdf').expanduser()

a=6
b=5
fig, axs = pl.subplots(ncols=a, nrows=b, sharex=True, sharey=True, figsize=(15,6))
for dim, ax in zip(range(30), axs.flat):
        sns.distplot(parameters_log[:, dim], color=colors[dim], ax=ax, kde_kws={"lw": 0.5, "label": "variable: " + pysb_sampled_parameter_names[dim]}, bins=np.linspace(-4,1,50))
        ax.legend(loc='upper left', fontsize = 6)
        
        
#pl.savefig(fn, bbox_inches='tight')

