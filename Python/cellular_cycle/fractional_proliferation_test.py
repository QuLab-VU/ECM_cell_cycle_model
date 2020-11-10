# -*- coding: utf-8 -*-
"""
Created on Mon Aug 17 10:03:17 2020

@author: Grim.M
"""


path=r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code"
import sys
sys.path.insert(0, path + "\cellular_cycle") #add a pathway to the model

from runAKTERK_2 import runAKTERK
from retrieve_cellState_dw import retrieve_cellState
from plot_frac_div_q_ap1 import all_events
import numpy as np
import pandas as pd
from pysb.simulator import BngSimulator
from pysb.simulator import ScipyOdeSimulator
import matplotlib.pyplot as plt
import pylab as pl
import fracprolif as fp


init_cell = 8
init = init_cell
run_time= 120
omega=200
from Austin_proliferation_model_population_3 import model
model_population=model
from Austin_proliferation_model_concentration_3 import model
model_concentration=model

###################################


condition= '0.6 kPa DMSO'

pysb_sampled_parameter_names = ["V_s2skp","V_Fp","k_M","k_EMY","k_CDF","k_CDS","k_ED","K_S","K_MY","K_DF"]

pysb_sampled_parameter_names = ["V_s2skp","V_Fp","k_M","k_EMY","k_CDF","k_CDS","k_ED","K_S","K_MY","K_DF","d_CD","k_P1", "k_P2", "K_CD", "k_E", "K_M", "K_CE", "V_2E2F", "K_2E2F", "k_RE", "k_b"]


#parameter_vector = np.array([4.00110767e-04, 8.50439757e-03, 4.49244718e-03, 7.62739688e-02, 5.72752269e-02, 5.91004292e-03, 4.58566275e-01, 1.03895467e-01, 4.17757245e-04, 2.99759014e-02]) #check_quiescence

#parameter_vector = np.array([2.88095027e-02, 1.08653936e-04, 1.06190007e-03, 5.91535448e-02, 9.63052769e-04, 1.12765757e-04, 3.35453630e-02, 1.18449277e-04, 1.24550538e-01, 2.00076132e-03, 1.28792094e-01, 6.23484593e-03, 8.39769019e-02, 2.02874490e-02, 2.16102188e+00, 6.79877479e-04, 1.18330213e-03, 8.82600428e+00, 3.06767094e-04, 1.46608884e-04, 3.84589838e-02])

#parameter_vector = np.array([1.57794006e-02, 1.39278870e-03, 4.12809683e+00, 1.57724527e-03, 2.75334197e-01, 4.35474584e-03, 3.11591117e-03, 2.61573101e-01, 1.50878679e+00, 2.15940417e-01, 3.60595546e-01, 9.69991577e-02, 5.06821907e-02, 9.53477326e-03, 1.39697557e+00, 2.87024812e-04, 5.77630385e-04, 5.25778682e+00, 7.90494665e-04, 2.46037958e-01, 2.96088698e-02])

#parameter_vector = np.array([6.24703998e-03, 2.91053884e-02, 2.67601660e-01, 4.93563624e-04, 1.31082069e-01, 1.02348783e-02, 2.37892797e-04, 2.54393252e+00, 7.94080174e-01, 1.53741761e-03, 1.28451604e+00, 1.90780716e+00, 1.28684268e-01, 4.11319145e+00, 2.26691878e+00, 7.93580063e-03, 9.71416614e-04, 7.60307733e+00, 1.56310211e-02, 1.46211690e-04, 8.18644759e-03])

#after quiescent cell checking
#parameter_vector = np.array([8.85976504e-02, 3.34310621e-04, 8.13978039e+00, 3.96451559e+00, 2.37123607e+00, 2.65445804e-02, 2.20803318e+00, 4.26345643e+00, 3.02795773e-04, 9.35211594e+00, 1.46166716e+00, 6.81092439e+00, 2.03478846e-03, 1.70050765e-04, 2.57055792e+00, 8.27820344e+00, 1.92676187e-02, 5.74259865e+00, 8.03687914e-03, 2.64191854e-01, 4.76005540e-02]) #A

#parameter_vector = np.array([2.98470799e-02, 1.35792858e-01, 5.02069869e-01, 1.02498893e+00, 5.56769488e-01, 4.05491396e-03, 2.51495467e-03, 6.38434298e+00, 5.72584015e-02, 3.12471250e-03, 3.52223372e-03, 7.75525775e-04, 3.29290719e-01, 3.14627969e-02, 1.45356753e+00, 1.98845982e-04, 8.52356219e-03, 6.62128006e+00, 2.54615009e-04, 6.37370204e-04, 5.12835520e-02]) #B

#parameter_vector = np.array([1.39053450e-02, 2.02916787e-02, 9.82700658e-01, 6.74432358e-01, 1.39102146e-01, 1.24143139e-03, 3.07128989e-03, 2.12392683e+00, 1.93583408e-01, 1.23987327e-03, 1.35932103e-03, 7.74671661e-04, 6.11135888e-01, 3.69925277e-01, 9.16830465e-01, 3.46726549e-04, 5.73744770e-03, 4.00964819e+00, 5.36818197e-04, 3.84433956e-04, 4.92648527e-02]) #C

       
parameter_vector = np.array([2.56569681e-04, 3.32574208e+00, 3.85314581e+00, 3.49511021e+00, 4.06194451e-04, 2.74040574e-03, 1.32515877e-03, 2.12705099e-01, 8.65564214e-02, 8.63985580e-01, 3.99766315e+00, 2.94834405e-04, 1.80765381e-01, 5.52830805e-01, 7.90070207e-02, 1.81336116e-01, 7.17570453e-04, 1.31294440e+00, 2.47376238e-01, 5.79824482e-04, 3.87295191e-02])


param_dict = {pname: pvalue for pname, pvalue in zip(pysb_sampled_parameter_names, parameter_vector)}  
for pname, pvalue in param_dict.items():  
    #Change model parameter values to current location in parameter space
    model_population.parameters[pname].value = (pvalue)
    model_concentration.parameters[pname].value = (pvalue)
    
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


#test des param
tmp_t = pl.linspace(0, 500, 10000)

#A simple simulation should work
sim = BngSimulator(model_population, tspan=tmp_t, verbose=True)
traj = sim.run()
plt.plot(traj.observables['obsCdc20'])


sim = ScipyOdeSimulator(model_concentration, tspan=tmp_t, verbose=True)
traj = sim.run()
plt.plot(traj.observables['obsCdc20'])



observed_cells_time = 30
init = 30

(cells, cell_count) = runAKTERK(init, run_time, observed_cells_time, model_concentration, model_population)

#####
##simulation w Austin parameters
#
#init = 5
#run_time= 120.3
#omega=200
#from Austin_proliferation_model_population_3 import model
#model_population=model
#from Austin_proliferation_model_concentration_3 import model
#model_concentration=model
#
#import mat73
#data_dict = mat73.loadmat(path+'/v1/params_8_14_20.mat')
#new_param=data_dict['p']
#del new_param['GF']
#del new_param['Stress']
#
#for key in new_param.keys():
#    model_population.parameters[key].value = new_param[key]
#    model_concentration.parameters[key].value = new_param[key]
#       
#model_population.parameters['FAK'].value = (new_param['FAK']*200)
#model_population.parameters['ERK'].value = (new_param['ERK']*200)
#model_population.parameters['AKT'].value = (new_param['AKT']*200)
#
#model_population.parameters['s'].value = (0.41)
#model_concentration.parameters['s'].value = (0.41)
#

#####################################################################################################################
#Model


IMT_eoe_counts = retrieve_cellState(run_time, init, cells)

all_events_simulated = all_events(run_time, IMT_eoe_counts)

all_events_simulated_2 = all_events_simulated[all_events_simulated['Birthtime']<30]


import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path
fn = Path(r'C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\report\figures\fractional_proliferation\frac_div_sim_data_06kpaDMSO_over_30.pdf').expanduser()

grid = sns.JointGrid(x='Birthtime', y='Lifespan', data=all_events_simulated)
g = grid.plot_joint(sns.scatterplot, hue='event', data=all_events_simulated)
#sns.kdeplot(all_events.loc[all_events['event']=='mitosis', 'Birthtime'], ax=g.ax_marg_x, legend=False)
sns.kdeplot(all_events_simulated.loc[all_events_simulated['event']=='mitosis', 'Lifespan'], ax=g.ax_marg_y, vertical=True, legend=False, bw=4)
plt.ylim(0, 110)

pl.savefig(fn, bbox_inches='tight')



#############frac prolif model test

IMT_eoe_counts = retrieve_cellState(run_time, init, cells)
#ii=55

all_events_simulated = all_events(run_time, IMT_eoe_counts)

mitotic_lifespans = all_events_simulated['Lifespan'][all_events_simulated['event']=='mitosis']

censored_lifespans = all_events_simulated['Lifespan'][all_events_simulated['event']=='End_of_Expt']


m = fp.mle(mitotic_lifespans, censored_lifespans)
m['x']

v = fp.q_rates(m)

card_total = np.log2(cell_count[:,1])

d0=5
q0=0
Time_h = cell_count[:,0]


#half_imt, untreated_division_rate, a = fp.death_rate_without_delay(Time_h, card_total, v) 

number_of_deaths = len(all_events_simulated) - (all_events_simulated['event']=='apoptosis').value_counts()[0]

cell_time_total = sum(all_events_simulated['Lifespan'])

a = np.array([number_of_deaths / cell_time_total])

y = np.log2(cell_count[:,1])-np.log2(cell_count[:,1])[0]

plt.plot(Time_h, y)


model_total = np.array([np.log2(fp.split_tot(t,  v['d'], v['d'], v['q'], v['q'], a, a, a, a, 1, 0, 0))[0] for t in Time_h])

model_fq = np.array([fp.split_fq(t,  v['d'], v['d'], v['q'], v['q'], a, a, a, a, 1, 0, 0)[0] for t in Time_h]) * model_total

model_fd = np.array([fp.split_fd(t,  v['d'], v['d'], v['q'], v['q'], a, a, a, a, 1, 0, 0)[0] for t in Time_h]) * model_total


fn = Path(r'C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\report\figures\fractional_proliferation\frac_div_model_sim_data_06kpaDMSO_over_30.pdf').expanduser()

pl.savefig(fn, bbox_inches='tight')

plt.scatter(Time_h, y, s=40, facecolors='none', edgecolors='black')
plt.plot(Time_h, model_total, c="green")
plt.plot(Time_h, model_fq, c="red")
plt.plot(Time_h, model_fd, c="blue")
pl.savefig(fn, bbox_inches='tight')


###############################################################################################################"
#Exp

IMT_eoe_counts = retrieve_cellState(run_time, init, cells)
#ii=55

all_events_simulated = all_events(run_time, IMT_eoe_counts)

mitotic_lifespans = all_events_simulated['Lifespan'][all_events_simulated['event']=='mitosis']

censored_lifespans = all_events_simulated['Lifespan'][all_events_simulated['event']=='End_of_Expt']


m = fp.mle(mitotic_lifespans, censored_lifespans)
m['x']

v = fp.q_rates(m)

card_total = np.log2(cell_count[:,1])

d0=5
q0=0
Time_h = cell_count[:,0]


#half_imt, untreated_division_rate, a = fp.death_rate_without_delay(Time_h, card_total, v) 

number_of_deaths = len(all_events_simulated) - (all_events_simulated['event']=='apoptosis').value_counts()[0]

cell_time_total = sum(all_events_simulated['Lifespan'])

a = np.array([number_of_deaths / cell_time_total])

y = np.log2(cell_count[:,1])-np.log2(cell_count[:,1])[0]

plt.plot(Time_h, y)


model_total = np.array([np.log2(fp.split_tot(t,  v['d'], v['d'], v['q'], v['q'], a, a, a, a, 1, 0, 0))[0] for t in Time_h])

model_fq = np.array([fp.split_fq(t,  v['d'], v['d'], v['q'], v['q'], a, a, a, a, 1, 0, 0)[0] for t in Time_h]) * model_total

model_fd = np.array([fp.split_fd(t,  v['d'], v['d'], v['q'], v['q'], a, a, a, a, 1, 0, 0)[0] for t in Time_h]) * model_total


fn = Path(r'C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\report\figures\fractional_proliferation\frac_div_model_sim_data_06kpaDMSO_over_30.pdf').expanduser()

pl.savefig(fn, bbox_inches='tight')

plt.scatter(Time_h, y, s=40, facecolors='none', edgecolors='black')
plt.plot(Time_h, model_total, c="green")
plt.plot(Time_h, model_fq, c="red")
plt.plot(Time_h, model_fd, c="blue")
pl.savefig(fn, bbox_inches='tight')