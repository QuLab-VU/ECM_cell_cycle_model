# -*- coding: utf-8 -*-
"""
Created on Fri Aug 14 18:11:34 2020

@author: Grim.M
"""

import seaborn as sns
from matplotlib import pyplot as plt
from SMs_G1s_IMTs_simulation_BngSimulator_reduced import IMTs

from Austin_proliferation_model_population import model
model_population=model

from Austin_proliferation_model_concentration import model
model_concentration=model


pysb_sampled_parameter_names = ["V_s2skp","V_Fp","k_M","k_EMY","k_CDF","k_CDS","k_ED","K_S","K_MY","K_DF","d_CD","k_P1", "k_P2", "K_CD", "k_E", "K_M", "K_CE", "V_2E2F", "K_2E2F", "k_RE", "k_b"]



logps={}
sampled_params_dict={}

for chain in range(10):   
    logps_temp = np.array([])
    params_temp = np.array([])
    
    logps_temp = vars()['cellular_cycle_dreamzs_5chain_logps_chain_' + str(chain) + "_036500"]
    params_temp = vars()['cellular_cycle_dreamzs_5chain_sampled_params_chain_' + str(chain) + "_036500"]          
    logps['cellular_cycle_dreamzs_5chain_logps_chain_' + str(chain)] = logps_temp
    sampled_params_dict['cellular_cycle_dreamzs_5chain_sampled_params_chain_' + str(chain)] = params_temp



sample_parameters=[]
logps_tot=[]

for (logps_chain, chain) in zip(logps.keys(), sampled_params_dict.keys()):
    sample_parameters+=[sampled_params_dict[chain]]
    logps_tot+=[logps[logps_chain].T[0]]

GR = Gelman_Rubin(sample_parameters)

old_samples = sample_parameters


#####################################################################
total_iterations = len(old_samples[0])
burnin = int(total_iterations / 2)
samples = np.concatenate([old_samples[i][burnin:, :] for i in range(10)])
logps_tot_p = np.concatenate([logps_tot[i][burnin:] for i in range(10)])

#removing bad values :

np.where(logps_tot_p >-10000)
np.where(logps_tot_p <-10000)

samples=samples[np.where(logps_tot_p >-10000)[0],:]
logps_tot_p=logps_tot_p[np.where(logps_tot_p >-10000)[0]]


from pathlib import Path
fn = Path(r'C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\report\figures\pydream\bis_log_quiescence_DMSO_only_21_parameters\logps_distribution.pdf').expanduser()
sns.distplot(logps_tot_p)
pl.savefig(fn, bbox_inches='tight')


samples=samples[np.where(logps_tot_p >-64)[0],:]
logps_tot_p=logps_tot_p[np.where(logps_tot_p >-64)[0]]

condition='O.6 kPa DMSO'

if 'DMSO' in condition:
    if '0.6 kPa' in condition :
        model_population.parameters['ERK'].value = 1.06970165050038 * omega
        model_population.parameters['FAK'].value = 0.337143388171177 * omega
        model_population.parameters['AKT'].value = 0.616658111482869 * omega
    elif '4 kPa' in condition :
        model_population.parameters['ERK'].value = 0.897485095798471 * omega
        model_population.parameters['FAK'].value = 0.594150207445117 * omega
        model_population.parameters['AKT'].value = 1.01975420186633 * omega
    else : #55 kPa
        model_population.parameters['ERK'].value = 1 * omega
        model_population.parameters['FAK'].value = 1 * omega
        model_population.parameters['AKT'].value = 1 * omega
else : #ERL treated
    if '0.6 kPa' in condition :
        model_population.parameters['ERK'].value = 0.186745227418267 * omega
        model_population.parameters['FAK'].value = 0.397060244114162 * omega
        model_population.parameters['AKT'].value = 0.245702711121333 * omega
    elif '4 kPa' in condition :
        model_population.parameters['ERK'].value = 0.211453830597783 * omega
        model_population.parameters['FAK'].value = 0.709172297576571 * omega
        model_population.parameters['AKT'].value = 0.598518145949655 * omega
    else : #55 kPa
        model_population.parameters['ERK'].value = 0.273618423688386 * omega
        model_population.parameters['FAK'].value = 1.075820621939 * omega
        model_population.parameters['AKT'].value = 0.585067267459499 * omega


tmp_t = np.linspace(0, 1000, 10000)

accepted=[]

for sample in samples:
    parameter_vector = 10**sample   
    param_dict = {pname: pvalue for pname, pvalue in zip(pysb_sampled_parameter_names, parameter_vector)}
    for pname, pvalue in param_dict.items():
        #Change model parameter values to current location in parameter space
        model_population.parameters[pname].value = (pvalue)
    (IMT_tmp, a, b) = IMTs(tmp_t, model_population)
    if sum(np.array(IMT_tmp)>100) > 0:
        print(parameter_vector)
        accepted+=[parameter_vector]
    print(len(accepted))


np.save(path + "/fit_data/pydream/bis_log_DMSO_only_quiescence_21_parameters_check_quiescence/accepted_100.npy", accepted)

accepted_2=[]

for sample in accepted:
    parameter_vector = 10**sample   
    param_dict = {pname: pvalue for pname, pvalue in zip(pysb_sampled_parameter_names, parameter_vector)}
    for pname, pvalue in param_dict.items():
        #Change model parameter values to current location in parameter space
        model_population.parameters[pname].value = (pvalue)
    (IMT_tmp, a, b) = IMTs(tmp_t, model_population)
    if sum(np.array(IMT_tmp)>120) > 0:
        print(parameter_vector)
        accepted_2+=[parameter_vector]
    print(len(accepted_2))