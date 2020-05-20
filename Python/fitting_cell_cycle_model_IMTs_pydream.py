# -*- coding: utf-8 -*-
"""
Created on Mon Mar  2 14:13:16 2020

@author: Jérémie Marlhens
Title : fitting Austin's model
"""

path="/home/jeremie"
#path=r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code"
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

import math

from pysb.simulator import BngSimulator

from SMs_G1s_simulation_OpenCL import IMTs
from rbf_mmd import nmmd2_rbf

'''
#Initialize PySB solver object for running simulations.  Simulation timespan should match experimental data.
tmp_t = np.linspace(0, 1000, 10000)
sim = BngSimulator(model, tspan=tmp_t, verbose=True)
simres = sim.run(initials=Initial)
'''

#Load experimental data to which the model will be fit here.  The "experimental data" in this case is the IMTs distribution

experimental_SMs={}
experimental_G1s={}
experimental_IMTs={}



df = pd.read_csv(path+'/IMT_times/DS8_SM_0_6_4_55.csv') #path to be changed
dg = pd.read_csv(path+'/IMT_times/ab03aa03proliferationrates.xlsx') #path to be changed


for col in df.columns :
    experimental_SMs[col]=df[col].dropna().tolist()


df = pd.read_csv(path+'/IMT_times/DS8_G1_0_6_4_55.csv')

for col in df.columns :
    experimental_G1s[col]=df[col].dropna().tolist()

for k in experimental_G1s.keys() :
    minimum=min(len(experimental_G1s[k]),len(experimental_SMs[k]))
    experimental_IMTs[k]=[experimental_G1s[k][i]+experimental_SMs[k][i] for i in range(minimum)]

####################################################################################
    
#Cost function

from Austin_proliferation_model_concentration import model
model_concentration=model

from Austin_proliferation_model_population import model
model_population=model

from SMs_G1s_simulation_OpenCL import IMTs
from rbf_mmd import nmmd2_rbf

import os
os.environ["PYOPENCL_CTX"] ="0"

import pylab as pl

pysb_sampled_parameter_names = ["V_s2skp","V_Fp","k_M","k_EMY","k_CDF","k_CDS","k_ED","K_S","K_MY","K_DF"]


def cost_function(parameter_vector):
    
    print("les paramètres sont ")
    print(parameter_vector)
    
    global experimental_IMTs, model_concentration, model_population
    
    param_dict = {pname: pvalue for pname, pvalue in zip(pysb_sampled_parameter_names, parameter_vector)}
    
    tmp_t=pl.linspace(0, 500, 2500)
    
    for pname, pvalue in param_dict.items():
        
        #Change model parameter values to current location in parameter space
        model_concentration.parameters[pname].value = (pvalue)
        model_population.parameters[pname].value = (pvalue)
        
    #Simulate experimentally IMTs
    n_mmd2_tot=[]
    simulated_G1s=[]
    simulated_SMs=[]
    omega=200
    
    for condition in experimental_IMTs.keys():
        print(condition)
        
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
        simulated_G1s=[]
        simulated_SMs=[]
        
        while len(simulated_G1s)< len(experimental_G1s[condition]) or len(simulated_SMs)< len(experimental_SMs[condition]):
            (a,b)=IMTs(tmp_t,model_concentration, model_population)
            print(a)
            print(b)
            if a==math.inf:
                print("stop?")
                return 10000
            else:
                simulated_G1s+=b
                simulated_SMs+=a
        n_mmd2_tot+=[nmmd2_rbf(np.array([simulated_G1s[0:len(experimental_G1s[condition])]]).transpose(), np.array([experimental_G1s[condition]]).transpose(), 1)]
        n_mmd2_tot+=[nmmd2_rbf(np.array([simulated_SMs[0:len(experimental_SMs[condition])]]).transpose(), np.array([experimental_SMs[condition]]).transpose(), 1)]
    print("ok")
    return sum(n_mmd2_tot)


############################################################################"
parameter_vector=np.array([1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1])


cost_function(parameter_vector)




#Set upper and lower limits for uniform prior to be 3 orders of magnitude above and below original parameter values.
lower_limits   = np.array(original_params)  #a changer

lower_limits = [0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0]

parameters_to_sample = SampledParam(uniform, loc=lower_limits, scale=2) #loc is the lower bdry and scale is the size of the distribution

sampled_parameter_names = [parameters_to_sample]

niterations = 1000
converged = False
total_iterations = niterations
nchains = 6




if __name__ == '__main__':

    #Run DREAM sampling.  Documentation of DREAM options is in Dream.py.
    sampled_params, log_ps = run_dream(sampled_parameter_names, cost_function, niterations=niterations, nchains=nchains, multitry=False, gamma_levels=4, adapt_gamma=True, history_thin=1, model_name='cellular_cycle_dreamzs_5chain', verbose=True)

    #Save sampling output (sampled parameter values and their corresponding logps).
    for chain in range(len(sampled_params)):
        np.save('cellular_cycle_dreamzs_5chain_sampled_params_chain_'+str(chain)+'_'+str(total_iterations), sampled_params[chain])
        np.save('cellular_cycle_dreamzs_5chain_logps_chain_'+str(chain)+'_'+str(total_iterations), log_ps[chain])

    # Check convergence and continue sampling if not converged

    GR = Gelman_Rubin(sampled_params)
    print('At iteration: ', total_iterations, ' GR = ', GR)
    np.savetxt('cellular_cycle_dreamzs_5chain_GelmanRubin_iteration_' + str(total_iterations) + '.txt', GR)

    old_samples = sampled_params
    if np.any(GR > 1.2):
        starts = [sampled_params[chain][-1, :] for chain in range(nchains)]
        while not converged:
            total_iterations += niterations

            sampled_params, log_ps = run_dream(sampled_parameter_names, likelihood, start=starts, niterations=niterations, nchains=nchains, multitry=False, gamma_levels=4, adapt_gamma=True, history_thin=1, model_name='cellular_cycle_dreamzs_5chain', verbose=True, restart=True)

            for chain in range(len(sampled_params)):
                np.save('cellular_cycle_dreamzs_5chain_sampled_params_chain_' + str(chain) + '_' + str(total_iterations), sampled_params[chain])
                np.save('cellular_cycle_dreamzs_5chain_logps_chain_' + str(chain) + '_' + str(total_iterations),log_ps[chain])

            old_samples = [np.concatenate((old_samples[chain], sampled_params[chain])) for chain in range(nchains)]
            GR = Gelman_Rubin(old_samples)
            print('At iteration: ', total_iterations, ' GR = ', GR)
            np.savetxt('cellular_cycle_dreamzs_5chain_GelmanRubin_iteration_' + str(total_iterations)+'.txt', GR)

            if np.all(GR < 1.2):
                converged = True

    try:
        # Plot output
        import seaborn as sns
        from matplotlib import pyplot as plt

        total_iterations = len(old_samples[0])
        burnin = total_iterations / 2
        samples = np.concatenate((old_samples[0][burnin:, :], old_samples[1][burnin:, :], old_samples[2][burnin:, :], old_samples[3][burnin:, :], old_samples[4][burnin:, :]))

        ndims = len(old_samples[0][0])
        colors = sns.color_palette(n_colors=ndims)
        for dim in range(ndims):
            fig = plt.figure()
            sns.distplot(samples[:, dim], color=colors[dim])
            fig.savefig('PyDREAM_example_cellular_cycle_dimension_' + str(dim))

    except ImportError:
        pass
else:
    run_kwargs = {'parameters':sampled_parameter_names, 'likelihood':cost_function, 'niterations':10000, 'nchains':nchains, 'multitry':False, 'gamma_levels':4, 'adapt_gamma':True, 'history_thin':1, 'model_name':'cellular_cycle_dreamzs_5chain', 'verbose':True}
