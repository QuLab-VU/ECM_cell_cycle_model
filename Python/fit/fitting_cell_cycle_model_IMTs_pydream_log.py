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

import os
os.environ["PYOPENCL_CTX"] ="0"

'''
#Initialize PySB solver object for running simulations.  Simulation timespan should match experimental data.
tmp_t = np.linspace(0, 1000, 10000)
sim = BngSimulator(model, tspan=tmp_t, verbose=True)
simres = sim.run(initials=Initial)
'''

#Load experimental data to which the model will be fit here.  The "experimental data" in this case is the IMTs distribution

####################################################################################
    
#Cost function

#from cost_function_BngSimulator_bis import cost_function

from cost_function_BngSimulator_bis_reduced import cost_function_reduced

###############################################################################

parameter_vector=np.array([1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1])
parameter_vector=np.log10(np.array([1.24091223, 0.31218308, 0.37257587, 1.8666591 , 1.4299013 , 0.02153797, 0.07613843, 1.72653966, 0.10895379, 0.02389209]))

parameter_vector=np.array([1.24091223, 0.31218308, 0.37257587, 1.8666591 , 1.4299013 , 0.02153797, 0.07613843, 1.72653966, 0.10895379, 0.02389209])


cost_function_reduced(parameter_vector)


cost_function(parameter_vector)


def cost_function_f(parameter_vector):
    parameter_vector=[10**i for i in parameter_vector]
    return -cost_function(parameter_vector)

#Set upper and lower limits for uniform prior to be 3 orders of magnitude above and below original parameter values.

lower_limits = [-4 ,-4 ,-4 ,-4 ,-4 ,-4 ,-4 ,-4 ,-4 ,-4]
parameters_to_sample = SampledParam(uniform, loc=lower_limits, scale=5) #loc is the lower bdry and scale is the size of the distribution
sampled_parameter_names = [parameters_to_sample]
niterations = 50
converged = False
total_iterations = niterations
nchains = 10


niterations = 600
total_iterations = 5000



cellular_cycle_dreamzs_5chain_sampled_params_chain_0_05000=np.load(path + "/fit_data/pydream/bis_log/cellular_cycle_dreamzs_5chain_sampled_params_chain_0_0-5000.npy")


old_samples = [np.load(path + "/fit_data/pydream/bis_log/cellular_cycle_dreamzs_5chain_sampled_params_chain_0_0-5000.npy"), np.load(path + "/fit_data/pydream/bis_log/cellular_cycle_dreamzs_5chain_sampled_params_chain_1_0-5000.npy"), np.load(path + "/fit_data/pydream/bis_log/cellular_cycle_dreamzs_5chain_sampled_params_chain_2_0-5000.npy"), np.load(path + "/fit_data/pydream/bis_log/cellular_cycle_dreamzs_5chain_sampled_params_chain_3_0-5000.npy"), np.load(path + "/fit_data/pydream/bis_log/cellular_cycle_dreamzs_5chain_sampled_params_chain_4_0-5000.npy"), np.load(path + "/fit_data/pydream/bis_log/cellular_cycle_dreamzs_5chain_sampled_params_chain_5_0-5000.npy"), np.load(path + "/fit_data/pydream/bis_log/cellular_cycle_dreamzs_5chain_sampled_params_chain_6_0-5000.npy"), np.load(path + "/fit_data/pydream/bis_log/cellular_cycle_dreamzs_5chain_sampled_params_chain_7_0-5000.npy"), np.load(path + "/fit_data/pydream/bis_log/cellular_cycle_dreamzs_5chain_sampled_params_chain_8_0-5000.npy"), np.load(path + "/fit_data/pydream/bis_log/cellular_cycle_dreamzs_5chain_sampled_params_chain_9_0-5000.npy")]


sampled_params = [np.load(path + "/fit_data/pydream/bis_log/cellular_cycle_dreamzs_5chain_sampled_params_chain_0_5000.npy"),np.load(path + "/fit_data/pydream/bis_log/cellular_cycle_dreamzs_5chain_sampled_params_chain_1_5000.npy"),np.load(path + "/fit_data/pydream/bis_log/cellular_cycle_dreamzs_5chain_sampled_params_chain_2_5000.npy"),np.load(path + "/fit_data/pydream/bis_log/cellular_cycle_dreamzs_5chain_sampled_params_chain_3_5000.npy"),np.load(path + "/fit_data/pydream/bis_log/cellular_cycle_dreamzs_5chain_sampled_params_chain_4_5000.npy"),np.load(path + "/fit_data/pydream/bis_log/cellular_cycle_dreamzs_5chain_sampled_params_chain_5_5000.npy"),np.load(path + "/fit_data/pydream/bis_log/cellular_cycle_dreamzs_5chain_sampled_params_chain_6_5000.npy"),np.load(path + "/fit_data/pydream/bis_log/cellular_cycle_dreamzs_5chain_sampled_params_chain_7_5000.npy"),np.load(path + "/fit_data/pydream/bis_log/cellular_cycle_dreamzs_5chain_sampled_params_chain_8_5000.npy"),np.load(path + "/fit_data/pydream/bis_log/cellular_cycle_dreamzs_5chain_sampled_params_chain_9_5000.npy")]


#sampled_params, log_ps = run_dream(sampled_parameter_names, cost_function_f, niterations=niterations, nchains=nchains, multitry=False, gamma_levels=4, adapt_gamma=True, history_thin=1, model_name='cellular_cycle_dreamzs_5chain', verbose=True)


#sampled_params = np.load(path + "/fit_data/pydream/bis_log/sample_parameters.npy")
GR=np.loadtxt(path + "/fit_data/pydream/bis_log/cellular_cycle_dreamzs_5chain_GelmanRubin_iteration_5000.txt")

#for chain in range(len(sampled_params)):
#    np.save(path + "/fit_data/pydream/bis_log/cellular_cycle_dreamzs_5chain_sampled_params_chain_" + str(chain) + '_' + str(total_iterations), sampled_params[chain])
#    np.save(path + "/fit_data/pydream/bis_log/cellular_cycle_dreamzs_5chain_logps_chain_" + str(chain) + '_' + str(total_iterations),log_ps[chain])
                
#old_samples = [np.concatenate((old_samples[chain], sampled_params[chain])) for chain in range(nchains)]
#GR = Gelman_Rubin(old_samples)
#print('At iteration: ', total_iterations, ' GR = ', GR)
#np.savetxt(path + "/fit_data/pydream/bis_log/cellular_cycle_dreamzs_5chain_GelmanRubin_iteration_" + str(total_iterations)+'.txt', GR)
            
            
            
if __name__ == '__main__':
    #Run DREAM sampling.  Documentation of DREAM options is in Dream.py.
    #sampled_params, log_ps = run_dream(sampled_parameter_names, cost_function_f, niterations=niterations, nchains=nchains, multitry=False, gamma_levels=4, adapt_gamma=True, history_thin=1, model_name = 'cellular_cycle_dreamzs_5chain', verbose=True)
    
    
    
    #Save sampling output (sampled parameter values and their corresponding logps).
#    for chain in range(len(sampled_params)):
#        np.save(path + "/fit_data/pydream/bis_log/cellular_cycle_dreamzs_5chain_sampled_params_chain_"+str(chain)+'_'+str(total_iterations), sampled_params[chain])
#        np.save(path + "/fit_data/pydream/bis_log/cellular_cycle_dreamzs_5chain_logps_chain_"+str(chain)+'_'+str(total_iterations), log_ps[chain])
        
    # Check convergence and continue sampling if not converged
    
#    GR = Gelman_Rubin(sampled_params)
#    print('At iteration: ', total_iterations, ' GR = ', GR)
#    np.savetxt(path + "/fit_data/pydream/bis_log/cellular_cycle_dreamzs_5chain_GelmanRubin_iteration_" + str(total_iterations) + '.txt', GR)
    
    #old_samples = sampled_params
    if np.any(GR > 1.2):
        starts = [sampled_params[chain][-1, :] for chain in range(nchains)]
        while not converged:
            total_iterations += niterations
            
            sampled_params, log_ps = run_dream(sampled_parameter_names, cost_function_f, start=starts, niterations=niterations, nchains=nchains, multitry=False, gamma_levels=4, adapt_gamma=True, history_thin=1, model_name='cellular_cycle_dreamzs_5chain', verbose=True, restart=True)
            
            for chain in range(len(sampled_params)):
                np.save(path + "/fit_data/pydream/bis_log/cellular_cycle_dreamzs_5chain_sampled_params_chain_" + str(chain) + '_' + str(total_iterations), sampled_params[chain])
                np.save(path + "/fit_data/pydream/bis_log/cellular_cycle_dreamzs_5chain_logps_chain_" + str(chain) + '_' + str(total_iterations),log_ps[chain])
                
            old_samples = [np.concatenate((old_samples[chain], sampled_params[chain])) for chain in range(nchains)]
            GR = Gelman_Rubin(old_samples)
            print('At iteration: ', total_iterations, ' GR = ', GR)
            np.savetxt(path + "/fit_data/pydream/bis_log/cellular_cycle_dreamzs_5chain_GelmanRubin_iteration_" + str(total_iterations)+'.txt', GR)
            
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
    run_kwargs = {'parameters':sampled_parameter_names, 'likelihood':cost_function_f, 'niterations':10000, 'nchains':nchains, 'multitry':False, 'gamma_levels':4, 'adapt_gamma':True, 'history_thin':1, 'model_name':'cellular_cycle_dreamzs_5chain', 'verbose':True}

####################################################################################
###fitting with DMSO only
#Cost function

from cost_function_BngSimulator_bis_reduced_DMSO_only import cost_function_reduced_DMSO_only
###############################################################################

parameter_vector=np.array([1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1])
parameter_vector=np.log10(np.array([1.24091223, 0.31218308, 0.37257587, 1.8666591 , 1.4299013 , 0.02153797, 0.07613843, 1.72653966, 0.10895379, 0.02389209]))

parameter_vector=np.array([1.24091223, 0.31218308, 0.37257587, 1.8666591 , 1.4299013 , 0.02153797, 0.07613843, 1.72653966, 0.10895379, 0.02389209])


cost_function_reduced_DMSO_only(parameter_vector)



def cost_function_f(parameter_vector):
    parameter_vector=[10**i for i in parameter_vector]
    return -cost_function_reduced_DMSO_only(parameter_vector)

#Set upper and lower limits for uniform prior to be 3 orders of magnitude above and below original parameter values.

lower_limits = [-4 ,-4 ,-4 ,-4 ,-4 ,-4 ,-4 ,-4 ,-4 ,-4]
parameters_to_sample = SampledParam(uniform, loc=lower_limits, scale=5) #loc is the lower bdry and scale is the size of the distribution
sampled_parameter_names = [parameters_to_sample]
niterations = 50
converged = False
total_iterations = niterations
nchains = 10


niterations = 500
total_iterations = 0


if __name__ == '__main__':
    #Run DREAM sampling.  Documentation of DREAM options is in Dream.py.
    sampled_params, log_ps = run_dream(sampled_parameter_names, cost_function_f, niterations=niterations, nchains=nchains, multitry=False, gamma_levels=4, adapt_gamma=True, history_thin=1, model_name = 'cellular_cycle_dreamzs_5chain', verbose=True)
    #Save sampling output (sampled parameter values and their corresponding logps).
    for chain in range(len(sampled_params)):
        np.save(path + "/fit_data/pydream/bis_log_DMSO_only/cellular_cycle_dreamzs_5chain_sampled_params_chain_"+str(chain)+'_'+str(total_iterations), sampled_params[chain])
        np.save(path + "/fit_data/pydream/bis_log_DMSO_only/cellular_cycle_dreamzs_5chain_logps_chain_"+str(chain)+'_'+str(total_iterations), log_ps[chain])
     #Check convergence and continue sampling if not converged
    GR = Gelman_Rubin(sampled_params)
    print('At iteration: ', total_iterations, ' GR = ', GR)
    np.savetxt(path + "/fit_data/pydream/bis_log_DMSO_only/cellular_cycle_dreamzs_5chain_GelmanRubin_iteration_" + str(total_iterations) + '.txt', GR)
    old_samples = sampled_params
    if np.any(GR > 1.2):
        starts = [sampled_params[chain][-1, :] for chain in range(nchains)]
        while not converged:
            total_iterations += niterations
            sampled_params, log_ps = run_dream(sampled_parameter_names, cost_function_f, start=starts, niterations=niterations, nchains=nchains, multitry=False, gamma_levels=4, adapt_gamma=True, history_thin=1, model_name='cellular_cycle_dreamzs_5chain', verbose=True, restart=True)
            for chain in range(len(sampled_params)):
                np.save(path + "/fit_data/pydream/bis_log_DMSO_only/cellular_cycle_dreamzs_5chain_sampled_params_chain_" + str(chain) + '_' + str(total_iterations), sampled_params[chain])
                np.save(path + "/fit_data/pydream/bis_log_DMSO_only/cellular_cycle_dreamzs_5chain_logps_chain_" + str(chain) + '_' + str(total_iterations),log_ps[chain])
            old_samples = [np.concatenate((old_samples[chain], sampled_params[chain])) for chain in range(nchains)]
            GR = Gelman_Rubin(old_samples)
            print('At iteration: ', total_iterations, ' GR = ', GR)
            np.savetxt(path + "/fit_data/pydream/bis_log_DMSO_only/cellular_cycle_dreamzs_5chain_GelmanRubin_iteration_" + str(total_iterations)+'.txt', GR)
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
    run_kwargs = {'parameters':sampled_parameter_names, 'likelihood':cost_function_f, 'niterations':10000, 'nchains':nchains, 'multitry':False, 'gamma_levels':4, 'adapt_gamma':True, 'history_thin':1, 'model_name':'cellular_cycle_dreamzs_5chain', 'verbose':True}



####################################################################################
###fitting with DMSO only quiescence
#Cost function

from cost_function_BngSimulator_bis_reduced_quiescence_DMSO_only import cost_function_reduced_quiescence_DMSO_only
###############################################################################

parameter_vector=np.array([1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1])
parameter_vector=np.log10(np.array([1.24091223, 0.31218308, 0.37257587, 1.8666591 , 1.4299013 , 0.02153797, 0.07613843, 1.72653966, 0.10895379, 0.02389209]))

parameter_vector=np.array([1.24091223, 0.31218308, 0.37257587, 1.8666591 , 1.4299013 , 0.02153797, 0.07613843, 1.72653966, 0.10895379, 0.02389209])

cost_function_reduced_quiescence_DMSO_only(parameter_vector)



def cost_function_f(parameter_vector):
    parameter_vector=[10**i for i in parameter_vector]
    return -cost_function_reduced_quiescence_DMSO_only(parameter_vector)

#Set upper and lower limits for uniform prior to be 3 orders of magnitude above and below original parameter values.

lower_limits = [-4 ,-4 ,-4 ,-4 ,-4 ,-4 ,-4 ,-4 ,-4 ,-4]
parameters_to_sample = SampledParam(uniform, loc=lower_limits, scale=5) #loc is the lower bdry and scale is the size of the distribution
sampled_parameter_names = [parameters_to_sample]
niterations = 50
converged = False
total_iterations = niterations
nchains = 10


niterations = 500
total_iterations = 0


if __name__ == '__main__':
    #Run DREAM sampling.  Documentation of DREAM options is in Dream.py.
    sampled_params, log_ps = run_dream(sampled_parameter_names, cost_function_f, niterations=niterations, nchains=nchains, multitry=False, gamma_levels=4, adapt_gamma=True, history_thin=1, model_name = 'cellular_cycle_dreamzs_5chain', verbose=True)
    #Save sampling output (sampled parameter values and their corresponding logps).
    for chain in range(len(sampled_params)):
        np.save(path + "/fit_data/pydream/bis_log_DMSO_only_quiescence_2/cellular_cycle_dreamzs_5chain_sampled_params_chain_"+str(chain)+'_'+str(total_iterations), sampled_params[chain])
        np.save(path + "/fit_data/pydream/bis_log_DMSO_only_quiescence_2/cellular_cycle_dreamzs_5chain_logps_chain_"+str(chain)+'_'+str(total_iterations), log_ps[chain])
     #Check convergence and continue sampling if not converged
    GR = Gelman_Rubin(sampled_params)
    print('At iteration: ', total_iterations, ' GR = ', GR)
    np.savetxt(path + "/fit_data/pydream/bis_log_DMSO_only_quiescence_2/cellular_cycle_dreamzs_5chain_GelmanRubin_iteration_" + str(total_iterations) + '.txt', GR)
    old_samples = sampled_params
    if np.any(GR > 1.2):
        starts = [sampled_params[chain][-1, :] for chain in range(nchains)]
        while not converged:
            total_iterations += niterations
            sampled_params, log_ps = run_dream(sampled_parameter_names, cost_function_f, start=starts, niterations=niterations, nchains=nchains, multitry=False, gamma_levels=4, adapt_gamma=True, history_thin=1, model_name='cellular_cycle_dreamzs_5chain', verbose=True, restart=True)
            for chain in range(len(sampled_params)):
                np.save(path + "/fit_data/pydream/bis_log_DMSO_only_quiescence_2/cellular_cycle_dreamzs_5chain_sampled_params_chain_" + str(chain) + '_' + str(total_iterations), sampled_params[chain])
                np.save(path + "/fit_data/pydream/bis_log_DMSO_only_quiescence_2/cellular_cycle_dreamzs_5chain_logps_chain_" + str(chain) + '_' + str(total_iterations),log_ps[chain])
            old_samples = [np.concatenate((old_samples[chain], sampled_params[chain])) for chain in range(nchains)]
            GR = Gelman_Rubin(old_samples)
            print('At iteration: ', total_iterations, ' GR = ', GR)
            np.savetxt(path + "/fit_data/pydream/bis_log_DMSO_only_quiescence_2/cellular_cycle_dreamzs_5chain_GelmanRubin_iteration_" + str(total_iterations)+'.txt', GR)
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
    run_kwargs = {'parameters':sampled_parameter_names, 'likelihood':cost_function_f, 'niterations':10000, 'nchains':nchains, 'multitry':False, 'gamma_levels':4, 'adapt_gamma':True, 'history_thin':1, 'model_name':'cellular_cycle_dreamzs_5chain', 'verbose':True}


####################################################################################
###fitting with DMSO only quiescence only
#Cost function

from cost_function_BngSimulator_bis_reduced_quiescence_only_DMSO_only import cost_function_reduced_quiescence_only_DMSO_only
###############################################################################

parameter_vector=np.array([1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1])
parameter_vector=np.log10(np.array([1.24091223, 0.31218308, 0.37257587, 1.8666591 , 1.4299013 , 0.02153797, 0.07613843, 1.72653966, 0.10895379, 0.02389209]))

parameter_vector=np.array([1.24091223, 0.31218308, 0.37257587, 1.8666591 , 1.4299013 , 0.02153797, 0.07613843, 1.72653966, 0.10895379, 0.02389209])

cost_function_reduced_quiescence_only_DMSO_only(parameter_vector)



def cost_function_f(parameter_vector):
    parameter_vector=[10**i for i in parameter_vector]
    return -cost_function_reduced_quiescence_only_DMSO_only(parameter_vector)

#Set upper and lower limits for uniform prior to be 3 orders of magnitude above and below original parameter values.

lower_limits = [-4 ,-4 ,-4 ,-4 ,-4 ,-4 ,-4 ,-4 ,-4 ,-4]
parameters_to_sample = SampledParam(uniform, loc=lower_limits, scale=5) #loc is the lower bdry and scale is the size of the distribution
sampled_parameter_names = [parameters_to_sample]
niterations = 50
converged = False
total_iterations = niterations
nchains = 10


niterations = 500
total_iterations = 0


if __name__ == '__main__':
    #Run DREAM sampling.  Documentation of DREAM options is in Dream.py.
    sampled_params, log_ps = run_dream(sampled_parameter_names, cost_function_f, niterations=niterations, nchains=nchains, multitry=False, gamma_levels=4, adapt_gamma=True, history_thin=1, model_name = 'cellular_cycle_dreamzs_5chain', verbose=True)
    #Save sampling output (sampled parameter values and their corresponding logps).
    for chain in range(len(sampled_params)):
        np.save(path + "/fit_data/pydream/bis_log_DMSO_only_quiescence_only/cellular_cycle_dreamzs_5chain_sampled_params_chain_"+str(chain)+'_'+str(total_iterations), sampled_params[chain])
        np.save(path + "/fit_data/pydream/bis_log_DMSO_only_quiescence_only/cellular_cycle_dreamzs_5chain_logps_chain_"+str(chain)+'_'+str(total_iterations), log_ps[chain])
     #Check convergence and continue sampling if not converged
    GR = Gelman_Rubin(sampled_params)
    print('At iteration: ', total_iterations, ' GR = ', GR)
    np.savetxt(path + "/fit_data/pydream/bis_log_DMSO_only_quiescence_only/cellular_cycle_dreamzs_5chain_GelmanRubin_iteration_" + str(total_iterations) + '.txt', GR)
    old_samples = sampled_params
    if np.any(GR > 1.2):
        starts = [sampled_params[chain][-1, :] for chain in range(nchains)]
        while not converged:
            total_iterations += niterations
            sampled_params, log_ps = run_dream(sampled_parameter_names, cost_function_f, start=starts, niterations=niterations, nchains=nchains, multitry=False, gamma_levels=4, adapt_gamma=True, history_thin=1, model_name='cellular_cycle_dreamzs_5chain', verbose=True, restart=True)
            for chain in range(len(sampled_params)):
                np.save(path + "/fit_data/pydream/bis_log_DMSO_only_quiescence_only/cellular_cycle_dreamzs_5chain_sampled_params_chain_" + str(chain) + '_' + str(total_iterations), sampled_params[chain])
                np.save(path + "/fit_data/pydream/bis_log_DMSO_only_quiescence_only/cellular_cycle_dreamzs_5chain_logps_chain_" + str(chain) + '_' + str(total_iterations),log_ps[chain])
            old_samples = [np.concatenate((old_samples[chain], sampled_params[chain])) for chain in range(nchains)]
            GR = Gelman_Rubin(old_samples)
            print('At iteration: ', total_iterations, ' GR = ', GR)
            np.savetxt(path + "/fit_data/pydream/bis_log_DMSO_only_quiescence_only/cellular_cycle_dreamzs_5chain_GelmanRubin_iteration_" + str(total_iterations)+'.txt', GR)
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
    run_kwargs = {'parameters':sampled_parameter_names, 'likelihood':cost_function_f, 'niterations':10000, 'nchains':nchains, 'multitry':False, 'gamma_levels':4, 'adapt_gamma':True, 'history_thin':1, 'model_name':'cellular_cycle_dreamzs_5chain', 'verbose':True}
    
####################################################################################
###fitting with DMSO only quiescence with more parameters
#Cost function

from cost_function_BngSimulator_bis_reduced_quiescence_DMSO_only_21_parameters import cost_function_reduced_quiescence_DMSO_only_21_parameters
###############################################################################

pysb_sampled_parameter_names = ["V_s2skp","V_Fp","k_M","k_EMY","k_CDF","k_CDS","k_ED","K_S","K_MY","K_DF","d_CD","k_P1", "k_P2", "K_CD", "k_E", "K_M", "K_CE", "V_2E2F", "K_2E2F", "k_RE", "k_b"]



parameter_vector=np.array([1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1, 1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1, 1])
#parameter_vector=np.log10(np.array([1.24091223, 0.31218308, 0.37257587, 1.8666591 , 1.4299013 , 0.02153797, 0.07613843, 1.72653966, 0.10895379, 0.02389209]))

#parameter_vector=np.array([1.24091223, 0.31218308, 0.37257587, 1.8666591 , 1.4299013 , 0.02153797, 0.07613843, 1.72653966, 0.10895379, 0.02389209])

cost_function_reduced_quiescence_DMSO_only_21_parameters(parameter_vector)



def cost_function_f(parameter_vector):
    parameter_vector=[10**i for i in parameter_vector]
    return -cost_function_reduced_quiescence_DMSO_only_21_parameters(parameter_vector)

#Set upper and lower limits for uniform prior to be 3 orders of magnitude above and below original parameter values.

lower_limits = [-4 ,-4 ,-4 ,-4 ,-4 ,-4 ,-4 ,-4 ,-4 ,-4,-4 ,-4 ,-4 ,-4 ,-4 ,-4 ,-4 ,-4 ,-4 ,-4, -4]
parameters_to_sample = SampledParam(uniform, loc=lower_limits, scale=5) #loc is the lower bdry and scale is the size of the distribution
sampled_parameter_names = [parameters_to_sample]
niterations = 5
converged = False
total_iterations = niterations
nchains = 10


niterations = 5
total_iterations = 0


if __name__ == '__main__':
    #Run DREAM sampling.  Documentation of DREAM options is in Dream.py.
    sampled_params, log_ps = run_dream(sampled_parameter_names, cost_function_f, niterations=niterations, nchains=nchains, multitry=False, gamma_levels=4, adapt_gamma=True, history_thin=1, model_name = 'cellular_cycle_dreamzs_5chain', verbose=True)
    #Save sampling output (sampled parameter values and their corresponding logps).
    for chain in range(len(sampled_params)):
        np.save(path + "/fit_data/pydream/bis_log_DMSO_only_quiescence_21_parameters/cellular_cycle_dreamzs_5chain_sampled_params_chain_"+str(chain)+'_'+str(total_iterations), sampled_params[chain])
        np.save(path + "/fit_data/pydream/bis_log_DMSO_only_quiescence_21_parameters/cellular_cycle_dreamzs_5chain_logps_chain_"+str(chain)+'_'+str(total_iterations), log_ps[chain])
     #Check convergence and continue sampling if not converged
    GR = Gelman_Rubin(sampled_params)
    print('At iteration: ', total_iterations, ' GR = ', GR)
    np.savetxt(path + "/fit_data/pydream/bis_log_DMSO_only_quiescence_21_parameters/cellular_cycle_dreamzs_5chain_GelmanRubin_iteration_" + str(total_iterations) + '.txt', GR)
    old_samples = sampled_params
    if np.any(GR > 1.2):
        starts = [sampled_params[chain][-1, :] for chain in range(nchains)]
        while not converged:
            total_iterations += niterations
            sampled_params, log_ps = run_dream(sampled_parameter_names, cost_function_f, start=starts, niterations=niterations, nchains=nchains, multitry=False, gamma_levels=4, adapt_gamma=True, history_thin=1, model_name='cellular_cycle_dreamzs_5chain', verbose=True, restart=True)
            for chain in range(len(sampled_params)):
                np.save(path + "/fit_data/pydream/bis_log_DMSO_only_quiescence_21_parameters/cellular_cycle_dreamzs_5chain_sampled_params_chain_" + str(chain) + '_' + str(total_iterations), sampled_params[chain])
                np.save(path + "/fit_data/pydream/bis_log_DMSO_only_quiescence_21_parameters/cellular_cycle_dreamzs_5chain_logps_chain_" + str(chain) + '_' + str(total_iterations),log_ps[chain])
            old_samples = [np.concatenate((old_samples[chain], sampled_params[chain])) for chain in range(nchains)]
            GR = Gelman_Rubin(old_samples)
            print('At iteration: ', total_iterations, ' GR = ', GR)
            np.savetxt(path + "/fit_data/pydream/bis_log_DMSO_only_quiescence_21_parameters/cellular_cycle_dreamzs_5chain_GelmanRubin_iteration_" + str(total_iterations)+'.txt', GR)
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
    run_kwargs = {'parameters':sampled_parameter_names, 'likelihood':cost_function_f, 'niterations':10000, 'nchains':nchains, 'multitry':False, 'gamma_levels':4, 'adapt_gamma':True, 'history_thin':1, 'model_name':'cellular_cycle_dreamzs_5chain', 'verbose':True}



####################################################################################
###fitting with DMSO only 0.3 kPa with fracprolif and model no3
#Cost function

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
from cost_function_BngSimulator_frac_prolif_quiescence_DMSO_only_21_parameters import cost_function_reduced_quiescence_DMSO_only_21_parameters

from Austin_proliferation_model_population_3 import model
model_population=model
###############################################################################

pysb_sampled_parameter_names = ["V_s2skp","V_Fp","k_M","k_EMY","k_CDF","k_CDS","k_ED","K_S","K_MY","K_DF","d_CD","k_P1", "k_P2", "K_CD", "k_E", "K_M", "K_CE", "V_2E2F", "K_2E2F", "k_RE", "k_b"]


parameter_vector_original=[]
for pname in pysb_sampled_parameter_names:
    #Change model parameter values to current location in parameter space
    parameter_vector_original+=[model_population.parameters[pname].value]

cost_2 = cost_function_reduced_quiescence_DMSO_only_21_parameters(parameter_vector_original)


parameter_vector=np.array([1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1, 1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1, 1])

cost = cost_function_reduced_quiescence_DMSO_only_21_parameters(parameter_vector)


#parameter_vector=np.log10(np.array([1.24091223, 0.31218308, 0.37257587, 1.8666591 , 1.4299013 , 0.02153797, 0.07613843, 1.72653966, 0.10895379, 0.02389209]))

#parameter_vector=np.array([1.24091223, 0.31218308, 0.37257587, 1.8666591 , 1.4299013 , 0.02153797, 0.07613843, 1.72653966, 0.10895379, 0.02389209])



def cost_function_f(parameter_vector):
    parameter_vector=[10**i for i in parameter_vector]
    return -cost_function_reduced_quiescence_DMSO_only_21_parameters(parameter_vector)



parameter_vector_2=pl.rand(21)*6-2


cost_2 = cost_function_f(parameter_vector_2)


#Set upper and lower limits for uniform prior to be 3 orders of magnitude above and below original parameter values.

lower_limits = [-4 ,-4 ,-4 ,-4 ,-4 ,-4 ,-4 ,-4 ,-4 ,-4,-4 ,-4 ,-4 ,-4 ,-4 ,-4 ,-4 ,-4 ,-4 ,-4, -4]
parameters_to_sample = SampledParam(uniform, loc=lower_limits, scale=5) #loc is the lower bdry and scale is the size of the distribution
sampled_parameter_names = [parameters_to_sample]
niterations = 10
converged = False
total_iterations = niterations
nchains = 10


niterations = 5
total_iterations = 0


if __name__ == '__main__':
    #Run DREAM sampling.  Documentation of DREAM options is in Dream.py.
    sampled_params, log_ps = run_dream(sampled_parameter_names, cost_function_f, niterations=niterations, nchains=nchains, multitry=False, gamma_levels=4, adapt_gamma=True, history_thin=1, model_name = 'cellular_cycle_dreamzs_5chain', verbose=True)
    #Save sampling output (sampled parameter values and their corresponding logps).
    for chain in range(len(sampled_params)):
        np.save(path + "/fit_data/pydream/frac_prolif_quiescence_DMSO_only_21_parameters/cellular_cycle_dreamzs_5chain_sampled_params_chain_"+str(chain)+'_'+str(total_iterations), sampled_params[chain])
        np.save(path + "/fit_data/pydream/frac_prolif_quiescence_DMSO_only_21_parameters/cellular_cycle_dreamzs_5chain_logps_chain_"+str(chain)+'_'+str(total_iterations), log_ps[chain])
     #Check convergence and continue sampling if not converged
    GR = Gelman_Rubin(sampled_params)
    print('At iteration: ', total_iterations, ' GR = ', GR)
    np.savetxt(path + "/fit_data/pydream/frac_prolif_quiescence_DMSO_only_21_parameters/cellular_cycle_dreamzs_5chain_GelmanRubin_iteration_" + str(total_iterations) + '.txt', GR)
    old_samples = sampled_params
    if np.any(GR > 1.2):
        starts = [sampled_params[chain][-1, :] for chain in range(nchains)]
        while not converged:
            total_iterations += niterations
            sampled_params, log_ps = run_dream(sampled_parameter_names, cost_function_f, start=starts, niterations=niterations, nchains=nchains, multitry=False, gamma_levels=4, adapt_gamma=True, history_thin=1, model_name='cellular_cycle_dreamzs_5chain', verbose=True, restart=True)
            for chain in range(len(sampled_params)):
                np.save(path + "/fit_data/pydream/frac_prolif_quiescence_DMSO_only_21_parameters/cellular_cycle_dreamzs_5chain_sampled_params_chain_" + str(chain) + '_' + str(total_iterations), sampled_params[chain])
                np.save(path + "/fit_data/pydream/frac_prolif_quiescence_DMSO_only_21_parameters/cellular_cycle_dreamzs_5chain_logps_chain_" + str(chain) + '_' + str(total_iterations),log_ps[chain])
            old_samples = [np.concatenate((old_samples[chain], sampled_params[chain])) for chain in range(nchains)]
            GR = Gelman_Rubin(old_samples)
            print('At iteration: ', total_iterations, ' GR = ', GR)
            np.savetxt(path + "/fit_data/pydream/frac_prolif_quiescence_DMSO_only_21_parameters/cellular_cycle_dreamzs_5chain_GelmanRubin_iteration_" + str(total_iterations)+'.txt', GR)
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
    run_kwargs = {'parameters':sampled_parameter_names, 'likelihood':cost_function_f, 'niterations':10000, 'nchains':nchains, 'multitry':False, 'gamma_levels':4, 'adapt_gamma':True, 'history_thin':1, 'model_name':'cellular_cycle_dreamzs_5chain', 'verbose':True}




####################################################################################
###fitting with DMSO only 0.3 kPa with fracprolif and model no3
#Cost function

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
from cost_function_BngSimulator_frac_death_rate_DMSO_only_switch_death_parameters import  cost_function_reduced_death_rate_DMSO_switch_death_parameters

from Austin_proliferation_model_population_3 import model
model_population=model


###############################################################################



def cost_function_f(parameter_vector):
    parameter_vector=[10**i for i in parameter_vector]
    return -cost_function_reduced_death_rate_DMSO_switch_death_parameters(parameter_vector)



parameter_vector = [0.1, 3, 0.4, 80, 1, 0.005, 0.003, 0.008, 1, 0.005]
cost = cost_function_reduced_death_rate_DMSO_switch_death_parameters(parameter_vector)

pysb_sampled_parameter_names = ["kf1","kf2","kb","kasbaxmbcl2","kdsbaxmbcl2","k1sbh3","k2sbh3","kdbh3","kasbh3bcl2","kdsbh3bcl2"]
parameter_vector_original=[]
for pname in pysb_sampled_parameter_names:
    #Change model parameter values to current location in parameter space
    parameter_vector_original+=[model_population.parameters[pname].value]
cost_2 = cost_function_reduced_death_rate_DMSO_switch_death_parameters(parameter_vector_original)

parameter_vector_2=pl.rand(10)*6-2
cost_2 = cost_function_f(parameter_vector_2)


#Set upper and lower limits for uniform prior to be 3 orders of magnitude above and below original parameter values.

lower_limits = [-4 ,-4 ,-4 ,-4 ,-4 ,-4 ,-4 ,-4 ,-4, -4]
parameters_to_sample = SampledParam(uniform, loc=lower_limits, scale=6) #loc is the lower bdry and scale is the size of the distribution
sampled_parameter_names = [parameters_to_sample]

niterations = 20
converged = False

nchains = 10

total_iterations = 0


if __name__ == '__main__':
    #Run DREAM sampling.  Documentation of DREAM options is in Dream.py.
    sampled_params, log_ps = run_dream(sampled_parameter_names, cost_function_f, niterations=niterations, nchains=nchains, multitry=False, gamma_levels=4, adapt_gamma=True, history_thin=1, model_name = 'cellular_cycle_dreamzs_5chain', verbose=True)
    #Save sampling output (sampled parameter values and their corresponding logps).
    for chain in range(len(sampled_params)):
        np.save(path + "/fit_data/pydream/frac_prolif_quiescence_DMSO_death_parameters/cellular_cycle_dreamzs_5chain_sampled_params_chain_"+str(chain)+'_'+str(total_iterations), sampled_params[chain])
        np.save(path + "/fit_data/pydream/frac_prolif_quiescence_DMSO_death_parameters/cellular_cycle_dreamzs_5chain_logps_chain_"+str(chain)+'_'+str(total_iterations), log_ps[chain])
     #Check convergence and continue sampling if not converged
    GR = Gelman_Rubin(sampled_params)
    print('At iteration: ', total_iterations, ' GR = ', GR)
    np.savetxt(path + "/fit_data/pydream/frac_prolif_quiescence_DMSO_death_parameters/cellular_cycle_dreamzs_5chain_GelmanRubin_iteration_" + str(total_iterations) + '.txt', GR)
    old_samples = sampled_params
    if np.any(GR > 1.2):
        
        while not converged:
            starts = [sampled_params[chain][-1, :] for chain in range(nchains)]
            total_iterations += niterations
            sampled_params, log_ps = run_dream(sampled_parameter_names, cost_function_f, start=starts, niterations=niterations, nchains=nchains, multitry=False, gamma_levels=4, adapt_gamma=True, history_thin=1, model_name='cellular_cycle_dreamzs_5chain', verbose=True, restart=True)
            for chain in range(len(sampled_params)):
                np.save(path + "/fit_data/pydream/frac_prolif_quiescence_DMSO_death_parameters/cellular_cycle_dreamzs_5chain_sampled_params_chain_" + str(chain) + '_' + str(total_iterations), sampled_params[chain])
                np.save(path + "/fit_data/pydream/frac_prolif_quiescence_DMSO_death_parameters/cellular_cycle_dreamzs_5chain_logps_chain_" + str(chain) + '_' + str(total_iterations),log_ps[chain])
            old_samples = [np.concatenate((old_samples[chain], sampled_params[chain])) for chain in range(nchains)]
            GR = Gelman_Rubin(old_samples)
            print('At iteration: ', total_iterations, ' GR = ', GR)
            np.savetxt(path + "/fit_data/pydream/frac_prolif_quiescence_DMSO_death_parameters/cellular_cycle_dreamzs_5chain_GelmanRubin_iteration_" + str(total_iterations)+'.txt', GR)
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
    run_kwargs = {'parameters':sampled_parameter_names, 'likelihood':cost_function_f, 'niterations':10000, 'nchains':nchains, 'multitry':False, 'gamma_levels':4, 'adapt_gamma':True, 'history_thin':1, 'model_name':'cellular_cycle_dreamzs_5chain', 'verbose':True}


####################################################################################
###fitting with DMSO only 0.3 kPa with fracprolif and model no3
#Cost function

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
from cost_function_BngSimulator_frac_prolif_quiescence_DMSO_only_30_parameters import  cost_function_reduced_quiescence_DMSO_only_30_parameters

from Austin_proliferation_model_population_3 import model
model_population=model

###############################################################################

parameter_vector = [4.800635907279269, 3.9774300061066046, 0.0387054504992875, 0.32499986371412365, 0.35777622204261267, 0.0009944953058487452, 2.237824499761186, 0.28043611679549096, 0.021259592285147864, 0.006260825171175394, 0.0009632938154961734, 0.40008970395372623, 0.0008261834973168124, 0.0004808385449294265, 1.0125180610037363, 0.0020144839403416154, 0.3903389858441723, 0.9663335772673012, 1.2920585358743464, 2.179996170318332, 0.08588331032846658, 0.04121752679283465, 0.0034152812637972276, 0.650689069719618, 0.738165002953985, 0.006021052129410297, 0.001183347044793872, 5.033805671111385, 14.46792006528087, 0.029428573646856398]

parameter_vector = [7.09698462592513, 1.4106440558775277, 9.344724006360165, 13.56943889467159, 0.062499924302699336, 0.011251847411602133, 37.07864288707931, 0.42969242570255894, 0.04448402040986224, 0.0003588886046907576, 0.0007497183316475982, 1.6173514368404565, 0.008236520864891912, 0.003753895712771778, 0.2520808713489205, 0.00022560992699410936, 2.580054373455581, 5.9504955699203546, 0.1572599137099059, 15.205498329912354, 0.05699574585320892, 0.003489022874309077, 0.0016770883398079064, 0.29843523615760786, 1.2325484452321356, 0.00237113573463495, 0.0010247951539084355, 5.585134808848889, 51.52911949628718, 0.006419876379536079]

def cost_function_f(parameter_vector):
    parameter_vector=[10**i for i in parameter_vector]
    return -cost_function_reduced_quiescence_DMSO_only_30_parameters(parameter_vector)


cost_2 = cost_function_f(parameter_vector)


pysb_sampled_parameter_names = ["V_s2skp","V_Fp","k_M","k_EMY","k_CDF","k_CDS","k_ED","K_S","K_MY","K_DF","d_CD","k_P1", "k_P2", "K_CD", "k_E", "K_M", "K_CE", "V_2E2F", "K_2E2F", "k_RE", "kf1", "kf2", "kb", "kasbaxmbcl2", "kdsbaxmbcl2", "k1sbh3", "k2sbh3", "kdbh3", "kasbh3bcl2", "kdsbh3bcl2"]
parameter_vector_original=[]
for pname in pysb_sampled_parameter_names:
    #Change model parameter values to current location in parameter space
    parameter_vector_original+=[model_population.parameters[pname].value]


cost_2 = cost_function_reduced_quiescence_DMSO_only_30_parameters(parameter_vector_original)



parameter_vector_2=pl.rand(30)*6-2
cost_2 = cost_function_f(parameter_vector_2)


#Set upper and lower limits for uniform prior to be 3 orders of magnitude above and below original parameter values.


pysb_sampled_parameter_names = ["V_s2skp","V_Fp","k_M","k_EMY","k_CDF","k_CDS","k_ED","K_S","K_MY","K_DF","d_CD","k_P1", "k_P2", "K_CD", "k_E", "K_M", "K_CE", "V_2E2F", "K_2E2F", "k_RE", "kf1", "kf2", "kb", "kasbaxmbcl2", "kdsbaxmbcl2", "k1sbh3", "k2sbh3", "kdbh3", "kasbh3bcl2", "kdsbh3bcl2"]

lower_limits = [-4 ,-4 ,-4 ,-4 ,-4 ,-4 ,-4 ,-4 ,-4, -4, -4 ,-4 ,-4 ,-4 ,-4 ,-4 ,-4 ,-4 ,-4, -4, -4 ,-4 ,-4 ,-4 ,-4 ,-4 ,-4 ,-4 ,-4, -4]
parameters_to_sample = SampledParam(uniform, loc=lower_limits, scale=6) #loc is the lower bdry and scale is the size of the distribution
sampled_parameter_names = [parameters_to_sample]

niterations = 4
converged = False

nchains = 10

total_iterations = 0


if __name__ == '__main__':
    #Run DREAM sampling.  Documentation of DREAM options is in Dream.py.
    sampled_params, log_ps = run_dream(sampled_parameter_names, cost_function_f, niterations=niterations, nchains=nchains, multitry=False, gamma_levels=1, adapt_gamma=True, history_thin=1, model_name = 'cellular_cycle_dreamzs_5chain', verbose=True)
    print('ok')
    #Save sampling output (sampled parameter values and their corresponding logps).
    for chain in range(len(sampled_params)):
        np.save(path + "/fit_data/pydream/frac_prolif_quiescence_DMSO_only_30_parameters/cellular_cycle_dreamzs_5chain_sampled_params_chain_"+str(chain)+'_'+str(total_iterations), sampled_params[chain])
        np.save(path + "/fit_data/pydream/frac_prolif_quiescence_DMSO_only_30_parameters/cellular_cycle_dreamzs_5chain_logps_chain_"+str(chain)+'_'+str(total_iterations), log_ps[chain])
    
    print('saved')
     #Check convergence and continue sampling if not converged
    GR = Gelman_Rubin(sampled_params)
    print('At iteration: ', total_iterations, ' GR = ', GR)
    np.savetxt(path + "/fit_data/pydream/frac_prolif_quiescence_DMSO_only_30_parameters/cellular_cycle_dreamzs_5chain_GelmanRubin_iteration_" + str(total_iterations) + '.txt', GR)
    old_samples = sampled_params
    if np.any(GR > 1.2):
        
        while not converged:
            print("new_start")
            starts = [sampled_params[chain][-1, :] for chain in range(nchains)]
            total_iterations += niterations
            sampled_params, log_ps = run_dream(sampled_parameter_names, cost_function_f, start=starts, niterations=niterations, nchains=nchains, multitry=False, gamma_levels=1, adapt_gamma=True, history_thin=1, model_name='cellular_cycle_dreamzs_5chain', verbose=True, restart=True)
            
            for chain in range(len(sampled_params)):
                np.save(path + "/fit_data/pydream/frac_prolif_quiescence_DMSO_only_30_parameters/cellular_cycle_dreamzs_5chain_sampled_params_chain_" + str(chain) + '_' + str(total_iterations), sampled_params[chain])
                np.save(path + "/fit_data/pydream/frac_prolif_quiescence_DMSO_only_30_parameters/cellular_cycle_dreamzs_5chain_logps_chain_" + str(chain) + '_' + str(total_iterations),log_ps[chain])
            print("saved")
            old_samples = [np.concatenate((old_samples[chain], sampled_params[chain])) for chain in range(nchains)]
            GR = Gelman_Rubin(old_samples)
            print('At iteration: ', total_iterations, ' GR = ', GR)
            np.savetxt(path + "/fit_data/pydream/frac_prolif_quiescence_DMSO_only_30_parameters/cellular_cycle_dreamzs_5chain_GelmanRubin_iteration_" + str(total_iterations)+'.txt', GR)
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
    run_kwargs = {'parameters':sampled_parameter_names, 'likelihood':cost_function_f, 'niterations':10000, 'nchains':nchains, 'multitry':False, 'gamma_levels':4, 'adapt_gamma':True, 'history_thin':1, 'model_name':'cellular_cycle_dreamzs_5chain', 'verbose':True}
    