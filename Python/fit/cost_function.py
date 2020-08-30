# -*- coding: utf-8 -*-
"""
Created on Wed May 20 10:42:26 2020

@author: Grim.M
"""

from Austin_proliferation_model_concentration import model
model_concentration=model

from Austin_proliferation_model_population import model
model_population=model

from pysb.simulator import BngSimulator

from SMs_G1s_simulation import IMTs
from rbf_mmd import nmmd2_rbf

import os
os.environ["PYOPENCL_CTX"] ="0"



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
    
    
    for condition in experimental_IMTs.keys():
        print(condition)
        
        if 'DMSO' in condition:            
            if '0.6 kPa' in condition :
                model_concentration.parameters['ERK'].value = 1.06970165050038
                model_concentration.parameters['FAK'].value = 0.337143388171177
                model_concentration.parameters['AKT'].value = 0.616658111482869
                
                model_population.parameters['ERK'].value = 1.06970165050038
                model_population.parameters['FAK'].value = 0.337143388171177
                model_population.parameters['AKT'].value = 0.616658111482869
            elif '4 kPa' in condition :
                model_concentration.parameters['ERK'].value = 0.897485095798471
                model_concentration.parameters['FAK'].value = 0.594150207445117
                model_concentration.parameters['AKT'].value = 1.01975420186633
                
                model_population.parameters['ERK'].value = 0.897485095798471
                model_population.parameters['FAK'].value = 0.594150207445117
                model_population.parameters['AKT'].value = 1.01975420186633
            else : #55 kPa
                model_concentration.parameters['ERK'].value = 1
                model_concentration.parameters['FAK'].value = 1
                model_concentration.parameters['AKT'].value = 1
                
                model_population.parameters['ERK'].value = 1
                model_population.parameters['FAK'].value = 1
                model_population.parameters['AKT'].value = 1
        else : #ERL treated
            if '0.6 kPa' in condition :
                model_concentration.parameters['ERK'].value = 0.186745227418267
                model_concentration.parameters['FAK'].value = 0.397060244114162
                model_concentration.parameters['AKT'].value = 0.245702711121333
                
                model_population.parameters['ERK'].value = 0.186745227418267
                model_population.parameters['FAK'].value = 0.397060244114162
                model_population.parameters['AKT'].value = 0.245702711121333
            elif '4 kPa' in condition :
                model_concentration.parameters['ERK'].value = 0.211453830597783
                model_concentration.parameters['FAK'].value = 0.709172297576571
                model_concentration.parameters['AKT'].value = 0.598518145949655
                
                model_population.parameters['ERK'].value = 0.211453830597783
                model_population.parameters['FAK'].value = 0.709172297576571
                model_population.parameters['AKT'].value = 0.598518145949655
            else : #55 kPa
                model_concentration.parameters['ERK'].value = 0.273618423688386
                model_concentration.parameters['FAK'].value = 1.075820621939
                model_concentration.parameters['AKT'].value = 0.585067267459499
                
                model_population.parameters['ERK'].value = 0.273618423688386
                model_population.parameters['FAK'].value = 1.075820621939
                model_population.parameters['AKT'].value = 0.585067267459499
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