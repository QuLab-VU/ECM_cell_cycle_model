# -*- coding: utf-8 -*-
"""
Created on Wed Apr 29 22:34:14 2020

@author: Grim.M
"""
path="/home/jeremie"
#path=r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code"

import sys
sys.path.insert(0, path+"/cellular_cycle") #add a pathway to the model

from peakfinder import peakfinder
from pysb.simulator import ScipyOdeSimulator
import numpy as np
import pylab as pl

import sobol_seq

from Austin_proliferation_model_concentration import model

pysb_sampled_parameter_names = ["V_s2skp","V_Fp","k_M","k_EMY","k_CDF","k_CDS","k_ED","K_S","K_MY","K_DF"]


l=sobol_seq.i4_sobol_generate(10, 50000)*2
t = np.linspace(0, 200, 120)
accepted=[]
rejected=[]

for parameter_vector in l :
    param_dict = {pname: pvalue for pname, pvalue in zip(pysb_sampled_parameter_names, parameter_vector)}
    for pname, pvalue in param_dict.items():
        #Change model parameter values to current location in parameter space
        model.parameters[pname].value = (pvalue)
    simres = ScipyOdeSimulator(model, t).run() 
    yout = simres.all
    a, maxima = peakfinder(yout['obsCdc20'][40:120],[],0.2,includeEndpoints=1)   #find the peak of cdc20 - last protein to peak
    b, maxima_b = peakfinder(yout['obsCycE'][40:120],[],0.2,includeEndpoints=1)   #find the peak of cdc20 - last protein to peak
    if len(a)>1 and len(b)>1:
        accepted+=[parameter_vector]
    else :
        rejected+=[parameter_vector]
        print(len(accepted)+len(rejected))


np.save("/home/jeremie/sobol_seq_set1_again_accepted.npy", accepted)
np.save("/home/jeremie/sobol_seq_set1_again_rejected.npy", rejected)








l=sobol_seq.i4_sobol_generate(10, 25000, skip=50000)*2
t = np.linspace(0, 200, 120)
accepted=[]
rejected=[]

for parameter_vector in l :
    param_dict = {pname: pvalue for pname, pvalue in zip(pysb_sampled_parameter_names, parameter_vector)}
    for pname, pvalue in param_dict.items():
        #Change model parameter values to current location in parameter space
        model.parameters[pname].value = (pvalue)
    simres = ScipyOdeSimulator(model, t).run() 
    yout = simres.all
    a, maxima = peakfinder(yout['obsCdc20'][40:120],[],0.2,includeEndpoints=1)   #find the peak of cdc20 - last protein to peak
    b, maxima_b = peakfinder(yout['obsCycE'][40:120],[],0.2,includeEndpoints=1)   #find the peak of cdc20 - last protein to peak
    if len(a)>1 and len(b)>1:
        accepted+=[parameter_vector]
    else :
        rejected+=[parameter_vector]



np.save(r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code\fit_data\sobol_seq_set2_again_accepted.npy", accepted)
np.save(r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code\fit_data\sobol_seq_set2_again_rejected.npy", rejected)







l=sobol_seq.i4_sobol_generate(10, 25000, skip=75000)*2
t = np.linspace(0, 200, 120)
accepted=[]
rejected=[]

for parameter_vector in l :
    param_dict = {pname: pvalue for pname, pvalue in zip(pysb_sampled_parameter_names, parameter_vector)}
    for pname, pvalue in param_dict.items():
        #Change model parameter values to current location in parameter space
        model.parameters[pname].value = (pvalue)
    simres = ScipyOdeSimulator(model, t).run() 
    yout = simres.all
    a, maxima = peakfinder(yout['obsCdc20'][40:120],[],0.2,includeEndpoints=1)   #find the peak of cdc20 - last protein to peak
    b, maxima_b = peakfinder(yout['obsCycE'][40:120],[],0.2,includeEndpoints=1)   #find the peak of cdc20 - last protein to peak
    if len(a)>1 and len(b)>1:
        accepted+=[parameter_vector]
    else :
        rejected+=[parameter_vector]
    print(len(accepted)+len(rejected))



np.save(r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code\fit_data\sobol_seq_set3_again_accepted.npy", accepted)
np.save(r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code\fit_data\sobol_seq_set3_again_rejected.npy", rejected)






l=sobol_seq.i4_sobol_generate(10, 50000, skip=100000)*2
t = np.linspace(0, 200, 120)
accepted=[]
rejected=[]

for parameter_vector in l :
    param_dict = {pname: pvalue for pname, pvalue in zip(pysb_sampled_parameter_names, parameter_vector)}
    for pname, pvalue in param_dict.items():
        #Change model parameter values to current location in parameter space
        model.parameters[pname].value = (pvalue)
    simres = ScipyOdeSimulator(model, t).run() 
    yout = simres.all
    a, maxima = peakfinder(yout['obsCdc20'][40:120],[],0.2,includeEndpoints=1)   #find the peak of cdc20 - last protein to peak
    b, maxima_b = peakfinder(yout['obsCycE'][40:120],[],0.2,includeEndpoints=1)   #find the peak of cdc20 - last protein to peak
    if len(a)>1 and len(b)>1:
        accepted+=[parameter_vector]
    else :
        rejected+=[parameter_vector]
    print(len(accepted)+len(rejected))



np.save("/home/jeremie/fit_data/sobol_seq_set4_again_accepted.npy", accepted)
np.save("/home/jeremie/fit_data/sobol_seq_set4_again_rejected.npy", rejected)




l=sobol_seq.i4_sobol_generate(10, 25000, skip=150000)*2
t = np.linspace(0, 200, 120)
accepted=[]
rejected=[]

for parameter_vector in l :
    param_dict = {pname: pvalue for pname, pvalue in zip(pysb_sampled_parameter_names, parameter_vector)}
    for pname, pvalue in param_dict.items():
        #Change model parameter values to current location in parameter space
        model.parameters[pname].value = (pvalue)
    simres = ScipyOdeSimulator(model, t).run() 
    yout = simres.all
    a, maxima = peakfinder(yout['obsCdc20'][40:120],[],0.2,includeEndpoints=1)   #find the peak of cdc20 - last protein to peak
    b, maxima_b = peakfinder(yout['obsCycE'][40:120],[],0.2,includeEndpoints=1)   #find the peak of cdc20 - last protein to peak
    if len(a)>1 and len(b)>1:
        accepted+=[parameter_vector]
    else :
        rejected+=[parameter_vector]
    print(len(accepted)+len(rejected))



np.save(r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code\fit_data\sobol_seq_set5_again_accepted.npy", accepted)
np.save(r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code\fit_data\sobol_seq_set5_again_rejected.npy", rejected)




l=sobol_seq.i4_sobol_generate(10, 50000, skip=175000)*2
t = np.linspace(0, 200, 120)
accepted=[]
rejected=[]

for parameter_vector in l :
    param_dict = {pname: pvalue for pname, pvalue in zip(pysb_sampled_parameter_names, parameter_vector)}
    for pname, pvalue in param_dict.items():
        #Change model parameter values to current location in parameter space
        model.parameters[pname].value = (pvalue)
    simres = ScipyOdeSimulator(model, t).run() 
    yout = simres.all
    a, maxima = peakfinder(yout['obsCdc20'][40:120],[],0.2,includeEndpoints=1)   #find the peak of cdc20 - last protein to peak
    b, maxima_b = peakfinder(yout['obsCycE'][40:120],[],0.2,includeEndpoints=1)   #find the peak of cdc20 - last protein to peak
    if len(a)>1 and len(b)>1:
        accepted+=[parameter_vector]
    else :
        rejected+=[parameter_vector]
    print(len(accepted)+len(rejected))



np.save("/home/jeremie/fit_data/sobol_seq_set6_again_accepted.npy", accepted)
np.save("/home/jeremie/fit_data/sobol_seq_set6_again_rejected.npy", rejected)
