# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 17:48:30 2020

@author: Grim.M
Title: runStochCell
"""

import sys
sys.path.insert(0, "../cellular_cycle") #add a pathway to the model
import os

from peakfinder import peakfinder
from pysb.simulator import ScipyOdeSimulator
from pysb.simulator import BngSimulator
import numpy as np
import pylab as pl
import random
import scipy.stats as ss



#model=model_concentration

def runStochcell(model_population, model_concentration, index, condition) :
    
    v = 0.05;     
    sys_size=200;
    run_time=500;    
    CI=0;
    t = np.linspace(0, 800, 800)
    
    #print("ode simulation starting")
    simres = ScipyOdeSimulator(model_concentration, t, verbose=False).run();
    yout = simres.all;
    
    simres_2 = BngSimulator(model_population, t, verbose=False).run();
    
    for obs in model_concentration.observables : 
        #pl.plot(yout[obs.name])
        #print(obs.name)
        #print(np.any(yout[obs.name]>15)) #concentration too high
        
        if np.any(yout[obs.name]>15) :
            return {}, False, 1
    
    """
    a, maxima = peakfinder(yout['obsCdc20'], [], 0.2);   #find the peak of cdc20 - last protein to peak
    b, maxima_b = peakfinder(yout['obsCycD'], [], 0.2);   #find the peak of cdc20 - last protein to peak
    
    variable_range={};
    
    if len(a)>=2 and len(b)>=2:
        for k in range(len(model_concentration.observables)) :
            variable_range[model_concentration.observables[k].name]=[yout[model_concentration.observables[k].name][i] for i in range(maxima[4],maxima[5])]
    """
   #Are they oscillating 
    
    T = t[1] - t[0]  # sampling interval 
    N = yout['obsCycD'][600:].size
    # 1/T = frequency
    f = np.linspace(0, 1 / T, N)
    
    fft_Cdc20 = np.fft.fft(yout['obsCdc20'][600:])
    fft_CycD = np.fft.fft(yout['obsCycD'][600:])
    fft_CycE = np.fft.fft(yout['obsCycE'][600:])
    
    f[(f>0.005)*(f<0.1)]
    [(f>0.005)*(f<0.1)]
    
    arg_freq_max_Cdc20 = np.argmax(np.abs(fft_Cdc20[(f>0.005)*(f<0.1)]) * 1 / N)
    arg_freq_max_CycE = np.argmax(np.abs(fft_CycE[(f>0.005)*(f<0.1)]) * 1 / N)
    arg_freq_max_CycD = np.argmax(np.abs(fft_CycD[(f>0.005)*(f<0.1)]) * 1 / N) #they should be the same, they have the same period
    
    check=False
    
    #print(check)
    
    if arg_freq_max_Cdc20>1 and arg_freq_max_CycD>1 and arg_freq_max_CycE>1 and np.abs(fft_Cdc20[(f>0.005)*(f<0.1)])[arg_freq_max_Cdc20]>np.abs(fft_Cdc20[(f>0.005)*(f<0.1)])[arg_freq_max_Cdc20-1] and np.abs(fft_Cdc20[(f>0.005)*(f<0.1)])[arg_freq_max_Cdc20]>np.abs(fft_Cdc20[(f>0.005)*(f<0.1)])[arg_freq_max_Cdc20+1] and np.abs(fft_CycD[(f>0.005)*(f<0.1)])[arg_freq_max_CycD]>np.abs(fft_CycD[(f>0.005)*(f<0.1)])[arg_freq_max_CycD-1] and np.abs(fft_CycD[(f>0.005)*(f<0.1)])[arg_freq_max_CycD]>np.abs(fft_CycD[(f>0.005)*(f<0.1)])[arg_freq_max_CycD+1] and np.abs(fft_CycE[(f>0.005)*(f<0.1)])[arg_freq_max_CycE]>np.abs(fft_CycE[(f>0.005)*(f<0.1)])[arg_freq_max_CycE-1] and np.abs(fft_CycE[(f>0.005)*(f<0.1)])[arg_freq_max_CycE]>np.abs(fft_CycE[(f>0.005)*(f<0.1)])[arg_freq_max_CycE+1] :
        
        #Not necessary here
        check=True
        #print("Check=")
        #print(check)

        
    (a,max_tmp_cdc) = peakfinder(simres_2.observables['obsCdc20'], 0.15*sys_size, 0.15*sys_size, 1, includeEndpoints=1);
    (b,max_tmp_cyce) = peakfinder(simres_2.observables['obsCycE'], 0.15*sys_size, 0.15*sys_size, 1, includeEndpoints=1);
    
    tmp_score=1
    #print(tmp_score)
    #print(len(max_tmp_cdc))
    #print(len(max_tmp_cdc))
    
    if len(max_tmp_cdc)>2 and len(max_tmp_cyce)>2 :
        ranks = ss.rankdata(np.concatenate([max_tmp_cdc, max_tmp_cyce]))
        max_tmp = np.concatenate([['obsCdc20' for i in range(len(max_tmp_cdc))],['obsCycE' for i in range(len(max_tmp_cyce))]])
        
        Z = [x for _,x in sorted(zip(ranks,max_tmp))]
        
        tmp_score=0
        
        for i in range(len(Z)-1):
            if Z[i]==Z[i+1]:
                tmp_score+=1
        tmp_score = tmp_score/len(Z)
        
        
        variable_range={}
        
        
        if len(a)>=2 and len(b)>=2:
            for k in range(len(model_concentration.observables)) :
                variable_range[model_concentration.observables[k].name]=[yout[model_concentration.observables[k].name][i] for i in range(max_tmp_cdc[4],max_tmp_cdc[5])]
        
        
    #variable_range = np.transpose(np.array(variable_range))
    
    #print("tmp_score is")
    #print(tmp_score)
    #print(path)
    
#    if check==True and tmp_score<0.15:
#        print(path + "/fit_data/sobol_seq/cost_function_model_5/biochemical_data/ODE_timecourse_CycE_index_" + str(index) + "_" + condition.replace(' ','_').replace('.','0'))
#        np.save(path + "/fit_data/sobol_seq/cost_function_model_5/biochemical_data/ODE_timecourse_CycE_index_" + str(index) + "_" + condition.replace(' ','_').replace('.','0'), simres.observables['obsCycE'])
#        np.save(path + "/fit_data/sobol_seq/cost_function_model_5/biochemical_data/ODE_timecourse_CycD_index_" + str(index) + "_" + condition.replace(' ','_').replace('.','0'), simres.observables['obsCycD'])
#        np.save(path + "/fit_data/sobol_seq/cost_function_model_5/biochemical_data/ODE_timecourse_Cdc20_index_" + str(index) + "_" + condition.replace(' ','_').replace('.','0'), simres.observables['obsCdc20'])
#        
#        np.save(path + "/fit_data/sobol_seq/cost_function_model_5/biochemical_data/ssa_timecourse_CycE_index_" + str(index) + "_" + condition.replace(' ','_').replace('.','0'), simres_2.observables['obsCycE'])
#        np.save(path + "/fit_data/sobol_seq/cost_function_model_5/biochemical_data/ssa_timecourse_CycD_index_" + str(index) + "_" + condition.replace(' ','_').replace('.','0'), simres_2.observables['obsCycD'])
#        np.save(path + "/fit_data/sobol_seq/cost_function_model_5/biochemical_data/ssa_timecourse_Cdc20_index_" + str(index) + "_" + condition.replace(' ','_').replace('.','0'), simres_2.observables['obsCdc20'])
#        
#        print(path + "/fit_data/pso/high_and_low_stiffness_cost/biochemical_data/ODE_timecourse_CycE_index_" + str(index) + "_" + condition.replace(' ','_').replace('.','0'))
#        np.save(path + "/fit_data/pso/high_and_low_stiffness_cost/biochemical_data/ODE_timecourse_CycE_index_" + str(index) + "_" + condition.replace(' ','_').replace('.','0'), simres.observables['obsCycE'])
#        np.save(path + "/fit_data/pso/high_and_low_stiffness_cost/biochemical_data/ODE_timecourse_CycD_index_" + str(index) + "_" + condition.replace(' ','_').replace('.','0'), simres.observables['obsCycD'])
#        np.save(path + "/fit_data/pso/high_and_low_stiffness_cost/biochemical_data/ODE_timecourse_Cdc20_index_" + str(index) + "_" + condition.replace(' ','_').replace('.','0'), simres.observables['obsCdc20'])
#        
#        np.save(path + "/fit_data/pso/high_and_low_stiffness_cost/biochemical_data/ssa_timecourse_CycE_index_" + str(index) + "_" + condition.replace(' ','_').replace('.','0'), simres_2.observables['obsCycE'])
#        np.save(path + "/fit_data/pso/high_and_low_stiffness_cost/biochemical_data/ssa_timecourse_CycD_index_" + str(index) + "_" + condition.replace(' ','_').replace('.','0'), simres_2.observables['obsCycD'])
#        np.save(path + "/fit_data/pso/high_and_low_stiffness_cost/biochemical_data/ssa_timecourse_Cdc20_index_" + str(index) + "_" + condition.replace(' ','_').replace('.','0'), simres_2.observables['obsCdc20'])
#        
    
    #print(tmp_score)
    #print(check)
    
    return variable_range, check, tmp_score
    
    
def runStochcell_check(ssa, ode, index) : #doesn't include check for too hgih concentration! 
    check=False
    tmp_score=1

    v = 0.05;     
    sys_size=200;
    run_time=500;    
    CI=0;
    
    t = np.linspace(0, 800, 800)
    
   #Are they oscillating 
    
    T = t[1] - t[0]  # sampling interval 
    N = ode['obsCycD'][600:].size
    # 1/T = frequency
    f = np.linspace(0, 1 / T, N)
    
    fft_Cdc20 = np.fft.fft(ode['obsCdc20'][600:])
    fft_CycD = np.fft.fft(ode['obsCycD'][600:])
    fft_CycE = np.fft.fft(ode['obsCycE'][600:])
    
    f[(f>0.005)*(f<0.1)]
    
    arg_freq_max_Cdc20 = np.argmax(np.abs(fft_Cdc20[(f>0.005)*(f<0.1)]) * 1 / N)
    arg_freq_max_CycE = np.argmax(np.abs(fft_CycE[(f>0.005)*(f<0.1)]) * 1 / N)
    arg_freq_max_CycD = np.argmax(np.abs(fft_CycD[(f>0.005)*(f<0.1)]) * 1 / N) #they should be the same, they have the same period
    
    
    if arg_freq_max_Cdc20>1 and arg_freq_max_CycD>1 and arg_freq_max_CycE>1 and np.abs(fft_Cdc20[(f>0.005)])[arg_freq_max_Cdc20] > np.abs(fft_Cdc20[(f>0.005)])[arg_freq_max_Cdc20-1] and np.abs(fft_Cdc20[(f>0.005)])[arg_freq_max_Cdc20]>np.abs(fft_Cdc20[(f>0.005)])[arg_freq_max_Cdc20+1] and np.abs(fft_CycD[(f>0.005)])[arg_freq_max_CycD]>np.abs(fft_CycD[(f>0.005)])[arg_freq_max_CycD-1] and np.abs(fft_CycD[(f>0.005)])[arg_freq_max_CycD]>np.abs(fft_CycD[(f>0.005)])[arg_freq_max_CycD+1] and np.abs(fft_CycE[(f>0.005)])[arg_freq_max_CycE]>np.abs(fft_CycE[(f>0.005)])[arg_freq_max_CycE-1] and np.abs(fft_CycE[(f>0.005)])[arg_freq_max_CycE]>np.abs(fft_CycE[(f>0.005)])[arg_freq_max_CycE+1] :
        check=True
        
    (a,max_tmp_cdc) = peakfinder(ssa['obsCdc20'], 0.15*sys_size, 0.15*sys_size, 1, includeEndpoints=1);
    (a,max_tmp_cyce) = peakfinder(ssa['obsCycE'], 0.15*sys_size, 0.15*sys_size, 1, includeEndpoints=1);
        
    if len(max_tmp_cdc)>2 and len(max_tmp_cyce)>2:
        ranks = ss.rankdata(np.concatenate([max_tmp_cdc, max_tmp_cyce]))
        max_tmp = np.concatenate([['obsCdc20' for i in range(len(max_tmp_cdc))],['obsCycE' for i in range(len(max_tmp_cyce))]])
        
        Z = [x for _,x in sorted(zip(ranks,max_tmp))]
        
        tmp_score=0
        
        for i in range(len(Z)-1):
            if Z[i]==Z[i+1]:
                tmp_score+=1
        tmp_score = tmp_score/len(Z)
        
    return check, tmp_score
    

