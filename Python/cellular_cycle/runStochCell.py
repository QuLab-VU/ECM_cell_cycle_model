# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 17:48:30 2020

@author: Grim.M
Title: runStochCell
"""

import sys
sys.path.insert(0, "../cellular_cycle") # add a pathway to the model

# import necessary packages
from peakfinder import peakfinder
from pysb.simulator import ScipyOdeSimulator
from pysb.simulator import BngSimulator
import numpy as np
import scipy.stats as ss

def runStochcell(model, index, condition, thr):
    """
    Run a stochastic simulation and perform multiple checks on the simulation results.

    Args:
        model: the PySB model.
        index: the index of the sets of tested parameters.
        condition: treated with ERL or untreated (DMSO) and the stiffness of the media (to be removed?, no use here).
        thr: the threshold value for the peakfinder algorithm.

    Returns:
        variable_range: a dictionary that contains the names of the observables and their value ranges.
        check: a boolean value that indicates if the simulation is successful.
        tmp_score: score for the alternance between Cdc20 and Cyclin E
    """
    
    # Run deterministic simulation with ScipyOdeSimulator
    t_ode = np.linspace(0, 1000, 1000)  # Time points for ODE simulation
    simres = ScipyOdeSimulator(model, t_ode, verbose=False).run()
    yout = simres.all  # get simulation results
    
    # Run stochastic simulation with BngSimulator
    t = np.linspace(0, 1000, 1000)  # Time points for stochastic simulation
    simres_2 = BngSimulator(model, t, verbose=False).run()
    
    
    ######################################
    # Check 1: The system doesn't explode
    for obs in model.observables:
        print(np.any(yout[obs.name]/model.parameters['omega'] > 15))  # Check if any concentration exceeds the threshold
        
        # If any concentration exceeds the threshold, return False and error code 1
        if np.any(yout[obs.name]/model.parameters['omega'] > 15):
            return {}, False, 1
    
    ######################################
    # Check 2: Are they oscillating 
    T = t_ode[1] - t_ode[0]  # Sampling interval 
    N = yout['obsCycD'][600:].size  # Number of sample points for FFT
    f = np.linspace(0, 1 / T, N)  # Frequency axis for FFT
    
    fft_Cdc20 = np.fft.fft(yout['obsCdc20'][600:]/model.parameters['omega'])
    fft_CycD = np.fft.fft(yout['obsCycD'][600:]/model.parameters['omega'])
    fft_CycE = np.fft.fft(yout['obsCycE'][600:]/model.parameters['omega'])
    
    # Find frequency corresponding to max amplitude of FFT
    arg_freq_max_Cdc20 = np.argmax(np.abs(fft_Cdc20[(f > 0.005) * (f < 0.1)]) * 1 / N)
    arg_freq_max_CycE = np.argmax(np.abs(fft_CycE[(f > 0.005) * (f < 0.1)]) * 1 / N)
    arg_freq_max_CycD = np.argmax(np.abs(fft_CycD[(f > 0.005) * (f < 0.1)]) * 1 / N)  # should be same as arg_freq_max_Cdc20
    
    
    ######################################
    # Check 3: If frequency of oscillations is in a biological range
    check = False
    if arg_freq_max_Cdc20>1 and arg_freq_max_CycD>1 and arg_freq_max_CycE>1 and np.abs(fft_Cdc20[(f>0.005)*(f<0.1)])[arg_freq_max_Cdc20]>np.abs(fft_Cdc20[(f>0.005)*(f<0.1)])[arg_freq_max_Cdc20-1] and np.abs(fft_Cdc20[(f>0.005)*(f<0.1)])[arg_freq_max_Cdc20]>np.abs(fft_Cdc20[(f>0.005)*(f<0.1)])[arg_freq_max_Cdc20+1] and np.abs(fft_CycD[(f>0.005)*(f<0.1)])[arg_freq_max_CycD]>np.abs(fft_CycD[(f>0.005)*(f<0.1)])[arg_freq_max_CycD-1] and np.abs(fft_CycD[(f>0.005)*(f<0.1)])[arg_freq_max_CycD]>np.abs(fft_CycD[(f>0.005)*(f<0.1)])[arg_freq_max_CycD+1] and np.abs(fft_CycE[(f>0.005)*(f<0.1)])[arg_freq_max_CycE]>np.abs(fft_CycE[(f>0.005)*(f<0.1)])[arg_freq_max_CycE-1] and np.abs(fft_CycE[(f>0.005)*(f<0.1)])[arg_freq_max_CycE]>np.abs(fft_CycE[(f>0.005)*(f<0.1)])[arg_freq_max_CycE+1] :
        check=True
    
    
    ######################################
    # Check 4: Biological significance of the cycle, cdc20 and Cyclin E should alternate
    sys_size = 200  # System size
    # Find the maximum of the oscillation peaks in Cdc20 and CyeE
    (a,max_tmp_cdc) = peakfinder(simres_2.observables['obsCdc20'], thr*sys_size, thr*sys_size, 1, includeEndpoints=1);
    (b,max_tmp_cyce) = peakfinder(simres_2.observables['obsCycE'], thr*sys_size, thr*sys_size, 1, includeEndpoints=1);
    
    tmp_score=1 #score on how good they alternate, it may not be exactly one because of the stochasticity of the network, and peakfinder mught miss one peak
    variable_range={}
    print(max_tmp_cdc)
    if len(max_tmp_cdc)>2 and len(max_tmp_cyce)>2 :
        ranks = ss.rankdata(np.concatenate([max_tmp_cdc, max_tmp_cyce]))
        max_tmp = np.concatenate([['obsCdc20' for i in range(len(max_tmp_cdc))],['obsCycE' for i in range(len(max_tmp_cyce))]])
        Z = [x for _,x in sorted(zip(ranks,max_tmp))]
        tmp_score=0
        
        for i in range(len(Z)-1):
            if Z[i]==Z[i+1]:
                tmp_score+=1
        tmp_score = tmp_score/len(Z)
        
        if len(a)>=2 and len(b)>=2:
            for k in range(len(model.observables)) :
                variable_range[model.observables[k].name]=[yout[model.observables[k].name][i]/model.parameters['omega'] for i in range(max_tmp_cdc[-2], max_tmp_cdc[-1])]
    
    return variable_range, check, tmp_score

