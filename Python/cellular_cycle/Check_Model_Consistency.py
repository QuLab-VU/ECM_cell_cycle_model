# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 17:48:30 2020

@author: Grim.M
Title: runStochCell
"""

import sys
import os
sys.path.insert(0, os.path.dirname(__file__) + "/../cellular_cycle")

# import necessary packages
from peakfinder import peakfinder
from pysb.simulator import ScipyOdeSimulator
from pysb.simulator import BngSimulator
import numpy as np
import scipy.stats as ss

import matplotlib.pyplot as plt

def Check_Model_Consistency(model, cycling_species = ['obsCdc20', 'obsCycE', 'obsCycD', 'obsp21'], thr = [0.26, 0.18], return_plot = False):
    """
    Run a stochastic simulation and perform multiple checks on the simulation results.

    Args:
        model: the PySB model.
        thr: the threshold value for the peakfinder algorithm.

    Returns:
        variable_range: a dictionary that contains the names of the observables and their value ranges.
        check: a boolean value that indicates if the simulation is successful.
        tmp_score: score for the alternance between Cdc20 and Cyclin E
    """

    # # print Skp2_0 * omega
    # print("Cdc20_0 * omega: ")
    # print(model.parameters['Cdc20_0'].value * model.parameters['omega'].value)
    # print("\n")
    # # Print initial concentration for Cdc20
    # print("Cdc20_0: ")
    # print(model.initial_conditions[18])
    # # model.expressions['Cdc20_0_expr']

    # Run deterministic simulation with ScipyOdeSimulator
    t_ode = np.linspace(0, 1000, 1000)  # Time points for ODE simulation
    simres = ScipyOdeSimulator(model, t_ode, verbose=False).run()
    yout = simres.all  # get simulation results

    if return_plot:
        plt.figure()
        for species in cycling_species:
            plt.plot(t_ode, yout[species]/model.parameters['omega'].value, label=species)
        # plt.plot(t_ode, yout['obsCycE']/model.parameters['omega'].value, label='CycE')
        # plt.plot(t_ode, yout['obsCycA']/model.parameters['omega'].value, label='CycA')
        # plt.plot(t_ode, yout['obsCycB']/model.parameters['omega'].value, label='CycB')
        plt.legend()
        plt.xlabel('Time (min)')
        plt.ylabel('Concentration (nM)')
        plt.title('ODE simulation')
        plt.show()



    
    ######################################
    # Check 1: The system doesn't explode
    for obs in model.observables:
        # If any concentration exceeds the threshold, return False and error code 1
        if np.any(yout[obs.name]/model.parameters['omega'].value > 15):
            return {}, False, 1
    
    ######################################
    # Check 2: Are they oscillating 
    T = t_ode[1] - t_ode[0]  # Sampling interval 
    N = yout['obsCycD'][600:].size  # Number of sample points for FFT
    f = np.linspace(0, 1 / T, N)  # Frequency axis for FFT
    fft = {}
    for species in cycling_species:
        fft[species] = np.fft.fft(yout[species][600:]/model.parameters['omega'].value)

    # Find frequency corresponding to max amplitude of FFT
    arg_freq_max = {}
    for species in cycling_species:
        arg_freq_max[species] = np.argmax(np.abs(fft[species][(f > 0.005) * (f < 0.1)]) * 1 / N)

    # Plot the frequency spectrum as well as the arg with the corresponding period
    if return_plot:
        # add scatter plot for the max frequency, with x the frequency and y the amplitude
        plt.figure()
        for species in cycling_species:
            plt.plot(f[(f > 0.005) * (f < 0.1)], np.abs(fft[species][(f > 0.005) * (f < 0.1)]), label=species)
            plt.scatter(f[(f > 0.005) * (f < 0.1)][arg_freq_max[species]], np.abs(fft[species][(f > 0.005) * (f < 0.1)])[arg_freq_max[species]], s=100)
        plt.xlabel('Frequency ()')
        plt.ylabel('Amplitude')
        plt.legend()
        plt.title('Frequency spectrum')
        plt.show()

        # print('Frequency of oscillations in Cdc20: ' + str(f[arg_freq_max_Cdc20]) + ' Hz')
        # print('Frequency of oscillations in CycD: ' + str(f[arg_freq_max_CycD]) + ' Hz')
        # print('Frequency of oscillations in CycE: ' + str(f[arg_freq_max_CycE]) + ' Hz')

    ######################################
    # Check 3: If frequency of oscillations is in a biological range
    t = np.linspace(0, 1000, 1000)  # Time points for stochastic simulation
    simres_2 = BngSimulator(model, t, verbose=False).run()

    if return_plot:
        plt.figure()
        for species in cycling_species:
            plt.plot(t, simres_2.observables[species], label=species)
        # plt.plot(t, simres_2.observables['obsCycD'], label='CycD')
        # plt.plot(t, simres_2.observables['obsCycE'], label='CycE')
        # plt.plot(t, simres_2.observables['obsCycA'], label='CycA')
        # plt.plot(t, simres_2.observables['obsCycB'], label='CycB')
        plt.legend()
        plt.xlabel('Time (min)')
        plt.ylabel('Concentration (nM)')
        plt.title('SSA simulation')
        plt.show()

    check = False
    # All frequencies should be between 0.005 and 0.1 Hz

    for species in cycling_species:
        if arg_freq_max[species] > 1 and \
            np.abs(fft[species][(f > 0.005) * (f < 0.1)])[arg_freq_max[species]] > np.abs(fft[species][(f > 0.005) * (f < 0.1)])[arg_freq_max[species] - 1] and \
            np.abs(fft[species][(f > 0.005) * (f < 0.1)])[arg_freq_max[species]] > np.abs(fft[species][(f > 0.005) * (f < 0.1)])[arg_freq_max[species] + 1]:
            check = True
        else:
            check = False
            break

    # if arg_freq_max_Cdc20>1 and \
    #     arg_freq_max_CycD>1 and \
    #     arg_freq_max_CycE>1 and \
    #     np.abs(fft_Cdc20[(f>0.005)*(f<0.1)])[arg_freq_max_Cdc20]>np.abs(fft_Cdc20[(f>0.005)*(f<0.1)])[arg_freq_max_Cdc20-1] and \
    #     np.abs(fft_Cdc20[(f>0.005)*(f<0.1)])[arg_freq_max_Cdc20]>np.abs(fft_Cdc20[(f>0.005)*(f<0.1)])[arg_freq_max_Cdc20+1] and \
    #     np.abs(fft_CycD[(f>0.005)*(f<0.1)])[arg_freq_max_CycD]>np.abs(fft_CycD[(f>0.005)*(f<0.1)])[arg_freq_max_CycD-1] and \
    #     np.abs(fft_CycD[(f>0.005)*(f<0.1)])[arg_freq_max_CycD]>np.abs(fft_CycD[(f>0.005)*(f<0.1)])[arg_freq_max_CycD+1] and \
    #     np.abs(fft_CycE[(f>0.005)*(f<0.1)])[arg_freq_max_CycE]>np.abs(fft_CycE[(f>0.005)*(f<0.1)])[arg_freq_max_CycE-1] and \
    #     np.abs(fft_CycE[(f>0.005)*(f<0.1)])[arg_freq_max_CycE]>np.abs(fft_CycE[(f>0.005)*(f<0.1)])[arg_freq_max_CycE+1] :
    #     check=True
    
    
    ######################################
    # Check 4: Biological significance of the cycle, cdc20 and Cyclin E should alternate
    # Run stochastic simulation with BngSimulator

    sys_size = model.parameters['omega'].value
    # Find the maximum of the oscillation peaks in Cdc20 and CyeE
    (a,max_tmp_cdc) = peakfinder(simres_2.observables['obsCdc20'], thr[0]*sys_size, thr[0]*sys_size, 1, includeEndpoints=1);
    (b,max_tmp_cyce) = peakfinder(simres_2.observables['obsCycE'], thr[1]*sys_size, thr[1]*sys_size, 1, includeEndpoints=1);
   
    
    tmp_score=1 #score on how good they alternate, it may not be exactly one because of the stochasticity of the network, and peakfinder mught miss one peak
    variable_range={}
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
            #simres = ScipyOdeSimulator(model, t_ode, verbose=False).run()
            (a_ode,max_tmp_cdc_ode) = peakfinder(yout[model.observables['obsCdc20'].name], thr[0]*sys_size, thr[0]*sys_size, 1, includeEndpoints=1);
            (b_ode,max_tmp_cyce_ode) = peakfinder(yout[model.observables['obsCycE'].name], thr[1]*sys_size, thr[1]*sys_size, 1, includeEndpoints=1);
            
            
            for k in range(len(model.observables)) :
                variable_range[model.observables[k].name]=[yout[model.observables[k].name][i] for i in range(max_tmp_cdc_ode[-3], max_tmp_cdc_ode[-2])]

    # Plot on a first plot the trajectory of Cdc20 and CycE, with scatter plots on peaks of each with different colors. On a lower plots the scatter aligned on a line showing alternance between Cdc20 and CycE, and show the score
    if return_plot:
        # Top plot
        plt.figure()
        ax = plt.subplot(211)
        ax.plot(t, simres_2.observables['obsCdc20'], label='Cdc20')
        ax.plot(t, simres_2.observables['obsCycE'], label='CycE')
        ax.scatter(t[max_tmp_cdc], simres_2.observables['obsCdc20'][max_tmp_cdc], c='r')
        ax.scatter(t[max_tmp_cyce], simres_2.observables['obsCycE'][max_tmp_cyce], c='g')
        ax.legend()
        ax.set_xlabel('Time (min)')
        ax.set_ylabel('Concentration (nM)')
        ax.set_title('Alternance between Cdc20 and CycE')

        # Bottom plot
        ax = plt.subplot(212)
        ax.scatter(t[max_tmp_cdc], np.ones(len(max_tmp_cdc)), c='r', label='Cdc20', s=8)
        ax.scatter(t[max_tmp_cyce], np.ones(len(max_tmp_cyce)), c='g', label='CycE', s=8)
        # remove y ticks
        ax.set_yticks([])
        ax.set_xticks([])
        #xlim
        ax.set_xlim([t[0], t[-1]])
        # Remove the top and right axes from the data plot
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.set_title('Score: ' + str(tmp_score))
        plt.tight_layout()
        plt.show()

    return variable_range, check, tmp_score

###############################################################################
# Test Check_Model_Consistency

if __name__=="__main__":
    # import model
    from cell_cycle_death_model import model
    import pandas as pd

    # import pysb sampled parameters
    sampled_parameter_names = np.loadtxt(os.path.dirname(__file__) + '/../model_data/sampled_parameter_names.txt', dtype=str)

    # Importing the parameters:
    new_param = pd.read_csv(os.path.dirname(__file__) + '/../model_data/parameters.csv')
    del new_param['Stress']
    for key in new_param.keys():
        try:
            if model.parameters[key].value != new_param[key][0]:
                model.parameters[key].value = new_param[key]
        except:
            print(key)
    model.parameters['omega'].value = 200

    # Setting initial concentrations '/../model_data/initial_concentrations.csv'
    initial_concentrations = pd.read_csv(os.path.dirname(__file__) + '/../model_data/initial_concentrations.csv')
    for key in initial_concentrations.keys():
        try:
            model.parameters[key].value = initial_concentrations[key][0] / model.parameters['omega'].value
        except:
            print(key)

    # # Import ECM_conditions.csv
    ECM_conditions = pd.read_csv(os.path.dirname(__file__) + '/../model_data/ECM_conditions_2.csv')
    condition = "4 kPa erlotinib"
    # Set the species AKT, ERK, FAK with species and value col
    for input in ["AKT", "ERK", "FAK"]:
        model.parameters[input].value = ECM_conditions[(ECM_conditions["condition"] == condition) & (ECM_conditions["species"] == input)].value

    variable_range, check, tmp_score = Check_Model_Consistency(model, return_plot=True)

    print(variable_range)
    print(check)
    print(tmp_score)

    condition = "4 kPa DMSO"
    # Set the species AKT, ERK, FAK with species and value col
    for input in ["AKT", "ERK", "FAK"]:
        model.parameters[input].value = ECM_conditions[(ECM_conditions["condition"] == condition) & (ECM_conditions["species"] == input)].value

    variable_range, check, tmp_score = Check_Model_Consistency(model, return_plot=True)
    print(condition)
    print(check)

    condition = "55 kPa erlotinib"
    # Set the species AKT, ERK, FAK with species and value col
    for input in ["AKT", "ERK", "FAK"]:
        model.parameters[input].value = ECM_conditions[(ECM_conditions["condition"] == condition) & (ECM_conditions["species"] == input)].value

    variable_range, check, tmp_score = Check_Model_Consistency(model, return_plot=True)
    print(condition)
    print(check)

    
    condition = "55 kPa DMSO"
    # Set the species AKT, ERK, FAK with species and value col
    for input in ["AKT", "ERK", "FAK"]:
        model.parameters[input].value = ECM_conditions[(ECM_conditions["condition"] == condition) & (ECM_conditions["species"] == input)].value

    variable_range, check, tmp_score = Check_Model_Consistency(model, return_plot=True)
    print(condition)
    print(check)


    condition = "0.6 kPa erlotinib"
    # Set the species AKT, ERK, FAK with species and value col
    for input in ["AKT", "ERK", "FAK"]:
        model.parameters[input].value = ECM_conditions[(ECM_conditions["condition"] == condition) & (ECM_conditions["species"] == input)].value

    variable_range, check, tmp_score = Check_Model_Consistency(model, return_plot=True)
    print(condition)
    print(check)


    condition = "0.6 kPa DMSO"
    # Set the species AKT, ERK, FAK with species and value col
    for input in ["AKT", "ERK", "FAK"]:
        model.parameters[input].value = ECM_conditions[(ECM_conditions["condition"] == condition) & (ECM_conditions["species"] == input)].value

    variable_range, check, tmp_score = Check_Model_Consistency(model, return_plot=True)
    print(condition)
    print(check)
