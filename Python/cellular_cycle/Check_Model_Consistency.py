# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 17:48:30 2020

@author: Grim.M
Title: runStochCell
"""
# import necessary packages
from Python.cellular_cycle.peakfinder import peakfinder
from pysb.simulator import ScipyOdeSimulator
from pysb.simulator import BngSimulator
import numpy as np
import scipy.stats as ss
import matplotlib.pyplot as plt


def check_model_consistency(model, t=None, cycling_species=None, thr=None,
                            return_plot=False, save_path=None):
    """
    Run a stochastic simulation and perform multiple checks on the simulation results.

    Args:
        model: the PySB model.
        thr: the threshold value for the peakfinder algorithm.
        cycling_species: the species that are expected to oscillate.
        return_plot: a boolean value that indicates if the plots should be returned.
        t: the time points for the simulation.
        save_path: the path where the plots should be saved.

    Returns:
        variable_range: a dictionary that contains the names of the observables and their value ranges.
        check: a boolean value that indicates if the simulation is successful.
        tmp_score: score for the alternance between Cdc20 and Cyclin E
    """

    if thr is None:
        thr = [0.26, 0.18]
    if cycling_species is None:
        cycling_species = ['obsCdc20', 'obsCycE', 'obsCycD', 'obsp21']
    if t is None:
        t = 1000
    if save_path is None:
        save_path = ''

    # Run deterministic simulation with ScipyOdeSimulator
    if cycling_species is None:
        cycling_species = ['obsCdc20', 'obsCycE', 'obsCycD', 'obsp21']
    t_ode = np.linspace(0, t, t)  # Time points for ODE simulation
    simres = ScipyOdeSimulator(model, t_ode, verbose=False, compiler='python').run()
    yout = simres.all  # get simulation results

    if return_plot:
        plt.figure()
        for species in cycling_species:
            plt.plot(t_ode, yout[species] / model.parameters['omega'].value, label=species)
        plt.legend()
        plt.xlabel('Time (min)')
        plt.ylabel('Concentration (nM)')
        plt.title('ODE simulation')
        plt.show()

    ######################################
    # Check 1: The system doesn't explode
    for obs in model.observables:
        # If any concentration exceeds the threshold, return False and error code 1
        if np.any(yout[obs.name] / model.parameters['omega'].value > 15):
            return {}, False, 1

    ######################################
    # Check 2: Are they oscillating 
    T = t_ode[1] - t_ode[0]  # Sampling interval
    N = yout['obsCycD'][600:].size  # Number of sample points for FFT
    f = np.linspace(0, 1 / T, N)  # Frequency axis for FFT
    fft = {}
    for species in cycling_species:
        fft[species] = np.fft.fft(yout[species][600:] / model.parameters['omega'].value)

    # Find frequency corresponding to max amplitude of FFT
    arg_freq_max = {}
    for species in cycling_species:
        arg_freq_max[species] = np.argmax(np.abs(fft[species][(f > 0.005) * (f < 0.1)]) * 1 / N)

    # Plot the frequency spectrum as well as the arg with the corresponding period
    if return_plot:
        # Create figure with 2 rows and n columns (n = number of species)
        n_species = len(cycling_species)
        fig, axs = plt.subplots(2, n_species, figsize=(5 * n_species, 8))

        # Define a color for each species
        colors = plt.cm.Set2(np.linspace(0, 1, n_species))  # Using Set2 colormap

        for idx, species in enumerate(cycling_species):
            # Top row: time series
            axs[0, idx].plot(t_ode,
                             yout[species] / model.parameters['omega'].value,
                             color=colors[idx])
            axs[0, idx].set_title(species)
            axs[0, idx].set_xlabel('Time (min)')
            if idx == 0:  # Only add y-label for first column
                axs[0, idx].set_ylabel('Concentration (nM)')

            # Bottom row: frequency spectrum
            axs[1, idx].bar(f[(f > 0.005) * (f < 0.1)],
                            np.abs(fft[species][(f > 0.005) * (f < 0.1)]),
                            width=0.001,
                            alpha=0.7,
                            color=colors[idx])

            # Add scatter for max frequency
            axs[1, idx].scatter(f[(f > 0.005) * (f < 0.1)][arg_freq_max[species]],
                                np.abs(fft[species][(f > 0.005) * (f < 0.1)])[arg_freq_max[species]],
                                s=100,
                                color='red',
                                zorder=3)

            axs[1, idx].set_xlabel('frequency')
            if idx == 0:  # Only add y-label for first column
                axs[1, idx].set_ylabel('Amplitude')

            # remove x and y ticks
            axs[1, idx].set_yticks([])
            axs[1, idx].set_xticks([])
            axs[0, idx].set_xticks([])
            axs[0, idx].set_yticks([])

        # Adjust layout to prevent overlap
        plt.tight_layout()
        plt.savefig(save_path + 'time_and_frequency.svg')
        plt.show()


    ######################################
    # Check 3: If frequency of oscillations is in a biological range
    t = np.linspace(0, 1000, 1000)  # Time points for stochastic simulation
    simres_2 = BngSimulator(model, t, verbose=False).run()
    check_parameters = False

    # All frequencies should be between 0.005 and 0.1 Hz
    for species in cycling_species:
        if arg_freq_max[species] > 1 and \
                np.abs(fft[species][(f > 0.005) * (f < 0.1)])[arg_freq_max[species]] > \
                np.abs(fft[species][(f > 0.005) * (f < 0.1)])[arg_freq_max[species] - 1] and \
                np.abs(fft[species][(f > 0.005) * (f < 0.1)])[arg_freq_max[species]] > \
                np.abs(fft[species][(f > 0.005) * (f < 0.1)])[arg_freq_max[species] + 1]:
            check_parameters = True
        else:
            check_parameters = False
            break

    ######################################
    # Check 4: Biological significance of the cycle, cdc20 and Cyclin E should alternate
    # Run stochastic simulation with BngSimulator

    sys_size = model.parameters['omega'].value
    # Find the maximum of the oscillation peaks in Cdc20 and CyeE
    (a, max_tmp_cdc) = peakfinder(simres_2.observables['obsCdc20'], thr[0] * sys_size, thr[0] * sys_size, 1,
                                  includeEndpoints=1)
    (b, max_tmp_cyce) = peakfinder(simres_2.observables['obsCycE'], thr[1] * sys_size, thr[1] * sys_size, 1,
                                   includeEndpoints=1)

    score = 1  # score on how good they alternate, it may not be exactly one because of the stochasticity of the network, and peakfinder mught miss one peak
    variable_range = {}
    if len(max_tmp_cdc) >= 3 and len(max_tmp_cyce) >= 3:
        ranks = ss.rankdata(np.concatenate([max_tmp_cdc, max_tmp_cyce]))
        max_tmp = np.concatenate(
            [['obsCdc20' for _ in range(len(max_tmp_cdc))], ['obsCycE' for _ in range(len(max_tmp_cyce))]])
        Z = [x for _, x in sorted(zip(ranks, max_tmp))]
        score = 0

        for i in range(len(Z) - 1):
            if Z[i] == Z[i + 1]:
                score += 1
        score = score / len(Z)

        if len(a) >= 3 and len(b) >= 3:
            # simres = ScipyOdeSimulator(model, t_ode, verbose=False).run()
            (a_ode, max_tmp_cdc_ode) = peakfinder(yout[model.observables['obsCdc20'].name], thr[0] * sys_size,
                                                  thr[0] * sys_size, 1, includeEndpoints=1)
            (b_ode, max_tmp_cyce_ode) = peakfinder(yout[model.observables['obsCycE'].name], thr[1] * sys_size,
                                                   thr[1] * sys_size, 1, includeEndpoints=1)

            if len(max_tmp_cdc_ode) >= 3:
                for k in range(len(model.observables)):
                    variable_range[model.observables[k].name] = [yout[model.observables[k].name][i] for i in
                                                                 range(max_tmp_cdc_ode[-3], max_tmp_cdc_ode[-2])]
            else:
                for k in range(len(model.observables)):
                    variable_range[model.observables[k].name] = [yout[model.observables[k].name][i] for i in
                                                                 range(900, 999)]

        else:
            for k in range(len(model.observables)):
                variable_range[model.observables[k].name] = [yout[model.observables[k].name][i] for i in
                                                             range(900, 999)]
    # Plot on a first plot the trajectory of Cdc20 and CycE, with scatter plots on peaks of each with different colors. On a lower plots the scatter aligned on a line showing alternance between Cdc20 and CycE, and show the score
    else:
        for k in range(len(model.observables)):
            variable_range[model.observables[k].name] = [yout[model.observables[k].name][i] for i in range(900, 999)]


    if return_plot:
        # Top plot
        plt.figure()
        ax = plt.subplot(211)
        ax.plot(t, simres_2.observables['obsCdc20'], label='Cdc20', c='r')
        ax.plot(t, simres_2.observables['obsCycE'], label='CycE', c='g')
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
        # xlim
        ax.set_xlim([t[0], t[-1]])
        # Remove the top and right axes from the data plot
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.set_title('Score: ' + str(score))
        plt.tight_layout()
        plt.show()

    return variable_range, check_parameters, score


###############################################################################
# Test check_model_consistency

def set_condition_parameters(model, ecm_conditions, condition):
    """
    Set the AKT, ERK, and FAK parameters for a given condition.

    Args:
        model: The PySB model
        ecm_conditions: DataFrame containing the ECM conditions
        condition: String specifying the condition to set
    """
    for input_param in ["AKT", "ERK", "FAK"]:
        model.parameters[input_param].value = ecm_conditions[
            (ecm_conditions["condition"] == condition) &
            (ecm_conditions["species"] == input_param)
            ].value


def check_condition(model, condition, t=1000):
    """
    Run model consistency check for a given condition and print results.

    Args:
        model: The PySB model
        condition: String specifying the condition being checked

    Returns:
        tuple: (variable_range, check, tmp_score) from check_model_consistency
    """
    results = check_model_consistency(model, return_plot=True, t=t)
    print(condition)
    print(results[1])  # Print the 'check' value
    return results


if __name__ == "__main__":
    import os
    from Python.cellular_cycle.cell_cycle_death_model import model
    import pandas as pd

    # Import parameters
    sampled_parameter_names = np.loadtxt(
        os.path.dirname(__file__) + '/../model_data/sampled_parameter_names.txt',
        dtype=str
    )

    # Set model parameters
    new_param = pd.read_csv(os.path.dirname(__file__) + '/../model_data/parameters.csv')
    for key in new_param.keys():
        if model.parameters[key].value != new_param[key][0]:
            model.parameters[key].value = new_param[key]
    model.parameters['omega'].value = 200

    # Set initial concentrations
    initial_concentrations = pd.read_csv(
        os.path.dirname(__file__) + '/../model_data/initial_concentrations.csv'
    )
    for key in initial_concentrations.keys():
        try:
            model.parameters[key].value = initial_concentrations[key][0] / model.parameters['omega'].value
        except:
            print(key)

    # Import ECM conditions
    ECM_conditions = pd.read_csv(os.path.dirname(__file__) + '/../model_data/ECM_conditions_2.csv')

    # Test all conditions
    conditions = [
        "4 kPa erlotinib",
        "4 kPa DMSO",
        "55 kPa erlotinib",
        "55 kPa DMSO",
        "0.6 kPa erlotinib",
        "0.6 kPa DMSO"
    ]

    # Store first condition results for full output
    first_condition = conditions[0]
    set_condition_parameters(model, ECM_conditions, first_condition)
    variable_range, check, tmp_score = check_condition(model, first_condition)
    print(variable_range)
    print(check)
    print(tmp_score)

    # Process remaining conditions
    for condition in conditions[1:]:
        set_condition_parameters(model, ECM_conditions, condition)
        check_condition(model, condition, t=1000)

