# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 18:44:09 2023

@author: Grim.M
"""
from rbf_mmd import mmd_rbf

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors

import pandas as pd
import numpy as np

###############################################################################
path_cost = r'C:\Users\Grim.M\Ecole\Ens\Master\S2\Stage\Code\fit_data\cost_sensitivity_cell_population\cost'
path_cell_count = r'C:\Users\Grim.M\Ecole\Ens\Master\S2\Stage\Code\fit_data\cost_sensitivity_cell_population\cell_count'
path_distribution = r'C:\Users\Grim.M\Ecole\Ens\Master\S2\Stage\Code\fit_data\cost_sensitivity_cell_population\distribution'

cell_count_dict = {condition:{init:[] for init in range(1, 11)} for condition in conditions}
time_dict = {}
for key in ['G1', 'S-M', 'IMT']: time_dict[key] = {condition:{init:[] for init in range(1, 11)} for condition in conditions}
frames = []

condition = '0.6 kPa DMSO'
init_cell_number = 1
i = 0

for condition in conditions[:]:
    for init_cell_number in range(1,11):
        for i in range(10):
                                
            filename = '/parameters_particle_test_initial_cell_number_' + condition.replace(' ', '_') + '_' + str(init_cell_number) + '_initial_cells_iteration_'  + str(i) + '_pso.pkl'
            frac_prolif_parameters_temp = pd.read_pickle(path_cost + "\\" + filename)
            
            frac_prolif_parameters_temp['initial cells number'] = init_cell_number
            frac_prolif_parameters_temp['condition'] = condition
            frames += [frac_prolif_parameters_temp]
            
            
            cell_count_dict[condition][init_cell_number] += [np.load(path_cell_count + '/cell_count_particle_test_initial_cell_number_' + condition.replace(" ", "_" ) + '_' + str(init_cell_number) + '_initial_cells_iteration_' + str(i) + '_' + condition.replace(" ", "_" ) + '.npy')] 
            
        
            time_dict['G1'][condition][init_cell_number] += [np.load(path_distribution + '/G1_distribution_particle_test_initial_cell_number_' +  condition.replace(" ", "_" )  + '_' + str(init_cell_number) + '_initial_cells_iteration_' + str(i) + '_' + condition.replace(" ", "_" ) + '.npy')]
            time_dict['SM'][condition][init_cell_number] += [np.load(path_distribution + '/SM_distribution_particle_test_initial_cell_number_' +  condition.replace(" ", "_" )  + '_' + str(init_cell_number) + '_initial_cells_iteration_' + str(i) + '_' + condition.replace(" ", "_" ) + '.npy')]
            
            try :           
                all_events_tmp = pd.read_pickle(path_distribution + "\Events_particle_test_initial_cell_number_" + condition.replace(" ", "_" ) + '_' + str(init_cell_number) + '_initial_cells_iteration_' + str(i) + '_' + condition.replace(" ", "_" ) + ".pkl")
                time_dict['IMT'][condition][init_cell_number] += [all_events_tmp['Lifespan'][all_events_tmp['event']=='mitosis'].to_list()]
            except:
                time_dict['IMT'][condition][init_cell_number]+=[[]]

frac_prolif_parameters_pso_initial_cell = pd.concat(frames, sort=False)



def import_data_pso(dimensions,iterations, path_pso):
    path_cost = path_pso + '/cost'
    path_cell_count = path_pso + '/cell_count'
    path_distribution = path_pso + '/distribution'
    conditions = ['0.6 kPa DMSO', '0.6 kPa ERL', '4 kPa DMSO', '4 kPa ERL', '55 kPa DMSO', '55 kPa ERL']
    cell_count_dict = {'chain ' + str(chain+1):{'iteration ' + str(iteration+1):{condition:[] for condition in conditions} for iteration in range(iterations)} for chain in range(dimensions)}
    time_dict = {key:{'chain ' + str(chain+1):{'iteration ' + str(iteration+1):{condition:[] for condition in conditions} for iteration in range(iterations)} for chain in range(dimensions)} for key in ['G1', 'S-M', 'IMT']}
    frames = []
    
    time_dict_original = {key:{condition:[] for condition in conditions}for key in ['G1', 'S-M', 'IMT']}
    cell_count_dict_original = {}
    for condition in conditions:
    
        cell_count_dict_original[condition] = np.load(path_cell_count + "/../original_parameters/cell_count_particle_" + 'original_parameters' + "_iteration_" + 'original_parameters' + "_" + condition.replace(' ', '_') + ".npy") 
        for key in ['G1', 'S-M', 'IMT']:
            time_dict_original[key][condition] = np.load(path_distribution + "/../original_parameters/" + key.replace("-", "" ) + "_distribution_particle_original_parameters_iteration_original_parameters_" + condition.replace(' ', '_') + ".npy") 
    
    
    for iteration in range(1, 1 + iterations):
        for chain in range(1, 1 + dimensions):

            filename = '/parameters_particle_' + str(chain) + "_iteration_" + str(iteration) + "_pso.pkl"
            frac_prolif_parameters_temp = pd.read_pickle(path_cost + filename)
        
            frames += [frac_prolif_parameters_temp]
            for condition in conditions:
                try:
                    cell_count_dict['chain ' + str(chain)]['iteration ' +str(iteration)][condition] = np.load(path_cell_count + "/cell_count_particle_" + str(chain) + "_iteration_" + str(iteration) + "_" + condition.replace(' ', '_') + ".npy") 
            
                    for key in ['G1', 'S-M', 'IMT']:
                        time_dict[key]['chain ' + str(chain)]['iteration ' +str(iteration)][condition] = np.load(path_distribution + "/" + key.replace("-", "" ) +"_distribution_particle_"  + str(chain) + "_iteration_" + str(iteration) + "_" + condition.replace(' ', '_') + ".npy") 
                except:
                    pass
    filename = '/../original_parameters/parameters_particle_' + 'original_parameters' + "_iteration_" + 'original_parameters' + "_pso.pkl"
    frac_prolif_parameters_original = pd.read_pickle(path_cost + filename)
    
            

    frac_prolif_parameters_pso = pd.concat(frames)
    frac_prolif_parameters_pso[('cost', 'total cost likelihood')] = frac_prolif_parameters_pso['cost'].drop('- log prior', axis=1).sum(axis=1, skipna=False)    
    frac_prolif_parameters_pso.replace([np.inf, -np.inf], np.nan, inplace=True)
    
    frac_prolif_parameters_original[('cost', 'total cost likelihood')] = frac_prolif_parameters_original['cost'].drop('- log prior', axis=1).sum(axis=1, skipna=False)    
    frac_prolif_parameters_original.replace([np.inf, -np.inf], np.nan, inplace=True)
    
    def plot_dist(frac_prolif_parameters_pso, var):
        values = frac_prolif_parameters_pso[('cost', var)]
        sorted_values = np.sort(values[~np.isnan(values) & np.isfinite(values)])
        x = np.arange(len(sorted_values))
        plt.plot(x, sorted_values)
        plt.fill_between(x, sorted_values, min(sorted_values), alpha=0.5)
        plt.xlabel('Index')
        plt.ylabel(var)
        plt.xlim(0, len(frac_prolif_parameters_pso))
        
        
    plot_dist(frac_prolif_parameters_pso, '- log prior')
    plt.plot([0, iterations*dimensions],2*[frac_prolif_parameters_original[('cost', '- log prior')]], color='red')
    plt.show()
    
    plot_dist(frac_prolif_parameters_pso, 'total cost likelihood')
    plt.plot([0, iterations*dimensions],2*[frac_prolif_parameters_original[('cost', 'total cost likelihood')]], color='red')
    plt.show()
    
    return time_dict_original, cell_count_dict_original, frac_prolif_parameters_original, frac_prolif_parameters_pso, cell_count_dict, time_dict






###############################################################################

def figure_initial_cell_number_cost(different_cost, conditions, frac_prolif_parameters_pso_initial_cell):
    fig, axs = plt.subplots(nrows=len(conditions), ncols=len(different_cost), sharey=True, figsize=(5 * len(different_cost), 5 * len(conditions)))
    
    for i, condition in enumerate(conditions[:]):
        for j, cost_name in enumerate(different_cost):
            ax = axs[i, j]
            
            cost_tmp = frac_prolif_parameters_pso_initial_cell[frac_prolif_parameters_pso_initial_cell['condition'] == condition]
            
            ax.scatter(cost_tmp['initial cells number'], cost_tmp[('cost', cost_name + ' '+ condition)])
            
            result = cost_tmp.groupby(['initial cells number'], as_index=False).agg({('cost', cost_name + ' '+ condition):['mean','std']})
            ax.plot(result['initial cells number'], result[('cost', cost_name + ' '+ condition, 'mean')], c='red', label='mean')
            ax.plot(result['initial cells number'], result[('cost', cost_name + ' '+ condition, 'std')], c='orange', label='std')
            
#            if j == 0: 
#                ax.set_title(condition, fontsize=16)  # Set condition as title only for the first column
            if i == 0:
                ax.set_title(cost_name, fontsize=16)  # Set cost_name as title only for the first row
            
            ax.legend(loc='best')
            ax.tick_params(axis='both', labelsize=10)
            # ax.set_yscale('log')
        
        axs[i, 0].set_ylabel(condition, fontsize=16)
    
    plt.tight_layout()


def figure_initial_cell_number_cell_count(cell_count_dict, conditions, run_times, means, stdevs, cell_numbers = range(1,11)):
    colormap = cm.RdYlBu
    normalize = mcolors.Normalize(vmin=min(cell_numbers), vmax=max(cell_numbers))
    
    fig, axs = plt.subplots(nrows=3, ncols=2, sharey=False, figsize=(10,14))
    
    for ax,run_time, condition, num_condition in zip(axs.flatten(), run_times, conditions, range(len(conditions))):
        for init_cell_number in cell_numbers:
            for i in range(10):  
                
                t_sim = np.arange(0,run_time,1/3)
                
                cell_count_t = np.zeros(len(t_sim))
                time_tmp_m = 0        
                
                for time_tmp_p, cell_count_time_tmp in cell_count_dict[condition][init_cell_number][i]:
                    cell_count_t[(time_tmp_m <= t_sim) & (t_sim < time_tmp_p)] = cell_count_time_tmp
                    time_tmp_m = time_tmp_p
                
                
                DIP_sim = np.log2(cell_count_t/cell_count_t[0])
                
                ax.scatter(t_sim, DIP_sim, color = colormap(normalize(init_cell_number)), s=1)
        ax.plot(t_sim, means[num_condition][: len(t_sim)], color = 'black')
        ax.fill_between(t_sim, means[num_condition][: len(t_sim)] -stdevs[num_condition][: len(t_sim)], means[num_condition][: len(t_sim)] +stdevs[num_condition][: len(t_sim)], color='black', alpha = 0.1)
        ax.set_title(condition)
    
    # Create an axis for the colorbar
    cax = fig.add_axes([0.15, -0.05, 0.7, 0.02])  # [left, bottom, width, height]

    scalarmappaple = cm.ScalarMappable(norm=normalize, cmap=colormap)

    cb = fig.colorbar(scalarmappaple, cax=cax, orientation='horizontal')
    cb.set_label(label='initial cell number', fontsize=14)
    cb.ax.tick_params(labelsize=10)

    plt.tight_layout()
    plt.show()
    
    
def plot_g1_distributions(conditions, time_dict, cell_phase , cell_numbers, all_IMTs, colorcode='Initial Cell Number', sigma=1):
    if colorcode=='Cost':
        max_cost = 0
        costs = {}
        costs[cell_phase] = {}
        for condition in conditions: 
            costs[cell_phase][condition] = {}
            stiffness, _, treatment = condition.split(maxsplit=2)  # split the string into two parts
            stiffness = float(stiffness)  # convert the stiffness from string to float
            treatment = treatment.strip()  # remove any leading or trailing white spaces
            experimental = all_IMTs[(all_IMTs['Drug']==treatment) & (all_IMTs['Stiffness']==stiffness)][cell_phase + '_t']
            for cell_number in cell_numbers:
                costs[cell_phase][condition][cell_number] = []
                time_lists = time_dict[cell_phase][condition][cell_number]
                for time in time_lists:
                    cost_tmp = mmd_rbf(np.array([time]).T, np.array([experimental]).transpose(), sigma)
                    costs[cell_phase][condition][cell_number] += [cost_tmp]
                    max_cost = max(max_cost, cost_tmp)
        normalize = mcolors.Normalize(vmin=0, vmax=max_cost)
    else: 
        normalize = mcolors.Normalize(vmin=min(cell_numbers), vmax=max(cell_numbers))
        
    colormap = cm.RdYlBu
    fig, axs = plt.subplots(nrows=3, ncols=2, sharex=True, sharey=True, figsize=(10,14))
    
    g1_min = 0
    g1_max = 100
    bins = np.linspace(g1_min, g1_max, int((g1_max - g1_min)/2))
    
    for ax, condition in zip(axs.flatten(), conditions):
        for cell_number in cell_numbers:
            g1_lists = time_dict[cell_phase][condition][cell_number]
            for j, g1 in enumerate(g1_lists):
                if colorcode=='Cost': 
                    color_tmp = colormap(normalize(costs[cell_phase][condition][cell_number][j]))
                elif colorcode=='Initial Cell Number': 
                    color_tmp = colormap(normalize(cell_number))
    
                ax.hist(g1, bins=bins, alpha=0.4, histtype='step', color=color_tmp, label=f'Cell Number: {cell_number} - ' + cell_phase, density =True)

        ax.set_title(condition)
        
        stiffness, _, treatment = condition.split(maxsplit=2)  # split the string into two parts
        stiffness = float(stiffness)  # convert the stiffness from string to float
        treatment = treatment.strip()  # remove any leading or trailing white spaces
            
        experimental = all_IMTs[(all_IMTs['Drug']==treatment) & (all_IMTs['Stiffness']==stiffness)][cell_phase + '_t']
        ax.hist(experimental, bins=bins, alpha=0.6, histtype='step', color='black', density =True, linewidth=2)
        
        
    cax = fig.add_axes([1.05, 0.1, 0.02, 0.8])

    scalarmappaple = cm.ScalarMappable(norm=normalize, cmap=colormap)

    cb = fig.colorbar(scalarmappaple, cax=cax)
    cb.set_label(label=colorcode, fontsize=14)
    cb.ax.tick_params(labelsize=10)

    fig.suptitle(cell_phase + ' Time Distribution', fontsize=16)
    plt.tight_layout()
    plt.show()


def plot_distributions_cost(conditions, time_dict, cell_numbers, frac_prolif_parameters_pso_initial_cell, all_IMTs, sigma = 1, colorcode = 'Initial Cell Number'):
    if colorcode=='Cost':
        max_cost = 0
        costs = {}
        for cell_phase in ['G1', 'S-M', 'IMT']:
            costs[cell_phase] = {}
            for condition in conditions: 
                costs[cell_phase][condition] = {}
                stiffness, _, treatment = condition.split(maxsplit=2)  # split the string into two parts
                stiffness = float(stiffness)  # convert the stiffness from string to float
                treatment = treatment.strip()  # remove any leading or trailing white spaces
                experimental = all_IMTs[(all_IMTs['Drug']==treatment) & (all_IMTs['Stiffness']==stiffness)][cell_phase + '_t']
                for cell_number in cell_numbers:
                    costs[cell_phase][condition][cell_number] = []
                    time_lists = time_dict[cell_phase][condition][cell_number]
                    for time in time_lists:
                        cost_tmp = mmd_rbf(np.array([time]).T, np.array([experimental]).transpose(), sigma)
                        costs[cell_phase][condition][cell_number] += [cost_tmp]
                        max_cost = max(max_cost, cost_tmp)
            normalize = mcolors.Normalize(vmin=0, vmax=max_cost)
    else: 
        normalize = mcolors.Normalize(vmin=min(cell_numbers), vmax=max(cell_numbers))
    
    g1_min = sm_min = imt_min = 0
    g1_max = 80
    sm_max = 40
    imt_max = 100
    bins = {}

    bins['G1'] = np.linspace(g1_min, g1_max, int((g1_max - g1_min)/2))
    bins['S-M'] = np.linspace(sm_min, sm_max, int((sm_max - sm_min)/2))
    bins['IMT'] = np.linspace(imt_min, imt_max, int((imt_max - imt_min)/2))
    
    colormap = cm.RdYlBu
    fig, axs = plt.subplots(nrows=6, ncols=3, figsize=(12, 20), sharex=True, sharey=True)
    
    for i, condition in enumerate(conditions):
        for j, cell_phase in enumerate(['G1', 'S-M', 'IMT']):
            ax = axs[i, j]
    
            ax.set_title(cell_phase + f' - {condition}')
            ax.set_xlabel('Time (hours)')
            ax.set_ylabel('Density')
            
            stiffness, _, treatment = condition.split(maxsplit=2)
            stiffness = float(stiffness)
            treatment = treatment.strip()
                
            experimental = all_IMTs[(all_IMTs['Drug']==treatment) & (all_IMTs['Stiffness']==stiffness)][cell_phase + '_t']
            
            for cell_number in cell_numbers:
                time_lists = time_dict[cell_phase][condition][cell_number]
                
                for j, time in enumerate(time_lists):
                    if colorcode=='Cost': 
                        color_tmp = colormap(normalize(costs[cell_phase][condition][cell_number][j]))
                    elif colorcode=='Initial Cell Number': 
                        color_tmp = colormap(normalize(cell_number))
        
                    ax.hist(time, bins=bins[cell_phase], alpha=0.4, histtype='step', color=color_tmp, label=f'Cell Number: {cell_number}', density=True)
                    ax.set_title(cell_phase + f' - {condition}')
                
            ax.hist(experimental, bins=bins[cell_phase], alpha=0.6, histtype='step', color='black', density=True, linewidth=2)
            
    cax = fig.add_axes([1.02, 0.1, 0.02, 0.8])
    scalarmappaple = cm.ScalarMappable(norm=normalize, cmap=colormap)
    cb = fig.colorbar(scalarmappaple, cax=cax)
    cb.set_label(label=colorcode, fontsize=14)
    cb.ax.tick_params(labelsize=10)
    
    fig.text(0.15, 1.01, 'G1 Time Distributions', fontsize=16, ha='center')
    fig.text(0.5, 1.01, 'S-M Time Distributions', fontsize=16, ha='center')
    fig.text(0.85, 1.01, 'IMT Time Distributions', fontsize=16, ha='center')
    
    plt.tight_layout()
    plt.show()





def figure_pso_cell_count(cell_count_dict_original, frac_prolif_parameters_original, frac_prolif_parameters_pso, cell_count_dict, conditions, run_times, means, stdevs):
    costs_array = np.array([frac_prolif_parameters_pso[('cost', 'cost population dynamics ' + condition)].to_numpy() for condition in conditions])

    vmax = np.nanmax(costs_array)
    vmin = np.nanmin(costs_array)
    normalize = mcolors.Normalize(vmin=vmin, vmax=vmax)
    colormap = cm.RdYlBu   
    fig, axs = plt.subplots(nrows=3, ncols=2, sharey=False, figsize=(10,14))
    for ax,run_time, condition, num_condition in zip(axs.flatten(), run_times, conditions, range(len(conditions))):
        t_sim = np.arange(0,run_time,1/3)
        for chain, iteration, cost_tmp in zip(frac_prolif_parameters_pso['pso']['chain'], frac_prolif_parameters_pso['pso']['iteration'], frac_prolif_parameters_pso[('cost', 'cost population dynamics ' + condition)]): 
            cell_count_t = np.zeros(len(t_sim))
            time_tmp_m = 0        
            
            for time_tmp_p, cell_count_time_tmp in cell_count_dict_pso['chain '+str(chain)]['iteration ' + str(iteration)][condition]:
                cell_count_t[(time_tmp_m <= t_sim) & (t_sim < time_tmp_p)] = cell_count_time_tmp
                time_tmp_m = time_tmp_p
            
            DIP_sim = np.log2(cell_count_t/cell_count_t[0])
            ax.scatter(t_sim, DIP_sim, color = colormap(normalize(cost_tmp)), s=1)
            
        cell_count_t = np.zeros(len(t_sim))
        time_tmp_m = 0 
        for time_tmp_p, cell_count_time_tmp in cell_count_dict_original[condition]:
            cell_count_t[(time_tmp_m <= t_sim) & (t_sim < time_tmp_p)] = cell_count_time_tmp
            time_tmp_m = time_tmp_p
        
        DIP_sim = np.log2(cell_count_t/cell_count_t[0])
         
        ax.plot(t_sim, DIP_sim, color = colormap(normalize(frac_prolif_parameters_original[('cost', 'cost population dynamics ' + condition)])), linewidth = 5)
        ax.set_title(condition)
        
        ax.plot(t_sim, means[num_condition][: len(t_sim)], color = 'black')
        ax.fill_between(t_sim, means[num_condition][: len(t_sim)] -stdevs[num_condition][: len(t_sim)], means[num_condition][: len(t_sim)] +stdevs[num_condition][: len(t_sim)], color='black', alpha = 0.1)
        ax.set_title(condition)
    
    # Create an axis for the colorbar
    cax = fig.add_axes([0.15, -0.05, 0.7, 0.02])  # [left, bottom, width, height]

    scalarmappaple = cm.ScalarMappable(norm=normalize, cmap=colormap)

    cb = fig.colorbar(scalarmappaple, cax=cax, orientation='horizontal')
    cb.set_label(label='cost population dynamics', fontsize=14)
    cb.ax.tick_params(labelsize=10)

    plt.tight_layout()
    plt.show()
    


def plot_distributions_cost_pso(conditions, time_dict_pso, frac_prolif_parameters_pso, all_IMTs, frac_prolif_parameters_original, time_dict_original, times=['G1', 'S-M']):
    costs_array = np.array([frac_prolif_parameters_pso[('cost', 'cost ' + cell_phase + ' distributions ' + condition)].to_numpy() for condition in conditions for cell_phase in times])
    
    vmax = np.nanmax(costs_array)
    vmin = min(np.nanmin(costs_array), 0)
    normalize = mcolors.Normalize(vmin=vmin, vmax=vmax)

    g1_min = sm_min = imt_min = 0
    g1_max = 80
    sm_max = 40
    imt_max = 100
    bins = {}

    bins['G1'] = np.linspace(g1_min, g1_max, int((g1_max - g1_min)/2))
    bins['SM'] = np.linspace(sm_min, sm_max, int((sm_max - sm_min)/2))
    bins['IMT'] = np.linspace(imt_min, imt_max, int((imt_max - imt_min)/2))
    
    colormap = cm.RdYlBu
    fig, axs = plt.subplots(nrows=len(conditions), ncols=len(times), figsize=(len(times)*4, 20), sharex=True, sharey=True)
    
    for i, condition in enumerate(conditions):
        for j, cell_phase in enumerate(times):
            ax = axs[i, j]
    
            ax.set_title(cell_phase + f' - {condition}')
            ax.set_xlabel('Time (hours)')
            ax.set_ylabel('Density')
            
            stiffness, _, treatment = condition.split(maxsplit=2)
            stiffness = float(stiffness)
            treatment = treatment.strip()
                
            experimental = all_IMTs[(all_IMTs['Drug']==treatment) & (all_IMTs['Stiffness']==stiffness)][cell_phase + '_t']
            
            for chain, iteration, cost_tmp in zip(frac_prolif_parameters_pso['pso']['chain'], frac_prolif_parameters_pso['pso']['iteration'], frac_prolif_parameters_pso[('cost', 'cost ' + cell_phase + ' distributions ' + condition)]): 
                
                time = time_dict_pso[cell_phase]['chain '+str(chain)]['iteration ' + str(iteration)][condition]
                ax.hist(time, bins=bins[cell_phase], alpha=0.4, histtype='step', color=colormap(normalize(cost_tmp)), density=True)
                ax.set_title(cell_phase + f' - {condition}')
                
            ax.hist(time_dict_original[cell_phase][condition], bins=bins[cell_phase], alpha=0.6, histtype='step', color=colormap(normalize(frac_prolif_parameters_original[('cost', 'cost ' + cell_phase + ' distributions ' + condition)])), density=True, linewidth=2)
            ax.hist(experimental, bins=bins[cell_phase], alpha=0.6, histtype='step', color='black', density=True, linewidth=2)
            
    cax = fig.add_axes([1.02, 0.1, 0.02, 0.8])
    scalarmappaple = cm.ScalarMappable(norm=normalize, cmap=colormap)
    cb = fig.colorbar(scalarmappaple, cax=cax)
    cb.set_label(label='Cost mmd', fontsize=14)
    cb.ax.tick_params(labelsize=10)
    
    fig.text(0.15, 1.01, 'G1 Time Distributions', fontsize=16, ha='center')
    fig.text(0.5, 1.01, 'S-M Time Distributions', fontsize=16, ha='center')
    #fig.text(0.85, 1.01, 'IMT Time Distributions', fontsize=16, ha='center')
    
    plt.tight_layout()
    plt.show()

########################################################################

run_times = [120, 120, 120, 120, 60, 60] 
cell_numbers = range(1,11)


figure_initial_cell_number_cost([ 'cost population dynamics', 'cost S-M distributions', 'cost G1 distributions'], conditions,  frac_prolif_parameters_pso_initial_cell)

figure_initial_cell_number_cell_count(cell_count_dict, conditions, run_times, means, stdevs, cell_numbers)

cell_phase = 'G1'

plot_g1_distributions(conditions, time_dict, 'G1' , cell_numbers, all_IMTs, colorcode='Initial Cell Number', sigma=1)
plot_g1_distributions(conditions, time_dict, 'S-M' , cell_numbers, all_IMTs, colorcode='Initial Cell Number', sigma=1)
plot_g1_distributions(conditions, time_dict, 'IMT' , cell_numbers, all_IMTs, colorcode='Initial Cell Number', sigma=1)

plot_g1_distributions(conditions, time_dict, 'G1' , cell_numbers, all_IMTs, colorcode='Cost', sigma=1)
plot_g1_distributions(conditions, time_dict, 'S-M' , cell_numbers, all_IMTs, colorcode='Cost', sigma=1)
plot_g1_distributions(conditions, time_dict, 'IMT' , cell_numbers, all_IMTs, colorcode='Cost', sigma=1)


 
plot_distributions_cost(conditions, time_dict, cell_numbers, frac_prolif_parameters_pso_initial_cell, all_IMTs, sigma = 4, colorcode = 'Initial Cell Number')
plot_distributions_cost(conditions, time_dict, cell_numbers, frac_prolif_parameters_pso_initial_cell, all_IMTs, sigma = 4, colorcode = 'Cost')





dimensions = 35
iterations = 80

conditions = ['0.6 kPa DMSO', '0.6 kPa ERL', '4 kPa DMSO', '4 kPa ERL', '55 kPa DMSO', '55 kPa ERL']
run_times = [120, 120, 120, 120, 60, 60] 

dimensions = 35
iterations = 8

path_pso = r'C:\Users\Grim.M\Ecole\Ens\Master\S2\Stage\Code\fit_data\pso\model6'

time_dict_original, cell_count_dict_original, frac_prolif_parameters_original, frac_prolif_parameters_pso, cell_count_dict_pso, time_dict_pso = import_data_pso(dimensions, iterations, path_pso)

frac_prolif_parameters_pso.to_csv(r'C:\Users\Grim.M\Documents\git_repos\ECM_cell_cycle_model\Python\fit_data_parameters_pso_simulations_results.csv')

frac_prolif_parameters_pso[('cost', 'cost population dynamics')] = np.array(sum([frac_prolif_parameters_pso[('cost', 'cost population dynamics ' + condition)] for condition in conditions]))
 
frac_prolif_parameters_pso['cost'].columns

frac_prolif_parameters_pso = frac_prolif_parameters_pso.sort_values(by=[('cost', 'cost population dynamics')])
frac_prolif_parameters_pso = frac_prolif_parameters_pso.sort_values(by=[('cost', 'total cost likelihood')])

figure_pso_cell_count(cell_count_dict_original, frac_prolif_parameters_original, frac_prolif_parameters_pso.iloc[0:40], cell_count_dict_pso, conditions, run_times, means, stdevs)



times = ['G1', 'S-M']
frac_prolif_parameters_pso[('cost', 'cost distribution')] = np.array(sum([frac_prolif_parameters_pso[('cost', 'cost ' + cell_phase + ' distributions ' + condition)].to_numpy() for condition in conditions for cell_phase in times]))

frac_prolif_parameters_pso = frac_prolif_parameters_pso.sort_values(by=('cost', 'cost distribution'))
frac_prolif_parameters_pso = frac_prolif_parameters_pso.sort_values(by=[('cost', 'total cost likelihood')])


#for condition in conditions: frac_prolif_parameters_pso[('cost', 'cost IMT distributions '+condition)]=0

plot_distributions_cost_pso(conditions, time_dict_pso, frac_prolif_parameters_pso.iloc[0:40], all_IMTs, frac_prolif_parameters_original, time_dict_original)

