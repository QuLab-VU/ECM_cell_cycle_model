# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 18:44:09 2023

@author: Grim.M
"""

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors

import pandas as pd
import numpy as np
path_cost = r'C:\Users\Grim.M\Ecole\Ens\Master\S2\Stage\Code\fit_data\cost_sensitivity_cell_population\cost'
path_cell_count = r'C:\Users\Grim.M\Ecole\Ens\Master\S2\Stage\Code\fit_data\cost_sensitivity_cell_population\cell_count'
path_distribution = r'C:\Users\Grim.M\Ecole\Ens\Master\S2\Stage\Code\fit_data\cost_sensitivity_cell_population\distribution'

conditions = ['0.6 kPa DMSO', '0.6 kPa ERL', '4 kPa DMSO', '4 kPa ERL', '55 kPa DMSO', '55 kPa ERL']
cell_count_dict = {condition:{init:[] for init in range(1, 11)} for condition in conditions}
SM_dict = {condition:{init:[] for init in range(1, 11)} for condition in conditions}
G1_dict = {condition:{init:[] for init in range(1, 11)} for condition in conditions}
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
        
            G1_dict[condition][init_cell_number] += [np.load(path_distribution + '/G1_distribution_particle_test_initial_cell_number_' +  condition.replace(" ", "_" )  + '_' + str(init_cell_number) + '_initial_cells_iteration_' + str(i) + '_' + condition.replace(" ", "_" ) + '.npy')]
            SM_dict[condition][init_cell_number] += [np.load(path_distribution + '/SM_distribution_particle_test_initial_cell_number_' +  condition.replace(" ", "_" )  + '_' + str(init_cell_number) + '_initial_cells_iteration_' + str(i) + '_' + condition.replace(" ", "_" ) + '.npy')]


frac_prolif_parameters_pso_initial_cell = pd.concat(frames, sort=False)




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
    
    
def plot_g1_distributions(conditions, G1_dict, cell_phase , cell_numbers):
    colormap = cm.RdYlBu
    normalize = mcolors.Normalize(vmin=min(cell_numbers), vmax=max(cell_numbers))

    fig, axs = plt.subplots(nrows=3, ncols=2, sharex=True, sharey=True, figsize=(10,14))
    
    g1_min = 0
    g1_max = 100
    bins = np.linspace(g1_min, g1_max, 60)
    
    for ax, condition in zip(axs.flatten(), conditions):
        for cell_number in cell_numbers:
            g1_lists = G1_dict[condition][cell_number]
            for g1 in g1_lists:
                ax.hist(g1, bins=bins, alpha=0.4, histtype='step', color=colormap(normalize(cell_number)), label=f'Cell Number: {cell_number} - ' + cell_phase, density =True)

        ax.set_title(condition)
        
        stiffness, _, treatment = condition.split(maxsplit=2)  # split the string into two parts
        stiffness = float(stiffness)  # convert the stiffness from string to float
        treatment = treatment.strip()  # remove any leading or trailing white spaces
            
        experimental = all_IMTs[(all_IMTs['Drug']==treatment) & (all_IMTs['Stiffness']==stiffness)][cell_phase + '_t']
        ax.hist(g1, bins=bins, alpha=0.6, histtype='step', color='black', density =True, linewidth=3)
        
        
    cax = fig.add_axes([1.05, 0.1, 0.02, 0.8])

    scalarmappaple = cm.ScalarMappable(norm=normalize, cmap=colormap)

    cb = fig.colorbar(scalarmappaple, cax=cax)
    cb.set_label(label='Initial Cell Number', fontsize=14)
    cb.ax.tick_params(labelsize=10)

    fig.suptitle(cell_phase + ' Time Distribution', fontsize=16)
    plt.tight_layout()
    plt.show()



########################################################################

run_times = [120, 120, 120, 120, 60, 60] 
cell_numbers = range(1,11)


figure_initial_cell_number_cost([ 'cost population dynamics', 'cost S-M distributions', 'cost G1 distributions'], conditions,  frac_prolif_parameters_pso_initial_cell)

figure_initial_cell_number_cell_count(cell_count_dict, conditions, run_times, means, stdevs, cell_numbers)


plot_g1_distributions(conditions, G1_dict, 'G1' , cell_numbers)
plot_g1_distributions(conditions, SM_dict,  'S-M' , cell_numbers)

