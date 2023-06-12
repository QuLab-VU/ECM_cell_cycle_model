# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 18:44:09 2023

@author: Grim.M
"""


path_cost = r'C:\Users\Grim.M\Ecole\Ens\Master\S2\Stage\Code\fit_data\cost_sensitivity_cell_population\cost'

cell_count_dict = {condition:{init:[] for init in range(10)} for condition in conditions}
frames = []

condition = '0.6 kPa DMSO'
init_cell_number = 0
i = 0

for condition in conditions:
    for init_cell_number in range(10):
        for i in range(1):
            i=0
                    
            filename = '/parameters_particle_test_initial_cell_number_' + condition.replace(' ', '_') + '_' + str(init_cell_number) + '_initial_cells_iteration_'  + str(i) + '_pso.pkl'
            frac_prolif_parameters_temp = pd.read_pickle(path_cost + "\\" + filename)
            frames += [frac_prolif_parameters_temp]

frac_prolif_parameters_pso_initial_cell = pd.concat(frames, sort=False)



fig, axs = plt.subplots(nrows=3, ncols=2, sharey=True, figsize=(10,14))
    
for ax, condition in zip(axs.flatten(), conditions):
    initial_cell_nber = frac_prolif_parameters_pso_initial_cell['initial cell numbers '+condition]
    cost = frac_prolif_parameters_pso_initial_cell['Cost '+condition]
    ax.scatter(initial_cell_nber, cost)
    ax.set_title(condition)

    df = pd.concat([initial_cell_nber,cost], axis=1)
    result = df.groupby(['initial cell numbers ' + condition], as_index=False).agg({'Cost '+condition:['mean','std']})
    ax.plot(result['initial cell numbers '+condition], result['Cost '+condition,'mean'], c='red')
    ax.plot(result['initial cell numbers '+condition], result['Cost '+condition,'std'], c='orange')
