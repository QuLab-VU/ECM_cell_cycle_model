import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
path = os.getcwd()
from pso_figures import figure_pso_cell_count, plot_distributions_cost_pso, plot_dist
import re
from pathlib import Path
import glob


###############################################################################
# Importing data

def import_experimental_data(path):
    DATADIR = path + "/../../Data/Exp_Cell_Counts"
    filepaths = os.listdir(DATADIR)
    filenames = [f for f in filepaths if re.match(r'^DIP_DS8', f)]
    filepaths = [Path(DATADIR, f) for f in filenames]
    filepaths.sort()
    filepaths

    all_response = []
    for f in filepaths:
        all_response.append(np.loadtxt(f, delimiter=','))

    means = []
    for vals in range(len(all_response)):
        means.append(np.mean(all_response[vals], axis=1))

    stdevs = []
    for vals in range(len(all_response)):
        stdevs.append(np.std(all_response[vals], axis=1))

    # t = np.arange(0, 420 / 3, 1 / 3)

    IMT_DATADIR = Path(os.path.dirname(__file__) + "/../../Data/Exp_Intermitotic_Times")
    # IMT_DATADIR = Path(path + "/Data/Exp_Intermitotic_Times")

    print(IMT_DATADIR)
    file_location = os.path.join(IMT_DATADIR, '*imt_*.csv')
    filenames = glob.glob(file_location)
    filenames.sort()
    # print(filenames)


    IMTs = []
    for f in filenames:
        tmp_df = pd.read_csv(f)
        tmp_df = tmp_df.replace(to_replace=-1, value=np.nan)
        IMTs.append(tmp_df)

    ii = 0
    for df in IMTs:
        df["IMT"] = df["Division"] - df["Mitotic Frame"]
        df["G1"] = df["Green On"] - df["Mitotic Frame"]
        df["S-M"] = df["Division"] - df["Green On"]
        if ii == 2 or ii == 3 or ii == 6 or ii == 7:
            df["IMT_t"] = df["IMT"] * (32 / 60)
            df["G1_t"] = df["G1"] * (32 / 60)
            df["S-M_t"] = df["S-M"] * (32 / 60)
        else:
            df["IMT_t"] = df["IMT"] * (20 / 60)
            df["G1_t"] = df["G1"] * (20 / 60)
            df["S-M_t"] = df["S-M"] * (20 / 60)
        # df = df[(df["IMT_t"] ==27)].reset_index(drop=True)
        ii += 1

    stiffs = [0.2, 0.2, 0.6, 0.6, 4, 4, 55, 55]
    ii = 0
    for df in IMTs:
        tmp_stiff = np.ones(len(df)) * stiffs[ii]
        df.insert(0, "Stiffness", tmp_stiff)
        if ii % 2 == 0:
            df["Drug"] = "DMSO"
        else:
            df["Drug"] = "erlotinib"
        ii += 1

    all_IMTs = pd.concat(IMTs, ignore_index=True, sort=False)
    centers = np.arange(0.5, 100.5, 1)  # histogram bin centers

    all_IMTs = all_IMTs[all_IMTs['IMT_t'].notna()]

    all_IMTs = all_IMTs[all_IMTs['S-M_t'] > 0]
    all_IMTs = all_IMTs[all_IMTs['G1_t'] > 0]

    return means, stdevs, all_IMTs

##############################################################################


def import_data_pso(dimensions, iterations, conditions, path_pso):
    path_cost = path_pso + 'cost'
    path_cell_count = path_pso + 'cell_count'
    path_distribution = path_pso + 'distribution'
    cell_count_dict = {
        'chain ' + str(chain + 1): {'iteration ' + str(iteration + 1): {condition: [] for condition in conditions} for
                                    iteration in range(iterations)} for chain in range(dimensions)}
    time_dict = {key: {
        'chain ' + str(chain + 1): {'iteration ' + str(iteration + 1): {condition: [] for condition in conditions} for
                                    iteration in range(iterations)} for chain in range(dimensions)} for key in
                 ['G1', 'S-M', 'IMT']}
    frames = []

    ### Original parameters
    time_dict_original = {key: {condition: [] for condition in conditions} for key in ['G1', 'S-M', 'IMT']}
    cell_count_dict_original = {}
    for condition in conditions:
        cell_count_dict_original[condition] = np.load(
            path_cell_count + "/../original_parameters/cell_count_particle_" + 'original_parameters' + "_iteration_" + 'original_parameters' + "_" + condition.replace(
                ' ', '_') + ".npy")
        for key in ['G1', 'S-M', 'IMT']:
            time_dict_original[key][condition] = np.load(
                path_distribution + "/../original_parameters/" + key.replace("-",
                                                                             "") + "_distribution_particle_original_parameters_iteration_original_parameters_" + condition.replace(
                    ' ', '_') + ".npy")

    for iteration in range(1, 1 + iterations):
        for chain in range(1, 1 + dimensions):


            try:
                filename = '/parameters_particle_' + str(chain) + "_iteration_" + str(iteration) + "_pso.pkl"
                frac_prolif_parameters_temp = pd.read_pickle(path_cost + filename)
                frames += [frac_prolif_parameters_temp]
            except:
                pass
            for condition in conditions:
                try:
                    cell_count_dict['chain ' + str(chain)]['iteration ' + str(iteration)][condition] = np.load(
                        path_cell_count + "/cell_count_particle_" + str(chain) + "_iteration_" + str(
                            iteration) + "_" + condition.replace(' ', '_') + ".npy")

                    for key in ['G1', 'S-M', 'IMT']:
                        time_dict[key]['chain ' + str(chain)]['iteration ' + str(iteration)][condition] = np.load(
                            path_distribution + "/" + key.replace("-", "") + "_distribution_particle_" + str(
                                chain) + "_iteration_" + str(iteration) + "_" + condition.replace(' ', '_') + ".npy")
                except:
                    pass
    filename = '/../original_parameters/parameters_particle_' + 'original_parameters' + "_iteration_" + 'original_parameters' + "_pso.pkl"
    frac_prolif_parameters_original = pd.read_pickle(path_cost + filename)

    frac_prolif_parameters_pso = pd.concat(frames)
    frac_prolif_parameters_pso[('cost', 'total cost likelihood')] = frac_prolif_parameters_pso['cost'].drop(
        '- log prior', axis=1).sum(axis=1, skipna=False)
    frac_prolif_parameters_pso.replace([np.inf, -np.inf], np.nan, inplace=True)

    frac_prolif_parameters_original[('cost', 'total cost likelihood')] = frac_prolif_parameters_original['cost'].drop(
        '- log prior', axis=1).sum(axis=1, skipna=False)
    frac_prolif_parameters_original.replace([np.inf, -np.inf], np.nan, inplace=True)

    # def plot_dist(frac_prolif_parameters_pso, var):
    #     values = frac_prolif_parameters_pso[('cost', var)]
    #     sorted_values = np.sort(values[~np.isnan(values) & np.isfinite(values)])
    #     x = np.arange(len(sorted_values))
    #     plt.plot(x, sorted_values)
    #     plt.fill_between(x, sorted_values, min(sorted_values), alpha=0.5)
    #     plt.xlabel('Index')
    #     plt.ylabel(var)
    #     plt.xlim(0, len(frac_prolif_parameters_pso))

    # plot_dist(frac_prolif_parameters_pso, 'total cost likelihood')
    # plt.plot([0, iterations * dimensions],
    #          2 * [frac_prolif_parameters_original['cost']['total cost likelihood'][0]])
    # plt.show()
    # plot_dist(frac_prolif_parameters_pso, '- log prior')
    # plt.plot([0, iterations * dimensions], 2 * [frac_prolif_parameters_original['cost']['- log prior'][0]],
    #          color='red')
    # plt.show()

    return time_dict_original, cell_count_dict_original, frac_prolif_parameters_original, frac_prolif_parameters_pso, cell_count_dict, time_dict



if __name__ == '__main__':
    path_pso = path + "/../fit_data/raw_data/pso/"
    iterations = 62
    dimensions = 33
    conditions = ['0.6 kPa DMSO', '0.6 kPa erlotinib', '4 kPa DMSO', '4 kPa erlotinib', '55 kPa DMSO', '55 kPa erlotinib']

    time_dict_original, cell_count_dict_original, frac_prolif_parameters_original, frac_prolif_parameters_pso, cell_count_dict_pso, time_dict_pso = import_data_pso(dimensions, iterations, conditions, path_pso)

    plot_dist(frac_prolif_parameters_pso, '- log prior')
    plt.plot([0, iterations * dimensions], 2 * [frac_prolif_parameters_original[('cost', '- log prior')].values[0,0]], color='red')
    # save in /Users/jeremie/my_papers/ECM_cell_cycle/figures in svg
    # plt.savefig(path + "/../fit_data/pso/log_prior_pso.svg")
    plt.show()

    plot_dist(frac_prolif_parameters_pso, 'total cost likelihood')
    plt.plot([0, iterations * dimensions], 2 * [frac_prolif_parameters_original[('cost', 'total cost likelihood')].values[0,0]],
             color='red')
    plt.savefig("/Users/jeremie/my_papers/ECM_cell_cycle/figures/total_cost_likelihood_pso.svg")
    plt.show()

    # value of frac_prolif_parameters_original[('cost', 'total cost likelihood')]

    # save frac_prolif_parameters_pso in fit_data/pso
    frac_prolif_parameters_pso.to_pickle(path + "/../fit_data/pso/frac_prolif_parameters_pso.pkl")
    # Save original parameters in fit_data/pso
    frac_prolif_parameters_original.to_pickle(path + "/../fit_data/pso/frac_prolif_parameters_original.pkl")

    means, stdevs, all_IMTs = import_experimental_data(path)
    run_times = [120, 120, 120, 120, 60, 60]

    frac_prolif_parameters_pso[('cost', 'cost custom')] = frac_prolif_parameters_pso[('cost', 'cost population dynamics 0.6 kPa DMSO')] +frac_prolif_parameters_pso[('cost', 'cost population dynamics 55 kPa erlotinib')]
    # Sort frac_prolif_parameters_pso by best cost
    frac_prolif_parameters_pso = frac_prolif_parameters_pso.sort_values(('cost', 'total cost likelihood'), ascending=True)
    frac_prolif_parameters_pso = frac_prolif_parameters_pso.sort_values(('cost', 'cost population dynamics 0.6 kPa DMSO'), ascending=True)
    frac_prolif_parameters_pso = frac_prolif_parameters_pso.sort_values(('cost', 'cost custom'), ascending=True)
    frac_prolif_parameters_pso = frac_prolif_parameters_pso.sort_values(('cost', 'total cost likelihood'), ascending=True)

    figure_pso_cell_count(cell_count_dict_original, frac_prolif_parameters_original, frac_prolif_parameters_pso.iloc[0:40], cell_count_dict_pso, conditions, run_times, means, stdevs)
    times = ['G1', 'S-M']
    plot_distributions_cost_pso(conditions, time_dict_pso, frac_prolif_parameters_pso.iloc[0:40], all_IMTs, frac_prolif_parameters_original, time_dict_original, times)

    frac_prolif_parameters_pso_parameters = frac_prolif_parameters_pso['model parameters']

    # umap

    frac_prolif_parameters_pso['cost'].columns


    frac_prolif_parameters_pso[('cost', 'cost 0.6 kPa DMSO')] =frac_prolif_parameters_pso['cost']['cost G1 distributions 0.6 kPa DMSO'] + frac_prolif_parameters_pso['cost']['cost S-M distributions 0.6 kPa DMSO'] + frac_prolif_parameters_pso['cost']['cost population dynamics 0.6 kPa DMSO']

    frac_prolif_parameters_pso[('cost', 'cost 0.6 kPa erlotinib')] =frac_prolif_parameters_pso['cost']['cost G1 distributions 0.6 kPa erlotinib'] + frac_prolif_parameters_pso['cost']['cost S-M distributions 0.6 kPa erlotinib'] + frac_prolif_parameters_pso['cost']['cost population dynamics 0.6 kPa erlotinib']

    frac_prolif_parameters_pso[('cost', 'cost 0.6 kPa')] = frac_prolif_parameters_pso['cost']['cost 0.6 kPa DMSO'] + frac_prolif_parameters_pso['cost']['cost 0.6 kPa erlotinib']

    frac_prolif_parameters_pso[('cost', 'cost 55 kPa DMSO')] =frac_prolif_parameters_pso['cost']['cost G1 distributions 55 kPa DMSO'] + frac_prolif_parameters_pso['cost']['cost S-M distributions 55 kPa DMSO'] + frac_prolif_parameters_pso['cost']['cost population dynamics 55 kPa DMSO']

    frac_prolif_parameters_pso[('cost', 'cost 55 kPa erlotinib')] =frac_prolif_parameters_pso['cost']['cost G1 distributions 55 kPa erlotinib'] + frac_prolif_parameters_pso['cost']['cost S-M distributions 55 kPa erlotinib'] + frac_prolif_parameters_pso['cost']['cost population dynamics 55 kPa erlotinib']

    frac_prolif_parameters_pso[('cost', 'cost 55 kPa')] = frac_prolif_parameters_pso['cost']['cost 55 kPa DMSO'] + frac_prolif_parameters_pso['cost']['cost 55 kPa erlotinib']

    plt.plot(frac_prolif_parameters_pso['cost', 'cost 0.6 kPa'], frac_prolif_parameters_pso['cost']['cost 55 kPa'], 'o')
    plt.show()

    plt.plot(frac_prolif_parameters_pso['cost', 'cost population dynamics 0.6 kPa DMSO'], frac_prolif_parameters_pso['cost']['cost population dynamics 55 kPa DMSO'], 'o')
    plt.xscale('log')
    plt.yscale('log')
    plt.show()


    import umap

    X = frac_prolif_parameters_pso['model parameters'].to_numpy()

    n_components = 2
    umap_parameters = umap.UMAP(n_components=n_components)
    umap_embedding = umap_parameters.fit_transform(X)

    # add these as col in frac_prolif_parameters_pso, with dim 1 and dim2
    frac_prolif_parameters_pso[('umap', 'dim 1')] = umap_embedding[:, 0]
    frac_prolif_parameters_pso[('umap', 'dim 2')] = umap_embedding[:, 1]
    # KeyError: 'key of type tuple not found and not a MultiIndex'

    # scatter plot of umap, with the color being the likelihood
    fig, ax = plt.subplots()
    x = umap_embedding[:, 0]
    y = umap_embedding[:, 1]
    c = frac_prolif_parameters_pso[('cost', 'total cost likelihood')]
    ax.scatter(x, y, c=c, cmap='viridis', s=1)
    plt.show()


    import seaborn as sns

    color_map = "discrete"
    # color_map = "continuous"

    number_of_chains = frac_prolif_parameters_pso['pso']['chain'].max()

    colors = sns.color_palette(n_colors=number_of_chains)
    number_of_chains = frac_prolif_parameters_pso['pso']['chain'].max()
    number_of_iterations = frac_prolif_parameters_pso['pso']['iteration'].max()


    if color_map == "continuous":
        # vmap with min and max value
        vmin = frac_prolif_parameters_pso[('cost', 'total cost likelihood')].min()
        vmax = frac_prolif_parameters_pso[('cost', 'total cost likelihood')].max()
        # vmap with min and max value
        vmap = plt.cm.ScalarMappable(cmap='viridis', norm=plt.Normalize(vmin=vmin, vmax=vmax))
        # vmap with min and max value

    # sort frac_prolif_parameters_pso by iteration frac_prolif_parameters_pso['pso]['iteration']
    frac_prolif_parameters_pso = frac_prolif_parameters_pso.sort_index(level=['pso', 'iteration'])
    fig, ax = plt.subplots()

    colors = sns.color_palette(n_colors=number_of_chains)
    for chain in range(1, number_of_chains + 1):
        frac_prolif_parameters_pso_chain_tmp = frac_prolif_parameters_pso[frac_prolif_parameters_pso['pso']['chain'] == chain]
        # sort it
        frac_prolif_parameters_pso_chain_tmp = frac_prolif_parameters_pso_chain_tmp.sort_values(by=('pso', 'iteration'))
        # plot it
        plt.plot(frac_prolif_parameters_pso_chain_tmp['umap']['dim 1'], frac_prolif_parameters_pso_chain_tmp['umap']['dim 2'], linewidth=0.3, color='grey', linestyle='--')
        # plt.scatter(frac_prolif_parameters_pso_chain_tmp['umap']['dim 1'], frac_prolif_parameters_pso_chain_tmp['umap']['dim 2'], c=colors[chain], s=1)
        # <input>:7: UserWarning: *c* argument looks like a single numeric RGB or RGBA sequence, which should be avoided as value-mapping will have precedence in case its length matches with *x* & *y*.  Please use the *color* keyword-argument or provide a 2D array with a single row if you intend to specify the same RGB or RGBA value for all points
        # colors[chain] is a tuple, so it is interpreted as an RGB value, turn it into a list, remove first and last one
        if color_map == "discrete":
            color_tmp = [colors[chain - 1]]
            color_tmp_0 = [colors[chain - 1]]
            color_tmp_n = [colors[chain - 1]]
        else:
            color_tmp = frac_prolif_parameters_pso_chain_tmp[('cost', 'total cost likelihood')].tolist()[1:-1]
            color_tmp_0 = frac_prolif_parameters_pso_chain_tmp[('cost', 'total cost likelihood')].tolist()[0]
            color_tmp_n = frac_prolif_parameters_pso_chain_tmp[('cost', 'total cost likelihood')].tolist()[-1]
        # scatter as a square
        plt.scatter(frac_prolif_parameters_pso_chain_tmp['umap']['dim 1'].tolist()[1:-1], frac_prolif_parameters_pso_chain_tmp['umap']['dim 2'].tolist()[1:-1], c=color_tmp, s=1, marker='s')
        # scatter the first point as a circle
        plt.scatter(frac_prolif_parameters_pso_chain_tmp['umap']['dim 1'].tolist()[0], frac_prolif_parameters_pso_chain_tmp['umap']['dim 2'].tolist()[0], c=color_tmp_0, s=30, marker='o')
        # scatter the last point as a cross
        plt.scatter(frac_prolif_parameters_pso_chain_tmp['umap']['dim 1'].tolist()[-1], frac_prolif_parameters_pso_chain_tmp['umap']['dim 2'].tolist()[-1], c=color_tmp_n, s=30, marker='x')

    if color_map == "discrete":
        # plt.legend(['chain ' + str(i) for i in range(1, number_of_chains + 1)])
        pass
    else: # add a cbar
        plt.colorbar()

    plt.show()




    cost_columns = [('cost', 'cost 0.6 kPa'), ('cost', 'cost 55 kPa')]

    # Identify the Pareto set
    is_dominated = (frac_prolif_parameters_pso[cost_columns].values[:, None, :] < frac_prolif_parameters_pso[cost_columns].values[None, :, :]).all(axis=-1)
    is_pareto = ~is_dominated.any(axis=0)
    # nan are false
    is_pareto = is_pareto & ~frac_prolif_parameters_pso[cost_columns].isna().any(axis=1)
    frac_prolif_parameters_pso[('pareto', ' '.join([cost[1] for cost in cost_columns]))] = is_pareto


    plt.plot(frac_prolif_parameters_pso['cost', 'cost 0.6 kPa'], frac_prolif_parameters_pso['cost']['cost 55 kPa'], 'o')
    plt.plot(frac_prolif_parameters_pso['cost', 'cost 0.6 kPa'][is_pareto], frac_prolif_parameters_pso['cost']['cost 55 kPa'][is_pareto], 'o')
    plt.show()

    figure_pso_cell_count(cell_count_dict_original, frac_prolif_parameters_original, frac_prolif_parameters_pso[is_pareto], cell_count_dict_pso, conditions, run_times, means, stdevs)
    figure_pso_cell_count(cell_count_dict_original, frac_prolif_parameters_original, frac_prolif_parameters_pso[is_pareto][2:3], cell_count_dict_pso, conditions, run_times, means, stdevs)


    # ALL COLUMNS ('cost', cost) for cost in cost_columns
    cost_columns = [('cost', 'cost ' + stiffness + ' kPa ' + drug) for stiffness in ['0.6', '4', '55'] for drug in ['DMSO', 'erlotinib']]

    # Identify the Pareto set
    is_dominated = (frac_prolif_parameters_pso[cost_columns].values[:, None, :] < frac_prolif_parameters_pso[cost_columns].values[None, :, :]).all(axis=-1)
    is_pareto = ~is_dominated.any(axis=0)
    frac_prolif_parameters_pso[('pareto', ' '.join([cost[1] for cost in cost_columns]))] = is_pareto

    # select only the pareto set
    frac_prolif_parameters_pso[('pareto', ' '.join([cost[1] for cost in cost_columns]))]