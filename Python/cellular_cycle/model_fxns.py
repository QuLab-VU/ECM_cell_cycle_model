# functions for analyzing PySB model output
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import pandas as pd
import numpy as np

def plotObs(model, obs, times, tpos, title=None, ax=None, legend=True, **kwargs):
    """Plot an observable of a model"""
    if ax is None:
        ax = plt.gca()
    # plt.plot(tspan,model[obs][tspan], label=obs, lw=2)

    for o in obs:
        ax.plot(times[tpos], model[o][tpos], label=o)
    ax.set_xlabel('Time (h)')
    ax.set_ylabel('Concentration')
    if title is not None:
        ax.set_title(title)
    if legend:
        ax.legend(loc="upper right")
    return ax


# def minmax_scale(dat, keep_cols):  # original; shoudl remove keep_cols argument
#     '''
#     dat: pandas dataframe
#     keep_cols: list of columns to keep
#     '''
#     d_minmax_scaled = (dat - dat.min(axis=0)) / ((dat.max(axis=0) - dat.min(axis=0)))
#     d_minmax_scaled = d_minmax_scaled[keep_cols]
#     return d_minmax_scaled

def minmax_scale(dat):
    '''
    dat: pandas dataframe
    '''
    return (dat - dat.min(axis=0)) / ((dat.max(axis=0) - dat.min(axis=0)))


def stack_sliding_window(dat, window_size, idx):
    '''
    dat: pandas dataframe
    window_size: int
    idx: list of indices to use as center of sliding window
    '''
    half_window = int(window_size / 2)
    out = pd.DataFrame([np.array(dat[i - half_window : i + half_window].unstack()) for i in idx])
    return out

def get_model_observable_names(model_code_filepath):
    out = pd.read_csv(model_code_filepath, sep="\0", header=None)
    out = out[out[0].str.contains("Obs") == True]
    out = out[0].to_list()
    out = [x.split(" = ")[0] for x in out]
    out = [x.split("Observable('obs")[1] for x in out if "Observable('obs" in x]
    out = [x.split("', ")[0] for x in out if "', " in x]
    return out

def calc_rates(population_df):
    population_history = population_df
    parameters_sim = {}
    if len(population_history) > 0:
        number_of_deaths = len(population_history) - (population_history['Event'] == 'Apoptosis').value_counts()[0]
        cell_time_total = sum(population_history['Lifetime'])
        parameters_sim['a'] = number_of_deaths / cell_time_total
    else:
        parameters_sim['a'] = np.nan
        number_of_deaths = np.nan

    population_df = population_history[population_history['Starting Cell'] == False]
    mitotic_lifespans = population_df['Lifetime'][population_df['Event'] == 'Mitosis']
    censored_lifespans = population_df['Lifetime'][population_df['Event'] == 'End of experiment']

    if len(mitotic_lifespans) > 2:
        m = fp.mle(mitotic_lifespans, censored_lifespans)
        parameters_sim['mu'] = m['x'][0]
        parameters_sim['sigma'] = m['x'][1]
        parameters_sim['lambda'] = m['x'][2]
        parameters_sim.update(fp.q_rates(m))
    else:
        parameters_sim['d'] = np.nan
        parameters_sim['q'] = np.nan
        parameters_sim['mu'] = np.nan
        parameters_sim['sigma'] = np.nan
        parameters_sim['lambda'] = np.nan

    return parameters_sim, number_of_deaths, len(mitotic_lifespans)


def plot_cell_count(cell_count, run_time, title="Unknown condition"):
    t_sim = np.arange(0, run_time, 1 / 3)
    cell_count_t = np.zeros(len(t_sim))
    time_tmp_m = 0
    for time_tmp_p, cell_count_time_tmp in cell_count:
        cell_count_t[(time_tmp_m <= t_sim) & (t_sim < time_tmp_p)] = cell_count_time_tmp
        time_tmp_m = time_tmp_p
    pop_dubs = np.log2(cell_count_t / cell_count_t[0])
    fig, ax = plt.subplots()
    ax.scatter(t_sim, pop_dubs, color="black", s=1)
    ax.set_xlabel('Time (h)')
    ax.set_ylabel('Population doublings')
    ax.set_ylim(0,5)
    plt.text(0,4.5,title)
    plt.show()

    return fig, ax


def make_fig(cell_count, run_time, title="Unknown condition"):
    fig, ax = plt.subplots()

    ax.set_xlabel('Time (h)')
    ax.set_ylabel('Population doublings')
    ax.set_ylim(0, 5)
    plt.text(0, 4.5, title)

    return fig, ax

def add_plot(cell_count, run_time, fig, ax, title="Unknown condition", color="black"):
    t_sim = np.arange(0, run_time, 1 / 3)
    cell_count_t = np.zeros(len(t_sim))
    time_tmp_m = 0
    for time_tmp_p, cell_count_time_tmp in cell_count:
        cell_count_t[(time_tmp_m <= t_sim) & (t_sim < time_tmp_p)] = cell_count_time_tmp
        time_tmp_m = time_tmp_p
    pop_dubs = np.log2(cell_count_t / cell_count_t[0])
    ax.plot(t_sim, pop_dubs, color=color, lw=0.5, alpha=0.8)

    return

def make_cellcount_df(cell_count_dat):
    dtp = pd.DataFrame(cell_count_dat)
    dtp.columns = ['Time (h)','Cell count','condition']
    dtp['Population doublings'] = pd.Series([dtp['Cell count'][dtp['condition']==cond] for cond in dtp['condition'].unique()])

def plot_popgrowth(cell_count_df):
    sns.scatterplot(data=cell_count_df, x='Time (h)', y='Population doublings', hue='condition')


def calc_death(cells, all_events_simulated):
    EXPT_TIME = 120
    total_time = 0
    for n in range(len(cells['time'])):
        if cells['time'][n][-1] > EXPT_TIME:
            total_time = total_time + (EXPT_TIME - cells['time'][n][0])
        else:
            total_time = total_time + (cells['time'][n][-1] - cells['time'][n][0])
    try:
        num_deaths = len(all_events_simulated[all_events_simulated['event'] == 'apoptosis'])
    except:
        num_deaths = 0
    a_rate = num_deaths / total_time

    return [a_rate, num_deaths, total_time]


def run_treatment(inhib_name, inhib_conc):
    model.parameters[inhib_name].value = inhib_conc

    (check, population_df, cells, cell_count) = Cell_Population_Simulation(init, run_time, observed_cells_time, model,
                                                                           ECM_conditions, condition,
                                                                           sampled_parameter_names, thr=[0.26, 0.14],
                                                                           return_plot=False, bypass=True)
    # IMT_eoe_counts, all_prot = retrieve_cellState(run_time, cells, return_all_prot=True)
    # all_events_simulated = call_events(run_time, IMT_eoe_counts)

    return cells, cell_count, population_df


def fit_DIP(cell_count, run_time):
    t_sim = np.arange(0, run_time, 1 / 3)
    cell_count_t = np.zeros(len(t_sim))
    time_tmp_m = 20
    for time_tmp_p, cell_count_time_tmp in cell_count:
        cell_count_t[(time_tmp_m <= t_sim) & (t_sim < time_tmp_p)] = cell_count_time_tmp
        time_tmp_m = time_tmp_p
    pop_dubs = np.log2(cell_count_t / cell_count_t[0])
    dip_rate, b = np.polyfit(t_sim, pop_dubs, 1)

    return dip_rate


def plot_dip(concs, DIPs, title="Unknown Condition"):
    fig, ax = plt.subplots()

    ax.set_xlabel('Conc (um)')
    ax.set_ylabel('DIP')
    plt.xscale("log")
    ax.set_ylim(-0.01, 0.03)
    # plt.text(0,4.5,title)
    for conc in concs:
        tmp_c = [conc * 10e-7] * len(DIPs[conc])
        ax.scatter(tmp_c, DIPs[conc], color='black')

    return fig, ax


def plot_arate(concs, population_dfs, title="Unknown Condition"):
    fig1 = plt.figure()
    ax = fig1.add_subplot(1, 1, 1)

    ax.set_xlabel('Conc (um)')
    ax.set_ylabel('Apoptosis rate')
    ax.set_xscale("log")
    ax.set_xlim(1.0e-10, 1.0e-4)
    # ax.set_ylim(-0.01,0.03)
    # plt.text(0,4.5,title)
    for conc in concs:
        a_rates = []
        for num in range(3):
            tmp_rates, _, _ = calc_rates(population_dfs[conc][num])
            a_rates.append(tmp_rates['a'])

        tmp_c = [conc * 10e-7] * 3
        ax.scatter(tmp_c, a_rates, color='black')

    return fig, ax


def plot_drate(concs, population_dfs, title="Unknown Condition"):
    fig1 = plt.figure()
    ax = fig1.add_subplot(1, 1, 1)

    ax.set_xlabel('CycDi (M)')
    ax.set_ylabel('Division rate')
    ax.set_xscale("log")
    ax.set_xlim(1.0e-10, 1.0e-4)
    # ax.set_ylim(-0.01,0.03)
    # plt.text(0,4.5,title)
    for conc in concs:
        d_rates = []
        for num in range(3):
            tmp_rates, _, _ = calc_rates(population_dfs[conc][num])
            d_rates.append(tmp_rates['d'])

        tmp_c = [conc * 10e-7] * 3
        ax.scatter(tmp_c, d_rates, color='black')

    return fig, ax


def plot_qrate(concs, population_dfs, title="Unknown Condition"):
    fig, ax = plt.subplots()

    ax.set_xlabel('Conc (um)')
    ax.set_ylabel('Quiescence rate')
    plt.xscale("log")
    # ax.set_ylim(-0.01,0.03)
    # plt.text(0,4.5,title)
    for conc in concs:
        q_rates = []
        for num in range(3):
            tmp_rates = calc_rates(population_dfs[conc][num])
            q_rates.append(tmp_rates['q'])

        tmp_c = [conc * 10e-7] * 3
        ax.scatter(tmp_c, q_rates, color='black')

    return fig, ax


def plot_daratio(concs, population_dfs, title="Unknown Condition"):
    fig, ax = plt.subplots()

    ax.set_xlabel('Conc (um)')
    ax.set_ylabel('Division/Apoptosis rate')
    plt.xscale("log")
    # ax.set_ylim(-0.01,0.03)
    # plt.text(0,4.5,title)
    for conc in concs:
        d_rates = []
        a_rates = []
        for num in range(3):
            tmp_rates = calc_rates(population_dfs[conc][num])
            d_rates.append(tmp_rates['d'])
            a_rates.append(tmp_rates['a'])

            d_a = np.mean(d_rates) / np.mean(a_rates)

        tmp_c = [conc * 10e-7] * 1
        ax.scatter(tmp_c, d_a, color='black')

    return fig, ax


def plot_divdeath_totals(concs, population_dfs, title="Unknown Condition"):
    fig, ax = plt.subplots()

    ax.set_xlabel('Conc (um)')
    plt.xscale("log")

    for conc in concs:
        d_rates = []
        a_rates = []
        for num in range(3):
            tmp_rates, death_num, div_num = calc_rates(population_dfs[conc][num])
            tot_events = death_num + div_num
            a_rates.append(death_num / tot_events)
            d_rates.append(div_num / tot_events)

        tmp_c = [conc * 10e-7] * 3
        # ax.scatter(tmp_c,d_rates, color='tab:blue')
        ax.scatter(tmp_c, a_rates, color='tab:red')

    return fig, ax


def plot_numberdeaths(concs, population_dfs, title="Unknown Condition"):
    fig1 = plt.figure()
    ax = fig1.add_subplot(1, 1, 1)

    ax.set_xlabel('Conc (um)')
    ax.set_ylabel('Death Count')
    ax.set_xscale("log")
    ax.set_xlim(1.0e-10, 1.0e-4)
    # ax.set_ylim(-0.01,0.03)
    # plt.text(0,4.5,title)
    for conc in concs:
        a_rates = []
        deaths = []
        for num in range(3):
            tmp_rates, tot_deaths, _ = calc_rates(population_dfs[conc][num])
            a_rates.append(tmp_rates['a'])
            deaths.append(tot_deaths)

        tmp_c = [conc * 10e-7] * 3
        ax.scatter(tmp_c, deaths, color='black')

    return fig, ax


def plot_death_avgdiv(concs, population_dfs, title="Unknown Condition"):
    fig1 = plt.figure()
    ax = fig1.add_subplot(1, 1, 1)

    ax.set_xlabel('Conc (um)')
    ax.set_ylabel('Death Count/Avg. Division Time')
    ax.set_xscale("log")
    ax.set_xlim(1.0e-10, 1.0e-4)
    # ax.set_ylim(-0.01,0.03)
    # plt.text(0,4.5,title)
    n = 0
    for conc in concs:
        a_rates = []
        deaths = []
        deathperdiv = []
        for num in range(3):
            tmp_rates, tot_deaths, _ = calc_rates(population_dfs[conc][num])
            a_rates.append(tmp_rates['a'])
            deaths.append(tot_deaths)
            deathperdiv.append(tot_deaths / avg_div_len[n])

        tmp_c = [conc * 10e-7]
        ax.scatter(tmp_c, tot_deaths / avg_div_len[n], color='black')
        n += 1
    return fig, ax


def runIMTs(init, run_time, observed_cells_time, model, ECM_conditions, condition, sampled_parameter_names, IMTs_count,
            inhib):
    all_IMTs = []
    n = 0
    while len(all_IMTs) < IMTs_count:
        check, population_df, cells, cell_count = Cell_Population_Simulation(init, run_time, observed_cells_time, model,
                                                                             ECM_conditions, condition,
                                                                             sampled_parameter_names, thr=[0.26, 0.14],
                                                                             return_plot=False, bypass=True)
        all_events_simulated = call_events(run_time, IMT_eoe_counts)

        mits = all_events_simulated.loc[all_events_simulated['event'] == 'mitosis'].reset_index()
        mits_list = mits['Lifespan'].tolist()
        mits_list_r = [round(x, 1) for x in mits_list]

        IMT_tmp, s_m_tmp, g1_tmp = SMs_G1s_IMTs_extraction(cells)

        if n == 0:
            all_IMTs = mits_list_r
            all_g1s = g1_tmp
        else:
            all_IMTs.extend(mits_list_r)
            all_g1s.extend(g1_tmp)

        # filename_cell_count = Path(SAVE_DIR, f'cell_counts_{condition}_{n}_cycdi_0_4.csv')
        # np.savetxt(filename_cell_count, cell_count, delimiter=',')

        n += 1

    filename_IMT = Path(SAVE_DIR, f'IMTs_{condition}_{inhib}.csv')
    filename_g1 = Path(SAVE_DIR, f'G1s_{condition}_{inhib}.csv')
    np.savetxt(filename_IMT, all_IMTs, delimiter=',')
    np.savetxt(filename_g1, all_g1s, delimiter=',')

