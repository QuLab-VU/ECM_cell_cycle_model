# functions for analyzing PySB model output
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import pandas as pd
import numpy as np
from Cell_Population_Simulation import Cell_Population_Simulation
from call_cell_events import call_events
from retrieve_cellState import retrieve_cellState
from SMs_G1s_IMTs_extraction import SMs_G1s_IMTs_extraction
from pathlib import Path
from typing import List
import fracprolif as fp
from scipy.ndimage import gaussian_filter
from scipy.signal import find_peaks
from statistics import mean


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


# def minmax_scale(dat, keep_cols):  # original; should remove keep_cols argument
#     '''
#     dat: pandas dataframe
#     keep_cols: list of columns to keep
#     '''
#     d_minmax_scaled = (dat - dat.min(axis=0)) / ((dat.max(axis=0) - dat.min(axis=0)))
#     d_minmax_scaled = d_minmax_scaled[keep_cols]
#     return d_minmax_scaled

def minmax_scale(dat: pd.DataFrame) -> pd.DataFrame:
    return (dat - dat.min(axis=0)) / ((dat.max(axis=0) - dat.min(axis=0)))


def stack_sliding_window(dat: pd.DataFrame, window_size: int, idx: List[int]) -> pd.DataFrame:
    '''
    dat: pandas dataframe
    window_size: int
    idx: list of indices to use as center of sliding window
    '''
    half_window = int(window_size / 2)
    out = pd.DataFrame([np.array(dat[i - half_window : i + half_window].unstack()) for i in idx])
    return out

def get_model_observable_names(model_code_filepath: str) -> List[str]:
    out = pd.read_csv(model_code_filepath, sep="\0", header=None)
    out = out[out[0].str.contains("Obs") == True]
    out = out[0].to_list()
    out = [x.split(" = ")[0] for x in out]
    out = [x.split("Observable('obs")[1] for x in out if "Observable('obs" in x]
    out = [x.split("', ")[0] for x in out if "', " in x]
    return out

def calc_rates(population_df, fp=fp):
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


def run_treatment(model, inhib_name, inhib_conc):
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

    ax.set_xlabel('Conc (µm)')
    ax.set_ylabel('DIP')
    plt.xscale("log")
    ax.set_ylim(-0.01, 0.03)
    # plt.text(0,4.5,title)
    for conc in concs:
        tmp_c = [conc * 10e-7] * len(DIPs[conc])
        ax.scatter(tmp_c, DIPs[conc], color='black')

    return fig, ax


def plot_arate(concs, population_dfs, title="Unknown Condition"):
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)

    ax.set_xlabel('Conc (µm)')
    ax.set_ylabel('Apoptosis rate')
    ax.set_xscale("log")
    ax.set_xlim(1.0e-10, 1.0e-4)
    # ax.set_ylim(-0.01,0.03)
    # plt.text(0,4.5,title)
    n = len(population_dfs[0])
    for conc in concs:
        a_rates = []
        for num in range(n):
            tmp_rates, _, _ = calc_rates(population_dfs[conc][num])
            a_rates.append(tmp_rates['a'])

        tmp_c = [conc * 10e-7] * n
        ax.scatter(tmp_c, a_rates, color='black')

    return fig, ax



def plot_drate(concs, population_dfs, title="Unknown Condition"):
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)

    ax.set_xlabel('CDK4/6i (M)')
    ax.set_ylabel('Division rate')
    ax.set_xscale("log")
    ax.set_xlim(1.0e-10, 1.0e-4)
    # ax.set_ylim(-0.01,0.03)
    # plt.text(0,4.5,title)
    n = len(population_dfs[0])
    for conc in concs:
        d_rates = []
        for num in range(n):
            tmp_rates, _, _ = calc_rates(population_dfs[conc][num])
            d_rates.append(tmp_rates['d'])

        tmp_c = [conc * 10e-7] * n
        ax.scatter(tmp_c, d_rates, color='black')

    return fig, ax


def plot_qrate(concs, population_dfs, title="Unknown Condition"):
    fig, ax = plt.subplots()

    ax.set_xlabel('Conc (µm)')
    ax.set_ylabel('Quiescence rate')
    plt.xscale("log")
    # ax.set_ylim(-0.01,0.03)
    # plt.text(0,4.5,title)
    n = len(population_dfs[0])

    for conc in concs:
        q_rates = []
        for num in range(3):
            tmp_rates = calc_rates(population_dfs[conc][num])
            q_rates.append(tmp_rates['q'])

        tmp_c = [conc * 10e-7] * n
        ax.scatter(tmp_c, q_rates, color='black')

    return fig, ax


def plot_daratio(concs, population_dfs, title="Unknown Condition"):
    fig, ax = plt.subplots()

    ax.set_xlabel('Conc (µm)')
    ax.set_ylabel('Division/Apoptosis rate')
    plt.xscale("log")
    # ax.set_ylim(-0.01,0.03)
    # plt.text(0,4.5,title)
    n = len(population_dfs[0])
    for conc in concs:
        d_rates = []
        a_rates = []
        for num in range(n):
            tmp_rates = calc_rates(population_dfs[conc][num])
            d_rates.append(tmp_rates['d'])
            a_rates.append(tmp_rates['a'])

            d_a = np.mean(d_rates) / np.mean(a_rates)

        tmp_c = [conc * 10e-7] * n
        ax.scatter(tmp_c, d_a, color='black')

    return fig, ax


def plot_divdeath_totals(concs, population_dfs, title="Unknown Condition"):
    fig, ax = plt.subplots()

    ax.set_xlabel('Conc (µm)')
    plt.xscale("log")

    for conc in concs:
        d_rates = []
        a_rates = []
        n = len(population_dfs[0])
        for num in range(n):
            tmp_rates, death_num, div_num = calc_rates(population_dfs[conc][num])
            tot_events = death_num + div_num
            a_rates.append(death_num / tot_events)
            d_rates.append(div_num / tot_events)

        tmp_c = [conc * 10e-7] * n
        # ax.scatter(tmp_c,d_rates, color='tab:blue')
        ax.scatter(tmp_c, a_rates, color='tab:red')

    return fig, ax

def plot_numberdeaths(concs, population_dfs, title="Unknown Condition"):
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)

    ax.set_xlabel('Conc (µm)')
    ax.set_ylabel('Death Count')
    ax.set_xscale("log")
    ax.set_xlim(1.0e-10, 1.0e-4)
    # ax.set_ylim(-0.01,0.03)
    # plt.text(0,4.5,title)
    n = len(population_dfs[0])
    for conc in concs:
        a_rates = []
        deaths = []
        for num in range(n):
            tmp_rates, tot_deaths, _ = calc_rates(population_dfs[conc][num])
            a_rates.append(tmp_rates['a'])
            deaths.append(tot_deaths)

        tmp_c = [conc * 10e-7] * n
        ax.scatter(tmp_c, deaths, color='black')

    return fig, ax


def plot_death_avgdiv(concs, population_dfs, title="Unknown Condition"):
    fig1 = plt.figure()
    ax = fig1.add_subplot(1, 1, 1)

    ax.set_xlabel('Conc (µm)')
    ax.set_ylabel('Death Count/Avg. Division Time')
    ax.set_xscale("log")
    ax.set_xlim(1.0e-10, 1.0e-4)
    # ax.set_ylim(-0.01,0.03)
    # plt.text(0,4.5,title)
    n_meas = len(population_dfs[0])

    n = 0
    for conc in concs:
        a_rates = []
        deaths = []
        deathperdiv = []
        for num in range(n_meas):
            tmp_rates, tot_deaths, _ = calc_rates(population_dfs[conc][num])
            a_rates.append(tmp_rates['a'])
            deaths.append(tot_deaths)
            deathperdiv.append(tot_deaths / avg_div_len[n])

        tmp_c = [conc * 10e-7]
        ax.scatter(tmp_c, tot_deaths / avg_div_len[n], color='black')
        n += 1
    return fig, ax


def runIMTs(init, run_time, observed_cells_time, model, ECM_conditions, condition, 
            sampled_parameter_names, IMTs_count, inhib, SAVE_DIR, verbose=False):
    all_IMTs = []
    n = 0
    while len(all_IMTs) < IMTs_count:
        check, population_df, cells, cell_count = Cell_Population_Simulation(init, run_time, observed_cells_time, model,
                                                                             ECM_conditions, condition,
                                                                             sampled_parameter_names, thr=[0.26, 0.14],
                                                                             return_plot=False, bypass=True)
        IMT_eoe_counts, all_prot = retrieve_cellState(run_time, cells, return_all_prot=True, verbose=verbose)
        all_events_simulated = call_events(run_time, IMT_eoe_counts)

        mits = all_events_simulated.loc[all_events_simulated['event'] == 'mitosis'].reset_index()
        mits_list = mits['Lifespan'].tolist()
        mits_list_r = [round(x, 1) for x in mits_list]

        IMT_tmp, s_m_tmp, g1_tmp = SMs_G1s_IMTs_extraction(cells)

        if n == 0:
            all_IMTs = mits_list_r
            all_g1s = g1_tmp
            all_sm = s_m_tmp
        else:
            all_IMTs.extend(mits_list_r)
            all_g1s.extend(g1_tmp)
            all_sm.extend(s_m_tmp)
        if verbose:
            print(f'runIMT round {n} has {len(all_IMTs)} IMT values')
        n += 1

    filename_IMT = SAVE_DIR / f'IMTs_{condition}_{inhib}.csv'
    filename_g1 = SAVE_DIR / f'G1s_{condition}_{inhib}.csv'
    filename_sm = SAVE_DIR / f'SMs_{condition}_{inhib}.csv'
    print(f'Saving {filename_IMT}, {filename_g1}, and {filename_sm}...')
    np.savetxt(filename_IMT, all_IMTs, delimiter=',')
    np.savetxt(filename_g1, all_g1s, delimiter=',')
    np.savetxt(filename_sm, all_sm, delimiter=',')


def load_sim_output(filenames: list) -> pd.DataFrame:
    # Function copied from 05_UMAP_of_stoch_sims.ipynb 20240821
    sim_data = pd.DataFrame()
    dfs = []
    for i, f in enumerate(filenames):
        # Extract base file name without extension
        base_name = Path(f).parts[-1].replace('.csv','')

        df = pd.read_csv(f, sep=",", header=None, usecols = range(0,19))
        df.columns = ["skp2", "p21", "pp21", "cdh1_act", "cdh1_inact", "myc", "cycD", "cycD-p21", 
        "Rb", "pRb", "E2F", "E2F-Rb", "cycE", "cycE-p21", "cycA", "cycA-p21",
        "cycB", "cycB-p21", "cdc20"]
        df = df.drop(columns=['cdh1_inact'])
        df.index = pd.MultiIndex.from_arrays(
            [[base_name] * len(df), range(1, len(df) + 1)],
            names=['File Name', 'Time index']
        )
        dfs.append(df)
    sim_data = pd.concat(dfs)
    return sim_data

def calc_minmax(filenames):
    # Function copied from 05_UMAP_of_stoch_sims.ipynb 20240821
    dat = load_sim_output(filenames)
    min_prots = dat.min(axis=0).to_list()
    max_prots = dat.max(axis=0).tolist()
    return max_prots, min_prots

def load_norm_data(filenames, filenames_ode, max_prots, min_prots, time_all, init_dt, dt, window, time_step, plot_phase):
    # Function copied from 05_UMAP_of_stoch_sims.ipynb 20240821
    if plot_phase:
        plot_phase(filenames, filenames_ode)
        plt.savefig(Path(SAVE_DIR, f"phaseplot_{num}.png"), bbox_inches='tight')
        plot_phase_norm(filenames, filenames_ode, max_prots, min_prots)
        plt.savefig(Path(SAVE_DIR, f"phaseplot_norm_{num}.png"), bbox_inches='tight')
    
    win_len = int(window/dt)
    all_prot_tmp = np.empty((0,win_len*18))     #careful, this could end up being massive 
    all_prot_tmp_ode = np.empty((0,win_len*18))
    all_time_tmp = np.empty((0,win_len))
    all_time_at_tmp = np.empty((0,1))
    range_ind = np.arange(int(win_len/2),int(len(time_all)-(int(win_len/2))),int(time_step/dt))
        
    for (f, f_ode) in zip(filenames, filenames_ode):
        prot_tmp = pd.read_csv(f, sep=",", header=None, usecols = range(0,19))
        prot_tmp.columns = ["skp2", "p21", "pp21", "cdh1_act", "cdh1_inact", "myc", "cycD", "cycD-p21", 
        "Rb", "pRb", "E2F", "E2F-Rb", "cycE", "cycE-p21", "cycA", "cycA-p21",
        "cycB", "cycB-p21", "cdc20"]
        prot_tmp = prot_tmp.drop(columns=['cdh1_inact'])
        #prot_tmp = prot_tmp/200
        for ii in range(18):
            prot_tmp.iloc[:,ii] = (prot_tmp.iloc[:,ii] - min_prots[ii])/(max_prots[ii] - min_prots[ii])

        prot_re_tmp = np.empty((len(range_ind),win_len*18))
        time_re_tmp = np.empty((len(range_ind),win_len))
        time_at_tmp = np.empty((len(range_ind),1))

        ii=0
        for n in range(int(win_len/2),int(len(time_all)-(int(win_len/2))),int(time_step/dt)):
            tmp = prot_tmp.iloc[n-int(win_len/2):n+int(win_len/2)]
            tmp_t = time_all['time'].iloc[n-int(win_len/2):n+int(win_len/2)]
            tmp=tmp.to_numpy()
            tmp=np.reshape(np.transpose(tmp),win_len*18)
            tmp_at = time_all['time'].iloc[n]

            prot_re_tmp[ii,:] = np.array([tmp])
            time_re_tmp[ii,:] = np.array([tmp_t])
            time_at_tmp[ii,:] = np.array([[tmp_at]])
            ii+=1

        all_prot_tmp = np.append(all_prot_tmp, prot_re_tmp, axis=0)
        all_time_tmp = np.append(all_time_tmp, time_re_tmp, axis=0)
        all_time_at_tmp = np.append(all_time_at_tmp, time_at_tmp, axis=0)
        
        prot_tmp_ode = pd.read_csv(f_ode, sep=",", header=None, usecols = range(0,19))
        prot_tmp_ode.columns = ["skp2", "p21", "pp21", "cdh1_act", "cdh1_inact", "myc", "cycD", "cycD-p21", 
        "Rb", "pRb", "E2F", "E2F-Rb", "cycE", "cycE-p21", "cycA", "cycA-p21",
        "cycB", "cycB-p21", "cdc20"]
        prot_tmp_ode = prot_tmp_ode.drop(columns=['cdh1_inact'])
        for ii in range(18):
            prot_tmp_ode.iloc[:,ii] = (prot_tmp_ode.iloc[:,ii] - min_prots[ii])/(max_prots[ii] - min_prots[ii])

        prot_re_tmp = np.empty((len(range_ind),win_len*18))


        ii=0
        for n in range(int(win_len/2),int(len(time_all)-(int(win_len/2))),int(time_step/dt)):
            tmp = prot_tmp_ode.iloc[n-int(win_len/2):n+int(win_len/2)]
            tmp=tmp.to_numpy()
            tmp=np.reshape(np.transpose(tmp),win_len*18)

            prot_re_tmp[ii,:] = np.array([tmp])
            time_re_tmp[ii,:] = np.array([tmp_t])
            time_at_tmp[ii,:] = np.array([[tmp_at]])
            ii+=1

        all_prot_tmp_ode = np.append(all_prot_tmp_ode, prot_re_tmp, axis=0)
        
    return all_prot_tmp, all_prot_tmp_ode, all_time_tmp, all_time_at_tmp
    
def get_traj(embed_dat, pos, trk_len):
    if pos > len(embed_dat)/trk_len:
        print('Position out of range')
        return np.nan   
    embed_traj = embed_dat[pos*trk_len:(pos+1)*trk_len]
    return embed_traj


def plot_traj(embed_dat, pos, trk_len, ax=None, color='gray', title=None):
    if ax is None:
        fig, ax = plt.subplots(figsize=(5,5))
    
    params = {'axes.titley': 1,
            'axes.titlepad': 5,
            'axes.labelsize': 16,
            'axes.edgecolor': "black",
            'axes.titlesize': 18,
            'xtick.labelsize': 14,
            'ytick.labelsize': 14}
    plt.rcParams.update(params)

    embed_ode = get_traj(embed_dat, pos, trk_len)
    ax.set_title(title)
    ax.plot(embed_ode[:,0], embed_ode[:,1], color=color, lw=0.1)
    return ax


def restructProt(protDF, embedDF, timeDF, win_len=int(4/0.1)): ## win_len=int(window/dt)
    peaks_E, _ = find_peaks(protDF.iloc[:,int((win_len*11+win_len/2)-1)], prominence = 30/200, distance = 10)
    peaks_cdc, _ = find_peaks(protDF.iloc[:,int((win_len*17+win_len/2)-1)], prominence = 30/200, distance = 10)
    
    mit_win = pd.DataFrame()
    time_win = pd.DataFrame()
    umap_win = pd.DataFrame()
    num = 1
    
    for n in range(2,len(peaks_cdc)-3):   #make sure we have divisions at ends 
        tmp_mit = peaks_cdc[n]
        tmp_mit2 = peaks_cdc[n+1]
        tmp_mit3 = peaks_cdc[n+2]
        tmp_mit4 = peaks_cdc[n+3]
        
        tmp_g1_s = peaks_E[(tmp_mit < peaks_E) & (peaks_E < tmp_mit2)]
        tmp_g1_s_2 = peaks_E[(tmp_mit2 < peaks_E) & (peaks_E < tmp_mit3)]
        tmp_g1_s_3 = peaks_E[(tmp_mit3 < peaks_E) & (peaks_E < tmp_mit4)]    
             
        if (np.size(tmp_g1_s) > 0) & (np.size(tmp_g1_s_2)>0) & (np.size(tmp_g1_s_3)>0):
            ordered = [tmp_mit, tmp_g1_s[0], tmp_mit2, tmp_g1_s_2[0], tmp_mit3, tmp_g1_s_3[0], tmp_mit4]
            g1_len_tmp = np.array(timeDF.iloc[ordered[3],0]) - np.array(timeDF.iloc[[ordered[2]],0])
            tmp_prot = protDF.iloc[int(ordered[1]):int(ordered[4]),:]
            tmp_umap = embedDF.iloc[int(ordered[1]):int(ordered[4]),:]
            tmp_time = timeDF.iloc[int(ordered[1]):int(ordered[4]),:]
            tmp = [g1_len_tmp]*len(tmp_time)
            tmp_g1s = [int(ordered[1])]*len(tmp_time)
            tmp_mit = [int(ordered[2])]*len(tmp_time)
            tmp_g1s2 = [int(ordered[3])]*len(tmp_time)
            tmp_mit2 = [int(ordered[4])]*len(tmp_time)
            
            tmp_time.insert(0, 'mit2', tmp_mit2)
            tmp_time.insert(0, 'G1-S_2', tmp_g1s2)
            tmp_time.insert(0, 'mit', tmp_mit)
            tmp_time.insert(0, 'G1-S', tmp_g1s)
            tmp_time.insert(0, 'G1_len', tmp)
            tmp_time.insert(0, 'div', [num]*len(tmp_time))
            num += 1
            
            mit_win = pd.concat([mit_win,tmp_prot],ignore_index=True)
            time_win = pd.concat([time_win, tmp_time], ignore_index=True)
            umap_win = pd.concat([umap_win, tmp_umap], ignore_index=True)
            
    return mit_win, time_win, umap_win

def cosine_similarity(X, Y):
    '''
    Calculate the cosine similarity between two vectors X and Y.
    '''
    # Convert lists to numpy arrays if they are not already
    X = np.array(X)
    Y = np.array(Y)
    
    # Calculate the dot product
    dot_product = np.dot(X, Y)
    
    # Calculate the norm (magnitude) of each vector
    norm_X = np.linalg.norm(X)
    norm_Y = np.linalg.norm(Y)
    
    # Calculate cosine similarity
    similarity = dot_product / (norm_X * norm_Y)
    
    return similarity

def get_cos_similarities(traj, v_len, offset):
    '''
    Calculate the cosine similarities between the vectors in a trajectory.
    '''
    cos_similarities = []
    
    # Calculate the cosine similarity between each vector in the trajectory and the vector v_len steps ahead of it
    for i in range(len(traj) - v_len - offset):
        v1 = traj[i + v_len,:] - traj[i,:]
        v2 = traj[i + v_len + offset,:] - traj[i + offset,:]
        cos_similarities.append(cosine_similarity(v1,v2))
    
    return cos_similarities