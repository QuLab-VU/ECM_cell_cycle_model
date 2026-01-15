import os
import sys
import numpy as np
import pandas as pd
from pathlib import Path
from numpy.linalg import norm
import datetime
import math
import multiprocess as mp
import time
import pickle
from pysb.simulator import BngSimulator
from pysb.simulator import ScipyOdeSimulator
from scipy.signal import find_peaks
import pyarrow.feather as feather
sys.path.insert(0, "./cellular_cycle")
from cell_cycle_death_model import model

PARAM = pd.read_csv(os.getcwd() + '/model_data/parameters.csv')  # parameter file to use
INIT_CONCS = pd.read_csv(os.getcwd() + '/model_data/initial_concentrations.csv')  # initial concentrations for simulations
TREATMENT = 'LOW'  # set AKT, ERK, and FAK by setting 'ERL', 'DMSO', or 'LOW'
CYCDI_CONC = 0
POINTS = 5  # number to divide each condition into, sets the number of points to divide cell cycle
PROCESSES = 5  # number of pool processes to run at once
OVERWRITE = True  # Overwrite existing simulation files
VERBOSE = True
SET_JPARAM = False
PARAM_NUM = 0
integrator_options = {
    "rtol": 1e-12,
    "atol": 1e-12
}


# Where to save large simulation output files
SAVE_DIR = Path(Path.home(), 'ECM_mod_output')
if not SAVE_DIR.exists():
    SAVE_DIR.mkdir()


# def run_lyaps(kk, np=np, model=model, points=POINTS, ScipyOdeSimulator=ScipyOdeSimulator,
#               find_peaks=find_peaks, integrator_options=integrator_options, norm=norm, math=math,
#               verbose=VERBOSE, peaks_half2=peaks_half2, y_out=y_out):
    
#     if len(peaks_half2) > 3:
#         IMTs = [t_0[y] - t_0[x] for x, y in zip(peaks_half2, peaks_half2[1:])]
#         meanIMT = sum(IMTs) / len(IMTs)
#         variable_range = y_out[peaks_half[-3]:peaks_half[-2]]
#         tmp_loc = (round(len(variable_range) / (points)) * (kk))
#     else:
#         meanIMT = 0
#         variable_range = y_out[575000:600000]
#         tmp_loc = (round(len(variable_range) / (points)) * (kk))
#     if verbose:
#         print(meanIMT)


#     Initial_1 = variable_range[tmp_loc]
#     Initial_1 = [val + ep for val in Initial_1]

#     prots = ['Skp2_0', 'p21_0', 'p21p_0', 'cdh1act_0', 'cdh1inact_0', 'Myc_0', 'CycD_0', 'CycDp21_0', 'pRb_0',
#              'pRbp_0', 'E2F_0', 'pRbE2F_0', 'CycE_0', 'CycEp21_0', 'CycA_0', 'CycAp21_0', 'CycB_0', 'CycBp21_0',
#              'Cdc20_0']
#     prot_init = {prots[i]: Initial_1[i] for i in range(len(prots))}

#     for key in prot_init.keys():
#         model.parameters[key].value = prot_init[key] / model.parameters['omega'].value

#     sim_ode1 = ScipyOdeSimulator(model, t, verbose=False, compiler='python',
#                                  integrator_options=integrator_options).run()
#     y_out1 = sim_ode1.observables

#     z = [0] * (round((len(y_out1) / 2000)) + 1)
#     n = 0
#     for i in range(0, len(y_out1), 2000):
#         if len(peaks_half2) > 0:
#             window = y_out[peaks_half[-3]+i+tmp_loc-5000:peaks_half[-2]+i+tmp_loc+5000]
#         else:
#             window = y_out[575000+i+tmp_loc-5000:600000+i+tmp_loc+5000]

#         dists = [0] * len(window)
#         for j in range(0, len(window)):
#             dif = np.asarray(window[j].tolist())-np.asarray(y_out1[i].tolist())
#             dists[j] = norm(dif[0:19]/200)

#         z[n] = min(dists)
#         n = n + 1

#     logz = [math.log(num) for num in z]
#     all_logz[kk + 1] = logz

#     return all_logz, meanIMT


CURDATE = datetime.datetime.now().strftime('%Y%m%d')

t = np.linspace(0, 1000, 1000001)
# t = np.linspace(0, 20000, 200001)
# t = np.linspace(0, 50000, 500001)

model.parameters['AKT'].value = 1
model.parameters['ERK'].value = 1
model.parameters['FAK'].value = 1

# Importing the parameters:

for key in PARAM.keys():
    try:
        if model.parameters[key].value != PARAM[key][0]:
            model.parameters[key].value = PARAM[key]
    except:
        print(key)
model.parameters['omega'].value = 200

if SET_JPARAM:
    with open('fit_data/pso/frac_prolif_parameters_pso.pkl', 'rb') as f:
        frac_prolif_parameters_pso = pickle.load(f)

    frac_prolif_parameters_pso[('cost', 'custom cost populations')] = frac_prolif_parameters_pso[(
    'cost', 'cost population dynamics 0.6 kPa DMSO')] + frac_prolif_parameters_pso[(
    'cost', 'cost population dynamics 4 kPa erlotinib')] + frac_prolif_parameters_pso[(
    'cost', 'cost population dynamics 55 kPa erlotinib')]
    pso_sort_pops = frac_prolif_parameters_pso.sort_values(by=[('cost', 'custom cost populations')])
    new_params = pso_sort_pops['model parameters']
    # Importing the parameters:
    for key in new_params.keys():
        try:
            if model.parameters[key].value != new_params[key].iloc[PARAM_NUM]:
                model.parameters[key].value = new_params[key].iloc[PARAM_NUM]
        except:
            print(key)
    model.parameters['omega'].value = 200


for key in INIT_CONCS.keys():
    try:
        model.parameters[key].value = INIT_CONCS[key][0] / model.parameters['omega'].value
    except:
        print(key)

model.parameters['CycDi_0'].value = CYCDI_CONC

if TREATMENT == 'ERL':
    FAKs = list(np.linspace(0,0.4,51))
    tmp = list(np.linspace(0.4,0.6,51))
    FAKs.extend(tmp[1:])
    tmp = list(np.linspace(0.6,1,51))
    FAKs.extend(tmp[1:])

    AKTs = list(np.linspace(0,0.3,51))
    tmp = list(np.linspace(0.3,0.6,51))
    AKTs.extend(tmp[1:])
    tmp = [0.6]*50
    AKTs.extend(tmp)
    
    ERKs = [0.2]*151
    
elif TREATMENT == 'DMSO':
    FAKs = list(np.linspace(0,0.4,51))
    tmp = list(np.linspace(0.4,0.6,51))
    FAKs.extend(tmp[1:])
    tmp = list(np.linspace(0.6,1,51))
    FAKs.extend(tmp[1:])

    AKTs = list(np.linspace(0,0.6,51))
    tmp = list(np.linspace(0.6,1,51))
    AKTs.extend(tmp[1:])
    tmp = [1]*50
    AKTs.extend(tmp)
    
    ERKs = [1]*151
else:
    FAKs = list(np.linspace(0,0.4,51))
    tmp = list(np.linspace(0.4,0.6,51))
    FAKs.extend(tmp[1:])
    tmp = list(np.linspace(0.6,1,51))
    FAKs.extend(tmp[1:])

    AKTs = list(np.linspace(0,0.3,51))
    tmp = list(np.linspace(0.3,0.6,51))
    AKTs.extend(tmp[1:])
    tmp = [0.6]*50
    AKTs.extend(tmp)
    
    ERKs = [0.1]*151

conds = [(AKT, ERK, FAK) for AKT, ERK, FAK, in zip(AKTs, ERKs, FAKs)]
# conds = conds[0:1]

for cond in conds:
    model.parameters['AKT'].value = cond[0]
    model.parameters['ERK'].value = cond[1]
    model.parameters['FAK'].value = cond[2]
    ep = 0.001 * 200
    t = np.linspace(0, 1000, 1000001)
    t_0 = np.linspace(0, 2000, 2000001)
    all_logz = dict.fromkeys(list(range(1, POINTS + 1)))

    sim_ode = ScipyOdeSimulator(model, t_0, verbose=False, compiler='python',
                                integrator_options=integrator_options).run()
    y_out = sim_ode.observables

    peaks_cdc, _ = find_peaks(y_out['obsCdc20'], prominence=0.6, distance=10)
    # new_y = y_out[peaks_cdc[-11]:peaks_cdc[-2]]
    peaks_half = peaks_cdc[peaks_cdc < len(y_out)/2]
    peaks_half2 = peaks_cdc[peaks_cdc > len(y_out) / 2]


    def run_lyaps(kk, np=np, model=model, points=POINTS, ScipyOdeSimulator=ScipyOdeSimulator,
                  find_peaks=find_peaks, integrator_options=integrator_options, norm=norm, math=math,
                  verbose=VERBOSE, peaks_half2=peaks_half2, y_out=y_out):
        
        if len(peaks_half2) > 3:
            IMTs = [t_0[y] - t_0[x] for x, y in zip(peaks_half2, peaks_half2[1:])]
            meanIMT = sum(IMTs) / len(IMTs)
            variable_range = y_out[peaks_half[-3]:peaks_half[-2]]
            tmp_loc = (round(len(variable_range) / (points)) * (kk))
        else:
            meanIMT = 0
            variable_range = y_out[575000:600000]
            tmp_loc = (round(len(variable_range) / (points)) * (kk))
        if verbose:
            print(meanIMT)
    
    
        Initial_1 = variable_range[tmp_loc]
        Initial_1 = [val + ep for val in Initial_1]
    
        prots = ['Skp2_0', 'p21_0', 'p21p_0', 'cdh1act_0', 'cdh1inact_0', 'Myc_0', 'CycD_0', 'CycDp21_0', 'pRb_0',
                 'pRbp_0', 'E2F_0', 'pRbE2F_0', 'CycE_0', 'CycEp21_0', 'CycA_0', 'CycAp21_0', 'CycB_0', 'CycBp21_0',
                 'Cdc20_0']
        prot_init = {prots[i]: Initial_1[i] for i in range(len(prots))}
    
        for key in prot_init.keys():
            model.parameters[key].value = prot_init[key] / model.parameters['omega'].value
    
        sim_ode1 = ScipyOdeSimulator(model, t, verbose=False, compiler='python',
                                     integrator_options=integrator_options).run()
        y_out1 = sim_ode1.observables
    
        z = [0] * (round((len(y_out1) / 2000)) + 1)
        n = 0
        for i in range(0, len(y_out1), 2000):
            if len(peaks_half2) > 0:
                window = y_out[peaks_half[-3]+i+tmp_loc-5000:peaks_half[-2]+i+tmp_loc+5000]
            else:
                window = y_out[575000+i+tmp_loc-5000:600000+i+tmp_loc+5000]
    
            dists = [0] * len(window)
            for j in range(0, len(window)):
                dif = np.asarray(window[j].tolist())-np.asarray(y_out1[i].tolist())
                dists[j] = norm(dif[0:19]/200)
    
            z[n] = min(dists)
            n = n + 1
    
        logz = [math.log(num) for num in z]
        all_logz[kk + 1] = logz
    
        return all_logz, meanIMT

    kk = list(range(POINTS))
    start = time.localtime()
    print(f'Start time is {time.strftime("%Y-%m-%d %H:%M:%S",start)}')
    results = []
    with mp.Pool(processes = PROCESSES) as p:        
        results = p.map(run_lyaps, kk)
    end = time.localtime()
    print(f'End time is {time.strftime("%Y-%m-%d %H:%M:%S",end)}')
    print(f'Process took {np.round((time.mktime(end) - time.mktime(start))/60, 2)} minutes')
    print(cond)
    
    z_list = [num[0] for num in results]
    mean_IMT = [num[1] for num in results]

    df = pd.DataFrame(z_list)
    mean_df = pd.DataFrame(mean_IMT)
    
    cycdi = str(round(model.parameters['CycDi_0'].value % 1,2)).split('.')[1]

    if SET_JPARAM:
        filename_z = str(Path(SAVE_DIR, f'test_lyaps_{TREATMENT}_pso_{PARAM_NUM}_cycdi_0_{cycdi}.fea'))
        filename_IMT = str(Path(SAVE_DIR, f'test_mINT_{TREATMENT}_pso_{PARAM_NUM}_cycdi_0_{cycdi}.fea'))
    else:
        filename_z = str(Path(SAVE_DIR, f'test_lyaps_{TREATMENT}_cycdi_0_{cycdi}.fea'))
        filename_IMT = str(Path(SAVE_DIR, f'test_mINT_{TREATMENT}_cycdi_0_{cycdi}.fea'))
    
    # Handle z_list file
    if os.path.exists(filename_z):
        existing_df = feather.read_feather(filename_z)
        combined_df = pd.concat([existing_df, df], ignore_index=True)
    else:
        combined_df = df
    # print(filename_z)
    feather.write_feather(combined_df, filename_z)
    
    # Handle mean_IMT file
    if os.path.exists(filename_IMT):
        existing_mean_df = feather.read_feather(filename_IMT)
        combined_mean_df = pd.concat([existing_mean_df, mean_df], ignore_index=True)
    else:
        combined_mean_df = mean_df
    
    feather.write_feather(combined_mean_df, filename_IMT)


# if SET_JPARAM:
#     filename = str(Path(SAVE_DIR, f'lyaps_{TREATMENT}_pso_{PARAM_NUM}_cycdi_0_{cycdi}.fea')) # feather might have an issue with Path object
# else:
#     filename = str(Path(SAVE_DIR, f'lyaps_{TREATMENT}_cycdi_0_{cycdi}.fea'))
# feather.write_feather(df, filename)

# if SET_JPARAM:
#     filename = str(Path(SAVE_DIR, f'mINT_{TREATMENT}_pso_{PARAM_NUM}_cycdi_0_{cycdi}.fea'))
# else:
#     filename = str(Path(SAVE_DIR, f'mINT_{TREATMENT}_cycdi_0_{cycdi}.fea'))

# feather.write_feather(mean_df, filename)



