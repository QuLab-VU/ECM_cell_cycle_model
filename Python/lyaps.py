import os
import sys
import numpy as np
import pandas as pd
from pathlib import Path
from numpy.linalg import norm
import datetime
import math
# import multiprocess as mp
from multiprocessing import Pool
import lyaps
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
TREATMENT = 'DMSO'  # set AKT, ERK, and FAK by setting 'ERL' or 'DMSO'
CYCDI_CONC = 0
POINTS = 5  # number to divide each condition into, sets the number of points to divide cell cycle
PROCESSES = 5  # number of pool processes to run at once
OVERWRITE = True  # Overwrite existing simulation files
VERBOSE = True
SET_JPARAM = False
PARAM_NUM = 2
integrator_options = {
    "rtol": 1e-12,
    "atol": 1e-12
}

def run_lyaps(AKT, ERK, FAK, np=np, model=model, points=POINTS, ScipyOdeSimulator=ScipyOdeSimulator,
              find_peaks=find_peaks, integrator_options=integrator_options, norm=norm, math=math,
              verbose=VERBOSE):
    model.parameters['AKT'].value = AKT
    model.parameters['ERK'].value = ERK
    model.parameters['FAK'].value = FAK
    ep = 0.001 * 200
    t = np.linspace(0, 1000, 1000001)
    t_0 = np.linspace(0, 2000, 2000001)
    all_logz = dict.fromkeys(list(range(1, points + 1)))

    for kk in range(points):
        sim_ode = ScipyOdeSimulator(model, t_0, verbose=False, compiler='python',
                                    integrator_options=integrator_options).run()
        y_out = sim_ode.observables

        peaks_cdc, _ = find_peaks(y_out['obsCdc20'], prominence=0.6, distance=10)
        # new_y = y_out[peaks_cdc[-11]:peaks_cdc[-2]]
        peaks_half = peaks_cdc[peaks_cdc < len(y_out)/2]
        peaks_half2 = peaks_cdc[peaks_cdc > len(y_out) / 2]

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
        if verbose:
            print([AKT, ERK, FAK, kk])

    return all_logz, meanIMT

def run_lyap_one_input(np=np, model=model, points=POINTS, ScipyOdeSimulator=ScipyOdeSimulator,
              find_peaks=find_peaks, integrator_options=integrator_options, norm=norm, math=math,
              verbose=VERBOSE):

    for kk in range(points):
        sim_ode = ScipyOdeSimulator(model, t_0, verbose=False, compiler='python',
                                    integrator_options=integrator_options).run()
        y_out = sim_ode.observables

        peaks_cdc, _ = find_peaks(y_out['obsCdc20'], prominence=0.6, distance=10)
        # new_y = y_out[peaks_cdc[-11]:peaks_cdc[-2]]
        peaks_half = peaks_cdc[peaks_cdc < len(y_out)/2]
        peaks_half2 = peaks_cdc[peaks_cdc > len(y_out) / 2]

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
        if verbose:
            print([AKT, ERK, FAK, kk])

    return all_logz, meanIMT
