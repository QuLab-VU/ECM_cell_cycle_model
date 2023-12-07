import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import random
from pathlib import Path
from numpy.linalg import norm
import datetime
import math
import multiprocess as mp
import time
from pysb.simulator import BngSimulator
from pysb.simulator import ScipyOdeSimulator
from cell_cycle_death_model import model
from scipy.signal import find_peaks
import pickle

TREATMENT = 'ERL'  # set AKT, ERK, and FAK by setting 'ERL' or 'DMSO'
POINTS = 1  # number to divide each condition into, sets the number of points to divide cell cycle
PROCESSES = 12  # number of pool processes to run at once
integrator_options = {
    "rtol": 1e-12,
    "atol": 1e-12
}


def run_lyaps(AKT, ERK, FAK, np=np, model=model, points=POINTS, ScipyOdeSimulator=ScipyOdeSimulator,
              find_peaks=find_peaks,
              integrator_options=integrator_options, norm=norm, math=math):
    model.parameters['AKT'].value = AKT
    model.parameters['ERK'].value = ERK
    model.parameters['FAK'].value = FAK
    ep = 0.001 * 200
    t = np.linspace(0, 1000, 1000001)
    all_logz = dict.fromkeys(list(range(1, points + 1)))

    for kk in range(points):
        sim_ode = ScipyOdeSimulator(model, t, verbose=False, compiler='python',
                                    integrator_options=integrator_options).run()
        y_out = sim_ode.all

        peaks_cdc, _ = find_peaks(y_out['obsCdc20'], prominence=0.6, distance=10)
        new_y = y_out[peaks_cdc[-11]:peaks_cdc[-2]]

        variable_range = y_out[peaks_cdc[-11]:peaks_cdc[-10]]
        tmp_loc = (round(len(variable_range) / (kk + 1)) * (kk))

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
            dists = [0] * len(variable_range)
            for j in range(0, len(variable_range)):
                dif = [variable_range[j][num] / 200 - y_out1[i][num] / 200 for num in range(19)]
                dists[j] = norm(dif)

            z[n] = min(dists)
            n = n + 1

        logz = [math.log(num) for num in z]
        all_logz[kk + 1] = logz

    return all_logz


sys.path.insert(0, "../../Python/cellular_cycle") # add path to the model


# Where to save large simulation output files
# SAVE_DIR = Path(Path.home(),'ECM_mod_output')
# if not SAVE_DIR.exists():
#     SAVE_DIR.mkdir()

OVERWRITE = False # Overwrite existing simulation files

CURDATE = datetime.datetime.now().strftime('%Y%m%d')

t = np.linspace(0, 1000, 1000001)
# t = np.linspace(0, 20000, 200001)
# t = np.linspace(0, 50000, 500001)

model.parameters['AKT'].value = 1
model.parameters['ERK'].value = 1
model.parameters['FAK'].value = 1

# Importing the parameters:
new_param = pd.read_csv(Path(os.getcwd()+ '/Python/model_data/parameters.csv'))

for key in new_param.keys():
    try:
        if model.parameters[key].value != new_param[key][0]:
            model.parameters[key].value = new_param[key]
    except:
        print(key)
model.parameters['omega'].value = 200

initial_concentrations = pd.read_csv(os.getcwd() + '/Python/model_data/initial_concentrations.csv')
for key in initial_concentrations.keys():
    try:
        model.parameters[key].value = initial_concentrations[key][0] / model.parameters['omega'].value
    except:
        print(key)

if TREATMENT == 'ERL':
    FAKs = list(np.linspace(0.4,0.6,50))
    tmp = list(np.linspace(0.6,1,50))
    FAKs.extend(tmp)
    AKTs = list(np.linspace(0.3,0.6,50))
    tmp = [0.6]*50
    AKTs.extend(tmp)
    ERKs = [0.2]*100
else:
    FAKs = list(np.linspace(0.4,0.6,50))
    tmp = list(np.linspace(0.6,1,50))
    FAKs.extend(tmp)
    AKTs = list(np.linspace(0.6,1,50))
    tmp = [1]*50
    AKTs.extend(tmp)
    ERKs = [1]*100

conds = [(AKT, ERK, FAK) for AKT, ERK, FAK, in zip(AKTs, ERKs, FAKs)]

start = time.time()
with mp.Pool(processes = PROCESSES) as p:
    results = p.starmap(run_lyaps, [(0.3,0.2,0.4),(0.6,0.2,0.6)])
end = time.time()
print(end - start)

with open('lyaps.pkl', 'wb') as f:
    pickle.dump(results, f)
