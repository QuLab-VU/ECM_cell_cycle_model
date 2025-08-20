# ------------------ set BLAS threads BEFORE numpy/scipy imports ------------------
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")
# -------------------------------------------------------------------------------

import sys
import time
import math
import pickle
import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from numpy.linalg import norm
from scipy.signal import find_peaks
import pyarrow.feather as feather
import multiprocessing as mp  # stdlib is fine once we avoid heavy pickling

# --------------------------- user-configurable ----------------------------------
PARAM = pd.read_csv(os.path.join(os.getcwd(), 'model_data', 'parameters.csv'))
INIT_CONCS = pd.read_csv(os.path.join(os.getcwd(), 'model_data', 'initial_concentrations.csv'))

TREATMENT   = 'ERL'   # 'ERL', 'DMSO', 'LOW', or 'TEST'
CYCDI_CONC  = 0
POINTS      = 10
PROCESSES   = 5  # sensible default
OVERWRITE   = True
VERBOSE     = True
SET_JPARAM  = False
PARAM_NUM   = 0

integrator_options = {"rtol": 1e-12, "atol": 1e-12}

SAVE_DIR = Path(Path.home(), 'ECM_mod_output')
SAVE_DIR.mkdir(exist_ok=True)

# ------------------------- per-process globals ----------------------------------
_GLOBALS = {}

def _init_worker(param_csv_path, init_concs_csv_path, set_jparam, param_num, cycdi_conc, integrator_opts):
    """Runs once per worker process: import heavy modules & model, apply params, precompute arrays."""
    sys.path.insert(0, "./cellular_cycle")
    from pysb.simulator import ScipyOdeSimulator
    from cell_cycle_death_model import model as _model

    _GLOBALS['model'] = _model
    _GLOBALS['ScipyOdeSimulator'] = ScipyOdeSimulator
    _GLOBALS['find_peaks'] = find_peaks
    _GLOBALS['integrator_options'] = integrator_opts

    # precompute large time grids once per worker
    _GLOBALS['t']   = np.linspace(0, 1000, 1_000_001)
    _GLOBALS['t_0'] = np.linspace(0, 2000, 2_000_001)

    # apply params from CSVs
    PARAM_local = pd.read_csv(param_csv_path)
    INIT_local  = pd.read_csv(init_concs_csv_path)

    # baseline values (overwritten per task)
    _model.parameters['AKT'].value = 1
    _model.parameters['ERK'].value = 1
    _model.parameters['FAK'].value = 1

    # assign scalar values from CSVs
    for key in PARAM_local.keys():
        try:
            val = PARAM_local[key][0]
            if _model.parameters[key].value != val:
                _model.parameters[key].value = val
        except Exception:
            pass

    _model.parameters['omega'].value = 200

    if set_jparam:
        with open('fit_data/pso/frac_prolif_parameters_pso.pkl', 'rb') as f:
            frac_prolif = pickle.load(f)
        frac_prolif[('cost', 'custom cost populations')] = (
            frac_prolif[('cost', 'cost population dynamics 0.6 kPa DMSO')]
            + frac_prolif[('cost', 'cost population dynamics 4 kPa erlotinib')]
            + frac_prolif[('cost', 'cost population dynamics 55 kPa erlotinib')]
        )
        pso_sort = frac_prolif.sort_values(by=[('cost', 'custom cost populations')])
        new_params = pso_sort['model parameters']
        for key in new_params.keys():
            try:
                val = new_params[key].iloc[param_num]
                if _model.parameters[key].value != val:
                    _model.parameters[key].value = val
            except Exception:
                pass
        _model.parameters['omega'].value = 200

    # initial concentrations
    for key in INIT_local.keys():
        try:
            _model.parameters[key].value = INIT_local[key][0] / _model.parameters['omega'].value
        except Exception:
            pass

    _model.parameters['CycDi_0'].value = cycdi_conc


def run_lyaps(AKT, ERK, FAK, points=POINTS, verbose=VERBOSE):
    """Worker task for a single (AKT, ERK, FAK)."""
    model = _GLOBALS['model']
    ScipyOdeSimulator = _GLOBALS['ScipyOdeSimulator']
    find_peaks = _GLOBALS['find_peaks']
    integrator_options = _GLOBALS['integrator_options']
    t = _GLOBALS['t']
    t_0 = _GLOBALS['t_0']

    # set per-condition parameters
    model.parameters['AKT'].value = AKT
    model.parameters['ERK'].value = ERK
    model.parameters['FAK'].value = FAK

    ep = 0.001 * 200
    all_logz = dict.fromkeys(range(1, points + 1))

    # long pre-run to find cycle window
    sim_ode = ScipyOdeSimulator(model, t_0, verbose=False, compiler='python',
                                integrator_options=integrator_options).run()
    y_out = sim_ode.observables  # structured array; y_out['obsCdc20'] works

    peaks_cdc, _ = find_peaks(y_out['obsCdc20'], prominence=0.6, distance=10)
    peaks_half  = peaks_cdc[peaks_cdc < len(y_out) // 2]
    peaks_half2 = peaks_cdc[peaks_cdc > len(y_out) // 2]

    if len(peaks_half2) > 3:
        IMTs = [t_0[y] - t_0[x] for x, y in zip(peaks_half2, peaks_half2[1:])]
        meanIMT = sum(IMTs) / len(IMTs)
        variable_range = y_out[peaks_half[-3]:peaks_half[-2]]
    else:
        meanIMT = 0
        variable_range = y_out[575000:600000]

    slice_len = max(1, round(len(variable_range) / points))
    prots = [
        'Skp2_0','p21_0','p21p_0','cdh1act_0','cdh1inact_0','Myc_0','CycD_0','CycDp21_0','pRb_0',
        'pRbp_0','E2F_0','pRbE2F_0','CycE_0','CycEp21_0','CycA_0','CycAp21_0','CycB_0','CycBp21_0','Cdc20_0'
    ]

    # run POINTS short sims from different initial points along the cycle window
    for kk in range(points):
        tmp_loc = min(len(variable_range)-1, slice_len * kk)

        Initial_1 = variable_range[tmp_loc]
        Initial_1 = [val + ep for val in Initial_1]

        prot_init = {prots[i]: Initial_1[i] for i in range(len(prots))}
        for key in prot_init.keys():
            model.parameters[key].value = prot_init[key] / model.parameters['omega'].value

        sim_ode1 = ScipyOdeSimulator(model, t, verbose=False, compiler='python',
                                     integrator_options=integrator_options).run()
        y_out1 = sim_ode1.observables

        # compute z (min distance to moving window) at coarse steps
        step = 2000
        z = [0] * (len(y_out1) // step + 1)
        n = 0
        for i in range(0, len(y_out1), step):
            if len(peaks_half2) > 0:
                start = peaks_half[-3] + i + tmp_loc - 5000
                stop  = peaks_half[-2] + i + tmp_loc + 5000
            else:
                start = 575000 + i + tmp_loc - 5000
                stop  = 600000 + i + tmp_loc + 5000

            # slice is auto-clipped by numpy even if start<0 or stop>len; safe
            window = y_out[max(0, start):min(len(y_out), stop)]

            # fallback when window is empty (edge case)
            if len(window) == 0:
                z[n] = 1e-12  # tiny positive to avoid log(0)
                n += 1
                continue

            # compute norms to current point (first 19 dims), avoid .tolist() churn
            target = np.asarray(y_out1[i].tolist())[:19] / 200.0
            # loop for safety with structured array
            best = float('inf')
            for row in window:
                dif = np.asarray(row.tolist())[:19] / 200.0 - target
                d = norm(dif)
                if d < best:
                    best = d
            z[n] = max(best, 1e-300)  # avoid log(0)
            n += 1

        logz = [math.log(num) for num in z]
        all_logz[kk + 1] = logz

        if verbose:
            print([AKT, ERK, FAK, kk, meanIMT])

    return all_logz, meanIMT


def build_conditions(treatment):
    if treatment == 'ERL':
        FAKs = list(np.linspace(0,0.4,51)); tmp = list(np.linspace(0.4,0.6,51)); FAKs.extend(tmp[1:])
        tmp = list(np.linspace(0.6,1,51)); FAKs.extend(tmp[1:])
        AKTs = list(np.linspace(0,0.3,51)); tmp = list(np.linspace(0.3,0.6,51)); AKTs.extend(tmp[1:])
        AKTs.extend([0.6]*50)
        ERKs = [0.2]*151
    elif treatment == 'DMSO':
        FAKs = list(np.linspace(0,0.4,51)); tmp = list(np.linspace(0.4,0.6,51)); FAKs.extend(tmp[1:])
        tmp = list(np.linspace(0.6,1,51)); FAKs.extend(tmp[1:])
        AKTs = list(np.linspace(0,0.6,51)); tmp = list(np.linspace(0.6,1,51)); AKTs.extend(tmp[1:])
        AKTs.extend([1]*50)
        ERKs = [1]*151
    elif treatment == 'TEST':
        FAKs = list(np.linspace(0,1,11))
        AKTs = list(np.linspace(0,1,11))
        ERKs = [1]*11
    else:  # LOW
        FAKs = list(np.linspace(0,0.4,51)); tmp=list(np.linspace(0.4,0.6,51)); FAKs.extend(tmp[1:])
        tmp = list(np.linspace(0.6,1,51)); FAKs.extend(tmp[1:])
        AKTs = list(np.linspace(0,0.3,51)); tmp=list(np.linspace(0.3,0.6,51)); AKTs.extend(tmp[1:])
        AKTs.extend([0.6]*50)
        ERKs = [0.1]*151
    return list(zip(AKTs, ERKs, FAKs))


if __name__ == '__main__':
    # REQUIRED for Windows/macOS spawn; also safer on Linux (avoids forking MKL state)
    try:
        mp.set_start_method('spawn', force=True)
    except RuntimeError:
        pass

    conds = build_conditions(TREATMENT)

    start = time.time()
    print(f"Start: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}  |  {len(conds)} conditions")

    ctx = mp.get_context('spawn')
    with ctx.Pool(
        processes=PROCESSES,
        initializer=_init_worker,
        initargs=(
            os.path.join(os.getcwd(), 'model_data', 'parameters.csv'),
            os.path.join(os.getcwd(), 'model_data', 'initial_concentrations.csv'),
            SET_JPARAM,
            PARAM_NUM,
            CYCDI_CONC,
            integrator_options,
        ),
        maxtasksperchild=1,  # keeps long runs stable
    ) as pool:
        # small chunksize to stream results and keep children responsive
        results = list(pool.starmap(run_lyaps, conds, chunksize=1))

    print(f"End:   {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Elapsed minutes: {round((time.time() - start)/60, 2)}")

    # collect & write
    z_list   = [num[0] for num in results]
    mean_IMT = [num[1] for num in results]

    df = pd.DataFrame(z_list)
    mean_df = pd.DataFrame(mean_IMT)

    cycdi_str = (
        f"{(float(CYCDI_CONC) % 1.0):.2f}".split(".")[1]
        if np.isfinite(float(CYCDI_CONC))
        else "x"
    )

    out1 = Path(SAVE_DIR, f'lyaps_{TREATMENT}{"_pso_"+str(PARAM_NUM) if SET_JPARAM else ""}_cycdi_0_{cycdi_str}.fea')
    out2 = Path(SAVE_DIR, f'mINT_{TREATMENT}{"_pso_"+str(PARAM_NUM) if SET_JPARAM else ""}_cycdi_0_{cycdi_str}.fea')

    feather.write_feather(df, str(out1))
    feather.write_feather(mean_df, str(out2))
    print("Saved:", out1)
    print("Saved:", out2)
