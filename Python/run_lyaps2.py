import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from numpy.linalg import norm
import math
from pysb.simulator import ScipyOdeSimulator
from scipy.signal import find_peaks
sys.path.insert(0, "./cellular_cycle")
from cell_cycle_death_model import model
from pysb.util import alias_model_components
alias_model_components()
from scipy.stats import linregress
import umap
import operator
import shutil


def get_fig_nrows_ncols(ndims):
    if not isinstance(ndims, int) or ndims <= 0:
        raise ValueError("'ndims' must be a positive integer")
    r1 = round(math.sqrt(ndims))  # round() returns an int
    r2 = math.ceil(math.sqrt(ndims))  # math.ceil() also returns an int
    while r1 * r2 >= ndims:
        r1 -= 1
        r2 += 1
    ncols = min(r1 + 1, r2 - 1)  # the smaller of the two integers is the # of columns
    nrows = math.ceil(ndims / ncols)
    return nrows, ncols


def find_closest_index(tspan, target):
    best_idx = 0
    min_diff = abs(tspan[0] - target)

    for i in range(1, len(tspan)):  # Start from index 1 to avoid redundant checks
        diff = abs(tspan[i] - target)
        if diff < min_diff:
            best_idx, min_diff = i, diff
        else:
            break  # Exit early if difference starts increasing

    return best_idx


def finite_diff_cycle(arr, t):
    """
    Finite difference for one-period template (shape: (M,d)).
    Handles non-uniform time spacing and enforces periodic wrap.
    Parameters
    ----------
    arr : ndarray, shape (M,d)
        One-period trajectory samples.
    t : ndarray, shape (M,)
        Corresponding times or phases (need not be uniform).
    Returns
    -------
    darr : ndarray, shape (M,d)
        Estimated derivatives at each point.
    """
    arr = np.asarray(arr, dtype=float)
    t = np.asarray(t, dtype=float)
    M, d = arr.shape
    darr = np.zeros_like(arr)
    # interior points: np.gradient handles non-uniform t
    for j in range(d):
        darr[:, j] = np.gradient(arr[:, j], t)
    # periodic wrap: replace first/last rows with central diff across boundary
    dt0 = (t[1] - t[0]) + (t[-1] - t[-2])  # average step near edges
    darr[0] = (arr[1] - arr[-1]) / (t[1] - t[-1] + (t[-1] - t[0]))
    darr[-1] = (arr[0] - arr[-2]) / (t[0] - t[-2] + (t[-1] - t[0]))

    return darr


def perturb_orthogonal_nonneg(x, ref_cycle, ref_times, magnitude=1e-3, safety=0.9, seed=None, rng=None):
    """
    Add a perturbation to x in a direction orthogonal to the cycle tangent, with preferred absolute size 'magnitude',
    while ensuring all components remain non-negative.
    Parameters
    ----------
    x : (d,) array
        Current state (must be non-negative).
    ref_cycle : (M,d) array
        One-period template trajectory.
    ref_times : (d,) array
        Times along ref_cycle.
    magnitude : float
        Preferred absolute size of perturbation (Euclidean norm).
    safety : float
        Shrink the maximum step slightly to avoid hitting 0 exactly.
    seed: int, optional
        Random seed (for reproducibility).
    rng : np.random.Generator, optional
        Random generator for fallback if needed.
    Returns
    -------
    x_new : (d,) array
        Perturbed, nonnegative state.
    alpha  : float
        Actual step size used (≤ magnitude).
    """
    x = np.asarray(x, float).ravel()
    # nearest point on template (phase alignment)
    j = int(np.argmin(np.linalg.norm(ref_cycle - x, axis=1)))
    # get tangent vectors along ref_cycle
    ref_dcycle = finite_diff_cycle(ref_cycle, ref_times)
    v = ref_dcycle[j]
    vhat = v / (np.linalg.norm(v) + 1e-12)
    # build a transverse direction
    u = x - ref_cycle[j]
    u = u - np.dot(u, vhat) * vhat
    if np.linalg.norm(u) < 1e-12:
        if rng is None:
            rng = np.random.default_rng(seed)
        w = rng.normal(size=x.shape)
        u = w - np.dot(w, vhat) * vhat
    uhat = u / (np.linalg.norm(u) + 1e-12)
    # maximum admissible alpha so that x + alpha * uhat >= 0
    neg_dir = uhat < 0
    if np.any(neg_dir):
        alpha_max = np.min(x[neg_dir] / (-uhat[neg_dir] + 1e-18))
    else:
        alpha_max = np.inf
    # choose actual step
    if np.isfinite(alpha_max):
        alpha = min(magnitude, safety * max(alpha_max, 0.0))
    else:
        alpha = magnitude
    # apply perturbation
    x_new = x + alpha * uhat
    x_new[x_new < 0] = 0.0  # clamp tiny negatives
    return x_new, float(alpha)


def create_umap_projection(prot_conc, downsample=None, n_neighbors=50, min_dist=0.01, n_components=2, random_state=42):
    proj = umap.UMAP(n_neighbors=n_neighbors, min_dist=min_dist, n_components=n_components, random_state=random_state)
    if downsample is not None:
        keep_pts = list(np.arange(0, len(prot_conc), downsample))
        prot_conc = prot_conc[keep_pts]
    reducer = proj.fit(prot_conc)
    u = reducer.transform(prot_conc)
    fig_umap = plt.figure(constrained_layout=True)
    plt.plot(u[:, 0], u[:, 1], lw=2)
    plt.xlabel('UMAP 1')
    plt.ylabel('UMAP 2')
    return reducer, fig_umap

# cwd = os.getcwd()
# PARAM = pd.read_csv(os.path.join(cwd, 'model_data', 'parameters.csv'))  # parameter file to use
# INIT_CONCS = pd.read_csv(os.path.join(cwd, 'model_data', 'initial_concentrations.csv'))  # initial concentrations
# TREATMENT = 'DMSO'  # set AKT, ERK, and FAK by setting 'ERL', 'DMSO', or 'LOW'
# CYCDI_CONC = 0
# POINTS = 10  # number to divide each condition into, sets the number of points to divide cell cycle
# PROCESSES = 6  # number of pool processes to run at once
# OVERWRITE = False  # Overwrite existing simulation files
# VERBOSE = True
# SET_JPARAM = False
# PARAM_NUM = 0  # index of parameter set to use from sets of fit parameters
# integrator_options = {"rtol": 1e-12, "atol": 1e-12}

def run_lyaps(AKT, ERK, FAK, np=np, model=model, points=10, ScipyOdeSimulator=ScipyOdeSimulator,
              find_peaks=find_peaks, integrator_options=None, norm=norm, math=math, overwrite=True, verbose=False):

    if integrator_options is None:
        integrator_options = {"rtol": 1e-12, "atol": 1e-12}

    # Set AKT, ERK, and FAK levels
    model.parameters['AKT'].value = AKT
    model.parameters['ERK'].value = ERK
    model.parameters['FAK'].value = FAK

    # Create simulator
    t_0 = np.linspace(0, 2000, int(2e6) + 1)
    sim = ScipyOdeSimulator(model, t_0, verbose=verbose, compiler='cython', integrator_options=integrator_options)

    # Protein species
    prots = [Skp2(), p21(b=None, S='u'), p21(b=None, S='p'), cdh1(A='yes'), cdh1(A='no'), Myc(), CycD(b=None),
             CycD(b=1) % p21(b=1, S='u'), pRb(b=None, S='u'), pRb(b=None, S='p'), E2F(b=None),
             E2F(b=1) % pRb(b=1, S='u'), CycE(b=None), CycE(b=1) % p21(b=1, S='u'), CycA(b=None),
             CycA(b=1) % p21(b=1, S='u'), CycB(b=None), CycB(b=1) % p21(b=1, S='u'), Cdc20()]

    observables = model.observables[:len(prots)]

    # Constraints on Cdh1, CycA, CycB, CycE, and Cdc20
    par = model.parameters
    prots_constraints = [  # sum(obs[i]) op constraint
        ((3, 4), operator.eq, (par['cdh1act_0'].value + par['cdh1inact_0'].value) * par['omega'].value),
        ((14, 15), operator.lt, par['CycA_tot'].value * par['omega'].value),
        ((16, 17), operator.lt, par['CycB_tot'].value * par['omega'].value),
        ((12, 13), operator.lt, par['CycE_tot'].value * par['omega'].value),
        ((18,), operator.lt, par['Cdc20_tot'].value * par['omega'].value)
    ]

    # Reference simulation
    output_ref = sim.run()
    initials = output_ref.species[-1]
    y_ref = output_ref.observables
    peaks_cdc, _ = find_peaks(y_ref['obsCdc20'], prominence=0.6, distance=10)
    peaks_ref = peaks_cdc[peaks_cdc > len(y_ref) / 2]  # peaks from last half of the simulation

    #####
    # plt.figure(constrained_layout=True)
    # plt.plot(t_0, y_ref['obsCdc20'], lw=2, label='obsCdc20 (ref)')
    # plt.xlim(left=1800, right=2000)
    # plt.show()
    # quit()
    #####

    # A more general method for identifying the repeating peaks in the oscillation, i.e., when dealing with
    #  complex oscillations can't assume the cycle length is the distance between two adjacent peaks
    # window = list(np.arange(peaks_ref[-4], peaks_ref[-2] + 1))  # window of points to compare against
    # y_cycle = y_ref[peaks_ref[-3]:peaks_ref[-2] + 1]  # concentrations for one full cycle
    i_start = -2
    ref_dt = t_0[peaks_ref[i_start]] - t_0[peaks_ref[i_start-1]]
    i_peak = []
    # print('ref_dt:', ref_dt)
    for i in np.arange(i_start - 1, i_start - len(peaks_ref) + 2, -1):
        # print('%d, %d:' % (i, i-1), end=' ')
        # print('%g, %g,' % (t_0[peaks_ref[i]], t_0[peaks_ref[i-1]]), end=' ')
        # print('%g' % (t_0[peaks_ref[i]] - t_0[peaks_ref[i-1]]), end=' ')
        # plt.axvline(t_0[peaks_ref[i]], color='r', ls='--')
        if np.isclose(ref_dt, t_0[peaks_ref[i]] - t_0[peaks_ref[i-1]], atol=0.01, rtol=1e-3):
            i_peak.append(i)
        #     print('True')
        # else:
        #     print('False')
        if len(i_peak) == 2:
            break
    # plt.show()
    if len(i_peak) != 2:
        raise Exception("Couldn't find two full cycles in the reference time course.")
    y_cycle = y_ref[peaks_ref[i_peak[0]]:peaks_ref[i_start] + 1]  # concentrations for one full cycle
    window = list(np.arange(peaks_ref[i_peak[1]], peaks_ref[i_start] + 1))  # window of points to compare against

    # quit()

    # FIGURES

    # 1. UMAP projection
    print('creating UMAP projection...', end='')
    conc_for_umap = np.array([y_cycle[obs.name] for obs in observables]).T
    reducer, fig_umap = create_umap_projection(conc_for_umap, downsample=100)
    fig_umap.suptitle('UMAP Embedding of the Cell Cycle')
    ax_umap = fig_umap.gca()
    x_mid, y_mid = (10, 4)  # point that define quadrants in UMAP (for placing annotations)
    print('DONE')

    # 2. Reference Cell Cycle (red)
    fig_window = plt.figure(constrained_layout=True)
    ax_window = fig_window.gca()
    ax_window.set_title('Reference Timecourse')
    ax_window.plot(t_0, y_ref['obsCdc20'], lw=2, label='obsCdc20 (ref)')
    ax_window.plot(t_0[window], y_ref['obsCdc20'][window], 'r', lw=3)
    ax_window.set_xlabel('time')
    ax_window.set_ylabel('concentration')
    ax_window.legend(loc='best')

    # 3. Peak comparisons
    nrows, ncols = get_fig_nrows_ncols(points)
    figsize = (6.4 / 3 * (ncols + 1), 4.8 * 0.4 * nrows)
    fig_tc, axs_tc = plt.subplots(nrows, ncols, sharex='all', sharey='all', constrained_layout=True,
                                  figsize=(4 / 3 * figsize[0], 4 / 3 * figsize[1]))  # increase width to include legend
    fig_tc.suptitle('Comparisons Between Reference and Perturbed Cycles')
    fig_tc.supxlabel('time')
    fig_tc.supylabel('concentration')
    axs_tc = axs_tc.flatten()

    # 4. Calculating Lyapunov Exponents
    fig_dist, axs_dist = plt.subplots(nrows, ncols, sharex='all', sharey='all', constrained_layout=True,
                                      figsize=figsize)
    fig_dist.suptitle('Lyapunov Exponents')
    fig_dist.supxlabel('time')
    fig_dist.supylabel('ln distance')
    axs_dist = axs_dist.flatten()

    # 5. Perturbed Cdc20 timecourses
    fig_obsCdc20, axs_obsCdc20 = plt.subplots(nrows, ncols, sharex='all', sharey='all', constrained_layout=True,
                                              figsize=np.array(figsize) * 2)
    fig_obsCdc20.suptitle('Perturbed Cdc20 Timecourses', fontsize=20)
    fig_obsCdc20.supxlabel('time', fontsize=20)
    fig_obsCdc20.supylabel('concentration', fontsize=20)
    axs_obsCdc20 = axs_obsCdc20.flatten()

    # 6. Perturbed p21-p timecourses
    fig_obsp21p, axs_obsp21p = plt.subplots(nrows, ncols, sharex='all', sharey='all', constrained_layout=True,
                                            figsize=np.array(figsize) * 2)
    fig_obsp21p.suptitle('Perturbed p21-p Timecourses vs. Reference', fontsize=20)
    fig_obsp21p.supxlabel('time', fontsize=20)
    fig_obsp21p.supylabel('concentration', fontsize=20)
    axs_obsp21p = axs_obsp21p.flatten()

    # 7. Perturbation distances along the cell cycle
    t_cycle = t_0[:len(y_cycle)]
    dt = 0.01 * t_cycle[-1]  # 1% of the cell cycle length
    fig_perturb = plt.figure(constrained_layout=True)
    ax_perturb = fig_perturb.gca()
    ax_perturb.plot(t_cycle, y_cycle['obsCdc20'], lw=2, label='obsCdc20')
    ax_perturb.set_title('Perturbation Distances Along the Cell Cycle')
    ax_perturb.set_xlabel('time')
    ax_perturb.set_ylabel('concentration')
    ax_perturb.legend(loc='best')

    figures = {
        fig_umap: 'umap',
        fig_window: 'window',
        fig_tc: 'tc',
        fig_dist: 'dist',
        fig_obsCdc20: 'obsCdc20',
        fig_obsp21p: 'obsp21p',
        fig_perturb: 'perturb',
    }

    print('POINTS:', points)

    slopes = []
    for kk in range(points):
        print('========\nPoint %d\n========' % kk)

        # plot reference timecourse we are comparing to
        ref_idx = round(len(y_cycle) / points) * kk
        window_kk = window[ref_idx:ref_idx + len(y_cycle)]
        t_0_window_kk = t_0[ref_idx:ref_idx + len(y_cycle)]
        axs_tc[kk].plot(t_0_window_kk, y_ref['obsCdc20'][window_kk], '0.75', lw=6, alpha=1, label='ref')

        outside_kk = [window[:ref_idx], window[ref_idx + len(y_cycle):]]
        t_0_outside_kk = [t_0[:ref_idx], t_0[ref_idx + len(y_cycle):len(window)]]
        for t_0_out, out in zip(t_0_outside_kk, outside_kk):
            axs_tc[kk].plot(t_0_out, y_ref['obsCdc20'][out], 'k--', lw=1, label='ref')

        # Get perturbed starting point
        t_idxs = [ref_idx, find_closest_index(t_cycle, t_0[ref_idx] + dt)]
        ax_perturb.plot(t_0[t_idxs], y_cycle['obsCdc20'][t_idxs], 'r-', lw=4)
        ep = float(norm([y_cycle[obs.name][t_idxs[0]] - y_cycle[obs.name][t_idxs[1]] for obs in observables]))
        ref_cycle = np.array([y_cycle[obs.name] for obs in observables]).T
        x_ref = [y_cycle[obs.name][t_idxs[0]] for obs in observables]
        x_new, alpha = perturb_orthogonal_nonneg(x_ref, ref_cycle, t_cycle, magnitude=ep, seed=42)

        # Impose constraints on species concentrations
        print('ep: %g, alpha: %g' % (ep, alpha))
        '''
        for obs, x0, xr, xn in zip(observables, output_ref.species[0], x_ref, x_new):
            print('%s | %g | %g | %g ' % (obs.name, x0, xr, xn))
        print()
        for pair in [(3, 4), (12, 13), (14, 15), (16, 17)]:
            print('%s + %s' % (observables[pair[0]].name, observables[pair[1]].name))
            print('%g | %g | %g' % (output_ref.species[0][pair[0]] + output_ref.species[0][pair[1]],
                                    x_ref[pair[0]] + x_ref[pair[1]], x_new[pair[0]] + x_new[pair[1]]))
        '''
        print('\nCheck constraints')
        for pc in prots_constraints:
            obs_idxs = pc[0]
            op = pc[1]
            constraint = pc[2]
            obs_sum = sum([x_new[i] for i in obs_idxs])
            print('-----')
            print('observables:', obs_idxs, [observables[i].name for i in obs_idxs])
            print('vals:', [float(x_new[i]) for i in obs_idxs])
            print('sum: %g' % obs_sum)
            print('constraint: %g' % constraint)
            print('Met?', end=' ')
            if not op(obs_sum, constraint):
                print('No')
                # print('fracs:', np.array([float(x_new[i]) for i in obs_idxs]) / obs_sum)
                for i in obs_idxs:
                    x_new[i] *= constraint / obs_sum
                '''
                sum_ref = x_ref[obs_idxs[0]] + x_ref[obs_idxs[1]]
                sum_new = x_new[obs_idxs[0]] + x_new[obs_idxs[1]]
                print('%g | %g + %g = %g | %g + %g = %g (%g, %g)' %
                      (constraint, x_ref[obs_idxs[0]], x_ref[obs_idxs[1]], sum_ref,
                       x_new[obs_idxs[0]], x_new[obs_idxs[1]], sum_new,
                       x_new[obs_idxs[0]] / sum_new, x_new[obs_idxs[1]] / sum_new))
                '''
                obs_sum = sum([x_new[i] for i in obs_idxs])
                print('new vals:', [float(x_new[i]) for i in obs_idxs])
                print('new sum: %g' % obs_sum)
            else:
                print('Yes')
        print()

        # Add perturbed point to UMAP plot
        u = reducer.transform(np.array([[x] for x in x_new]).T)
        ax_umap.plot(u[:, 0], u[:, 1], marker='x', ms=8, color='r')
        dx, dy = (8, 10)
        x_offset = -dx if u[:, 0] > x_mid else dx
        y_offset = -dy if y_mid < u[:, 1] < 12 else dy if u[:, 1] < y_mid else -1
        ha = "center" if u[:, 1] < 0 else "left" if x_offset > 0 else "right"
        ax_umap.annotate(
            rf'(t: {t_0[t_idxs[0]]:.1f}, $\Delta$: {alpha:.1f})',
            xy=(u[:, 0], u[:, 1]),  # point in data coords
            xytext=(x_offset, y_offset),  # offset from the point
            textcoords="offset points",  # interpret xytext in points (screen pixels)
            ha=ha, va="bottom"  # text alignment relative to the offset
        )

        # perturbation simulation
        initials[[i for i in range(len(prots))]] = x_new
        y_out = sim.run(initials=initials).observables
        peaks_cdc, _ = find_peaks(y_out['obsCdc20'], prominence=0.6, distance=10)

        # perturbed Cdc20 timecourses
        axs_obsCdc20[kk].plot(t_0, y_out['obsCdc20'], lw=2, label='Point %d' % kk)
        axs_obsCdc20[kk].tick_params(axis='both', which='major', labelsize=18)
        axs_obsCdc20[kk].legend(loc='upper right', fontsize=18, borderpad=0, frameon=False)

        # perturbed p21-p timecourses
        axs_obsp21p[kk].plot(t_0, y_ref['obsp21p'], lw=2, label='ref')
        axs_obsp21p[kk].plot(t_0, y_out['obsp21p'], lw=2, label='perturb')
        axs_obsp21p[kk].set_title('Point %d' % kk, fontsize=20)
        axs_obsp21p[kk].tick_params(axis='both', which='major', labelsize=18)
        axs_obsp21p[kk].legend(loc='upper right', fontsize=18, borderpad=0, frameon=False, labelspacing=0)

        # loop over peaks in perturbed simulation
        times = []
        log_dists = []
        for ii, peak in enumerate(peaks_cdc[:i_peak[1] - i_start]):
            start_idx = peak + ref_idx
            axs_obsCdc20[kk].plot(t_0[start_idx], 0, '2r', ms=8)  # starting points for comparisons

            # calculate Euclidian distance
            dist = [y_out[obs.name][start_idx:start_idx + len(y_cycle)] - \
                    y_ref[obs.name][window_kk] for obs in observables]
            log_dist = math.log(norm(dist))
            times.append(t_0[start_idx])
            log_dists.append(log_dist)

            # Plot comparisons between perturbed and reference cycles for this peak
            axs_tc[kk].plot(t_0[ref_idx:ref_idx + len(y_cycle)],
                            y_out['obsCdc20'][start_idx:start_idx + len(y_cycle)], lw=1, label=f'{ii}')
            axs_tc[kk].set_title('Point %d' % kk, fontsize=12)

            # For the last peak, plot the portion of the cell cycle compared against the reference
            if ii == len(peaks_cdc[:i_peak[1] - i_start]) - 1:
                axs_obsCdc20[kk].plot(t_0[start_idx:start_idx + len(y_cycle)],
                                      y_out['obsCdc20'][start_idx:start_idx + len(y_cycle)], 'r', lw=3)

        # Set figure legend for peak comparisons plot
        if kk == (ncols - 1):
            handles, labels = axs_tc[kk].get_legend_handles_labels()
            handles = [(handles[0], handles[1])] + list(handles[3:])
            labels = list(dict.fromkeys(labels))
            n_peaks = len(peaks_cdc[:i_peak[1] - i_start])
            leg_tc = fig_tc.legend(handles=handles, labels=labels, loc='upper left', title='Peak', fontsize=10,
                                       ncols=math.ceil(n_peaks / 30), bbox_to_anchor=(0.75, 1))
            leg_tc.get_title().set_fontweight('bold')
            leg_tc.get_title().set_fontsize(8)
            fig_tc.get_layout_engine().set(rect=(0, 0, 0.75, 1))  # (left, bottom, width, height)

        # filter ln(dists)
        if i_peak[0] - i_peak[1] > 1:  # indicates complex oscillations
            tt = []
            ld = []
            for i in range(1, len(log_dists)):
                if log_dists[i] > log_dists[i - 1]:
                    tt.append(times[i - 1])
                    ld.append(log_dists[i - 1])
            times = tt
            log_dists = ld

        # fit a regression line and calculate Lyaponov exponent (slope)
        max_idx = len(log_dists)
        r2val = [0]
        while r2val[-1] < 0.95 and max_idx >= 3:
            fit_line = linregress(times[:max_idx], log_dists[:max_idx])
            r2val.append(fit_line.rvalue ** 2)
            max_idx -= 1
        if max_idx > 3:
            max_idx += 1
        else:
            max_idx = r2val.index(max(r2val))
            max_idx = len(times) - max_idx + 1
            fit_line = linregress(times[:max_idx], log_dists[:max_idx])
        slopes.append(fit_line.slope)
        times = np.array(times)

        # plot log(dists) vs. time
        axs_dist[kk].plot(times, log_dists, 'o-', lw=2)
        axs_dist[kk].plot(times[:max_idx], fit_line.slope * times[:max_idx] + fit_line.intercept, '--r', lw=2,
                          label=f'{fit_line.slope:.3g}\n({fit_line.rvalue ** 2:.3g})')
        axs_dist[kk].set_title('Point %d' % kk, fontsize=12)
        leg_dist = axs_dist[kk].legend(loc='upper right', title=r'slope (R$^2$)', frameon=False, borderpad=0)
        leg_dist.get_title().set_fontweight('bold')
        leg_dist._legend_box.sep = 2  # points between the title row and the entries (default ~5–6)

        # remove x- and y-labels for inner plots of peak comparisons and Lyapunov exponents plots
        row = kk // ncols
        col = kk % ncols
        if row < nrows - 1:
            axs_tc[kk].set_xlabel('')
            axs_dist[kk].set_xlabel('')
        if col > 0:
            axs_tc[kk].set_ylabel('')
            axs_dist[kk].set_ylabel('')

    # adjust y-axis upper lim for Cdc20 and p21-p timecourse plots
    for fig in [fig_obsCdc20, fig_obsp21p]:
        fig.canvas.draw()
        axs = fig.get_axes()
        _, top = axs[0].get_ylim()
        for kk in range(points):
            axs[kk].set_ylim(top=top * 1.1)

    # save figures
    outdir = os.path.join(os.getcwd(), 'LYAP_FIGURES', 'DMSO' if ERK == 1 else 'ERL' if ERK == 0.2 else 'LOW',
                          'AKT_%g_FAK_%g_ERK_%g' % (AKT, FAK, ERK))
    p = Path(outdir)
    if p.exists():
        if overwrite or not Path(os.path.join(outdir, 'lyap.npy')).exists():
            shutil.rmtree(p)  # removes the directory and all contents
    p.mkdir(parents=True, exist_ok=True)

    for fig in figures.keys():
        fig.savefig(os.path.join(outdir, figures[fig]))

    # save and return Lyapunov exponents for all points along the cell cycle
    slopes = np.array(slopes)
    np.save(os.path.join(outdir, 'lyap'), slopes)

    return slopes


if __name__ == '__main__':

    # Calculate Lyapunov exponents for all conditions
    for TREATMENT in ['DMSO', 'ERL', 'LOW']:

        # FAK: 0 -> 0.4 -> 0.6 -> 1 for all conditions
        FAKs = np.concatenate((
            np.linspace(0, 0.4, 51),
            np.linspace(0.4, 0.6, 51)[1:],
            np.linspace(0.6, 1, 51)[1:]
        ))

        # AKT: 0 -> 0.6 -> 1 -> 1 for 'DMSO', 0 -> 0.3 -> 0.6 -> 0.6 otherwise
        AKTs = \
            np.concatenate((
                np.linspace(0, 0.6, 51),
                np.linspace(0.6, 1, 51)[1:],
                [1] * 50)) if TREATMENT == 'DMSO' else \
                np.concatenate((
                    np.linspace(0, 0.3, 51),
                    np.linspace(0.3, 0.6, 51)[1:],
                    [0.6] * 50))

        # ERK = 1 for 'DMSO', 0.2 for 'ERL', 0.1 for 'LOW'
        ERKs = [1] * 151 if TREATMENT == 'DMSO' else [0.2] * 151 if TREATMENT == 'ERL' else [0.1] * 151

        # loop over all conditions
        print('=== %s ===' % TREATMENT)
        for AKT, ERK, FAK in zip(AKTs, ERKs, FAKs):
            # [(0.792, 1, 0.496)]:  [(0.0, 1, 0.0)]:  [(0.736, 1, 0.468)]:  zip(AKTs, ERKs, FAKs):
            # [(0.744, 1, 0.472)]:
            print('AKT: %g, ERK: %g, FAK: %g' % (AKT, ERK, FAK))
            lyap_exps = run_lyaps(AKT, ERK, FAK, points=9, overwrite=False)
            print('Lyapunov exponents:', lyap_exps)
            print('Mean value: %g\n' % np.mean(lyap_exps))

    plt.show()
