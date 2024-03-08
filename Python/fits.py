# functions for fitting and analyzing cell cycle time distributions
import math
from math import log10, floor
from scipy.stats import exponnorm
import scipy.special as sse
import scipy.optimize as opt
import numpy as np
import pandas as pd
from matplotlib import rc


def calculate_aic(LL, num_params):
    aic = 2*num_params - 2*LL
    return aic


def calculate_bic(LL, num_params, samp_size):
    bic = num_params*math.log(samp_size) - 2*LL
    return bic


def round_sig(x, sig=2):
    return round(x, sig - int(floor(log10(abs(x)))) - 1)


def emg_fit_func(x, l, s, m):
    return 0.5*l*np.exp(0.5*l*(2*m+l*s*s-2*x))*sse.erfc((m+l*s*s-x)/(np.sqrt(2)*s)) # exponentially modified gaussian


def emg_fit_func_2(x, l, s, m, l_2, s_2, m_2, w):
    return w * (0.5 * l * np.exp(0.5 * l * (2 * m + l * s * s - 2 * x)) * sse.erfc(
        (m + l * s * s - x) / (np.sqrt(2) * s))) + (1 - w) * (
                0.5 * l_2 * np.exp(0.5 * l_2 * (2 * m_2 + l_2 * s_2 * s_2 - 2 * x)) * sse.erfc(
            (m_2 + l_2 * s_2 * s_2 - x) / (np.sqrt(2) * s_2)))  # exponential gaussian


def logl_emgv2(log_par, data):
    l, s, m = np.exp(log_par)
    eps = 10e-12
    return np.log(list(emg_fit_func(data, l, s, m) + eps)).sum()


def logl_emg2v2(log_par, data):
    l, s, m, l_2, s_2, m_2, w = np.exp(log_par)
    eps = 10e-9
    temg2 = list(emg_fit_func_2(data, l, s, m, l_2, s_2, m_2, w))
    emg2 = [x + eps for x in temg2]
    return np.log(emg2).sum()


def calc_aicc(AIC, num_params, samp_size):
    aicc = AIC + (2 * num_params * (num_params + 1) / (samp_size - num_params - 1))
    return aicc


def calc_aicp(AICs):
    aicp = [math.exp((min(AICs) - AIC) / 2) for AIC in AICs]
    return aicp


def emg_mle_v2(data, log_par, bnds):
    res = opt.minimize(
        fun=lambda log_params, data: -logl_emgv2(log_params, data),
        x0=np.array(log_par), args=(data,), bounds=bnds, method='L-BFGS-B')
    log_p = res.x
    params = [np.exp(y) for y in log_p]
    return log_p, params


def emg2_mle_v2(data, log_par2, bnds):
    res = opt.minimize(
        fun=lambda log_params, data: -logl_emg2v2(log_params, data),
        x0=np.array(log_par2), args=(data,), bounds=bnds, method='L-BFGS-B')
    log_p = res.x
    params = [np.exp(y) for y in log_p]
    return log_p, params


def crit(log_par, G1_times):
    data = np.array(G1_times)
    l, s, m, l_2, s_2, m_2, w = np.exp(log_par)

    LL2 = logl_emg2v2(log_par, G1_times)
    negLL2 = -LL2
    return negLL2


def smooth(y, box_pts):
    box = np.ones(box_pts) / box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth


def fit_EMGMM(dat, init, bounds, centers, cond, ax=None, loc=4):
    """
    dat:
    """
    G1_times = dat[(dat['Stiffness']==cond['stiff']) & (dat['Drug']==cond['drug'])]['G1_t']
    G1_times = G1_times[~np.isnan(G1_times)]
    data = G1_times[~np.isnan(G1_times)]

    h, edges =np.histogram(data, bins=range(0,101), density=True)
    tmp = smooth(h,7)
    out = ax.stairs(h, edges, fill=True, color="tab:red", alpha=0.2)
    params = exponnorm.fit(data)
    LL = exponnorm.logpdf(G1_times,*params).sum()
    AIC = calculate_aic(LL, 3)
    AICc = calc_aicc(AIC, 3, len(G1_times))

    popt, pcov = opt.curve_fit(emg_fit_func_2, centers, tmp, p0=init, bounds=bounds, maxfev=50000)

    tmp = [math.log(x) for x in popt]
    bounds = ((-5, 0),
                 (0, 80),
                 (0.5, 3),
                 (-5, 0),
                 (0, 80),
                 (3, 4),
                 (-10, -0.00000001))

    results = opt.minimize(crit, tmp, G1_times, method = 'Nelder-Mead', bounds=bounds, options={'maxiter':10000})


    LL2 = logl_emg2v2(results.x, G1_times)
    negLL2 = -LL2
    AIC2 = calculate_aic(LL2, 7)
    AICc2 = calc_aicc(AIC2, 7, len(G1_times))

    AICs = [AICc, AICc2]
    AICp = calc_aicp(AICs)
    AICp_r = [round_sig(AICp[0], 3), round_sig(AICp[1], 3)]


    x = np.arange(0,100,0.1)
    y = emg_fit_func_2(x,*np.exp(results.x))

    y1 = emg_fit_func(x, np.exp(results.x[0]), np.exp(results.x[1]), np.exp(results.x[2]))*np.exp(results.x[6])
    y2 = emg_fit_func(x, np.exp(results.x[3]), np.exp(results.x[4]), np.exp(results.x[5]))*(1-np.exp(results.x[6]))
    y3 = exponnorm.pdf(x, *params)
    ax.plot(x,y, color='black')
    ax.plot(x,y1, color='black', linestyle='dashed')
    ax.plot(x,y2, color='black', linestyle='dashed')
    ax.plot(x,y3, color='gray')
    h, edges =np.histogram(G1_times, bins=range(0,101), density=True)
    ax.stairs(h, edges, fill=True, color="tab:red", alpha=0.2)

    LL2 = logl_emg2v2(results.x, G1_times)
    AIC2 = calculate_aic(LL2, 7)
    AICc2 = calc_aicc(AIC2, 7, len(G1_times))

    AICs = [AICc, AICc2]
    AICp = calc_aicp(AICs)
    AICp_r = [round_sig(AICp[0], 3), round_sig(AICp[1], 3)]

    # rc('text', usetex=True) #  # must have LaTeX installed
    # ax.text(loc, 0.11, r'$\underline{AICp}$', fontsize=12)

    # # if no LaTeX, use:
    rc('text', usetex=False) #  # must have LaTeX installed
    ax.text(loc, 0.11, r'${AICp}$', fontsize=12)
    ax.text(loc, 0.11, '_____', fontsize=12)

    rc('text', usetex=False)
    ax.text(loc, 0.097, f'EMG: {AICp_r[0]}', fontsize=12, color='gray')
    ax.text(loc, 0.087, f'EMG2: {AICp_r[1]}', fontsize=12)

    out = pd.Series(data={'LL':LL, 'AIC':AIC, 'AICc':AICc, 'params':params,
                         'LL2':LL2, 'AIC2':AIC2, 'AICc2':AICc2, 'params2':results.x, 'AICp':AICp_r})

    return(out)
