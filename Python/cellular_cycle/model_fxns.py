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
