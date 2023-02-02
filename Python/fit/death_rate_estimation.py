# -*- coding: utf-8 -*-
"""
Created on Fri Jan 27 17:37:34 2023

@author: Grim.M
"""

###############################################################################

import sys
import os
import numpy as np
import pandas as pd
import glob
from pathlib import Path
import re
path= r'C:\Users\Grim.M\Documents\git_repos\ECM_cell_cycle_model'


sys.path.insert(0, path + "/Python/cellular_cycle")
sys.path.insert(0, path + "/Python/fit")

import fracprolif as fp
import matplotlib.pyplot as plt

import math

###############################################################################


IMT_DATADIR = Path(path + "/Data/Exp_Intermitotic_Times")

print(IMT_DATADIR)
file_location = os.path.join(IMT_DATADIR,'*imt_*.csv')
filenames = glob.glob(file_location)
filenames.sort()
#print(filenames)


IMTs = []
for f in filenames:
    tmp_df = pd.read_csv(f)
    tmp_df = tmp_df.replace(to_replace=-1, value=np.nan)
    IMTs.append(tmp_df)

ii=0
for df in IMTs:
    df["IMT"] = df["Division"] - df["Mitotic Frame"]
    df["G1"] = df["Green On"] - df["Mitotic Frame"]
    df["S-M"] = df["Division"] - df["Green On"]
    if ii == 2 or ii == 3 or ii == 6 or ii == 7:
        df["IMT_t"] = df["IMT"] * (32/60)
        df["G1_t"] = df["G1"] * (32/60)
        df["S-M_t"] = df["S-M"] * (32/60)
    else:
        df["IMT_t"] = df["IMT"] * (20/60)
        df["G1_t"] = df["G1"] * (20/60)
        df["S-M_t"] = df["S-M"] * (20/60)
    #df = df[(df["IMT_t"] ==27)].reset_index(drop=True)
    ii+=1
    
stiffs = [0.2, 0.2, 0.6, 0.6, 4, 4, 55, 55]
ii=0
for df in IMTs:
    tmp_stiff = np.ones(len(df))*stiffs[ii]
    df.insert(0, "Stiffness", tmp_stiff)
    if ii % 2 == 0:
        df["Drug"] = "DMSO"
    else:
        df["Drug"] = "ERL"    
    ii+=1

all_IMTs = pd.concat(IMTs, ignore_index=True, sort=False)
centers = np.arange(0.5, 100.5, 1)  # histogram bin centers 


all_IMTs = all_IMTs[all_IMTs['IMT_t'].notna()]

###############################################################################

DATADIR = Path(path + "/Data/Exp_Cell_Counts")

filepaths = os.listdir(DATADIR)
filenames = [f for f in filepaths if re.match(r'^DIP_DS8', f)]
filepaths = [Path(DATADIR, f) for f in filenames]
filepaths.sort()
filepaths

all_response = []
for f in filepaths:
    all_response.append(np.loadtxt(f,delimiter=','))

means=[]
for vals in range(len(all_response)):
    means.append(np.mean(all_response[vals],axis=1))
    
stdevs=[]
for vals in range(len(all_response)):
    stdevs.append(np.std(all_response[vals],axis=1))
  


run_times = [140, 140, 140, 140, 121, 121] 
conditions = ['0.2 kPa DMSO', '0.2 kPa ERL', '4 kPa DMSO', '4 kPa ERL', '55 kPa DMSO', '55 kPa ERL']
growth_rates = {}
growth_rates_end = {}
growth_rates_start = {}
growth_rates_intercept = {}
for mean, condition, run_time in zip(means, conditions, run_times):
    t = np.arange(0, run_time, 1/3)
    print(len(mean))
    print(len(t))
    
    growth_rate_intercept = LinearRegression()
    growth_rate_intercept.fit(t[:,np.newaxis], mean)

    growth_rate = LinearRegression(fit_intercept=False)
    growth_rate.fit(t[:,np.newaxis], mean)

    growth_rate_end = LinearRegression(fit_intercept=False)
    growth_rate_end.fit(t[int(len(t)/2):,np.newaxis], mean[int(len(t)/2):]-mean[int(len(t)/2)])

    growth_rate_start = LinearRegression(fit_intercept=False)
    growth_rate_start.fit(t[:int(len(t)/2),np.newaxis], mean[:int(len(t)/2)])


#    growth_rate, _, _, _ = np.linalg.lstsq(t[:,np.newaxis], mean[:])
#    growth_rate_start, _, _, _ = np.linalg.lstsq(t[:int(len(t)/2),np.newaxis], mean[:int(len(t)/2)])
#    growth_rate_end, _, _, _ = np.linalg.lstsq(t[int(len(t)/2):,np.newaxis], mean[int(len(t)/2):]-mean[int(len(t)/2)])

    plt.plot(t, mean, 'bo')
    plt.plot(t, growth_rate.coef_[0]*t, 'r-')
    plt.show()
    growth_rates_intercept[condition]=growth_rate_intercept.coef_[0]
    growth_rates[condition]=growth_rate.coef_[0]
    growth_rates_end[condition]=growth_rate_end.coef_[0]
    growth_rates_start[condition]=growth_rate_start.coef_[0]


###############################################################################

estimate_death_rate = {}
estimate_death_rate_start = {}
estimate_death_rate_end = {}
estimate_death_rate_intercept = {}
for mean, condition, run_time in zip(means, conditions, run_times):
    t = np.arange(0, run_time, 1/3)
    stiffness_tmp = float(condition.split(" ")[0])
    treatment_tmp = condition.split(" ")[2]
    
    #
    mitotic_lifespans = all_IMTs[(all_IMTs['Stiffness'] == stiffness_tmp) & (all_IMTs['Drug'] == treatment_tmp)]['IMT_t'].to_numpy()
    mitotic_lifespans = mitotic_lifespans[mitotic_lifespans>0]
    m = fp.mle(mitotic_lifespans, [])
    rates_tmp = fp.q_rates(m)
    
    #plot IMTs distribution fits
    demg_fit = [fp.demg(x, mu=m['x'][0], sigma=m['x'][1], lambd=m['x'][2]) for x in np.linspace(0, int(max(mitotic_lifespans))+10, 200)]
    plt.plot(np.linspace(0, int(max(mitotic_lifespans))+10, 200), demg_fit)
    plt.hist(mitotic_lifespans, bins=30, density=True)
    plt.show()
    
    #Estimate death rate  # d - (d-a) = a
    estimate_death_rate_intercept[condition] = rates_tmp['d']/math.log(2) - growth_rates_intercept[condition] # d - (d-a) = a
    estimate_death_rate[condition] = rates_tmp['d']/math.log(2) - growth_rates[condition] # d - (d-a) = a
    estimate_death_rate_start[condition] = rates_tmp['d']/math.log(2) - growth_rates_start[condition]
    estimate_death_rate_end[condition] = rates_tmp['d']/math.log(2) - growth_rates_end[condition]
    
    #Plot population doublings with doubling rate from IMT fit
    plt.plot(t, mean, 'bo')
    plt.plot(t, rates_tmp['d']/math.log(2)*t, 'r-')
    plt.show()
#
