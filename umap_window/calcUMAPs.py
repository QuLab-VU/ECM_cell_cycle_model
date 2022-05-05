# -*- coding: utf-8 -*-
"""
Created on Wed May  4 17:52:00 2022

@author: Austin
"""

import pandas as pd
import numpy as np
import umap
import glob
import os
import joblib

#The filenames to make the umaps from
filenames = ['erl_4-9_1_1.txt', 'erl_4-9_1_2.txt', 'erl_4-9_1_3.txt',
             'erl_4-9_2_1.txt', 'erl_4-9_2_2.txt', 'erl_4-9_2_3.txt',
             'erl_4-9_3_1.txt', 'erl_4-9_3_2.txt', 'erl_4-9_3_3.txt',
             'dmso_4-9_1_1.txt', 'dmso_4-9_1_2.txt', 'dmso_4-9_1_3.txt',
             'dmso_4-9_2_1.txt', 'dmso_4-9_2_2.txt', 'dmso_4-9_2_3.txt',
             'dmso_4-9_3_1.txt', 'dmso_4-9_3_2.txt', 'dmso_4-9_3_3.txt']

#if you have a lot of data for umaps
"""
file_location = os.path.join('C:\\Users\Austin\Documents\\UMAP_projects\\nsclc_model\\22_02_14_runs\\', '*.txt')
filenames = glob.glob(file_location)
"""

#Important params to set the window size and step time distance
init_dt = 0.01                #initial dt from the raw model output
dt = 0.1                      #dt for umap window
window = 5                    #window size in hours
time_step = 1                 #time step between windows

#Load umap reducer object
reducer = joblib.load('reducer_05_02.sav')

#Load and resample time based on umap window
time_tmp = pd.read_csv('t_4-9_3_3.txt', sep=",", header=None)
time_tmp = time_tmp.T
time_tmp.columns = ["time"]
time_tmp2 = time_tmp.iloc[::int(dt//init_dt)].reset_index().drop(columns=['index'])

#create window and total protein arrays
win_len = int(window/dt)
range_ind = np.arange(int(win_len/2),int(len(time_tmp2)-(int(win_len/2))),int(time_step/dt))


for f in filenames:
    prot_tmp = pd.read_csv(f, sep=",", header=None, usecols = range(0,18))
    prot_tmp.columns = ["skp2", "p21", "pp21", "cdh1", "myc", "cycD", "cycD-p21", 
    "Rb", "pRb", "E2F", "E2F-Rb", "cycE", "cycE-p21", "cycA", "cycA-p21",
    "cycB", "cycB-p21", "cdc20"]
    prot_tmp2 = prot_tmp.iloc[::int(dt//init_dt),:].reset_index().drop(columns=['index'])
    
    prot_re_tmp = np.empty((len(range_ind),win_len*18))
    time_re_tmp = np.empty((len(range_ind),win_len))
    time_at_tmp = np.empty((len(range_ind),1))

    ii=0
    for n in range(int(win_len/2),int(len(time_tmp2)-(int(win_len/2))),int(time_step/dt)):
        tmp = prot_tmp2.iloc[n-int(win_len/2):n+int(win_len/2)]
        tmp_t = time_tmp2['time'].iloc[n-int(win_len/2):n+int(win_len/2)]
        tmp=tmp.to_numpy()
        tmp=np.reshape(np.transpose(tmp),win_len*18)
        tmp_at = time_tmp2['time'].iloc[n]
        
        prot_re_tmp[ii,:] = np.array([tmp])
        time_re_tmp[ii,:] = np.array([tmp_t])
        time_at_tmp[ii,:] = np.array([[tmp_at]])
        ii+=1
    
    embedding_all_prot = reducer.transform(prot_re_tmp)
    split_s = f.split('.')
    name = 'embedding_erl_'+split_s[0]+'.csv'
    np.savetxt('embedding_'+split_s[0]+'.csv', embedding_all_prot, delimiter=',', fmt='%f')
    np.savetxt('all_time_'+split_s[0]+'.csv', time_re_tmp, delimiter=',', fmt='%f')
    np.savetxt('all_time_at_'+split_s[0]+'.csv', time_at_tmp, delimiter=',', fmt='%f')
    np.savetxt('all_prot_'+split_s[0]+'.csv', prot_re_tmp, delimiter=',', fmt='%f')

