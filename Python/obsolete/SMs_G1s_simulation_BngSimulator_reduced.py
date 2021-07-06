# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 15:25:14 2020

@author: Jérémie Marlhens
Title: IMTs_simulation
"""

from pysb.simulator import BngSimulator

from pysb.simulator import ScipyOdeSimulator, StochKitSimulator

import numpy as np

import math

import os
os.environ["PYOPENCL_CTX"] ="0"

path="/home/jeremie"
#path=r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code"
import sys
sys.path.insert(0, path+"/cellular_cycle") #add a pathway to the model
from peakfinder import peakfinder


def IMTs(tmp_t, model_population) :
    #Import the model                
                    
    sim = BngSimulator(model_population, tspan=tmp_t, verbose=False);
    traj = sim.run();
    
    result_Cdc20 = np.array(traj.observables)['obsCdc20']
    result_CycE = np.array(traj.observables)['obsCycE']
    
    a, max_tmp_cdc = peakfinder(result_Cdc20,0.15*200,0.15*200, includeEndpoints=0)
    a, max_tmp_cycE = peakfinder(result_CycE,0.15*200,0.15*200, includeEndpoints=0)
    
    if len(max_tmp_cdc)<4 or len(max_tmp_cycE)<4:
            return ([],[])
    
    IMT_tmp=[];
    g1_tmp =[];
    s_m_tmp = [];
    
    if len(max_tmp_cdc)>0 and len(max_tmp_cycE)>0 :          
        div_ind = max_tmp_cdc;
        S_ind = max_tmp_cycE;
        if len(div_ind)>=2 :
            for jj in range(1,len(div_ind)) :
                IMT_tmp.append(tmp_t[div_ind[jj]]-tmp_t[div_ind[jj-1]])
                
        for jj in range(len(max_tmp_cycE)) :
            tmp = [x for x in (max_tmp_cdc) if x > max_tmp_cycE[jj]];
            if len(tmp)>0 :
                s_m_tmp.append(tmp_t[tmp[0]] -tmp_t[max_tmp_cycE[jj]])
    
        for jj in range(len(max_tmp_cdc)) :
            tmp = [x for x in (max_tmp_cycE) if x > max_tmp_cdc[jj]];
            if len(tmp)>0 :
                g1_tmp.append(tmp_t[tmp[0]] - tmp_t[max_tmp_cdc[jj]]);
                
    return s_m_tmp, g1_tmp
