# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 15:25:14 2020

@author: Jérémie Marlhens
Title: IMTs_simulation
Extract the G1, S to M, and IMT times from cells structure
"""

import numpy as np
import math
import os
os.environ["PYOPENCL_CTX"] ="0"
from peakfinder import peakfinder


def SMs_G1s_IMTs_extraction(cells, thr = 0.4, omega=200):
    #Import the model                
    
    IMT_tmp=[];
    g1_tmp =[];
    s_m_tmp = [];
        
    for cell in range(len(cells['molecules'])) :
        
        tmp_t = cells['time'][cell]
        
        (max_cdc, max_tmp_cdc) = peakfinder(cells['molecules'][cell]['obsCdc20'], thr*omega, thr*omega, 1, includeEndpoints=1); #peakfinder parameters, can be changed
        (max_cycE, max_tmp_cycE) = peakfinder(cells['molecules'][cell]['obsCycE'], thr*omega, thr*omega, 1, includeEndpoints=1);


#        plt.plot(cells['molecules'][cell]['obsCdc20'])
#        plt.scatter(max_tmp_cdc, max_cdc)
#
#        plt.plot(cells['molecules'][cell]['obsCycE'])
#        plt.scatter(max_tmp_cycE, max_cycE)

        
        if len(max_tmp_cdc)>0 and len(max_tmp_cycE)>0 :          
            div_ind = max_tmp_cdc;
            # S_ind = max_tmp_cycE;
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
        
    return IMT_tmp, s_m_tmp, g1_tmp
