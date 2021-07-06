# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 15:25:14 2020

@author: Jérémie Marlhens
Title: IMTs_simulation
"""

from pysb.simulator import BngSimulator

from pysb.simulator import OpenCLSSASimulator, ScipyOdeSimulator

import numpy as np

path="/home/jeremie"
#path=r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code"
import sys
sys.path.insert(0, path+"/cellular_cycle") #add a pathway to the model
from peakfinder import peakfinder

from runStochCell import runStochcell

def IMTs(tmp_t,model_concentration, model_population):
    #Import the model
    
    n_sim=150
    print(n_sim)
    
    print(model_concentration)

    [variable_range, initial_divided]=runStochcell(model_concentration,1)

    ini={}
    a=np.random.randint(len(variable_range['obsBCL2']))
    for var in list(variable_range.keys()) :
        ini[var] = [round(variable_range[var][a]*200)]

    Initial=[]
    for j in ini :
        Initial=Initial+[ini[j][-1]]

                
    sim = OpenCLSSASimulator(model_population, tspan=tmp_t, verbose=True)
    traj = sim.run(tspan=tmp_t, number_sim=n_sim, initials=np.array(Initial))
  
    
    '''
    #define initial molecule populations
    [variable_range, initial_divided] = runStochcell(1)
     
    Ini={}
    
    a=np.random.randint(len(variable_range))
    
    for var in list(variable_range.keys()) :
        Ini[var] = [round(variable_range[var][a]*200)]
                 
    initials=[]        
    for j in Ini :
        initials=initials+[Ini[j][-1]]
    
    simres = sim.run(initials=initials)
    '''
    
    
    result_Cdc20 = np.array(traj.observables)['obsCdc20']
    result_CycE = np.array(traj.observables)['obsCycE']
    
    max_tmp_cdc=[]
    max_tmp_cycE=[]
    
    for chain in range(n_sim):
        a, b=peakfinder(result_Cdc20[chain,:],0.15*200,0.15*200)
        max_tmp_cdc=max_tmp_cdc+[b]
        a,b=peakfinder(result_CycE[chain,:],0.15*200,0.15*200)
        max_tmp_cycE=max_tmp_cycE+[b]
            
    
    IMT_tmp=[];
    g1_tmp =[];
    s_m_tmp = [];
    
    for chain in range(n_sim):
        
        if len(max_tmp_cdc[chain])>0 and len(max_tmp_cycE[chain])>0 :          
            div_ind = max_tmp_cdc[chain];
            if len(div_ind)>=2 :
                for jj in range(1,len(div_ind)) :
                    IMT_tmp.append(tmp_t[div_ind[jj]]-tmp_t[div_ind[jj-1]])
       
    return IMT_tmp
    
    
    