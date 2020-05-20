# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 15:25:14 2020

@author: Jérémie Marlhens
Title: IMTs_simulation
"""

from runStochCell import runStochcell

from pysb.simulator import BngSimulator

from peakfinder import peakfinder

    


def IMTs(tmp_t,model):
    #Import the model
    sim = BngSimulator(model, tspan=tmp_t, verbose=True)

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
    
    simres=sim.run()         
    
    yout = simres.observables
            
    a, max_tmp_cdc=peakfinder(yout['obsCdc20'],0.15*200,0.15*200);      #cdc20
    a, max_tmp_cycE=peakfinder((yout['obsCycE']),0.15*200,0.15*200);  #cycE
    
    IMT_tmp=[];
    g1_tmp =[];
    s_m_tmp = [];
    
    if len(max_tmp_cdc)>0 and len(max_tmp_cycE)>0 :          
        div_ind = max_tmp_cdc;
        if len(div_ind)>=2 :
            for jj in range(1,len(div_ind)) :
                IMT_tmp.append(tmp_t[div_ind[jj]]-tmp_t[div_ind[jj-1]])
       
    return IMT_tmp
    
    
    