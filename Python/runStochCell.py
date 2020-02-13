# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 17:48:30 2020

@author: Grim.M
Title: runStochCell
"""

import sys
sys.path.insert(0, r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code\cellular_cycle") #add a pathway to the model
from Austin_proliferation_model_concentration import model
from peakfinder import peakfinder
from pysb.simulator import ScipyOdeSimulator
import numpy as np
import pylab as pl

init=3 


def runStochcell(init) :

    init_cell=init;
    v = 0.05;     
    sys_size=200;
    run_time=500;    
    CI=0;
    
    t = np.linspace(0, 1000, 3000)


    simres = ScipyOdeSimulator(model, t).run() 
    yout = simres.all
    
    a, maxima = peakfinder(yout['obsCdc20'],[],0.2);   #find the peak of cdc20 - last protein to peak

    initials={}
    for k in range(len(model.observables)) :
        initials[model.observables[k].name]=yout[model.observables[k].name][maxima[6]]*sys_size
        
        
    initial_divided={}
    for i in initials :
        initial_divided[i]=round(initials[i]) 
        
        
    variable_range={}
    for k in range(len(model.observables)) :
        variable_range[model.observables[k].name]=[yout[model.observables[k].name][i] for i in range(maxima[6],maxima[7])]
       
    #variable_range=np.transpose(np.array(variable_range))

    
    