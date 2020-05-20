# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 17:48:30 2020

@author: Grim.M
Title: runStochCell
"""
#path="/home/jeremie"
path=r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code"

import sys
sys.path.insert(0, path+"/cellular_cycle") #add a pathway to the model

from peakfinder import peakfinder
from pysb.simulator import ScipyOdeSimulator
import numpy as np
import pylab as pl




def runStochcell(model) :
    
    v = 0.05;     
    sys_size=200;
    run_time=500;    
    CI=0;
    
    t = np.linspace(0, 200, 400)
    
    
    simres = ScipyOdeSimulator(model, t).run() 
    yout = simres.all
    
    a, maxima = peakfinder(yout['obsCdc20'], [], 0.2, includeEndpoints=1);   #find the peak of cdc20 - last protein to peak
    
    variable_range={}
    
    if len(a)>=2:
        for k in range(len(model.observables)) :
            variable_range[model.observables[k].name]=[yout[model.observables[k].name][i] for i in range(maxima[0],maxima[1])]
    
    #variable_range=np.transpose(np.array(variable_range))
    
    return variable_range
    