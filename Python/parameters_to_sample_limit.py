# -*- coding: utf-8 -*-
"""
Created on Fri Apr 24 12:07:36 2020

@author: Grim.M
"""

t = np.linspace(0, 200, 400)


simres = ScipyOdeSimulator(model, t).run() 
yout = simres.all
    
a, maxima = peakfinder(yout['obsCdc20'],[],0.2);   #find the peak of cdc20 - last protein to peak

if len(maxima)>= 2 :
    parameter_space+=[parameter_vector]













