# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 11:52:01 2020

@author: Jérémie Marlhens
Title: Austin_proliferation_model_simulation
"""

import sys
sys.path.insert(0, r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code\cellular_cycle") #add a pathway to the model
from Austin_proliferation_model_concentration import model

from pysb.simulator import ScipyOdeSimulator
from pysb.simulator import BngSimulator

import pylab as pl

from pysb.export import export

'''
####Export the model to matlab
matlab_output = export(model, 'matlab')

with open('output_matlab.m', 'w') as f:
    f.write(matlab_output)

'''

####Deterministic simulation
t = pl.linspace(0, 1000, 3000)

m=ScipyOdeSimulator(model, t)

simres = ScipyOdeSimulator(model, t).run() 
yout = simres.all

observables=['obspRb','obspRbp','obsCycE','obsCycD']

for obs in range(4) :    
    pl.plot(t,yout[observables[obs]], label=observables[obs], lw=2)
pl.legend(loc=0)
pl.show()
    
    
    
####Stochastic simulation
tmp_t = pl.linspace(0, 5000, 100000)

sim = BngSimulator(model, tspan=t, verbose=True)
simres = sim.run()
yout = simres.observables

observables=['obspRb','obspRbp','obsCycE','obsCycD','obsCdc20']

for obs in range(5) :    
    pl.plot(tmp_t,yout[observables[obs]], label=observables[obs], lw=2)
pl.legend(loc=0)
pl.show()
    

simres2= sim.run(tspan = t[-1] + t, initials = simres.species[-1])
yout2 = simres2.observables

for obs in model.observables:    
    pl.plot(t[-1] + t,yout2[obs.name], label=obs.name, lw=2)
pl.legend(loc=0)
pl.show()
