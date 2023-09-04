from pysb import *
from pysb.simulator import ScipyOdeSimulator
import numpy as np
# Minimal model giving errors, 

Model()
Monomer('A')
Parameter('k_s', 1)
Parameter('k_d', 1)
Rule('synthesize_and_degrade', None | A(), k_s, k_d)

t = np.linspace(0, 10, 10)
sim = ScipyOdeSimulator(model, tspan=t)
simres = sim.run()
