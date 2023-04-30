# import the pysb module and all its methods and functions
from pysb import *
import numpy as np
from pysb.simulator import OpenCLSSASimulator, ScipyOdeSimulator
from pysb.simulator import BngSimulator
PYOPENCL_CTX='0'

# instantiate a model
Model()

# declare monomers
Monomer('C8', ['b'])
Monomer('Bid', ['b', 'S'], {'S': ['u', 't']})

# input the parameter values
Parameter('kf', 1.0e-07)
Parameter('kr', 1.0e-03)
Parameter('kc', 1.0)


Observable('obsC8', C8(b=None))
Observable('obsBid', Bid(b=None, S='u'))
Observable('obstBid', Bid(b=None, S='t'))


#############################################################
Expression('kg', kf/(kf+kc)) #doesn't work,around 1.0e-07
Expression('kg', (kf/kc)) #doesn't work = 1.0e-07
Expression('kg', kr/(obsBid+kc)) #works
Expression('kg', obstBid/(obsBid+kc))   #works
Expression('kg', (kr*kc)) #works
Expression('kg', (kr+kf)) #works
Expression('kg', (kc)) #works
Expression('kg', (kr)) #works
Expression('kg', (kf)) #works 1.0e-07
#############################################################

# now input the rules
Rule('C8_Bid_bind', C8(b=None) + Bid(b=None, S='u') | C8(b=1) % Bid(b=1, S='u'), kf, kr)
Rule('tBid_from_C8Bid', C8(b=1) % Bid(b=1, S='u') >> C8(b=None) + Bid(b=None, S='t'), kg) #expression is in this rule

model.monomers
model.parameters
model.rules


Parameter('C8_0', 1000)
Parameter('Bid_0', 10000)
Initial(C8(b=None), C8_0)
Initial(Bid(b=None, S='u'), Bid_0)



#Check if the simulation ODE works 
t = np.linspace(0, 200)
simres = ScipyOdeSimulator(model, tspan=t).run()
yout = simres.all

#Check bng
sim = BngSimulator(model, tspan=t, verbose=True)
simres = sim.run(method='ssa', verbose=True)
yout = simres.observables

#And check openCLSSASimulator
n_sim=4
sim = OpenCLSSASimulator(model, tspan=t, verbose=True)
traj = sim.run(tspan=t, number_sim=n_sim)
