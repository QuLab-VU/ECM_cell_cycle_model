# -*- coding: utf-8 -*-
"""
Created on Tue Nov 22 00:35:20 2022

@author: olesk
"""

model_population.parameters['ERK'].value = 0.2 * omega
model_population.parameters['FAK'].value = 0.4 * omega
model_population.parameters['AKT'].value = 0.4 * omega

model_concentration.parameters['ERK'].value = 0.2 * 1
model_concentration.parameters['FAK'].value = 0.4 * 1
model_concentration.parameters['AKT'].value = 0.4 * 1

t = np.linspace(0, 500, 50001)
simres = ScipyOdeSimulator(model_concentration, t, verbose=False).run()
yout = simres.all;
observables=['obsCycE', 'obsCdc20', 'obsBAXm']; #plot Cyclin E, Cdc20 and Cyclin A
for obs in range(3) :
    plt.plot(t[:],yout[observables[obs]][:], label=observables[obs], lw=2)
plt.legend(loc=1)
plt.ylim(0, 1)
plt.show()

t = np.linspace(0, 500, 50001)
simres = ScipyOdeSimulator(model_concentration, t, verbose=False).run(); #deterministic simulation
yout = simres.all;
observables=['obsCdc20', 'obsBAXm']; #plot Cyclin E, Cdc20 and Cyclin A

for obs in range(2) :
    plt.plot(t[:],yout[observables[obs]][:], label=observables[obs], lw=2)
plt.legend(loc=0)
plt.show()

fig, ax = plt.subplots()
bins = list(range(0,501,5))
death_time = []
t = np.linspace(0, 500, 50001)
sim = BngSimulator(model_population, tspan=t, verbose=False); #stochastic simulator
for n in range(50):
    simres = sim.run(method='ssa', verbose=False); #running the simulation; remove seed for truly random; seed=33 is used for similarity to deterministic sim
    yout = simres.observables;
    dead = yout['obsBAXm'][:]>25
    d_ind = [i for i, x in enumerate(dead) if x]
    if d_ind:
        ax.plot(t[:d_ind[0]],yout['obsBAXm'][:d_ind[0]], label=observables[obs], lw=0.5, c='grey', alpha=0.5)
        death_time.append(t[d_ind[0]])
    else:
        ax.plot(t[:],yout['obsBAXm'][:], label=observables[obs], lw=0.5, c='grey', alpha=0.5)

        
if death_time:
    plt.hist(death_time, bins=bins, density=True)
frac_dead = len(death_time)/50

