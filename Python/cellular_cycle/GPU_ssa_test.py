###############################################################################
# GPU SSA test


PYOPENCL_CTX = '0'

import numpy as np
from pysb.examples.earm_1_0 import model
from pysb.simulator import OpenCLSSASimulator, ScipyOdeSimulator

if __name__ == "__main__":
    n_sim = 100
    tspan = np.linspace(0, 5000, 101)
    name = 'tBid_total'
    sim = OpenCLSSASimulator(model, tspan=tspan, verbose=True)
    traj = sim.run(tspan=tspan, number_sim=n_sim)

    result = np.array(traj.observables)[name]

    x = np.array([tr[:] for tr in result]).T

from pysb.examples.schloegl import model

if __name__ == '__main__':
    n_sim = 20
    x_0 = 100

    obs_name = 'X_total'
    tspan = np.linspace(0, 100, 101)
    x_0 = int(x_0)

    savename = 'schloegl_{}'.format(int(x_0))
    model.parameters['X_0'].value = x_0
    # create simulator and run simulations
    simulator = OpenCLSSASimulator(model, verbose=True,
                                   initials=[[10, 5] for i in range(10)] + [[1, 2] for i in range(10)])

    traj = simulator.run(
        tspan,
        number_sim=n_sim
    )

    x = traj.dataframe[obs_name].unstack(0).values

# import the pysb module and all its methods and functions
from pysb import *
import numpy as np
from pysb.simulator import OpenCLSSASimulator, ScipyOdeSimulator

# instantiate a model
Model()

# declare monomers
Monomer('C8', ['b'])
Monomer('Bid', ['b', 'S'], {'S': ['u', 't']})

# input the parameter values
Parameter('kf', 1.0e-07)
Parameter('kr', 1.0e-03)
Parameter('kc', 1.0)

Expression('kg', kf / (kf + kc))

# now input the rules
Rule('C8_Bid_bind', C8(b=None) + Bid(b=None, S='u') | C8(b=1) % Bid(b=1, S='u'), kf, kr)
Rule('tBid_from_C8Bid', C8(b=1) % Bid(b=1, S='u') >> C8(b=None) + Bid(b=None, S='t'), kg)

model.monomers
model.parameters
model.rules

Observable('obsC8', C8(b=None))
Observable('obsBid', Bid(b=None, S='u'))
Observable('obstBid', Bid(b=None, S='t'))

Parameter('C8_0', 1000)
Parameter('Bid_0', 10000)
Initial(C8(b=None), C8_0)
Initial(Bid(b=None, S='u'), Bid_0)

from pysb.simulator import ScipyOdeSimulator

t = np.linspace(0, 20000)
simres = ScipyOdeSimulator(model, tspan=t).run()
yout = simres.all

n_sim = 4
sim = OpenCLSSASimulator(model, tspan=t, verbose=True)
traj = sim.run(tspan=t, number_sim=n_sim)

# import the pysb module and all its methods and functions
from pysb import *
import numpy as np
from pysb.simulator import OpenCLSSASimulator, ScipyOdeSimulator
from pysb.simulator import BngSimulator

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

# Expression('kg', (kf+kc)) #expression! OpenCLSSASimulator will work with an expression of this type
Expression('kg', kf / (kf + obsBid))  # expression! OpenCLSSASimulator won't work with an expression of this type

# now input the rules
Rule('C8_Bid_bind', C8(b=None) + Bid(b=None, S='u') | C8(b=1) % Bid(b=1, S='u'), kf, kr)
Rule('tBid_from_C8Bid', C8(b=1) % Bid(b=1, S='u') >> C8(b=None) + Bid(b=None, S='t'), kg)

model.monomers
model.parameters
model.rules

Parameter('C8_0', 1000)
Parameter('Bid_0', 10000)
Initial(C8(b=None), C8_0)
Initial(Bid(b=None, S='u'), Bid_0)

# Check if the simulation ODE works
t = np.linspace(0, 20000)
simres = ScipyOdeSimulator(model, tspan=t).run()
yout = simres.all

# Check bng
sim = BngSimulator(model, tspan=t, verbose=True)
simres = sim.run(method='ssa', verbose=True)
yout = simres.observables

# And check openCLSSASimulator
n_sim = 4
sim = OpenCLSSASimulator(model, tspan=t, verbose=True)
traj = sim.run(tspan=t, number_sim=n_sim)

###############################################################################
# test w cell cycle


from runStochCell import runStochcell
from peakfinder import peakfinder

from pysb.simulator import OpenCLSSASimulator  # from gpu_ssa branch in Pysb github repos

pysb_sampled_parameter_names = ["V_s2skp", "V_Fp", "k_M", "k_EMY", "k_CDF", "k_CDS", "k_ED", "K_S", "K_MY", "K_DF",
                                "d_CD", "k_P1", "k_P2", "K_CD", "k_E", "K_M", "K_CE", "V_2E2F", "K_2E2F", "k_RE", "kf1",
                                "kf2", "kb", "kasbaxmbcl2", "kdsbaxmbcl2", "k1sbh3", "k2sbh3", "kdbh3", "kasbh3bcl2",
                                "kdsbh3bcl2"]
parameter_vector_original = []
for pname in pysb_sampled_parameter_names:
    # Change model parameter values to current location in parameter space
    parameter_vector_original += [model_population.parameters[pname].value]

v = 0.05;
sys_size = 200;  # cell size
thr = 0.26
run_time = 120
observed_cells_time = run_time

idx = 1
condition = ''

variable_range, check, tmp_score = runStochcell(model_population, model_concentration, idx, condition, thr);

check_ode = check == False or np.any(np.array([variable_range[key] for key in
                                               variable_range.keys()]) < 0) or tmp_score > 0.2  # the model doesn't oscillate or the

cells = {}
cells['molecules'] = [];  # Number of molecules taken randomly at one point during the cell cycle
cells['time'] = [];  #
cells['division'] = [];  #
cells['living_cell'] = [];  # living cell?
cells['recording'] = [];

# plt.plot(variable_range['obsCdc20'])

init = 2

for i in range(init):
    cells['molecules'] = cells['molecules'] + [{}]
    a = np.random.randint(len(variable_range['obsSkp2']))
    for var in list(variable_range.keys()):
        cells['molecules'][i][var] = [round(variable_range[var][a] * 200)]
    cells['time'] = cells['time'] + [[0]]
    cells['division'] = cells['division'] + [0]
    cells['living_cell'] = cells['living_cell'] + [True]  # living cell = True
    cells['recording'] = cells['recording'] + [
        True]  # recording the cell = True (if it isn't dead and if no event didn't occur after observed_cells_time)

cell_count = np.array(
    [[0., len(cells['molecules'])]]);  # start the cell counter at time zero and with number of initial cells
current_time = np.transpose(np.concatenate([np.zeros((1, len(cells['molecules']))), np.ones(
    (1, len(cells['molecules'])))]));  # keeps track of the current time of each cell and whether it's dead or not

# index of cell we still have to simulate
index_cell_recording = [i for i, x in enumerate(cells['recording']) if x]

# index of cells that are below the 'observed time', in which we take into account new cell births
cell_below_observed_time = current_time[:, 0] < observed_cells_time
index_cell_below_observed_time = [i for i, x in enumerate(cell_below_observed_time) if x]

# Initial molecule numbers for the different cells to be simulated
Initials = []
for ii in index_cell_recording:
    Initial = []
    for j in cells['molecules'][ii]:
        Initial = Initial + [cells['molecules'][ii][j][-1]]
    Initials += [Initial]

cell_count_tmp = np.zeros(
    (1, 2));  # temporary cell count with current cell time and whether a cell was added, taken away, or not

T = 1
tmp_t = np.linspace(0, T, int(T * 10))

simulator = OpenCLSSASimulator(model_population, verbose=3, initials=Initials)
simres = simulator.run(tspan=tmp_t, number_sim=len(index_cell_recording));

simres.dataframe['obsBAX'].unstack(0).values


