# Overview of MATLAB code

This is an overview of how to execute simulation of the ECM-dependent cell cycle model. We describe two different modes of simulation:  

1) Keeping track of number of cells necesary for the population dynamics model---including the model of cell death
2) Running for 5000 hours tracking only cell cycle progression. 

Necessary MATLAB toolbox: Statistics and Machine Learning Toolbox

## Simulations for population dynamics model

This is a stochastic simulation that captures the dynamics of the molecular species in all cells during the simulation and tracks cell death and 
division events.

**Inputs**: AKT, ERK, FAK activity levels each in a range from 0 - 1 (1 represents the maximum value as determined experimentally) and the initial number
of cells.  
**Output**: dataframe [cells] containing the molecular species behavior, time, number of divisions, and death. Matrix [cell_count] containing
time of division or death (i.e. adding or subtracting a cell from the simulation) used to measure cell population response.

Using MATLAB function runStochCell.m define AKT, ERK, FAK, and initial cell count. This function returns a cell structure [cells] and [cell_count].
Example:
``` 
AKT = 0.6; ERK = 0.2; FAK = 0.6; init = 1; dt = 0.01;
[cells, cell_count] = runStochCell(AKT,ERK,FAK,init,dt)

```

Default input values are related to the 4kPa substrate stiffness, erlotinib-treated experimental values.

rescale_prot_traj.m makes uniform time steps from a stochastic run for UMAP projection.

set random number seed