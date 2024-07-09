# Data Generation

### Commands and description for generating all the necessary simulated data

These files load the necessary functions from `cell_cycle_death_model`, `model_fxns`, and `Cell_Population_Simulation` to produce the simulated data analyzed in `/Figures`.

The script generates a folder, `~/ECM_mod_output`, in the user's home directory to store files locally. 

### Script Overview

`01_GenerateCellPopulations.ipynb`

- Creates the folder `~/ECM_mod_output` to store all simulated data in.
- Runs a simulated experiment at each experimental condition (low, intermediate, and high stiffness; with and without erlotinib) with three replicates and saves the total cell counts to file.
- Generates simulated IMT and G1 timing and saves both to file.
- Data are used in generation of Figure 3.

`02_runODEs.ipynb`

- Generates the ODE solutions of the model at 1000 stiffness incremental inputs with and without erlotinib.
- Uses the package `multiprocess` for parallel computing; the number of cores to be used will be determined based on system information.
- ODE time courses and calculated IMTs are stored in `~/ECM_mod_output/ODE_output/`.
- The analyzed data are used when generating Figure 4.

`03_ECM_model_output.ipynb`

- Generates ODE and stochastic solutions at each of the six experimental conditions.
- The same script is used to analyze the effects Cyclin inhibitor on the cell cycle by setting the inhibitor concentration.
- The output is used to generate the initial UMAP object and used in Figures 5 and 6.

`04_Stoch_sims_for_UMAP.ipynb`

- Generates ODE and stochastic solutions for full range of stiffnesses for Lyapunov analysis (1000 each DMSO and ERL)
- Identified conditions of interest using Lyapunov analysis (Fig 4C)
- Save specific solutions to embed in UMAP space for Figure 4D.
- Example filename: ODE_sim_DMSO_stiff_074_orig.csv

`05_UMAP_of_stoch_sims.ipynb`

- Creates the UMAP object used in the projections and analysis shown in Figures 4-7.
- Loads the data produced by `03_ECM_model_output.ipynb` and performs a normalization and windowing step.
- Saves the UMAP object, the stochastic projections, the ODE projections, and data points used in the projections.

`06_calcODE_UMAPs.ipynb`

- 