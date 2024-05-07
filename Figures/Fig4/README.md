# Steps to reproduce Figure 4
### Summary of steps
1) `Python/04_Stoch_sims_for_UMAP.ipynb`: Perform model simulations of each -/+ erlotinib condition for matrix stiffnesses ranging from 0-99  
1) `Python/05_UMAP_of_stoch_sims.ipynb`: Generate UMAP embedded space from (normalized) simulated data (subsampled?) from all conditions  
1) `Figures/Fig4/run_lyaps.py`: Calculate Lyapunov exponent for each condition  
1) `Figures/Fig4/calcLyapunov.ipynb`: Calculate Lyapunov exponent for each condition  
1) `Figures/Fig4/Lyapunov_analysis.ipynb`: Produce plots for figure:  
    A) Plot of cell cycle trajectories of DMSO or erlotinib in UMAP embedded space colored by IMT  
    B) Plot of cell cycle trajectories of DMSO or erlotinib in UMAP embedded space colored by Lyapunov exponent  
    C) Plot of Lyapunov exponent values of each simulated condition  
    D) Plot of ODE simulation trajectories for conditions at specific Lyapunov exponent values  
    
    

