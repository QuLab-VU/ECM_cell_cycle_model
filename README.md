# ECM_cell_cycle_model

### Creating conda env
To create a conda environment to run umap-generating code, use the following command:  

`conda create --file Python/ECM_model_env.yaml`  
_(Assumes you are in the `ECM_cell_cycle_model` top directory.)_

  
Activate the new environment:

`conda activate ECM_model`  


To allow access in JupyterHub server:

`ipython kernel install --user --name=ECM_model`
