# ECM_cell_cycle_model

### Creating conda env
To create a conda environment to run umap-generating code, use the following command:

`conda create -n ECM_model -c conda-forge umap-learn python=3.8 seaborn matplotlib pandas joblib numpy bokeh jupyter jupyterlab pylablib cython`
  
Add the pysb module:
  
`conda install -c alubbock pysb`  

Add the sobol_seq module:

`pip install sobol-seq`

To allow access in JupyterHub server:

`ipython kernel install --user --name=ECM_model`
