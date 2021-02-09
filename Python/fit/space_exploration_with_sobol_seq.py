# -*- coding: utf-8 -*-
"""
Created on Fri Jan 29 11:43:11 2021

@author: Grim.M
"""



path="/home/jeremie"
#path=r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code"
import sys
sys.path.insert(0, path + "/cellular_cycle") #add a pathway to the model
sys.path.insert(0, path + "/fit") #add a pathway 

from pysb.integrate import Solver
import numpy as np
from scipy.stats import norm, uniform
import os
import inspect
import pandas as pd
import pylab as pl
import sobol_seq
import math
from math import inf
import seaborn as sns



from cost_function_BngSimulator_frac_prolif_quiescence_and_mmd_DMSO_only_30_parameters_fixed_initial_concentrations_mitosis import  cost_function_reduced_quiescence_DMSO_only_30_parameters


from Austin_proliferation_model_population_3 import model
model_population=model


###############################################################################


def cost_function_f(parameter_vector):
    parameter_vector=[10**i for i in parameter_vector]
    return -cost_function_reduced_quiescence_DMSO_only_30_parameters(parameter_vector)


def cost_function_f_2(parameter_vector):
    return -cost_function_reduced_quiescence_DMSO_only_30_parameters(parameter_vector)



pysb_sampled_parameter_names = ["V_s2skp","V_Fp","k_M","k_EMY","k_CDF","k_CDS","k_ED","K_S","K_MY","K_DF","d_CD","k_P1", "k_P2", "K_CD", "k_E", "K_M", "K_CE", "V_2E2F", "K_2E2F", "k_RE", "kf1", "kf2", "kb", "kasbaxmbcl2", "kdsbaxmbcl2", "k1sbh3", "k2sbh3", "kdbh3", "kasbh3bcl2", "kdsbh3bcl2"]

loc=[]
for parameter in pysb_sampled_parameter_names:
    print(np.log10(model.parameters[parameter].value))
    loc+=[np.log10(model.parameters[parameter].value)]


####################################################################################

l=sobol_seq.i4_sobol_generate(30, 20000, skip=0)
l2 = 4*l + (np.array(loc)-2)
cost=[]

for parameter_vector in l2 :
    c = cost_function_f(parameter_vector)
    print(cost)
    cost+=[c]

np.save(path + "/fit_data/cost_sobol_seq_0000_20000.npy", cost)


###############################


cost = np.load(path+"/fit_data/cost_sobol_seq_0000_20000.npy")


sum(cost>-inf)/len(cost)

max(cost)


cost_ = cost[cost>-inf]


l=sobol_seq.i4_sobol_generate(30, 20000, skip=0)
l2 = 4*l + (np.array(loc)-2)


parameters_log = l2[cost>-inf]
parameters = 10**parameters_log
all_parameters = 10**l2

parameters_log_not_working = l2[cost==-inf]
parameters_not_working = 10**parameters_log_not_working



from sklearn import decomposition
import umap as umap

#
#pip install llvmlite --ignore-installed 
#pip install umap_learn

colors = sns.color_palette(n_colors=10)

#####################################
#####################################


#pca
pca = decomposition.PCA(n_components=2)
pca.fit(parameters)
#pca.fit(all_parameters)
X_pca = pca.transform(parameters)
X_pca_2 = pca.transform(parameters_not_working)


###
pca = decomposition.PCA(n_components=2)
pca.fit(parameters_log)
#pca.fit(l2)
X_pca = pca.transform(parameters_log)
X_pca_2 = pca.transform(parameters_log_not_working)


fig = pl.figure(figsize=(3,3),dpi=200)
ax = fig.add_subplot(111)
sc = ax.scatter(X_pca[:,0], X_pca[:,1], c=cost_, cmap="seismic", alpha = 0.8,s = 1);
ax.scatter(X_pca_2[:,0], X_pca_2[:,1], c="black", alpha = 0.1,s = 0.1);

pl.colorbar(sc)

#pl.savefig(r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\report\figures\sobol_seq\sobol_seq_cost_function\PCA_working_parameters.pdf")

pl.savefig(r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\report\figures\sobol_seq\sobol_seq_cost_function\PCA_working_parameters_log.pdf")




#umap   
trans = umap.UMAP().fit(parameters)
#trans = umap.UMAP().fit(all_parameters)
X_umap = trans.transform(parameters)
X_umap_2 = trans.transform(parameters_not_working)



###
trans = umap.UMAP().fit(parameters_log)
#trans = umap.UMAP().fit(l2)
X_umap = trans.transform(parameters_log)
X_umap_2 = trans.transform(parameters_log_not_working)



fig = pl.figure(figsize=(4,3),dpi=200)
ax = fig.add_subplot(111)
sc = ax.scatter(X_umap[:,0], X_umap[:,1], c=cost_, cmap="seismic", alpha = 0.8,s = 1);

ax.scatter(X_umap_2[:,0], X_umap_2[:,1], c="black", alpha = 0.1,s = 0.1);

ax.scatter(test_pts[:,0], test_pts[:,1], c="black", alpha = 1,s = 3);

pl.colorbar(sc)

#pl.savefig(r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\report\figures\sobol_seq\sobol_seq_cost_function\Umap_all_parameters.pdf")

#pl.savefig(r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\report\figures\sobol_seq\sobol_seq_cost_function\Umap_working_parameters_test_points.pdf")


pl.savefig(r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\report\figures\sobol_seq\sobol_seq_cost_function\Umap_working_parameters_test_points.pdf")


#inverse umap ?

corners = np.array([
    [11.1, 14.8],  # 
    [11.1, 14.8],  # 
    [11.2, 15.5],  # 
    [11.2, 15.5],  # 
])
    

#corners = np.array([
#    [8, 4],  # 
#    [8, 5],  # 
#    [9, 5],  # 
#    [9, 4],  # 
#])
#    
#corners = np.array([
#    [4, 3],  # 
#    [4, 5],  # 
#    [6, 5],  # 
#    [6, 3],  # 
#])
#        
#corners = np.array([
#    [11, 0],  # 
#    [11, 6],  # 
#    [15, 6],  # 
#    [15, 0],  # 
#])
#    
#corners = np.array([
#    [12, 3],  # 
#    [12, 5],  # 
#    [14, 5],  # 
#    [14, 3],  # 
#])
#    
#corners = np.array([
#    [12, 7],  # 
#    [12, 8],  # 
#    [13, 8],  # 
#    [13, 7],  # 
#])
#    
#        
#corners = np.array([
#    [11, 4],  # 
#    [11, 9],  # 
#    [15, 9],  # 
#    [15, 4],  # 
#])
#    
#corners = np.array([
#    [13, 4],  # 
#    [13, 6],  # 
#    [15, 6],  # 
#    [15, 4],  # 
#])
#    
#corners = np.array([
#    [3, 2],  # 
#    [3, 4],  # 
#    [5, 4],  # 
#    [5, 2],  # 
#])

    
test_pts = np.array([
    (corners[0]*(1-x) + corners[1]*x)*(1-y) +
    (corners[3]*(1-x) + corners[2]*x)*y
    for y in np.linspace(0, 1, 10)
    for x in np.linspace(0, 1, 10)
])

    
inv_transformed_points = trans.inverse_transform(test_pts)


np.save(path + "/fit_data/inv_transformed_points_log.npy", inv_transformed_points)
np.save(path + "/fit_data/inv_transformed_points.npy", inv_transformed_points)


inv_transformed_points_log = np.load(path+"/fit_data/inv_transformed_points_log.npy")
inv_transformed_points = np.load(path+"/fit_data/inv_transformed_points.npy")


cost=[]

for parameter_vector in inv_transformed_points_log :
    print(parameter_vector)
    c = cost_function_f(parameter_vector)
    print(cost)
    cost+=[c]


cost=[]

for parameter_vector in inv_transformed_points :
    print(parameter_vector)
    c = cost_function_f_2(parameter_vector)
    print(cost)
    cost+=[c]



sum(np.array(cost)>-inf)/len(np.array(cost))
#22%

