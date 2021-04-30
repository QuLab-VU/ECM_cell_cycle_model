# -*- coding: utf-8 -*-
"""
Created on Tue Mar 30 15:17:10 2021

@author: Grim.M
"""

import numpy as np
import pandas as pd
import sobol_seq

##
from Austin_proliferation_model_concentration_4 import model
model_concentration=model
from Austin_proliferation_model_population_4 import model
model_population=model

#import the model and its parameters

path=r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code"

#frac_prolif_parameters = pd.read_pickle(path + r"/fit_data/sobol_seq/cost_function_model_4/cost_and_IMT_distribution/frac_prolif_parameters_all_parameters_50_initial_cells.pkl")

frac_prolif_parameters = pd.read_pickle(path + r"/fit_data/sobol_seq/cost_function_model_4/cost_and_IMT_distribution/frac_prolif_parameters_all_parameters_006_55_kPa.pkl")


IMT_distribution = pd.read_pickle(path + r"/fit_data/sobol_seq/cost_function_model_4/cost_and_IMT_distribution/IMT_distribution_all_parameters_50_initial_cells.pkl")

#import parameters and their costs

frac_prolif_parameters.sort_values(by=['cost'])

#sort the parameters according to their glocal cost
#frac_prolif_parameters.sort_values(by=['cost_mmd'])

##########

#ploting cost distribution

#setting the global cost
sum_column = frac_prolif_parameters["cost_006"] +frac_prolif_parameters["cost_frac_prolif_model_006"] + 0 * frac_prolif_parameters["cost_mmd_55"]/2

frac_prolif_parameters["cost_total"] = sum_column

#different plots:
sns.distplot(frac_prolif_parameters['cost_mmd_006'])
#sns.distplot(frac_prolif_parameters['cost_frac_prolif_model_006'])
sns.distplot(frac_prolif_parameters['cost_006'])

plt.scatter(frac_prolif_parameters['cost_mmd_55']/2,frac_prolif_parameters['cost_006'])

sns.distplot(frac_prolif_parameters['cost_total'])
sns.distplot(frac_prolif_parameters['cost_mmd_55']/2)
sns.distplot(frac_prolif_parameters['cost_006'])

np.mean(frac_prolif_parameters['cost_mmd_55']/3)
np.mean(frac_prolif_parameters['cost_006'])

##########
#generate th sobol sequence

pysb_sampled_parameter_names = ["V_s2skp","V_Fp","k_M","k_EMY","k_CDF","k_CDS","k_ED","K_S","K_MY","K_DF","d_CD","k_P1", "k_P2", "K_CD", "k_E", "K_M", "K_CE", "V_2E2F", "K_2E2F", "k_RE", "kf1", "kf2", "kb", "kasbaxmbcl2", "kdsbaxmbcl2", "k1sbh3", "k2sbh3", "kdbh3", "kasbh3bcl2", "kdsbh3bcl2"]

loc=[]

for parameter in pysb_sampled_parameter_names:
    print((model.parameters[parameter].value))
    loc+=[np.log10(model.parameters[parameter].value)]



l = sobol_seq.i4_sobol_generate(30, 120000, skip = 0)
l2 = 4*l + (np.array(loc)-2)

######################
#selecting the best parameters
#example: the ten best

frac_prolif_parameters.sort_values(by=['cost_total'])[0:100]

indices = [int(i)-1 for i in frac_prolif_parameters.sort_values(by=['cost_total'])["index"].tolist()[0:100]] #indices in the sobol sequence

parameters = [l2[i] for i in indices]


#########################
