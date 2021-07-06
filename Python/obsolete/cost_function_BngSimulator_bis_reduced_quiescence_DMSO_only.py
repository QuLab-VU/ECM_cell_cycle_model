# -*- coding: utf-8 -*-
"""
Created on Wed May 20 10:42:26 2020

@author: Grim.M
"""

path="/home/jeremie"

path=r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code"

import sys
sys.path.insert(0, path+"/cellular_cycle") #add a pathway to the model
sys.path.insert(0, path+"/fit") #add a pathway to the model

from Austin_proliferation_model_population import model
model_population=model

from pysb.simulator import BngSimulator

from SMs_G1s_IMTs_simulation_BngSimulator_reduced import IMTs
from rbf_mmd import nmmd2_rbf


import random

import os
os.environ["PYOPENCL_CTX"] ="0"

import pylab as pl
import pandas as pd
import math
import numpy as np

tmp_t=pl.linspace(0, 500, 2500)
sim = BngSimulator(model_population, tspan=tmp_t, verbose=True)
traj = sim.run()

###############################################################################


experimental_SMs={}
experimental_G1s={}
experimental_IMTs={}



df = pd.read_csv(path+'/IMT_times/DS8_SM_0_2_55.csv')

for col in df.columns :
    experimental_SMs[col]=df[col].dropna().tolist()


df = pd.read_csv(path+'/IMT_times/DS8_G1_0_2_55.csv')

for col in df.columns :
    experimental_G1s[col]=df[col].dropna().tolist()

for k in experimental_G1s.keys() :
    minimum=min(len(experimental_G1s[k]),len(experimental_SMs[k]))
    experimental_IMTs[k]=[experimental_G1s[k][i]+experimental_SMs[k][i] for i in range(minimum)]


df = pd.read_csv(path+'/cell_fates/ds8_06_DMSO_fates.csv')
lines = {}
for col in df :
    lines[col]=df[col].dropna().tolist()

q=m=o=d=0
lifetimes={}
lifetimes['6 kPa DMSO']=[]

for (death, mitotic_frame, end_frame, lifetime) in zip(lines['Death'], lines['Mitotic Frame'], lines['End Frame'], lines['Lifetime']):
    if death==0:
        print('life')
        if mitotic_frame!=0 and end_frame!=186:
            print('mtosis')
            m+=1
            lifetimes['6 kPa DMSO']+=[lifetime]
            
        elif  mitotic_frame==0 and end_frame==186:
            print('quiescent')
            q+=1
            
        else :
            print('ther')
            o+=1 
    else:
        d+=1

p1={}
n1={}
p1["0.2 kPa DMSO"]=q/(q+m+d)
p1["55 kPa DMSO"]=0
n1["0.2 kPa DMSO"]=q+m+d
n1["55 kPa DMSO"] = len(experimental_G1s["55 kPa DMSO"]) #no important cuse p1=0

pysb_sampled_parameter_names = ["V_s2skp","V_Fp","k_M","k_EMY","k_CDF","k_CDS","k_ED","K_S","K_MY","K_DF"]

###############################################################################
#tmp_t for different conditions


tmp_t={}

for condition in experimental_IMTs.keys():

    mean_g1=sum(experimental_G1s[condition])/len(experimental_G1s[condition])
    mean_sm=sum(experimental_SMs[condition])/len(experimental_SMs[condition])
    
    time = ((25*(mean_g1+mean_sm)*1.3))
    tmp_t[condition] = np.linspace(0, time, int(time*10))
    
    
##############################################################################
conditions=['0.2 kPa DMSO','55 kPa DMSO']


def cost_function_reduced_quiescence_DMSO_only(parameter_vector):    
    
    print("les paramètres sont ")
    print(parameter_vector)
    
    global experimental_G1s, experimental_SMs, model_concentration, model_population, tmp_t, p1, n1
    
    param_dict = {pname: pvalue for pname, pvalue in zip(pysb_sampled_parameter_names, parameter_vector)}
    
    
    for pname, pvalue in param_dict.items():
        
        #Change model parameter values to current location in parameter space
        model_population.parameters[pname].value = (pvalue)
        
    #Simulate experimentally IMTs
    n_mmd2_tot=[]
    simulated_G1s=[]
    simulated_SMs=[]
    omega=200
    
    for condition in conditions:
        print(condition)
        
        if 'DMSO' in condition:            
            if '0.2 kPa' in condition :
                #model_concentration.parameters['ERK'].value = 1.06970165050038 
                #model_concentration.parameters['FAK'].value = 0.337143388171177
                #model_concentration.parameters['AKT'].value = 0.616658111482869
                
                model_population.parameters['ERK'].value = 1.06970165050038 * omega
                model_population.parameters['FAK'].value = 0.337143388171177 * omega
                model_population.parameters['AKT'].value = 0.616658111482869 * omega
            elif '4 kPa' in condition :
                #model_concentration.parameters['ERK'].value = 0.897485095798471
                #model_concentration.parameters['FAK'].value = 0.594150207445117
                #model_concentration.parameters['AKT'].value = 1.01975420186633
                
                model_population.parameters['ERK'].value = 0.897485095798471 * omega
                model_population.parameters['FAK'].value = 0.594150207445117 * omega
                model_population.parameters['AKT'].value = 1.01975420186633 * omega
            else : #55 kPa
                #model_concentration.parameters['ERK'].value = 1
                #model_concentration.parameters['FAK'].value = 1
                #model_concentration.parameters['AKT'].value = 1
                
                model_population.parameters['ERK'].value = 1 * omega
                model_population.parameters['FAK'].value = 1 * omega
                model_population.parameters['AKT'].value = 1 * omega
        else : #ERL treated
            if '0.2 kPa' in condition :
                #model_concentration.parameters['ERK'].value = 0.186745227418267
                #model_concentration.parameters['FAK'].value = 0.397060244114162
                #model_concentration.parameters['AKT'].value = 0.245702711121333
                
                model_population.parameters['ERK'].value = 0.186745227418267 * omega
                model_population.parameters['FAK'].value = 0.397060244114162 * omega
                model_population.parameters['AKT'].value = 0.245702711121333 * omega
            elif '4 kPa' in condition :
                #model_concentration.parameters['ERK'].value = 0.211453830597783
                #model_concentration.parameters['FAK'].value = 0.709172297576571
                #model_concentration.parameters['AKT'].value = 0.598518145949655
                
                model_population.parameters['ERK'].value = 0.211453830597783 * omega
                model_population.parameters['FAK'].value = 0.709172297576571 * omega
                model_population.parameters['AKT'].value = 0.598518145949655 * omega
            else : #55 kPa
                #model_concentration.parameters['ERK'].value = 0.273618423688386
                #model_concentration.parameters['FAK'].value = 1.075820621939
                #model_concentration.parameters['AKT'].value = 0.585067267459499
                
                model_population.parameters['ERK'].value = 0.273618423688386 * omega
                model_population.parameters['FAK'].value = 1.075820621939 * omega
                model_population.parameters['AKT'].value = 0.585067267459499 * omega
        simulated_G1s=[]
        simulated_SMs=[]
        
        (IMT_tmp, a, b) = IMTs(tmp_t[condition], model_population)
        print(a)
        print(len(a))
        print(b)
        print(len(b))
        
        
        if b==[] or a==[] or len(a)<25 or len(b)<25: #same nber of times
            return 10000
        
        n2=len(IMT_tmp)
        p2=sum(np.array(IMT_tmp)>102.3)/n2
        
        a=random.sample(a, 25)
        b=random.sample(b, 25)
        g1=random.sample(experimental_G1s[condition], 25)
        sm=random.sample(experimental_SMs[condition], 25)
        
        cost_quiescence=80/(1+math.exp(-25*abs(p1[condition]-p2)+1))
        
        n_mmd2_tot+=[nmmd2_rbf(np.array([b]).transpose(), np.array([g1]).transpose(), 1)]
        n_mmd2_tot+=[nmmd2_rbf(np.array([a]).transpose(), np.array([sm]).transpose(), 1)]
    print("ok")
    print("cost_quiescence")
    return sum(n_mmd2_tot)+cost_quiescence

#cost_quiescence test
#deltap_range=np.linspace(0,1,100)
#cost_quiescence_range=[80/(1+math.exp(-30*(deltap))) for deltap in deltap_range ] #avant60 et -16
#pl.plot(deltap_range, cost_quiescence_range)
