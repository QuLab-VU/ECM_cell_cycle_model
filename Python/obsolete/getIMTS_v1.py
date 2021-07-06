# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 18:18:24 2020

@author: Jérémie Marlhens
Title: getIMTS_v1

"""
  
IMTs=[];
IMT_tmp=[];
g1 = [];
g1_tmp =[];
s_m = [];
s_m_tmp = [];

for ii in range(len(cells['molecules'])) :

    if len(cells['time'][ii])>1 :
       a, max_tmp_cdc=peakfinder(cells['molecules'][ii]['obsCdc20'],0.15*200,0.15*200);      #cdc20
       a, max_tmp_cycE=peakfinder((cells['molecules'][ii]['obsCycE']),0.15*200,0.15*200);  #cycE
       if len(max_tmp_cdc)>0 and len(max_tmp_cycE)>0 :          
            div_ind = max_tmp_cdc;
            if len(div_ind)>=2 :
                for jj in range(1,len(div_ind)) :
                    IMT_tmp.append(cells['time'][ii][div_ind[jj]]-cells['time'][ii][div_ind[jj-1]])
    
            for jj in range(len(max_tmp_cycE)) :
                tmp = [x for x in (max_tmp_cdc) if x > max_tmp_cycE[jj]];
                if len(tmp)>0 :
                    s_m_tmp.append(cells['time'][ii][tmp[0]] - cells['time'][ii][max_tmp_cycE[jj]])
        
            for jj in range(len(max_tmp_cdc)) :
                tmp = [x for x in (max_tmp_cycE) if x > max_tmp_cdc[jj]];
                if len(tmp)>0 :
                    g1_tmp.append(cells['time'][ii][tmp[0]] - cells['time'][ii][max_tmp_cdc[jj]]);


    
    IMTs=IMTs + IMT_tmp;
    g1 = g1 + g1_tmp;
    s_m = s_m + s_m_tmp;
    IMT_tmp=[];
    g1_tmp = [];
    s_m_tmp = [];


###Tests
    
import math
import numpy as np
from scipy import stats
import statsmodels.api as sm
import matplotlib.pyplot as plt
from statsmodels.distributions.mixture_rvs import mixture_rvs
import pylab as pl

#Kernel density estimation (ou estimation de densité par noyau)
IMTs
kde = sm.nonparametric.KDEUnivariate(IMTs)
kde.fit() # Estimate the densities


fig = plt.figure(figsize=(12, 5))
ax = fig.add_subplot(111)

# Plot the KDE as fitted using the default arguments
ax.plot(kde.support, kde.density, lw=3, label='KDE IMTs from PySB simulation')

kde = sm.nonparametric.KDEUnivariate(IMTs_Matlab)
kde.fit() # Es
ax.plot(kde.support, kde.density, lw=3, label='KDE IMTs from Matlab simulation')

ax.legend(loc='best')
ax.grid(True, zorder=-5)


import scipy.io
scipy.io.savemat(r'C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code\cellular_cycle\max_tmp_cycE.mat',  {'max_tmp_cycE':max_tmp_cycE})
scipy.io.savemat(r'C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code\cellular_cycle\max_tmp_cdc.mat',  {'max_tmp_cdc':max_tmp_cdc})
scipy.io.savemat(r'C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code\cellular_cycle\time.mat',  {'time':cells['time'][0]})
scipy.io.savemat(r'C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code\cellular_cycle\IMTs_PySB.mat',  {'IMTs_PySB':IMTs_PySB})



from scipy.io import loadmat
IMTs_Matlab = loadmat('C:/Users/Grim.M/Ecole/ens/Master/S2/Stage/Code/cellular_cycle/IMTs_Matlab.mat')
l=IMTs_Matlab['IMTs_Matlab'].tolist()
IMTs_Matlab=[item for sublist in l for item in sublist]

IMTs_PySB = loadmat('C:/Users/Grim.M/Ecole/ens/Master/S2/Stage/Code/cellular_cycle/IMTs_PySB.mat')
l=IMTs_PySB['IMTs_PySB'].tolist()
IMTs_PySB=l[0]



def RBS(x,x0,sigma) :
    observations , dimensions= np.shape(x0)
    s=0
    for i in range(observations):
        s+=math.exp(-np.linalg.norm(x-x0[i])**2/(2*sigma**2))
    return s/observations

sigma=2

x=np.linspace(0,70,700)
y_Matlab=[RBS(i,np.transpose(np.array([IMTs_Matlab])),sigma) for i in x]
y_PySB=[RBS(i,np.transpose(np.array([IMTs_PySB])),sigma) for i in x]


fig = pl.figure(figsize=(12, 5))
ax = fig.add_subplot(111)

# Plot the KDE as fitted using the default arguments
ax.plot(x, y_Matlab, lw=3, label='KDE IMTs from experiments - 55 kPa DMSO')

ax.plot(x, y_PySB, lw=3, label='KDE IMTs from PySB simulation')

ax.legend(loc='best')
ax.grid(True, zorder=-5)

plt.savefig(r'C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\team meeting presentation\data2.pdf')  


path=r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code"
import sys
sys.path.insert(0, path+"/fit") #add a pathway to the model


from rbf_mmd import nmmd2_rbf

length=min(len(IMTs_PySB), len(IMTs_Matlab))
nmmd2_rbf(np.array([IMTs_Matlab[0:length]]).transpose(), np.array([IMTs_PySB[0:length]]).transpose(), 1)
