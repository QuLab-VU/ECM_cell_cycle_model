# -*- coding: utf-8 -*-
"""
Created on Mon Aug 10 11:05:56 2020

@author: Grim.M
"""

path=r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code"
import sys
sys.path.insert(0, path + "\cellular_cycle") #add a pathway to the model
from peakfinder import peakfinder
import pandas as pd
import numpy as np


def retrieve_cellState(run_time, init_cell, cells):
        
    omega=200
    IMT_eoe_counts = pd.DataFrame(columns=['End_of_Expt', 'Death' ,'Birthtime','Lifespan'])
#    divSum=0;
#    for ii in range(len(cells['molecules'])):
#        divSum=divSum+cells['division'][ii];
    
    for ii in range(len(cells['living_cell'])):
        end_index=[]
        end_index = [i for i, x in enumerate(np.array(cells['time'][ii]) < run_time) if x] #time within the experiment
    
        if len(end_index)>0 :
            (a,max_tmp_cdc) = peakfinder(cells['molecules'][ii]['obsCdc20'], 0.2*omega, 0.2*omega,1);
            (a,max_tmp_cycE) = peakfinder(cells['molecules'][ii]['obsCycE'], 0.2*omega, 0.2*omega,1);
            div_ind = max_tmp_cdc
            if ii <= (init_cell-1) : #initial cells don't include first point...
                if len(div_ind)>=2 and cells['living_cell'][ii] :
                    
                    for jj in range(len(div_ind)) :
                        if jj<(len(div_ind)-1) :
                            tmp_counts = {'End_of_Expt' : False, 'Death' : False, 'Birthtime' : cells['time'][ii][div_ind[jj]], 'Lifespan' : (cells['time'][ii][div_ind[jj+1]] - cells['time'][ii][div_ind[jj]])}
                        else :
                            tmp_counts = {'End_of_Expt' : True, 'Death' : False, 'Birthtime' : cells['time'][ii][div_ind[jj]], 'Lifespan' : (run_time - cells['time'][ii][div_ind[jj]])}
                        IMT_eoe_counts = IMT_eoe_counts.append(tmp_counts, ignore_index=True)
                        print(ii)
                        print(tmp_counts)
        
                elif len(div_ind)==1 and cells['living_cell'][ii] :
                    tmp_counts = {'End_of_Expt' : True, 'Death' : False, 'Birthtime' : cells['time'][ii][div_ind[0]], 'Lifespan' : (run_time - cells['time'][ii][div_ind[0]])}
                    IMT_eoe_counts = IMT_eoe_counts.append(tmp_counts, ignore_index=True)
                    print(ii)
                    print(tmp_counts)
                        
                elif len(div_ind)==0 and cells['living_cell'][ii] :
                    tmp_counts = {'End_of_Expt' : True, 'Death' : False, 'Birthtime' : cells['time'][ii][0], 'Lifespan' : (run_time - cells['time'][ii][0])}
                    IMT_eoe_counts = IMT_eoe_counts.append(tmp_counts, ignore_index=True)
                    print(ii)
                    print(tmp_counts)
                    
                elif len(div_ind)>=2 and cells['living_cell'][ii]==False :
                    for jj in range(len(div_ind)) :
                        if jj<(len(div_ind)-1) :
                            tmp_counts = {'End_of_Expt' : False, 'Death' : False, 'Birthtime' : cells['time'][ii][div_ind[jj]], 'Lifespan' : (cells['time'][ii][div_ind[jj+1]] - cells['time'][ii][div_ind[jj]])}
                        else :
                            tmp_counts = {'End_of_Expt' : False, 'Death' : True, 'Birthtime' : cells['time'][ii][div_ind[jj]], 'Lifespan' : (cells['time'][ii][-1] - cells['time'][ii][div_ind[jj]])}
                        IMT_eoe_counts = IMT_eoe_counts.append(tmp_counts, ignore_index=True)                
                        print(ii)
                        print(tmp_counts)
                elif len(div_ind)==1 and cells['living_cell'][ii]==False :
                    tmp_counts = {'End_of_Expt' : False, 'Death' : True, 'Birthtime' : cells['time'][ii][div_ind[0]], 'Lifespan' : (cells['time'][ii][-1] - cells['time'][ii][div_ind[0]])}
                    
                    IMT_eoe_counts = IMT_eoe_counts.append(tmp_counts, ignore_index=True)
                    print(ii)
                    print(tmp_counts)
    #            else #isempty(div_ind) && cells{ii,5}==0 :
    #                 tmp_counts{ii+1,1}=false;
    #                 tmp_counts{ii+1,2}=true;
    #                 tmp_counts{ii+1,3}=cells{ii,3}(1);
    #                 tmp_counts{ii+1,4}=cells{ii,3}(end)-cells{ii,3}(1);
    #                    IMT_eoe_counts = IMT_eoe_counts.append(tmp_counts, ignore_index=True)
            else:
                if len(div_ind)>=1 and cells['living_cell'][ii]:
                    for jj in range(len(div_ind)):
                        if jj==0:
                            tmp_counts = {'End_of_Expt' : False, 'Death' : False, 'Birthtime' : cells['time'][ii][0], 'Lifespan' : (cells['time'][ii][div_ind[0]] - cells['time'][ii][0])}
    
                        elif jj==(len(div_ind)-1):
    
                            tmp_counts = {'End_of_Expt' : True, 'Death' : False, 'Birthtime' : cells['time'][ii][div_ind[jj]], 'Lifespan' : run_time - cells['time'][ii][div_ind[jj]]}
                        
                        else:
                            tmp_counts = {'End_of_Expt' : False, 'Death' : False, 'Birthtime' : cells['time'][ii][div_ind[jj]], 'Lifespan' : (cells['time'][ii][div_ind[jj+1]] - cells['time'][ii][div_ind[jj]])}
                        IMT_eoe_counts = IMT_eoe_counts.append(tmp_counts, ignore_index=True)
                        print(ii)
                        print(tmp_counts)
                         
                elif len(div_ind)==1 and cells['living_cell'][ii]:               
                    tmp_counts = {'End_of_Expt' : False, 'Death' : False, 'Birthtime' : cells['time'][ii][0], 'Lifespan' : (cells['time'][ii][div_ind[0]] - cells['time'][ii][0])}
                    IMT_eoe_counts = IMT_eoe_counts.append(tmp_counts, ignore_index=True)
    
                    print(ii)
                    print(tmp_counts)
                    
                    tmp_counts = {'End_of_Expt' : True, 'Death' : False, 'Birthtime' : cells['time'][ii][div_ind[div_ind[0]]], 'Lifespan' : (run_time - cells['time'][ii][div_ind[0]])}
                    IMT_eoe_counts = IMT_eoe_counts.append(tmp_counts, ignore_index=True)
                    print(ii)
                    print(tmp_counts)
                    
                elif len(div_ind)==0 and cells['living_cell'][ii] :
                    tmp_counts = {'End_of_Expt' : True, 'Death' : False, 'Birthtime' : cells['time'][ii][0], 'Lifespan' : (run_time - cells['time'][ii][0])}
                    IMT_eoe_counts = IMT_eoe_counts.append(tmp_counts, ignore_index=True)
                    print(ii)
                    print(tmp_counts)
                    
                elif len(div_ind)>=2 and cells['living_cell'][ii]==False:
                    for jj in range(len(div_ind)):
                        if jj==0 :
                            tmp_counts = {'End_of_Expt' : False, 'Death' : False, 'Birthtime' : cells['time'][ii][0], 'Lifespan' : (cells['time'][ii][div_ind[0]] - cells['time'][ii][0])}
                        elif jj==(len(div_ind)-1):
                            tmp_counts = {'End_of_Expt' : False, 'Death' : True, 'Birthtime' : cells['time'][ii][div_ind[jj]], 'Lifespan' : (cells['time'][ii][-1] - cells['time'][ii][div_ind[jj]])}
                        else:
                            tmp_counts = {'End_of_Expt' : False, 'Death' : False, 'Birthtime' : cells['time'][ii][div_ind[jj]], 'Lifespan' : (cells['time'][ii][div_ind[jj+1]] - cells['time'][ii][div_ind[jj]])}
                        IMT_eoe_counts = IMT_eoe_counts.append(tmp_counts, ignore_index=True)
                    print(ii)
                    print(tmp_counts)
    
                elif len(div_ind)==1 and cells['living_cell'][ii]==False : 
    
                    tmp_counts = {'End_of_Expt' : False, 'Death' : False, 'Birthtime' : cells['time'][ii][0], 'Lifespan' : (cells['time'][ii][div_ind[0]] - cells['time'][ii][0])}
                    IMT_eoe_counts = IMT_eoe_counts.append(tmp_counts, ignore_index=True)
                    print(ii)
                    print(tmp_counts)
                    
                    tmp_counts = {'End_of_Expt' : False, 'Death' : True, 'Birthtime' : cells['time'][ii][div_ind[0]], 'Lifespan' : (cells['time'][ii][-1] - cells['time'][ii][div_ind[0]])}
                    IMT_eoe_counts = IMT_eoe_counts.append(tmp_counts, ignore_index=True)
                    print(ii)
                    print(tmp_counts)
                    
                    
                else: #len(div_ind)==0 and cells['living_cell'][ii]==False
                    tmp_counts = {'End_of_Expt' : False, 'Death' : True, 'Birthtime' : cells['time'][ii][0], 'Lifespan' : (cells['time'][ii][-1] - cells['time'][ii][0])}
                    IMT_eoe_counts = IMT_eoe_counts.append(tmp_counts, ignore_index=True)          
                    print(ii)
                    print(tmp_counts)
                    
    IMT_eoe_counts = IMT_eoe_counts[IMT_eoe_counts['Lifespan']>0]
    return IMT_eoe_counts
