# -*- coding: utf-8 -*-
"""
Created on Mon Aug 10 11:05:56 2020

@author: Grim.M
"""

path="/home/jeremie"
#path=r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code"


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
        #depends on when did we stop to record
        end_index=[]
        end_index = [i for i, x in enumerate(np.array(cells['time'][ii]) < run_time) if x] #time within the experiment
    
        if len(end_index)>0 :
            if len(cells['molecules'][ii]['obsCdc20'])>20:
                (a,max_tmp_cdc) = peakfinder(cells['molecules'][ii]['obsCdc20'], 0.15*omega, 0.15*omega, 1, includeEndpoints=1);
                (a,max_tmp_cycE) = peakfinder(cells['molecules'][ii]['obsCycE'], 0.15*omega, 0.15*omega, 1, includeEndpoints=1); #those values are chosen (threshold)
            else : #the course is very short so there can't be a division
                max_tmp_cdc=[]
            div_ind = max_tmp_cdc
            
            if cells['time'][ii][-1]<run_time : #if we didn't record the time until the end of the experiment
                if cells['time'][ii][0]>0 : #start at a dividing point and ends at a dividing point or death

                    if cells['living_cell'][ii]==False :#dead
                        tmp_counts = {'End_of_Expt' : False, 'Death' : True, 'Birthtime' : cells['time'][ii][0], 'Lifespan' : (cells['time'][ii][-1] - cells['time'][ii][0])}

                        IMT_eoe_counts = IMT_eoe_counts.append(tmp_counts, ignore_index=True)
                    elif len(div_ind) == 2: #Should divide once or die
                        tmp_counts = {'End_of_Expt' : False, 'Death' : False, 'Birthtime' : cells['time'][ii][0], 'Lifespan' : (cells['time'][ii][div_ind[1]] - cells['time'][ii][div_ind[0]])}
                        IMT_eoe_counts = IMT_eoe_counts.append(tmp_counts, ignore_index=True)


                else: #start at the beginning and ends as a death, a division, or quiescence, take into account only if two division, otherwise no information
                    #print(cells['time'][ii][div_ind[0]])
                    if len(div_ind) == 2 :
                        tmp_counts = {'End_of_Expt' : False, 'Death' : False, 'Birthtime' : cells['time'][ii][div_ind[0]], 'Lifespan' : (cells['time'][ii][div_ind[1]] - cells['time'][ii][div_ind[0]])}
                        
                        
                        IMT_eoe_counts = IMT_eoe_counts.append(tmp_counts, ignore_index=True)

            else:    
                if ii <= (init_cell-1) : #initial cells don't include first point...
                    if len(div_ind)>=2 and cells['living_cell'][ii] :
                        
                        for jj in range(len(div_ind)) :
                            if jj<(len(div_ind)-1) :
                                tmp_counts = {'End_of_Expt' : False, 'Death' : False, 'Birthtime' : cells['time'][ii][div_ind[jj]], 'Lifespan' : (cells['time'][ii][div_ind[jj+1]] - cells['time'][ii][div_ind[jj]])}
                            else :
                                tmp_counts = {'End_of_Expt' : True, 'Death' : False, 'Birthtime' : cells['time'][ii][div_ind[jj]], 'Lifespan' : (run_time - cells['time'][ii][div_ind[jj]])}
                            IMT_eoe_counts = IMT_eoe_counts.append(tmp_counts, ignore_index=True)
            
                    elif len(div_ind)==1 and cells['living_cell'][ii] : #one division
                        tmp_counts = {'End_of_Expt' : True, 'Death' : False, 'Birthtime' : cells['time'][ii][div_ind[0]], 'Lifespan' : (run_time - cells['time'][ii][div_ind[0]])}
                        IMT_eoe_counts = IMT_eoe_counts.append(tmp_counts, ignore_index=True)

                            
                    elif len(div_ind)==0 and cells['living_cell'][ii] :
                        tmp_counts = {'End_of_Expt' : True, 'Death' : False, 'Birthtime' : cells['time'][ii][0], 'Lifespan' : (run_time - cells['time'][ii][0])}
                        IMT_eoe_counts = IMT_eoe_counts.append(tmp_counts, ignore_index=True)

                        
                    elif len(div_ind)>=2 and cells['living_cell'][ii]==False :
                        for jj in range(len(div_ind)) :
                            if jj<(len(div_ind)-1) :
                                tmp_counts = {'End_of_Expt' : False, 'Death' : False, 'Birthtime' : cells['time'][ii][div_ind[jj]], 'Lifespan' : (cells['time'][ii][div_ind[jj+1]] - cells['time'][ii][div_ind[jj]])}
                            else :
                                tmp_counts = {'End_of_Expt' : False, 'Death' : True, 'Birthtime' : cells['time'][ii][div_ind[jj]], 'Lifespan' : (cells['time'][ii][-1] - cells['time'][ii][div_ind[jj]])}
                            IMT_eoe_counts = IMT_eoe_counts.append(tmp_counts, ignore_index=True)                

                    elif len(div_ind)==1 and cells['living_cell'][ii]==False :
                        tmp_counts = {'End_of_Expt' : False, 'Death' : True, 'Birthtime' : cells['time'][ii][div_ind[0]], 'Lifespan' : (cells['time'][ii][-1] - cells['time'][ii][div_ind[0]])}
                        
                        IMT_eoe_counts = IMT_eoe_counts.append(tmp_counts, ignore_index=True)
#                        print(ii)
#                        print(tmp_counts)
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

                             
                    elif len(div_ind)==1 and cells['living_cell'][ii]:               
                        tmp_counts = {'End_of_Expt' : False, 'Death' : False, 'Birthtime' : cells['time'][ii][0], 'Lifespan' : (cells['time'][ii][div_ind[0]] - cells['time'][ii][0])}
                        IMT_eoe_counts = IMT_eoe_counts.append(tmp_counts, ignore_index=True)
        

                        tmp_counts = {'End_of_Expt' : True, 'Death' : False, 'Birthtime' : cells['time'][ii][div_ind[div_ind[0]]], 'Lifespan' : (run_time - cells['time'][ii][div_ind[0]])}
                        IMT_eoe_counts = IMT_eoe_counts.append(tmp_counts, ignore_index=True)

                    elif len(div_ind)==0 and cells['living_cell'][ii] :
                        tmp_counts = {'End_of_Expt' : True, 'Death' : False, 'Birthtime' : cells['time'][ii][0], 'Lifespan' : (run_time - cells['time'][ii][0])}
                        IMT_eoe_counts = IMT_eoe_counts.append(tmp_counts, ignore_index=True)
  
                    elif len(div_ind)>=2 and cells['living_cell'][ii]==False:
                        for jj in range(len(div_ind)):
                            if jj==0 :
                                tmp_counts = {'End_of_Expt' : False, 'Death' : False, 'Birthtime' : cells['time'][ii][0], 'Lifespan' : (cells['time'][ii][div_ind[0]] - cells['time'][ii][0])}
                            elif jj==(len(div_ind)-1):
                                tmp_counts = {'End_of_Expt' : False, 'Death' : True, 'Birthtime' : cells['time'][ii][div_ind[jj]], 'Lifespan' : (cells['time'][ii][-1] - cells['time'][ii][div_ind[jj]])}
                            else:
                                tmp_counts = {'End_of_Expt' : False, 'Death' : False, 'Birthtime' : cells['time'][ii][div_ind[jj]], 'Lifespan' : (cells['time'][ii][div_ind[jj+1]] - cells['time'][ii][div_ind[jj]])}
                            IMT_eoe_counts = IMT_eoe_counts.append(tmp_counts, ignore_index=True)

                    elif len(div_ind)==1 and cells['living_cell'][ii]==False : 
        
                        tmp_counts = {'End_of_Expt' : False, 'Death' : False, 'Birthtime' : cells['time'][ii][0], 'Lifespan' : (cells['time'][ii][div_ind[0]] - cells['time'][ii][0])}
                        IMT_eoe_counts = IMT_eoe_counts.append(tmp_counts, ignore_index=True)

                        tmp_counts = {'End_of_Expt' : False, 'Death' : True, 'Birthtime' : cells['time'][ii][div_ind[0]], 'Lifespan' : (cells['time'][ii][-1] - cells['time'][ii][div_ind[0]])}
                        IMT_eoe_counts = IMT_eoe_counts.append(tmp_counts, ignore_index=True)
                    else: #len(div_ind)==0 and cells['living_cell'][ii]==False
                        tmp_counts = {'End_of_Expt' : False, 'Death' : True, 'Birthtime' : cells['time'][ii][0], 'Lifespan' : (cells['time'][ii][-1] - cells['time'][ii][0])}
                        IMT_eoe_counts = IMT_eoe_counts.append(tmp_counts, ignore_index=True)          

    IMT_eoe_counts = IMT_eoe_counts[IMT_eoe_counts['Lifespan']>0]
    return IMT_eoe_counts
