# -*- coding: utf-8 -*-
"""
Created on Mon Aug 10 11:05:56 2020

@author: Grim.M

Edited on Sat Aug 20 16:32:00 2022
@author: Darren Tyson
"""

from peakfinder import peakfinder
import pandas as pd
import numpy as np
from tempfile import TemporaryDirectory

def retrieve_cellState(run_time, init_cell, cells, return_all_prot=False):
    all_prot = np.empty((0,len(cells['molecules'][0])+3)) #cell ind, time, relative time, protein time courses
    index_cell=0 #cell identification nber
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
                try:
                    # choose threshold values
                    (a,max_tmp_cdc) = peakfinder(cells['molecules'][ii]['obsCdc20'], 0.15*omega, 0.15*omega, 1, includeEndpoints=1);
                    (a,max_tmp_cycE) = peakfinder(cells['molecules'][ii]['obsCycE'], 0.15*omega, 0.15*omega, 1, includeEndpoints=1); 
                except:
                    np.save(TemporaryDirectory("pydream_frac_prolif_quiescence_obscycE.npy"), 
                            cells['molecules'][ii]['obsCycE'])
                    np.save(TemporaryDirectory("pydream_frac_prolif_quiescence_obscdc20.npy"), 
                            cells['molecules'][ii]['obsCdc20'])
                    
            else : #the course is very short so there can't be a division
                max_tmp_cdc=[]
            
            div_ind = max_tmp_cdc
            
            if cells['time'][ii][-1]<run_time : # if we didn't record the time until the end of the experiment
                if cells['time'][ii][0]>0 : # start at a mitosis (cdc20 peak) and ends at next peak or death (over 

                    if cells['living_cell'][ii]==False :#dead
                        tmp_counts = {'End_of_Expt' : False, 'Death' : True, 'Birthtime' : 
                                      cells['time'][ii][0], 'Lifespan' : (cells['time'][ii][-1] - cells['time'][ii][0])}

                        IMT_eoe_counts = IMT_eoe_counts.append(tmp_counts, ignore_index=True)
                    elif len(div_ind) == 2: #Should divide once or die
                        tmp_counts = {'End_of_Expt' : False, 
                                      'Death' : False, 
                                      'Birthtime' : cells['time'][ii][0], 
                                      'Lifespan' : (cells['time'][ii][div_ind[1]] - cells['time'][ii][div_ind[0]])}
                        if return_all_prot:
#                            print(cells['time'][ii][0])
#                            print(cells['time'][ii][div_ind[0]])
                            cell_info_tmp = np.array(
                                [[index_cell for i in range( div_ind[1] -  div_ind[0])]] + 
                                [(np.array(cells['time'][ii][div_ind[0]: div_ind[1]]) - cells['time'][ii][div_ind[0]]).tolist()] + 
                                [cells['time'][ii][div_ind[0]: div_ind[1]]] + 
                                [cells['molecules'][ii][obs][div_ind[0]: div_ind[1]] for obs in cells['molecules'][0].keys()]).T
                            all_prot = np.append(all_prot, cell_info_tmp , axis=0)
                            index_cell+=1
                            
                        IMT_eoe_counts = IMT_eoe_counts.append(tmp_counts, ignore_index=True)
                else: 
                # start at the beginning of the experiment and ends as a death, a division, or quiescence, 
                # take into account only if two mitoses (cdc20 peaks), otherwise no cell fate registered                    
                    #print(cells['time'][ii][div_ind[0]])
                    if len(div_ind) == 2 :
                        tmp_counts = {
                            'End_of_Expt' : False, 
                            'Death' : False, 
                            'Birthtime' : cells['time'][ii][div_ind[0]], 
                            'Lifespan' : (cells['time'][ii][div_ind[1]] - cells['time'][ii][div_ind[0]])}
                        
                        if return_all_prot:
                            cell_info_tmp = np.array(
                                [[index_cell for i in range( div_ind[1] -  div_ind[0])]] + 
                                [(np.array(cells['time'][ii][div_ind[0]: div_ind[1]]) - cells['time'][ii][div_ind[0]]).tolist()] + 
                                [cells['time'][ii][div_ind[0]: div_ind[1]]] + 
                                [cells['molecules'][ii][obs][div_ind[0]: div_ind[1]] for obs in cells['molecules'][0].keys()]).T
                            all_prot = np.append(all_prot, cell_info_tmp , axis=0)
                            index_cell+=1
                        
                        IMT_eoe_counts = IMT_eoe_counts.append(tmp_counts, ignore_index=True)

            else:    
                if ii <= (init_cell-1) : # initial cells should not include first point...
                    if len(div_ind)>=2 and cells['living_cell'][ii] :
                        
                        for jj in range(len(div_ind)) :
                            if jj<(len(div_ind)-1) :
                                tmp_counts = {
                                    'End_of_Expt' : False, 
                                    'Death' : False, 
                                    'Birthtime' : cells['time'][ii][div_ind[jj]], 
                                    'Lifespan' : (cells['time'][ii][div_ind[jj+1]] - cells['time'][ii][div_ind[jj]])}
                                if return_all_prot:
                                    cell_info_tmp = np.array(
                                        [[index_cell for i in range( div_ind[jj+1] -  div_ind[jj])]] + 
                                        [(np.array(cells['time'][ii][div_ind[jj]: div_ind[jj+1]]) - 
                                          cells['time'][ii][div_ind[jj]]).tolist()] + 
                                        [cells['time'][ii][div_ind[jj]: div_ind[jj+1]]] + 
                                        [cells['molecules'][ii][obs][
                                            div_ind[jj]: div_ind[jj+1]] for obs in cells['molecules'][0].keys()]).T
                                    all_prot = np.append(all_prot, cell_info_tmp , axis=0)
                                    index_cell+=1

                                
                            else :
                                tmp_counts = {'End_of_Expt' : True, 
                                              'Death' : False, 
                                              'Birthtime' : cells['time'][ii][div_ind[jj]], 
                                              'Lifespan' : (run_time - cells['time'][ii][div_ind[jj]])}
                            IMT_eoe_counts = IMT_eoe_counts.append(tmp_counts, ignore_index=True)
            
                    elif len(div_ind)==1 and cells['living_cell'][ii] : # one division
                        tmp_counts = {'End_of_Expt' : True, 
                                      'Death' : False, 
                                      'Birthtime' : cells['time'][ii][div_ind[0]], 
                                      'Lifespan' : (run_time - cells['time'][ii][div_ind[0]])}
                        IMT_eoe_counts = IMT_eoe_counts.append(tmp_counts, ignore_index=True)

                            
                    elif len(div_ind)==0 and cells['living_cell'][ii] :
                        tmp_counts = {'End_of_Expt' : True, 
                                      'Death' : False, 
                                      'Birthtime' : cells['time'][ii][0], 
                                      'Lifespan' : (run_time - cells['time'][ii][0])}
                        IMT_eoe_counts = IMT_eoe_counts.append(tmp_counts, ignore_index=True)

                        
                    elif len(div_ind)>=2 and cells['living_cell'][ii]==False :
                        for jj in range(len(div_ind)) :
                            if jj<(len(div_ind)-1) :
                                tmp_counts = {'End_of_Expt' : False, 
                                              'Death' : False, 
                                              'Birthtime' : cells['time'][ii][div_ind[jj]], 
                                              'Lifespan' : (cells['time'][ii][div_ind[jj+1]] - cells['time'][ii][div_ind[jj]])}
                                if return_all_prot:
                                    cell_info_tmp = np.array(
                                        [[index_cell for i in range( div_ind[jj+1] -  div_ind[jj])]] + 
                                        [(np.array(cells['time'][ii][div_ind[jj]: div_ind[jj+1]]) - 
                                          cells['time'][ii][div_ind[jj]]).tolist()] + 
                                        [cells['time'][ii][div_ind[jj]: div_ind[jj+1]]] + 
                                        [cells['molecules'][ii][obs][
                                            div_ind[jj]: div_ind[jj+1]] for obs in cells['molecules'][0].keys()]).T
                                    all_prot = np.append(all_prot, cell_info_tmp , axis=0)
                                    index_cell+=1
                            else :
                                tmp_counts = {'End_of_Expt' : False, 
                                              'Death' : True, 
                                              'Birthtime' : cells['time'][ii][div_ind[jj]], 
                                              'Lifespan' : (cells['time'][ii][-1] - cells['time'][ii][div_ind[jj]])}
                            IMT_eoe_counts = IMT_eoe_counts.append(tmp_counts, ignore_index=True)                

                    elif len(div_ind)==1 and cells['living_cell'][ii]==False :
                        tmp_counts = {'End_of_Expt' : False, 
                                      'Death' : True, 
                                      'Birthtime' : cells['time'][ii][div_ind[0]], 
                                      'Lifespan' : (cells['time'][ii][-1] - cells['time'][ii][div_ind[0]])}
                        
                        IMT_eoe_counts = IMT_eoe_counts.append(tmp_counts, ignore_index=True)
                else:
                    if len(div_ind)>=1 and cells['living_cell'][ii]:
                        for jj in range(len(div_ind)):
                                   
                            if jj==(len(div_ind)-1):
        
                                tmp_counts = {'End_of_Expt' : True, 
                                              'Death' : False, 
                                              'Birthtime' : cells['time'][ii][div_ind[jj]], 
                                              'Lifespan' : run_time - cells['time'][ii][div_ind[jj]]}
                            
                            else:
                                tmp_counts = {'End_of_Expt' : False, 
                                              'Death' : False, 
                                              'Birthtime' : cells['time'][ii][div_ind[jj]], 
                                              'Lifespan' : (cells['time'][ii][div_ind[jj+1]] - cells['time'][ii][div_ind[jj]])}
                                if return_all_prot:
                                    cell_info_tmp = np.array(
                                        [[index_cell for i in range( div_ind[jj+1] -  div_ind[jj])]] + 
                                        [(np.array(cells['time'][ii][div_ind[jj]: div_ind[jj+1]]) - 
                                          cells['time'][ii][div_ind[jj]]).tolist()] + 
                                        [cells['time'][ii][div_ind[jj]: div_ind[jj+1]]] + 
                                        [cells['molecules'][ii][obs][
                                            div_ind[jj]: div_ind[jj+1]] for obs in cells['molecules'][0].keys()]).T
                                    all_prot = np.append(all_prot, cell_info_tmp , axis=0)
                                    index_cell+=1
                                    
                            IMT_eoe_counts = IMT_eoe_counts.append(tmp_counts, ignore_index=True)
  
                    elif len(div_ind)>=2 and cells['living_cell'][ii]==False:
                        for jj in range(len(div_ind)):

                            if jj==(len(div_ind)-1):
                                tmp_counts = {'End_of_Expt' : False, 
                                              'Death' : True, 
                                              'Birthtime' : cells['time'][ii][div_ind[jj]], 
                                              'Lifespan' : (cells['time'][ii][-1] - cells['time'][ii][div_ind[jj]])}
                            else:
                                tmp_counts = {'End_of_Expt' : False, 
                                              'Death' : False, 
                                              'Birthtime' : cells['time'][ii][div_ind[jj]], 
                                              'Lifespan' : (cells['time'][ii][div_ind[jj+1]] - cells['time'][ii][div_ind[jj]])}
                                if return_all_prot:
                                    cell_info_tmp = np.array(
                                        [[index_cell for i in range( div_ind[jj+1] -  div_ind[jj])]] + 
                                        [(np.array(cells['time'][ii][div_ind[jj]: div_ind[jj+1]]) - 
                                          cells['time'][ii][div_ind[jj]]).tolist()] + 
                                        [cells['time'][ii][div_ind[jj]: div_ind[jj+1]]] + 
                                        [cells['molecules'][ii][obs][
                                            div_ind[jj]: div_ind[jj+1]] for obs in cells['molecules'][0].keys()]).T
                                    all_prot = np.append(all_prot, cell_info_tmp , axis=0)
                                    index_cell+=1
                                    
                            IMT_eoe_counts = IMT_eoe_counts.append(tmp_counts, ignore_index=True)

                    elif len(div_ind)==1 and cells['living_cell'][ii]==False : 
        
                        tmp_counts = {'End_of_Expt' : False, 
                                      'Death' : False, 
                                      'Birthtime' : cells['time'][ii][0], 
                                      'Lifespan' : (cells['time'][ii][div_ind[0]] - cells['time'][ii][0])}
                        if return_all_prot:
                            cell_info_tmp = np.array(
                                [[index_cell for i in range( div_ind[0] -  0)]] + 
                                [(np.array(cells['time'][ii][0: div_ind[0]]) - 
                                  cells['time'][ii][0]).tolist()] + 
                                [cells['time'][ii][0: div_ind[0]]] + 
                                [cells['molecules'][ii][obs][0: div_ind[0]] for obs in cells['molecules'][0].keys()]).T
                            all_prot = np.append(all_prot, cell_info_tmp , axis=0)
                            index_cell+=1
                        
                        
                        IMT_eoe_counts = IMT_eoe_counts.append(tmp_counts, ignore_index=True)

                        tmp_counts = {'End_of_Expt' : False, 
                                      'Death' : True, 
                                      'Birthtime' : cells['time'][ii][div_ind[0]], 
                                      'Lifespan' : (cells['time'][ii][-1] - cells['time'][ii][div_ind[0]])}
                        IMT_eoe_counts = IMT_eoe_counts.append(tmp_counts, ignore_index=True)

    IMT_eoe_counts = IMT_eoe_counts[IMT_eoe_counts['Lifespan']>0]
    
    if return_all_prot:
        return IMT_eoe_counts, pd.DataFrame(all_prot, columns=['index', 'relative time', 'experiment time'] + 
                                            [obs.replace('obs','') for obs in cells['molecules'][0].keys()])
    else :
        return IMT_eoe_counts

