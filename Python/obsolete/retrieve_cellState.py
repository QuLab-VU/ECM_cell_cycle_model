# -*- coding: utf-8 -*-
"""
@author: Jérémie Mahrlens
@author: Darren Tyson

Created by Jérémie and edited by Darren
"""

from peakfinder import peakfinder
import pandas as pd
import numpy as np
from tempfile import TemporaryDirectory

def retrieve_cellState(run_time, cells, return_all_prot=False, cells_starting_at_mitosis=True, thr = 0.4, omega=200, verbose=False):
    if verbose:
    	print(f'retieving cell state of {len(cells["living_cell"])} cells')
    all_prot = np.empty((0,len(cells['molecules'][0])+6)) #cell ind, time, relative time, protein time courses
    index_cell=0 #cell identification nber
    IMT_eoe_counts = pd.DataFrame(columns=['Starting cell', 'End_of_Expt', 'Death' ,'Birthtime','Lifespan'])
    
    for ii in range(len(cells['living_cell'])):

        #depends on when did we stop to record
        end_index=[]
        end_index = [i for i, x in enumerate(np.array(cells['time'][ii]) < run_time) if x] #time within the experiment
    
        if len(end_index)>0 :
            if len(cells['molecules'][ii]['obsCdc20'])>20:
                try:
                    # choose threshold values
                    
                    (a,max_tmp_cdc) = peakfinder(cells['molecules'][ii]['obsCdc20'], thr*omega, thr*omega, 1, includeEndpoints=1); #peakfinder parameters, can be changed
                    (a,max_tmp_cycE) = peakfinder(cells['molecules'][ii]['obsCycE'], thr*omega, thr*omega, 1, includeEndpoints=1);
                          
                    max_tmp_cdc = [tmp for tmp in max_tmp_cdc if tmp<=end_index[-1]]
                    max_tmp_cycE = [tmp for tmp in max_tmp_cycE if tmp<=end_index[-1]]
                    
                except:
                    np.save(TemporaryDirectory("pydream_frac_prolif_quiescence_obscycE.npy"), 
                            cells['molecules'][ii]['obsCycE'])
                    np.save(TemporaryDirectory("pydream_frac_prolif_quiescence_obscdc20.npy"), 
                            cells['molecules'][ii]['obsCdc20'])
                    
            else : #the course is very short so there can't be a division
                max_tmp_cdc = []
            
            div_ind = max_tmp_cdc
            
            if cells['time'][ii][-1]<run_time : # if the simulation didn't reach the end of the experiment (i.e. it ends with mitosis or death)
                if cells['time'][ii][0]==0 and len(div_ind) > 0 and div_ind[0]>0: # start at a mitosis (cdc20 peak) and ends at next peak or death 
                    tmp_counts = {'Starting cell' : True, 
                                  'End_of_Expt' : False, 
                                  'Death' : False, 
                                  'Birthtime' : cells['time'][ii][0], 
                                  'Lifespan' : (cells['time'][ii][div_ind[0]] - cells['time'][ii][0])} # add first event = from start of experiment to first mitosis
                    IMT_eoe_counts = pd.concat([IMT_eoe_counts, pd.DataFrame(tmp_counts, index=[0])], ignore_index=True)
                    if return_all_prot:
                        cell_info_tmp = np.array(
                           [[True for i in range(div_ind[0] -  0)]] +
                            [[False for i in range(div_ind[0] -  0)]] +
                            [[False for i in range(div_ind[0] -  0)]] +
                            [[index_cell for i in range(div_ind[0] -  0)]] + 
                            [(np.array(cells['time'][ii][0: div_ind[0]]) - cells['time'][ii][0]).tolist()] + 
                            [cells['time'][ii][0: div_ind[0]]] + 
                            [cells['molecules'][ii][obs][0: div_ind[0]] for obs in cells['molecules'][0].keys()]).T
                        all_prot = np.append(all_prot, cell_info_tmp , axis=0)
                        index_cell+=1
                
                for jj in range(len(div_ind)-1) :
                    #print(jj)
                    tmp_counts = {
                    'Starting cell' : False, #cells['time'][ii][div_ind[jj]]==0
                    'End_of_Expt' : False, 
                    'Death' : False, 
                    'Birthtime' : cells['time'][ii][div_ind[jj]], 
                    'Lifespan' : (cells['time'][ii][div_ind[jj+1]] - cells['time'][ii][div_ind[jj]])}
                    IMT_eoe_counts = pd.concat([IMT_eoe_counts, pd.DataFrame(tmp_counts, index=[0])], ignore_index=True)
                    
                    if return_all_prot:
                        cell_info_tmp = np.array(
                            [[True for i in range(div_ind[jj+1] -  div_ind[jj])]] + 
                            [[False for i in range(div_ind[jj+1] -  div_ind[jj])]] + 
                            [[False for i in range(div_ind[jj+1] -  div_ind[jj])]] + 
                            [[index_cell for i in range(div_ind[jj+1] -  div_ind[jj])]] + 
                            [(np.array(cells['time'][ii][div_ind[jj]: div_ind[jj+1]]) - cells['time'][ii][div_ind[jj]]).tolist()] + 
                            [cells['time'][ii][div_ind[jj]: div_ind[jj+1]]] + 
                            [cells['molecules'][ii][obs][div_ind[jj]: div_ind[jj+1]] for obs in cells['molecules'][0].keys()]).T
                        all_prot = np.append(all_prot, cell_info_tmp , axis=0)
                        index_cell+=1
                
                if cells['living_cell'][ii]==False :# ends with death
                    start = 0 if len(div_ind)==0 else div_ind[-1]
                    tmp_counts = {'Starting cell' : False, 
                                  'End_of_Expt' : False, 
                                  'Death' : True, 
                                  'Birthtime' : cells['time'][ii][start], 
                                  'Lifespan' : (cells['time'][ii][-1] - cells['time'][ii][start])}
                    IMT_eoe_counts = pd.concat([IMT_eoe_counts, pd.DataFrame(tmp_counts, index=[0])], ignore_index=True)
                    
                    if return_all_prot:
                        cell_info_tmp = np.array(
                            [[False for i in range(end_index[-1] -  start + 1)]] + 
                            [[False for i in range(end_index[-1] -  start + 1)]] + 
                            [[True for i in range(end_index[-1] -  start + 1)]] + 
                            [[index_cell for i in range(end_index[-1] - start + 1)]] + 
                            [(np.array(cells['time'][ii][start:]) - cells['time'][ii][start]).tolist()] + 
                            [cells['time'][ii][start:]] + 
                            [cells['molecules'][ii][obs][start:] for obs in cells['molecules'][0].keys()]).T
                        all_prot = np.append(all_prot, cell_info_tmp , axis=0)
                        index_cell+=1

            else: # record until the end of experiment
                if cells['time'][ii][0]==0 and (len(div_ind)==0 or (div_ind[0]==0 and len(div_ind)==1)): # No mitosis at all during the whole experiment
                    tmp_counts = {'Starting cell' : True, 
                                  'End_of_Expt' : True, 
                                  'Death' : False, 
                                  'Birthtime' : cells['time'][ii][0], 
                                  'Lifespan' : run_time} # add first event = from start of experiment to first mitosis
                    IMT_eoe_counts = pd.concat([IMT_eoe_counts, pd.DataFrame(tmp_counts, index=[0])], ignore_index=True)
                    
                else:
                    if cells['time'][ii][0]==0 and div_ind[0]>0: # doesn't start at a mitosis (cdc20 peak) and ends at next peak or death 
                        tmp_counts = {'Starting cell' : True, 
                        			  'End_of_Expt' : False, 
                        			  'Death' : False, 
                        			  'Birthtime' : cells['time'][ii][0], 
                        			  'Lifespan' : (cells['time'][ii][div_ind[0]] - cells['time'][ii][0])} # add first event = from start of experiment to first mitosis
                        IMT_eoe_counts = pd.concat([IMT_eoe_counts, pd.DataFrame(tmp_counts, index=[0])], ignore_index=True)
                        
                    
                    for jj in range(len(div_ind)-1) :
                        
                        tmp_counts = {'Starting cell' : False,
                                      'End_of_Expt' : False, 
                                      'Death' : False, 
                                      'Birthtime' : cells['time'][ii][div_ind[jj]], 
                                      'Lifespan' : (cells['time'][ii][div_ind[jj+1]] - cells['time'][ii][div_ind[jj]])}
                        IMT_eoe_counts = pd.concat([IMT_eoe_counts, pd.DataFrame(tmp_counts, index=[0])], ignore_index=True)
                        if return_all_prot:
                            cell_info_tmp = np.array(
                                [[True if jj==0 else False for i in range(div_ind[jj+1] -  div_ind[jj])]] +
                                [[False for i in range(div_ind[jj+1] -  div_ind[jj])]] + 
                                [[False for i in range(div_ind[jj+1] -  div_ind[jj])]] + 
                                [[index_cell for i in range(div_ind[jj+1] -  div_ind[jj])]] + 
                                [(np.array(cells['time'][ii][div_ind[jj]: div_ind[jj+1]]) - cells['time'][ii][div_ind[jj]]).tolist()] + 
                                [cells['time'][ii][div_ind[jj]: div_ind[jj+1]]] + 
                                [cells['molecules'][ii][obs][div_ind[jj]: div_ind[jj+1]] for obs in cells['molecules'][0].keys()]).T
                            all_prot = np.append(all_prot, cell_info_tmp , axis=0)
                            index_cell+=1                        
                            
                    if cells['living_cell'][ii]==True :# end of experiment
                        start = 0 if len(div_ind)==0 else div_ind[-1]
                        tmp_counts = {'Starting cell' : False, 
                                      'End_of_Expt' : True, 
                                      'Death' : False, 
                                      'Birthtime' : cells['time'][ii][start], 
                                      'Lifespan' : (cells['time'][ii][end_index[-1]] - cells['time'][ii][start])}
                        IMT_eoe_counts = pd.concat([IMT_eoe_counts, pd.DataFrame(tmp_counts, index=[0])], ignore_index=True)
                        if return_all_prot:
                            cell_info_tmp = np.array(
                                [[False for i in range(end_index[-1] -  start)]] +
                                [[True for i in range(end_index[-1] -  start)]] + 
                                [[False for i in range(end_index[-1] -  start)]] + 
                                [[index_cell for i in range(end_index[-1] -  start)]] + 
                                [(np.array(cells['time'][ii][start:end_index[-1]]) - cells['time'][ii][start]).tolist()] + 
                                [cells['time'][ii][start:end_index[-1]]] + 
                                [cells['molecules'][ii][obs][start:end_index[-1]] for obs in cells['molecules'][0].keys()]).T
                            all_prot = np.append(all_prot, cell_info_tmp , axis=0)
                            index_cell+=1
                            
                    else: #ends with apoptosis
                        start = 0 if len(div_ind)==0 else div_ind[-1]
                        tmp_counts = {'Starting cell' : False, 'End_of_Expt' : False, 'Death' : True, 'Birthtime' : cells['time'][ii][start], 'Lifespan' : (cells['time'][ii][end_index[-1]] - cells['time'][ii][start])}
                        IMT_eoe_counts = pd.concat([IMT_eoe_counts, pd.DataFrame(tmp_counts, index=[0])], ignore_index=True)
                        if return_all_prot:
                            cell_info_tmp = np.array(
                                [[False for i in range(end_index[-1] - start)]] +
                                [[False for i in range(end_index[-1] - start)]] + 
                                [[True for i in range(end_index[-1] - start)]] + 
                                [[index_cell for i in range(end_index[-1] - start)]] + 
                                [(np.array(cells['time'][ii][start:end_index[-1]]) - cells['time'][ii][start]).tolist()] + 
                                [cells['time'][ii][start:end_index[-1]]] + 
                                [cells['molecules'][ii][obs][start:end_index[-1]] for obs in cells['molecules'][0].keys()]).T
                            all_prot = np.append(all_prot, cell_info_tmp , axis=0)
                            index_cell+=1
    
    
    IMT_eoe_counts = IMT_eoe_counts[IMT_eoe_counts['Lifespan']>0]
    
    if return_all_prot:
        return IMT_eoe_counts, pd.DataFrame(all_prot, 
                                            columns=['Starting cell', 'End of experiment', 'Apoptosis', 'index', 'relative time', 'experiment time'] + 
                                            [obs.replace('obs','') for obs in cells['molecules'][0].keys()])
    else :
        return IMT_eoe_counts


if __name__=="__main__":
    run_time = 120
    IMT_eoe_counts = pd.DataFrame(columns=['Starting cell', 'End_of_Expt', 'Death' ,'Birthtime','Lifespan'])
    tmp_counts = {'Starting cell' : True, 'End_of_Expt' : False, 'Death' : False, 'Birthtime' : 0, 'Lifespan' : 10}
    IMT_eoe_counts = pd.concat([IMT_eoe_counts, pd.DataFrame(tmp_counts, index=[0])], ignore_index=True)
    tmp_counts = {'Starting cell' : False, 'End_of_Expt' : False, 'Death' : False, 'Birthtime' : 10, 'Lifespan' : 10}
    IMT_eoe_counts = pd.concat([IMT_eoe_counts, pd.DataFrame(tmp_counts, index=[0])], ignore_index=True)

    print(IMT_eoe_counts)


