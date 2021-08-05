# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 18:19:31 2020

@author: Grim.M
Title: runAKTERK
"""

path="/home/jeremie"

#path=r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code"


import sys
#sys.path.insert(0, path + "\cellular_cycle") #add a pathway to the model

sys.path.insert(0, __file__ + "/../cellular_cycle")
from pysb import * 
import numpy as np
import pylab as pl

from runStochCell import runStochcell

from pysb.simulator import BngSimulator
from pysb.simulator import StochKitSimulator

from peakfinder import peakfinder

omega=200

#######################################################################

def runAKTERK(cells, init, run_time, observed_cells_time, model_concentration, model_population):
    
    
    
    pysb_sampled_parameter_names = ["V_s2skp","V_Fp","k_M","k_EMY","k_CDF","k_CDS","k_ED","K_S","K_MY","K_DF","d_CD","k_P1", "k_P2", "K_CD", "k_E", "K_M", "K_CE", "V_2E2F", "K_2E2F", "k_RE", "kf1", "kf2", "kb", "kasbaxmbcl2", "kdsbaxmbcl2", "k1sbh3", "k2sbh3", "kdbh3", "kasbh3bcl2", "kdsbh3bcl2"]
    parameter_vector_original=[]
    for pname in pysb_sampled_parameter_names:
        #Change model parameter values to current location in parameter space
        parameter_vector_original+=[model_population.parameters[pname].value]
    
    v = 0.05;     
    sys_size=200; #cell size
    
    
    
    cell_count=np.array([[0., len(cells['molecules'])]]);   #start the cell counter at time zero and with number of initial cells
    current_time=np.transpose(np.concatenate([np.zeros((1,len(cells['molecules']))), np.ones((1,len(cells['molecules'])))])); #keeps track of the current time of each cell and whether it's dead or not 
    
    while any(cells['recording']) : #run the code only if there are cells that are recorded
        print('while loop is starting')
        
        cell_count_tmp=np.zeros((1,2)); #temporary cell count with current cell time and whether a cell was added, taken away, or not
        
        print(cell_count_tmp)
        for ii in range(len(cells['molecules'])) :
            print(current_time)
            print("new cell")
            
            #If the time is between 0 to observed_cells_time
            
            if current_time[ii,0]<observed_cells_time and cells['living_cell'][ii] :#if there's time left and the cell has not died
                print('before')
                
                T = 20
                
                tmp_t = pl.linspace(0, T, int(T*10))
                
                sim = BngSimulator(model_population, tspan=tmp_t, verbose=False);
                
                if len(model_population.species)>0 :
                    Initial=[]
                    for j in cells['molecules'][ii] :
                        Initial=Initial+[cells['molecules'][ii][j][-1]]
                    simres = sim.run(initials=Initial, stop_if = 'obsBH3>200');
                    print(parameter_vector_original)
                    print("simulation done")
                    
                else :
                    simres = sim.run(verbose=False, stop_if = 'obsBH3>200');
                    print(parameter_vector_original)
                    print("simulation done")
                    
                tmp_t = simres.tout[0]
                    
                yout = simres.observables
                
                
                for obs in list(cells['molecules'][ii].keys()) :
                    cells['molecules'][ii][obs]=cells['molecules'][ii][obs][:-1]+yout[obs].tolist(); #evolution of different species
                    
                    
            #if type(cells['time'][ii])!=list : #if it's an array
            #    cells['time'][ii]=cells['time'][ii].tolist()
                
                
                cells['time'][ii] = cells['time'][ii] + [x+cells['time'][ii][-1] for x in tmp_t[1:].tolist()]; #time step
                current_ind_length=len(cells['molecules'][ii]['obsCdc20']);
                
                a, max_tmp_cdc = peakfinder(cells['molecules'][ii]['obsCdc20'], 0.15*omega, 0.15*omega, 1, includeEndpoints=1);      #cdc20
                a, max_tmp_cycE = peakfinder((cells['molecules'][ii]['obsCycE']), 0.15*omega, 0.15*omega, 1, includeEndpoints=1);  #cycE
            #min_tmp=peakfinder(cells{ii,2}(:,14),[],30,-1);   #cdc20
                
                if len([i for i,x in enumerate(cells['molecules'][ii]['obsBAXm']) if x > 25])>20 or np.any(np.array(cells['molecules'][ii]['obsBAXm'])>=40) or np.any(np.array(cells['molecules'][ii]['obsBH3'])>=200) :#dying cell   #baxm >25 then apoptosis for 20 time steps or BAXm>=40 (it means it has started with more than 40 molecules of BAXm)
                    cells['living_cell'][ii]=False; #the cell dies
                    cells['recording'][ii]=False; #the record stops
                    death=[i for i,x in enumerate(cells['molecules'][ii]['obsBAXm']) if x > 20];
                    current_time[ii,0]=cells['time'][ii][death[0]]; #death at time t
                    
                    for obs in list(cells['molecules'][ii].keys()) :
                        cells['molecules'][ii][obs]=cells['molecules'][ii][obs][0:death[0]+1]
                    
                    
                    cells['time'][ii]=cells['time'][ii][0:death[0]+1];
                    
                    current_time[ii,1]=0; #cell is death
                    
                    
                    cell_count_tmp = np.insert(cell_count_tmp, len(cell_count_tmp), [cells['time'][ii][death[0]], -1], axis=0); #one cell to be removed to the count
                    
                elif len(max_tmp_cdc) > cells['division'][ii] and len([i for i,x in enumerate(cells['molecules'][ii]['obsBAXm']) if x > 25])<20 : #dividing cell

                    for index in range(cells['division'][ii], len(max_tmp_cdc)) :

                        div_ind = (max_tmp_cdc[index]);
                        print('index')
                        
                        if index==0 : #first cell, there is no need to check if CycE has cycled (it can have cycled before recording), it will be checked between two mitosis
                            cells['division'][ii] = cells['division'][ii]+1;               
                            current_time[ii,0]=cells['time'][ii][-1]; #no change because this is the time and cell_count_tmp for the cell
                            cell_count_tmp = np.insert(cell_count_tmp, len(cell_count_tmp), [cells['time'][ii][div_ind], 1], axis=0);
                        #add a new cell : Change to the time and concentrations of the previous cell but at time = div_ind
                            tmp_num=len(cells['time']);
                            
                            if cells['time'][ii][div_ind] < observed_cells_time : #record the new cell only if it appears befor time of observation
                                
                                initial_new_cell={}
                                for obs in list(cells['molecules'][ii].keys()) :
                                    initial_new_cell[obs]=[cells['molecules'][ii][obs][div_ind]]; #evolution of different species but at time=div_ind
                        
                                cells['molecules'] = cells['molecules'] + [initial_new_cell]; #to change
                                cells['time'] = cells['time']+[[cells['time'][ii][div_ind]]];
                                cells['division'] = cells['division']+[1];
                                cells['living_cell'] =cells['living_cell']+[True];
                                
                                cells['recording']=cells['recording']+[True];
                                
                                current_time = np.concatenate((current_time, np.array([[cells['time'][ii][div_ind], 1]])), axis=0); #changed to the real time of the division
                            
                            
                        else :                                
                            if np.any(np.logical_and(max_tmp_cdc[index-1] < np.array(max_tmp_cycE) ,np.array(max_tmp_cycE) < max_tmp_cdc[index])) : #just checking if Cyc2 is actully oscillating between Cdc20
                            
                                #cells['time'][ii]=cells['time'][ii][0:div_ind];
                                cells['division'][ii] = cells['division'][ii]+1; #it divides so let keep this
                                
                                current_time[ii,0]=cells['time'][ii][-1]; #no change because this is the time and cell_count_tmp for the cell
                                
                                cell_count_tmp = np.insert(cell_count_tmp, len(cell_count_tmp), [cells['time'][ii][div_ind], 1], axis=0);
                                    
                            #add a new cell : Change to the time and concentrations of the previous cell but at time = div_ind
                            
                                tmp_num=len(cells['time']);
                                
                                if cells['time'][ii][div_ind] < observed_cells_time: #record the new cell only if it appears after the
                                    
                                    initial_new_cell={}
                                    for obs in list(cells['molecules'][ii].keys()) :
                                        initial_new_cell[obs]=[cells['molecules'][ii][obs][div_ind]]; #evolution of different species but at time=div_ind
                            
                                    cells['molecules'] = cells['molecules'] + [initial_new_cell]; #to change
                                    cells['time'] = cells['time']+[[cells['time'][ii][div_ind]]];
                                    cells['division'] = cells['division']+[1];
                                    cells['living_cell'] = cells['living_cell']+[True];
                                    
                                    if cells['time'][ii][div_ind]==run_time:
                                        cells['recording']=cells['recording']+[True];
                                    else : #mitosis happens right on the end time (almost impossible)
                                        cells['recording']=cells['recording']+[False];
                                        
                                    current_time=np.concatenate((current_time, np.array([[cells['time'][ii][div_ind], 1]])), axis=0); #changed to the real time of the division


                else : #size(max_tmp_cdc,1)==cells{ii,4} && length(find(cells{ii,2}(:,23)>25))<20 #quiescent cell
                    current_time[ii,0]=cells['time'][ii][-1];
                    cell_count_tmp = np.insert(cell_count_tmp, len(cell_count_tmp), [cells['time'][ii][-1], 0], axis=0);
            
            if current_time[ii,0]>=run_time : #if the cell reached the end of the experiment
                cells['recording'][ii]=False; #Stop the recording
                
                
            #if the time is > to observed_cells_time, we don't take into account new cells, only the fate of the cells that are already born.
            
            elif current_time[ii,0]<run_time and current_time[ii,0]>=observed_cells_time and cells['recording'][ii]:#if there's time left and the cell has not died/or no fate after observed_cells_time
                print('after')
                T = 41 - (current_time[ii,0]>(run_time-40))*(40 - run_time + current_time[ii,0])
                
                tmp_t = pl.linspace(0, T, int(T*10))
                
                sim = BngSimulator(model_population, tspan=tmp_t, verbose=False);
                if len(model_population.species)>0 :
                    Initial=[]
                    for j in cells['molecules'][ii] :
                        Initial=Initial+[cells['molecules'][ii][j][-1]]
                    simres = sim.run(initials=Initial, stop_if = 'obsBH3>200');
                    print(parameter_vector_original)
                    
                else :
                    simres = sim.run(verbose=False, stop_if = 'obsBH3>200');
                    print(parameter_vector_original)
                
                tmp_t = simres.tout[0]
                yout = simres.observables
                
                for obs in list(cells['molecules'][ii].keys()) :
                    cells['molecules'][ii][obs]=cells['molecules'][ii][obs][:-1]+yout[obs].tolist(); #evolution of different species
                    
                    
            #if type(cells['time'][ii])!=list : #if it's an array
            #    cells['time'][ii]=cells['time'][ii].tolist()
                
                
                cells['time'][ii] = cells['time'][ii] + [x+cells['time'][ii][-1] for x in tmp_t[1:].tolist()]; #time step
                current_ind_length=len(cells['molecules'][ii]['obsCdc20']);
                
                a, max_tmp_cdc = peakfinder(cells['molecules'][ii]['obsCdc20'],0.15*200,0.15*200);      #cdc20
                a, max_tmp_cycE = peakfinder(cells['molecules'][ii]['obsCycE'],0.15*200,0.15*200);  #cycE
            #min_tmp=peakfinder(cells{ii,2}(:,14),[],30,-1);   #cdc20
                
                if len([i for i,x in enumerate(cells['molecules'][ii]['obsBAXm']) if x > 25])>20 or np.any(np.array(cells['molecules'][ii]['obsBAXm'])>=40) or np.any(np.array(cells['molecules'][ii]['obsBH3'])>=200) :#dying cell   #baxm >20 then apoptosis
                    print('death')
                    cells['living_cell'][ii]=False; #the cell dies
                    cells['recording'][ii]=False; #Stop the recording
                    death=[i for i,x in enumerate(cells['molecules'][ii]['obsBAXm']) if x > 20];
                    current_time[ii,0]=cells['time'][ii][death[0]]; #death at time t
                    for obs in list(cells['molecules'][ii].keys()) :
                        cells['molecules'][ii][obs]=cells['molecules'][ii][obs][0:death[0]+1]
                    
                    
                    cells['time'][ii]=cells['time'][ii][0:death[0]+1];
                    
                    current_time[ii,1]=0; #cell is death
                    
                    cell_count_tmp = np.insert(cell_count_tmp, len(cell_count_tmp), [cells['time'][ii][death[0]], -1], axis=0); #one cell to be removed to the count
                    
                elif len(max_tmp_cdc) > cells['division'][ii] and len([i for i,x in enumerate(cells['molecules'][ii]['obsBAXm']) if x > 25])<20 : #dividing cell
                    
                    
                    for index in range(cells['division'][ii]-1, len(max_tmp_cdc)-1) :
                        cells['recording'][ii]=False

                        div_ind = (max_tmp_cdc[index]);
                        print('index')
                        
                        if index==0: #first cell, there is no need to check id CycE has cycled, it will be checked between two mitosis
                            cells['division'][ii] = cells['division'][ii]+1;               
                            current_time[ii,0]=cells['time'][ii][-1]; #no change because this is the time and cell_count_tmp for the cell
                            cell_count_tmp = np.insert(cell_count_tmp, len(cell_count_tmp), [cells['time'][ii][div_ind], 1], axis=0);

                        #add a new cell : Change to the time and concentrations of the previous cell but at time = div_ind
                            tmp_num=len(cells['time']);
                            if cells['time'][ii][div_ind] < observed_cells_time: #record the new cell only if it appears after the
                                initial_new_cell={}
                                for obs in list(cells['molecules'][ii].keys()) :
                                    initial_new_cell[obs]=[cells['molecules'][ii][obs][div_ind]]; #evolution of different species but at time=div_ind
                        
                                cells['molecules']=cells['molecules'] + [initial_new_cell]; #to change
                                cells['time']=cells['time']+[[cells['time'][ii][div_ind]]];
                                cells['division']=cells['division']+[1];
                                cells['living_cell']=cells['living_cell']+[True];
                                
                                cells['recording']=cells['recording']+[False];
                                
                                current_time=np.concatenate((current_time, np.array([[cells['time'][ii][div_ind], 1]])), axis=0); #changed to the real time of the division
                            
                        else :                                
                            if np.any(np.logical_and(max_tmp_cdc[index-1] < np.array(max_tmp_cycE) ,np.array(max_tmp_cycE) < max_tmp_cdc[index])) : #if cycE oscillates between two CDC20 peaks
                            
                                #cells['time'][ii]=cells['time'][ii][0:div_ind];
                                cells['division'][ii] = cells['division'][ii]+1; #it divides so let keep this
                                
                                current_time[ii,0]=cells['time'][ii][-1]; #no change because this is the time and cell_count_tmp for the cell
                                
                                cell_count_tmp = np.insert(cell_count_tmp, len(cell_count_tmp), [cells['time'][ii][div_ind], 1], axis=0);

                                
                            #add a new cell : Change to the time and concentrations of the previous cell but at time = div_ind
                            
                                tmp_num=len(cells['time']);
                                
                                if cells['time'][ii][div_ind] < observed_cells_time: #record the new cell only if it appears after the
                                    
                                    initial_new_cell={}
                                    for obs in list(cells['molecules'][ii].keys()) :
                                        initial_new_cell[obs]=[cells['molecules'][ii][obs][div_ind]]; #evolution of different species but at time=div_ind
                            
                                    cells['molecules']=cells['molecules'] + [initial_new_cell]; #to change
                                    cells['time']=cells['time']+[[cells['time'][ii][div_ind]]];
                                    cells['division']=cells['division']+[1];
                                    cells['living_cell']=cells['living_cell']+[True];
                                    cells['recording']=cells['recording']+[False];
                                        
                                    current_time=np.concatenate((current_time, np.array([[cells['time'][ii][div_ind], 1]])), axis=0); #changed to the real time of the division
                    
                    
                else : #size(max_tmp_cdc,1)==cells{ii,4} && length(find(cells{ii,2}(:,23)>25))<20 #quiescent cell$
                    print('quiescence')
                    current_time[ii,0]=cells['time'][ii][-1];
                    cell_count_tmp = np.insert(cell_count_tmp, len(cell_count_tmp), [cells['time'][ii][-1], 0], axis=0);
                    if current_time[ii,0]>=run_time : #if the cell reached the end of the experiment
                        cells['recording'][ii]=False; #Stop the recording
                        
                    
            else :
                cell_count_tmp = np.insert(cell_count_tmp, len(cell_count_tmp), [cells['time'][ii][-1], 0], axis=0);
                
            combo = [ii, cells['time'][ii][-1]];
            
        try :
            for ii in range(len(cell_count_tmp)) : #updating cell_count (cell population count)
            
                if cell_count_tmp[ii,1]==1 : #new cell
                    ins_ind=[i for i,x in enumerate(cell_count[:,0]) if x < cell_count_tmp[ii,0]];
                    
                    if ins_ind==[]: #the event happens at time=0
                        
                        current_cell_num=cell_count[0,1];
                        #cell_count=np.insert(cell_count, 0, cell_count_tmp[ii,:], axis=0)
                        cell_count[0,1]=current_cell_num;
                        cell_count[0:,1]=cell_count[0:,1]+1;
                        
                    else :
                        current_cell_num = cell_count[ins_ind[-1],1];
                        cell_count = np.insert(cell_count, ins_ind[-1]+1, cell_count_tmp[ii,:], axis=0)
                        cell_count[ins_ind[-1]+1,1] = current_cell_num;
                        cell_count[ins_ind[-1]+1:,1] = cell_count[ins_ind[-1]+1:,1]+1;
                
                
                elif cell_count_tmp[ii,1]==-1 : #apoptosis
                    ins_ind=[i for i,x in enumerate(cell_count[:,0]) if x < cell_count_tmp[ii,0]]
                    if ins_ind==[] : #the event happens at time=0
                        current_cell_num = cell_count[0,1]
                        cell_count=np.insert(cell_count, 0+1, cell_count_tmp[ii,:], axis=0)
                        cell_count[0+1,1]=current_cell_num;
                        cell_count[0+1:,1]=cell_count[0+1:,1]-1; 
                     
                    else :
                        current_cell_num=cell_count[ins_ind[-1],1];
                        cell_count=np.insert(cell_count, ins_ind[-1]+1, cell_count_tmp[ii,:], axis=0)
                        cell_count[ins_ind[-1]+1,1]=current_cell_num;
                        cell_count[ins_ind[-1]+1:,1]=cell_count[ins_ind[-1]+1:,1]-1;
                    
                elif cell_count_tmp[ii,1]==0 : #quiescnce
                    ins_ind=[i for i,x in enumerate(cell_count[:,0]) if x < cell_count_tmp[ii,0]]
                    
                    if ins_ind==[] :
                        current_cell_num = cell_count[0,1]
                        cell_count=np.insert(cell_count, 0+1, cell_count_tmp[ii,:], axis=0)
                        cell_count[0+1,1]=current_cell_num;
                        
                    else :
                        current_cell_num=cell_count[ins_ind[-1],1];
                        cell_count=np.insert(cell_count, ins_ind[-1]+1, cell_count_tmp[ii,:], axis=0)
                        cell_count[ins_ind[-1]+1,1]=current_cell_num;
                    
            u, indices = np.unique(cell_count[:,0], return_index=True)
            cell_count=cell_count[indices,:];
        #disp(cell_count(end,:)) 
        except:
            np.save(path + "/fit_data/pydream/frac_prolif_quiescence_DMSO_only_30_parameters/cell_count_tmp.npy", cell_count_tmp)
            np.save(path + "/fit_data/pydream/frac_prolif_quiescence_DMSO_only_30_parameters/cell_count_tmp.npy", cell_count)
            np.save(path + "/fit_data/pydream/frac_prolif_quiescence_DMSO_only_30_parameters/ii.npy", ii)
            
    return (cells, cell_count)
