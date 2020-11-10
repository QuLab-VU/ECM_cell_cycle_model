# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 18:19:31 2020

@author: Grim.M
Title: runAKTERK
"""

path="/home/jeremie"

#path=r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code"


import sys
sys.path.insert(0, path + "\cellular_cycle") #add a pathway to the model

from pysb import * 
import numpy as np
import pylab as pl

from runStochCell import runStochcell

from pysb.simulator import BngSimulator
from pysb.simulator import StochKitSimulator

from peakfinder import peakfinder

#######################################################################

def runAKTERK(init, run_time, observed_cells_time, model_concentration, model_population):
    
    
    v = 0.05;     
    sys_size=200; #cell size
    
    variable_range = runStochcell(model_concentration);
    
    if variable_range=={} or np.any(np.array([variable_range[key] for key in variable_range.keys()])<0): #the model doesn't oscillate o the model allows negativ concentrations (possible for the apoptosis module)
        return (0,0)
    
    
    
    cells={}
    cells['molecules']=[]; #Number of molecules taken randomly at one point during the cell cycle
    cells['time']=[]; #
    cells['division']=[]; #
    cells['living_cell']=[]; #living cell?
    cells['recording']=[];
    
    for i in range(init) :
        cells['molecules']=cells['molecules']+[{}]
        a=np.random.randint(len(variable_range['obsSkp2']))
        for var in list(variable_range.keys()) :
            cells['molecules'][i][var] = [round(variable_range[var][a]*200)]
        cells['time']=cells['time'] + [[0]]
        cells['division']=cells['division'] + [0]
        cells['living_cell']=cells['living_cell'] + [True] #living cell = True
        cells['recording']=cells['recording'] + [True] #recording the cell = True (if it isn't dead and if no event didn't occur after observed_cells_time)
        
    cell_count=np.array([[0., len(cells['molecules'])]]);   #start the cell counter at time zero and with number of initial cells
    current_time=np.transpose(np.concatenate([np.zeros((1,len(cells['molecules']))), np.ones((1,len(cells['molecules'])))])); #keeps track of the current time of each cell and whether it's dead or not 
    
    while any(cells['recording']) : #run the code only if there are cells that are recorded
        print('while loop is starting')
        
        cell_count_tmp=np.zeros((len(current_time),2)); #temporary cell count with current cell time and whether a cell was added, taken away, or not
        
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
                    simres = sim.run(initials=Initial, stop_if = 'obsBAXm>45');
                    print("simulation done")
                    
                else :
                    simres = sim.run(verbose=False, stop_if = 'obsBAXm>45');
                    print("simulation done")
                    
                tmp_t = simres.tout[0]
                    
                yout = simres.observables
                
                
                for obs in list(cells['molecules'][ii].keys()) :
                    cells['molecules'][ii][obs]=cells['molecules'][ii][obs][:-1]+yout[obs].tolist(); #evolution of different species
                    
                    
            #if type(cells['time'][ii])!=list : #if it's an array
            #    cells['time'][ii]=cells['time'][ii].tolist()
                
                
                cells['time'][ii] = cells['time'][ii] + [x+cells['time'][ii][-1] for x in tmp_t[1:].tolist()]; #time step
                current_ind_length=len(cells['molecules'][ii]['obsCdc20']);
                
                a, max_tmp_cdc=peakfinder(cells['molecules'][ii]['obsCdc20'],0.15*200,0.15*200);      #cdc20
                a, max_tmp_cycE=peakfinder((cells['molecules'][ii]['obsCycE']),0.15*200,0.15*200);  #cycE
            #min_tmp=peakfinder(cells{ii,2}(:,14),[],30,-1);   #cdc20
                
                if len([i for i,x in enumerate(cells['molecules'][ii]['obsBAXm']) if x > 25])>20 or np.any(np.array(cells['molecules'][ii]['obsBAXm'])>=45) :#dying cell   #baxm >25 then apoptosis for 20 time steps or BAXm>=45 (it means it has started with more than 45 molecules of BAXm)
                    cells['living_cell'][ii]=False; #the cell dies
                    cells['recording'][ii]=False; #the record stops
                    death=[i for i,x in enumerate(cells['molecules'][ii]['obsBAXm']) if x > 25];
                    current_time[ii,0]=cells['time'][ii][death[0]]; #death at time t
                    for obs in list(cells['molecules'][ii].keys()) :
                        cells['molecules'][ii][obs]=cells['molecules'][ii][obs][0:death[0]+1]
                    
                    
                    cells['time'][ii]=cells['time'][ii][0:death[0]+1];
                    
                    current_time[ii,1]=0; #cell is death
                    
                    cell_count_tmp[ii]=[cells['time'][ii][death[0]], -1];              #one cell to be removed to the count
                elif len(max_tmp_cdc) > cells['division'][ii] and len([i for i,x in enumerate(cells['molecules'][ii]['obsBAXm']) if x > 25])<20 : #dividing cell
#                 if isempty(max_tmp_cycE)
#                      current_time(ii)=cells{ii,3}(end);
#                      cell_count_tmp(ii,:)=[cells{ii,3}(end), 0];
                    if len(max_tmp_cdc)>cells['division'][ii] :# && isempty(max_tmp_cycE)==false#make sure cycE has cycled
                        for index in range(cells['division'][ii], len(max_tmp_cdc)) :
    
                            div_ind = (max_tmp_cdc[index]);
                            if len(max_tmp_cdc)<=cells['division'][ii] :
                                current_time[ii,0]=cells['time'][ii][-1];
                                cell_count_tmp[ii,:]=[cells['time'][ii][-1], 0];
                                
                            elif cells['division'][ii] == len(max_tmp_cycE) :
                                current_time[ii,0]=cells['time'][ii][-1];
                                cell_count_tmp[ii,:]=[cells['time'][ii][-1], 0];
                                
                            else :
                                #for obs in list(cells['molecules'][ii].keys()) :
                                    #cells['molecules'][ii][obs]=cells['molecules'][ii][obs][0:div_ind]; #evolution of different species
                                #keep the course for just one cell
                                
                                #cells['time'][ii]=cells['time'][ii][0:div_ind];
                                cells['division'][ii] = cells['division'][ii]+1; #it divides so let keep this
                                
                                current_time[ii,0]=cells['time'][ii][-1]; #no change because this is the time and cell_count_tmp for the cell
                                cell_count_tmp[ii]=[cells['time'][ii][div_ind], 1];
                                
                            #add a new cell : Change to the time and concentrations of the previous cell but at time = div_ind
                            
                                tmp_num=len(cells['time']);
                                
                                if cells['time'][ii][div_ind] < observed_cells_time: #record the new cell only if it appears after the
                                    
                                    initial_new_cell={}
                                    for obs in list(cells['molecules'][ii].keys()) :
                                        initial_new_cell[obs]=[cells['molecules'][ii][obs][div_ind]]; #evolution of different species but at time=div_ind
                            
                                    cells['molecules']=cells['molecules'] + [initial_new_cell]; #to change
                                    cells['time']=cells['time']+[[cells['time'][ii][div_ind]]];
                                    cells['division']=cells['division']+[0];
                                    cells['living_cell']=cells['living_cell']+[True];
                                    cells['recording']=cells['recording']+[True];
    
                                    current_time=np.concatenate((current_time, np.array([[cells['time'][ii][div_ind], 1]])), axis=0); #changed to the real time of the division
                    else : #else ?
                        current_time[ii,0]=cells['time'][ii][-1];
                        cell_count_tmp[ii]=[cells['time'][ii][-1], 0];
                        
                else : #size(max_tmp_cdc,1)==cells{ii,4} && length(find(cells{ii,2}(:,23)>25))<20 #quiescent cell
                    current_time[ii,0]=cells['time'][ii][-1];
                    cell_count_tmp[ii]=[cells['time'][ii][-1], 0];
            
            
            #if the time is > to observed_cells_time, we don't take into account new cells, only the fate of the cells that are already born.
            
            elif current_time[ii,0]<run_time and current_time[ii,0]>=observed_cells_time and cells['recording'][ii]:#if there's time left and the cell has not died/or no fate after observed_cells_time
                print('after')
                T = 4 - (current_time[ii,0]>(run_time-40))*(40 - run_time + current_time[ii,0])
                
                tmp_t = pl.linspace(0, T, int(T*10))
                
                sim = BngSimulator(model_population, tspan=tmp_t, verbose=False);
                if len(model_population.species)>0 :
                    Initial=[]
                    for j in cells['molecules'][ii] :
                        Initial=Initial+[cells['molecules'][ii][j][-1]]
                    simres = sim.run(initials=Initial, stop_if = 'obsBAXm>45');
                    
                else :
                    simres = sim.run(verbose=False, stop_if = 'obsBAXm>45');
                
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
                
                if len([i for i,x in enumerate(cells['molecules'][ii]['obsBAXm']) if x > 25])>20 or np.any(np.array(cells['molecules'][ii]['obsBAXm'])>=45) :#dying cell   #baxm >20 then apoptosis
                    print('death')
                    cells['living_cell'][ii]=False; #the cell dies
                    cells['recording'][ii]=False; #Stop the recording
                    death=[i for i,x in enumerate(cells['molecules'][ii]['obsBAXm']) if x > 25];
                    current_time[ii,0]=cells['time'][ii][death[0]]; #death at time t
                    for obs in list(cells['molecules'][ii].keys()) :
                        cells['molecules'][ii][obs]=cells['molecules'][ii][obs][0:death[0]+1]
                    
                    
                    cells['time'][ii]=cells['time'][ii][0:death[0]+1];
                    
                    current_time[ii,1]=0; #cell is death
                    
                    cell_count_tmp[ii]=[cells['time'][ii][death[0]], -1];              #one cell to be removed to the count
                    
                elif len(max_tmp_cdc)>cells['division'][ii] and len([i for i,x in enumerate(cells['molecules'][ii]['obsBAXm']) if x > 25])<20 : #dividing cell
#                 if isempty(max_tmp_cycE)
#                      current_time(ii)=cells{ii,3}(end);
#                      cell_count_tmp(ii,:)=[cells{ii,3}(end), 0];
                    print('dividing')
                    if len(max_tmp_cdc)>cells['division'][ii] :# && isempty(max_tmp_cycE)==false#make sure cycE has cycled

                        
                        div_ind = max_tmp_cdc[cells['division'][ii]]

                        if cells['division'][ii] == len(max_tmp_cycE) : #Cyc E doesn't oscillate
                            current_time[ii,0]=cells['time'][ii][-1];
                            cell_count_tmp[ii,:]=[cells['time'][ii][-1], 0];
                            if T<40 : #if the cell reached the end of the experiment
                                cells['recording'][ii]=False; #Stop the recording
                        
                            
                        else :

                            #cells['time'][ii]=cells['time'][ii][0:div_ind];
                            cells['division'][ii]=cells['division'][ii]+1; #it divides so let keep this
                            
                            
                            cell_count_tmp[ii]=[cells['time'][ii][-1], 1];
                            current_time[ii,0]=cells['time'][ii][-1];
                            
                            cells['recording'][ii]=False; #Stop the recording
                            #no adding new cell

                    else : #else, oscillations don't work as the system should be ~ quiescence
                        current_time[ii,0]=cells['time'][ii][-1];
                        cell_count_tmp[ii]=[cells['time'][ii][-1], 0];
                        if T<40 : #if the cell reached the end of the experiment
                            cells['recording'][ii]=False; #Stop the recording
                        
                        
                else : #size(max_tmp_cdc,1)==cells{ii,4} && length(find(cells{ii,2}(:,23)>25))<20 #quiescent cell$
                    print('quiescence')
                    current_time[ii,0]=cells['time'][ii][-1];
                    cell_count_tmp[ii]=[cells['time'][ii][-1], 0];
                    if T<40 : #if the cell reached the end of the experiment
                        cells['recording'][ii]=False; #Stop the recording
                        
                    
            else :
                cell_count_tmp[ii]=[cells['time'][ii][-1], 0];
                
            combo = [ii, cells['time'][ii][-1]];
            
        for ii in range(len(cell_count_tmp)) : #updating cell_count (cell population count)
        
            if cell_count_tmp[ii,1]==1 : #new cell
                ins_ind=[i for i,x in enumerate(cell_count[:,0]) if x < cell_count_tmp[ii,0]];
                if ins_ind==[]: #the event happens at time=0
                    current_cell_num=cell_count[0,1];
                    cell_count=np.insert(cell_count, ins_ind[-1]+1, cell_count_tmp[ii,:], axis=0)
                    cell_count[0,1]=current_cell_num;
                    cell_count[0:,1]=cell_count[0+1:,1]+1;
                else :
                    current_cell_num=cell_count[ins_ind[-1],1];
                    cell_count=np.insert(cell_count, ins_ind[-1]+1, cell_count_tmp[ii,:], axis=0)
                    cell_count[ins_ind[-1]+1,1]=current_cell_num;
                    cell_count[ins_ind[-1]+1:,1]=cell_count[ins_ind[-1]+1:,1]+1;
            
            
            
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
    
    return (cells, cell_count)
