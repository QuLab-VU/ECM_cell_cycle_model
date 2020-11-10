# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 18:19:31 2020

@author: Grim.M
Title: runAKTERK
"""

path=r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code"
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

def runAKTERK(init, run_time, model_concentration, model_population):
    tmp_t = pl.linspace(0, 10, 1000)
    sim = BngSimulator(model_population, tspan=tmp_t, verbose=False);
    simres = sim.run();
    v = 0.05;     
    sys_size=200; #cell size
    
    variable_range = runStochcell(model_concentration);
    
    cells={}
    cells['molecules']=[]; #Number of molecules taken randomly at one point during the cell cycle
    cells['time']=[]; #
    cells['division']=[]; #
    cells['living_cell']=[]; #living cell?
    
    
    for i in range(init) :
        cells['molecules']=cells['molecules']+[{}]
        a=np.random.randint(len(variable_range['obsSkp2']))
        for var in list(variable_range.keys()) :
            cells['molecules'][i][var] = [round(variable_range[var][a]*200)]
        cells['time']=cells['time'] + [[0]]
        cells['division']=cells['division'] + [0]
        cells['living_cell']=cells['living_cell'] + [True] #living cell = True
        
    cell_count=np.array([[0., len(cells['molecules'])]]);   #start the cell counter at time zero and with number of initial cells
    current_time=np.transpose(np.concatenate([np.zeros((1,len(cells['molecules']))), np.ones((1,len(cells['molecules'])))])); #keeps track of the current time of each cell and whether it's dead or not 
    while any([a and b for (a,b) in zip(list(map(bool,current_time[:,1])), (current_time[:,0]<run_time).tolist())]) : #run the code only if there are live cells under the run time
        cell_count_tmp=np.zeros((len(current_time),2)); #temporary cell count with current cell time and whether a cell was added, taken away, or not
        for ii in range(len(cells['molecules'])) :
            print(current_time)
            if current_time[ii,0]<run_time and cells['living_cell'][ii] :#if there's time left and the cell has not died
                tmp_t = pl.linspace(0, 10, 1000)
                
                sim = BngSimulator(model_population, tspan=tmp_t, verbose=False);
                if len(model_population.species)>0 :
                    Initial=[]
                    for j in cells['molecules'][ii] :
                        Initial=Initial+[cells['molecules'][ii][j][-1]]
                    simres = sim.run(initials=Initial);
                    
                else :
                    simres = sim.run(verbose=False);
                
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
                
                if len([i for i,x in enumerate(cells['molecules'][ii]['obsBAXm']) if x > 25])>20 :#dying cell   #baxm >20 then apoptosis
                    cells['living_cell'][ii]=False; #the cell dies
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
                    if len(max_tmp_cdc)>cells['division'][ii] :# && isempty(max_tmp_cycE)==false#make sure cycE has cycled
                        div_ind = (max_tmp_cdc[cells['division'][ii]]);
                        if len(max_tmp_cdc)<=cells['division'][ii] :
                            current_time[ii,0]=cells['time'][ii][-1];
                            cell_count_tmp[ii,:]=[cells['time'][ii][-1], 0];
                            
                        elif cells['division'][ii] == len(max_tmp_cycE) :
                            current_time[ii,0]=cells['time'][ii][-1];
                            cell_count_tmp[ii,:]=[cells['time'][ii][-1], 0];
                            
                        else :
                            for obs in list(cells['molecules'][ii].keys()) :
                                cells['molecules'][ii][obs]=cells['molecules'][ii][obs][0:div_ind]; #evolution of different species
                            cells['time'][ii]=cells['time'][ii][0:div_ind];
                            cells['division'][ii]=cells['division'][ii]+1;
                            
                            current_time[ii,0]=cells['time'][ii][-1];
                            cell_count_tmp[ii]=[cells['time'][ii][-1], 1];
                        #add a new cell
                            tmp_num=len(cells['time']);
                            
                            initial_new_cell={}
                            for obs in list(cells['molecules'][ii].keys()) :
                               initial_new_cell[obs]=[cells['molecules'][ii][obs][-1]]; #evolution of different species
                    
                            cells['molecules']=cells['molecules'] + [initial_new_cell]; #to change
                            cells['time']=cells['time']+[[current_time[ii,0]]];
                            cells['division']=cells['division']+[0];
                            cells['living_cell']=cells['living_cell']+[True];
                            current_time=np.concatenate((current_time, np.array([[current_time[ii,0], 1]])), axis=0);
                    else :            
                        current_time[ii,0]=cells['time'][ii][-1];
                        cell_count_tmp[ii]=[cells['time'][ii][-1], 0];
                        
                else : #size(max_tmp_cdc,1)==cells{ii,4} && length(find(cells{ii,2}(:,23)>25))<20 #quiescent cell$
                    current_time[ii,0]=cells['time'][ii][-1];
                    cell_count_tmp[ii]=[cells['time'][ii][-1], 0];
                    
            else :
                cell_count_tmp[ii]=[cells['time'][ii][-1], 0];
                
            combo = [ii, cells['time'][ii][-1]];
            
        for ii in range(len(cell_count_tmp)) :
        
            if cell_count_tmp[ii,1]==1 :
                ins_ind=[i for i,x in enumerate(cell_count[:,0]) if x < cell_count_tmp[ii,0]]
                current_cell_num=cell_count[ins_ind[-1],1];
                cell_count=np.insert(cell_count, ins_ind[-1]+1, cell_count_tmp[ii,:], axis=0)
                cell_count[ins_ind[-1]+1,1]=current_cell_num;
                cell_count[ins_ind[-1]+1:,1]=cell_count[ins_ind[-1]+1:,1]+1;
            
            elif cell_count_tmp[ii,1]==-1 :
                ins_ind=[i for i,x in enumerate(cell_count[:,0]) if x < cell_count_tmp[ii,0]]
                current_cell_num=cell_count[ins_ind[-1],1];
                cell_count=np.insert(cell_count, ins_ind[-1]+1, cell_count_tmp[ii,:], axis=0)
                cell_count[ins_ind[-1]+1,1]=current_cell_num;
                cell_count[ins_ind[-1]+1:,1]=cell_count[ins_ind[-1]+1:,1]-1;
            
            elif cell_count_tmp[ii,1]==0 :
                ins_ind=[i for i,x in enumerate(cell_count[:,0]) if x < cell_count_tmp[ii,0]]
                current_cell_num=cell_count[ins_ind[-1],1];
                cell_count=np.insert(cell_count, ins_ind[-1]+1, cell_count_tmp[ii,:], axis=0)
                cell_count[ins_ind[-1]+1,1]=current_cell_num;
                
        u, indices = np.unique(cell_count[:,0], return_index=True)
        cell_count=cell_count[indices,:];
    #disp(cell_count(end,:)) 
    
    return (cells, cell_count)

