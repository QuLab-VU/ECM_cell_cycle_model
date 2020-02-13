# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 18:19:31 2020

@author: Grim.M
Title: runAKTERK
"""

from pysb import * 

import numpy as np
import pylab as pl

import sys
sys.path.insert(0, r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code\cellular_cycle") #add a pathway to the model

from Austin_proliferation_model_population import model
from pysb.simulator import BngSimulator
from peakfinder import peakfinder


v = 0.05;     
sys_size=200; #cell size
run_time=120;
 
##Cells:

cells={}


cells['molecules']=[]; #Number of molecules taken randomly at one point during the cell cycle
cells['time']=[]; #
cells['division']=[]; #
cells['living_cell']=[]; #living cell?


for i in range(init) :
    cells['molecules']=cells['molecules']+[{}]
    a=np.random.randint(len(variable_range))
    for var in list(variable_range.keys()) :
        cells['molecules'][i][var] = [round(variable_range[var][a]*200)]
    cells['time']=cells['time'] + [[0]]
    cells['division']=cells['division'] + [0]
    cells['living_cell']=cells['living_cell'] + [True] #living cell = True
    
 
##
cell_count=np.array([[0, len(cells['molecules'])]]);   #start the cell counter at time zero and with number of initial cells
current_time=np.ones((len(cells['molecules']),2)); #keeps track of the current time of each cell and whether it's dead or not
while (np.isin(1, current_time[:,1]) and np.any(current_time[:,0]<run_time)) : #run the code only if there are live cells under the run time
    cell_count_tmp=np.zeros((len(current_time),2)); #temporary cell count with current cell time and whether a cell was added, taken away, or not
    for ii in range(len(cells['molecules'])) :
       if current_time[ii,0]<120 and cells['living_cell'][ii] :#if there's time left and the cell has not died
           '''
            if run_time-cells{ii,3}(end)<30
                p.run_time=run_time-cells{ii,3}(end);
            else
                p.run_time=30;
            '''
            
            tmp_t = pl.linspace(0, 500, 50000)
            
            ##Ajouter concentrations initiales
            #Model()
            from Austin_proliferation_model_population import model

            
            sim = BngSimulator(model, tspan=tmp_t, verbose=True)
            simres = sim.run()
            yout = simres.observables
            
            yout=yout[10002:50000]
            tmp_t=pl.linspace(0,400, 39998)
            
            for obs in list(cells['molecules'][ii].keys()) :
                cells['molecules'][ii][obs]=cells['molecules'][ii][obs]+yout[obs].tolist(); #evolution of different species
            
            cells['time'][ii]=cells['time'][ii] + tmp_t[1:-1].tolist(); #time step
            current_ind_length=len(cells['molecules'][ii]['obsCdc20']);
 
            a, max_tmp_cdc=peakfinder(cells['molecules'][ii]['obsCdc20'],0.15*200,0.15*200);      #cdc20
            a, max_tmp_cycE=peakfinder((cells['molecules'][ii]['obsCycE']),0.15*200,0.15*200);  #cycE
            #min_tmp=peakfinder(cells{ii,2}(:,14),[],30,-1);   #cdc20
 
            if len([i for i,x in enumerate(cells['molecules'][ii]['obsBAXm']) if x > 25])>20 #dying cell   #baxm >20 then apoptosis
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
                        cells['time'][ii]=cells['time'][ii][0:div_ind];
                        cells['division'][ii]=cells['division'][ii]+1;
                        current_time[ii,0]=cells['time'][ii][-1];
                        cell_count_tmp[ii]=[cells['time'][ii][-1], 1];
                        #add a new cell
                        tmp_num=len(cells['time']);
                        cells['molecules']=cells['molecules'] + [initial_divided]; #to change
                        cells['time']=cells['time']+[current_time[ii,0]];
                        cells['division']=cells['division']+[0];
                        cells['living_cell']=cells['living_cell']+[True];
                        current_time=np.concatenate((current_time, np.array([[current_time[ii,0], 1]])), axis=0);
                else :            
                    current_time[ii,0]=cells['time'][ii][-1];
                    cell_count_tmp[ii]=[cells['time'][ii][-1], 0];

            else : #size(max_tmp_cdc,1)==cells{ii,4} && length(find(cells{ii,2}(:,23)>25))<20 #quiescent cell
                current_time[ii,0]=cells['time'][ii][-1];
                cell_count_tmp[ii]=[cells['time'][ii][-1], 0];
            
       else :
           cell_count_tmp[ii]=[cells['time'][ii][-1], 0];
       
       combo = [ii, cells['time'][ii][-1]];
       print(combo)
    
    for ii in range(len(cell_count_tmp)) :
        
        if cell_count_tmp[ii,1]==1 :
            ins_ind=find(cell_count_tmp[ii,0]>cell_count[:,0]); [i for i,x in enumerate(cells['molecules'][ii]['obsBAXm']) if x > cell_count_tmp[ii,0]]
            current_cell_num=cell_count(ins_ind(end),2);
            cell_count=insertrows(cell_count,cell_count_tmp(ii,:),ins_ind(end));
            cell_count(ins_ind(end)+1,2)=current_cell_num;
            cell_count(ins_ind(end)+1:end,2)=cell_count(ins_ind(end)+1:end,2)+1;
            
        elif cell_count_tmp[ii,1]==-1 :
            ins_ind=find(cell_count_tmp(ii,1)>cell_count(:,1));
            current_cell_num=cell_count(ins_ind(end),2);
            cell_count=insertrows(cell_count,cell_count_tmp(ii,:),ins_ind(end));
            cell_count(ins_ind(end)+1,2)=current_cell_num;
            cell_count(ins_ind(end)+1:end,2)=cell_count(ins_ind(end)+1:end,2)-1;
            
        elif cell_count_tmp[ii,1]==0 :
            ins_ind=find(cell_count_tmp(ii,1)>cell_count(:,1));
            current_cell_num=cell_count(ins_ind(end),2);
            cell_count=insertrows(cell_count,cell_count_tmp(ii,:),ins_ind(end));
            cell_count(ins_ind(end)+1,2)=current_cell_num;
        
         
     
    [~,ia]=unique(cell_count(:,1));
    cell_count=cell_count(ia,:);
    #disp(cell_count(end,:))    
     
