# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 18:19:31 2020

@author: Grim.M
Title: runAKTERK
"""

import sys
sys.path.insert(0, __file__ + "/../cellular_cycle")
from pysb import * 
import numpy as np

from Check_Model_Consistency import CheckModelConsistency
from pysb.simulator import BngSimulator
from peakfinder import peakfinder


#######################################################################

def Cell_Population_Simulation(init, run_time, observed_cells_time, model, ECM_conditions, condition, sampled_parameter_names, thr = 0.26):
    """
    init: number of initial cells
    run_time: time of the experiment
    observed_cells_time: time of observation
    model: pysb model
    ECM_conditions: inputs for FAK, ERK, AKT signaling pathways related to the condition (stiffness and drug treatment)
    condition: condition (stiffness and drug treatment)
    cond_split_DMSO: condition for determinstic simulation, using DMSO treatment (initially, cells are untreated)
    """

    parameter_vector_original=[]
    for pname in sampled_parameter_names:
        #Change model parameter values to current location in parameter space
        parameter_vector_original+=[model.parameters[pname].value]

    cond_split_DMSO = condition.split(' ')
    condition_untreated = ' '.join(cond_split_DMSO)
    print(condition_untreated)
    #Set initial condition to DMSO treatment
    for signal_tmp in ['ERK', 'FAK', 'AKT'] :
        # ECM_conditions has col condition to be matched with the input splt DMSO, for each species from species col, and value col
        model.parameters[signal_tmp].value = ECM_conditions[(ECM_conditions["condition"] == condition_untreated) & (ECM_conditions["species"] == signal_tmp)]["value"]

    variable_range, check, tmp_score = CheckModelConsistency(model, thr = 0.26, return_plot = True)

    if check==False or np.any(np.array([variable_range[key] for key in variable_range.keys()])<0) or tmp_score>0.2: #the model doesn't oscillate or the model allows negativ concentrations (possible for the apoptosis module) or tmpscore>0.1
        print("check didn't pass");
        return (0,0)

    #Set the experiment to real condition
    for signal_tmp in ['ERK', 'FAK', 'AKT'] :
        model.parameters[signal_tmp].value = ECM_conditions[(ECM_conditions["condition"] == condition) & (ECM_conditions["species"] == signal_tmp)]["value"]
    
    cells={}
    cells['molecules']=[]; #Number of molecules taken randomly at one point during the cell cycle
    cells['time']=[]; #
    cells['division']=[]; #
    cells['living_cell']=[]; #living cell?
    cells['recording']=[];
    
    #plt.plot(variable_range['obsCdc20'])
    
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
    
    sim = BngSimulator(model_population, verbose=False);
    
    while any(cells['recording']) : #run the code only if there are cells that are recorded
        #print('while loop is starting')
        
        cell_count_tmp=np.zeros((1,2)); #temporary cell count with current cell time and whether a cell was added, taken away, or not
        
        #print(cell_count_tmp)
        for ii in range(len(cells['molecules'])) :
            #print(current_time)
            #print("new cell")
            
            #If the time is between 0 to observed_cells_time
            
            if current_time[ii,0]<observed_cells_time and cells['living_cell'][ii] :#if there's time left and the cell has not died
                #print('before')
                #print(ii)
                
                T = min(30, run_time+1-current_time[ii,0])  #80
                tmp_t = np.linspace(0, T, max(3,int(T*10)))
                
                
                
                if len(model_population.species)>0 :
                    Initial=[]
                    for j in cells['molecules'][ii] :
                        Initial=Initial+[cells['molecules'][ii][j][-1]]
                    simres = sim.run(initials=Initial, tspan = tmp_t);
                    #print(parameter_vector_original)
                    #print("simulation done")
                    
                else :
                    simres = sim.run(verbose=False);
                    #print(parameter_vector_original)
                    #print("simulation done")
                    
                tmp_t = simres.tout[0]
                    
                yout = simres.observables
                
                
                for obs in list(cells['molecules'][ii].keys()) :
                    cells['molecules'][ii][obs]=cells['molecules'][ii][obs][:-1]+yout[obs].tolist(); #evolution of different species
                    
                    
            #if type(cells['time'][ii])!=list : #if it's an array
            #    cells['time'][ii]=cells['time'][ii].tolist()
                
                
                cells['time'][ii] = cells['time'][ii] + [x+cells['time'][ii][-1] for x in tmp_t[1:].tolist()]; #time step
                current_ind_length=len(cells['molecules'][ii]['obsCdc20']);
                
                
                a, max_tmp_cdc = peakfinder(cells['molecules'][ii]['obsCdc20'], thr*model_population.parameters['omega'].value, thr*model_population.parameters['omega'].value, 1, includeEndpoints=1);      #cdc20
                
                a, max_tmp_cycE = peakfinder((cells['molecules'][ii]['obsCycE']), thr*model_population.parameters['omega'].value, thr*model_population.parameters['omega'].value, 1, includeEndpoints=1);  #cycE
                
                death = np.array([i for i,x in enumerate(cells['molecules'][ii]['obsBAXm']) if x > 20]) if len([i for i,x in enumerate(cells['molecules'][ii]['obsBAXm']) if x > 25])>20 else np.array([])
                
#                plt.plot(cells['molecules'][ii]['obsCycE'])
#                plt.plot(cells['molecules'][ii]['obsCdc20'])
###############################################################################
                #Cell apoptosis
                death_before_division = len(max_tmp_cdc)==0 or np.min(death) < max_tmp_cdc[min([cells['division'][ii], len(max_tmp_cdc)-1])] if len(death)>0 else False
                
                if death.size>0 and death_before_division:# or np.any(np.array(cells['molecules'][ii]['obsBH3'])>=200) :#dying cell   #baxm >25 then apoptosis for 20 time steps or BAXm>=40 (it means it has started with more than 40 molecules of BAXm) # or np.any(np.array(cells['molecules'][ii]['obsBAXm'])>=40)
                        
                    cells['living_cell'][ii]=False; #the cell dies
                    cells['recording'][ii]=False; #the record stops
                    
                    current_time[ii,0] = cells['time'][ii][death[0]]; #death at time t
                    for obs in list(cells['molecules'][ii].keys()) :
                        cells['molecules'][ii][obs]=cells['molecules'][ii][obs][0:death[0]+1]
                    cells['time'][ii]=cells['time'][ii][0:death[0]+1];
                    current_time[ii,1]=0; #cell is dead
                    
                    cell_count_tmp = np.insert(cell_count_tmp, len(cell_count_tmp), [cells['time'][ii][death[0]], -1], axis=0); #one cell to be removed to the count
                    #print('death 1')
                
###############################################################################
                #Dividing cell 
                
                elif len(max_tmp_cdc) > cells['division'][ii]:
                    
                    divisions_tmp = np.count_nonzero(np.min(death) > np.array(max_tmp_cdc)) if  death.size != 0 else len(max_tmp_cdc) #divisions before apoptosis happens
                    
                    #print('division 1')
                    for index in range(cells['division'][ii], divisions_tmp) :
                        
                        
                        div_ind = (max_tmp_cdc[index]);
                        #print('index')
                        #print(index)
                        
                        if index==0 : #first cell, there is no need to check if CycE has cycled (it can have cycled before recording), it will be checked between two mitosis
                            #print('division 1 0')
                            
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
                            #print('division 1 1')
                            if np.any(np.logical_and(max_tmp_cdc[index-1] < np.array(max_tmp_cycE) ,np.array(max_tmp_cycE) < max_tmp_cdc[index])) : #just checking if Cyc2 is actully oscillating between Cdc20
                            
                                #cells['time'][ii]=cells['time'][ii][0:div_ind];
                                cells['division'][ii] = cells['division'][ii]+1; #it divides so let keep this
                                
                                current_time[ii,0]=cells['time'][ii][-1]; #no change because this is the time and cell_count_tmp for the cell
                                
                                cell_count_tmp = np.insert(cell_count_tmp, len(cell_count_tmp), [cells['time'][ii][div_ind], 1], axis=0);
                                    
                            #add a new cell : Change to the time and concentrations of the previous cell but at time = div_ind
                            
                                tmp_num=len(cells['time']);
                                
                                if cells['time'][ii][div_ind] < observed_cells_time: #record the new cell only if the cdc20 appears after the cyc E peak
                                    
                                    initial_new_cell={}
                                    for obs in list(cells['molecules'][ii].keys()) :
                                        initial_new_cell[obs]=[cells['molecules'][ii][obs][div_ind]]; #evolution of different species but at time=div_ind
                            
                                    cells['molecules'] = cells['molecules'] + [initial_new_cell]; #to change
                                    cells['time'] = cells['time']+[[cells['time'][ii][div_ind]]];
                                    cells['division'] = cells['division']+[1];
                                    cells['living_cell'] = cells['living_cell']+[True];
                                    
                                    if cells['time'][ii][div_ind]==run_time:
                                        cells['recording']=cells['recording']+[False];
                                    else : #mitosis happens right on the end time (almost impossible)
                                        cells['recording']=cells['recording']+[True];
                                        
                                    current_time = np.concatenate((current_time, np.array([[cells['time'][ii][div_ind], 1]])), axis=0); #changed to the real time of the division
                            else: 
                                current_time[ii,0]=cells['time'][ii][-1]; #no change because this is the time and cell_count_tmp for the cell
                                
                    if death.size != 0: #apoptosis right after mitosis
                    
                        cells['living_cell'][ii]=False; #the cell dies
                        cells['recording'][ii]=False; #the record stops
                        
                        current_time[ii,0] = cells['time'][ii][death[0]]; #death at time t
                        for obs in list(cells['molecules'][ii].keys()) :
                            cells['molecules'][ii][obs]=cells['molecules'][ii][obs][0:death[0]+1]
                        cells['time'][ii]=cells['time'][ii][0:death[0]+1];
                        current_time[ii,1]=0; #cell is dead
                        
                        cell_count_tmp = np.insert(cell_count_tmp, len(cell_count_tmp), [cells['time'][ii][death[0]], -1], axis=0); #one cell to be removed to the count
                        #print('apoptosis after division 1')

###############################################################################
                        #quiescent cell

                else : #size(max_tmp_cdc,1)==cells{ii,4} && length(find(cells{ii,2}(:,23)>25))<20 #quiescent cell
                    current_time[ii,0]=cells['time'][ii][-1];
                    cell_count_tmp = np.insert(cell_count_tmp, len(cell_count_tmp), [cells['time'][ii][-1], 0], axis=0);
                    #print('quiescence 1')

###############################################################################


            if current_time[ii,0]>=run_time : #if the cell reached the end of the experiment
                #print('stop cell recording')
                cells['recording'][ii]=False; #Stop the recording
                
                ind_tmp = np.where(np.array(cells['time'][ii])<=run_time)[-1][-1]
                
                cells['time'][ii] = cells['time'][ii][0:ind_tmp+2]
                for obs in list(cells['molecules'][ii].keys()) : cells['molecules'][ii][obs]=cells['molecules'][ii][obs][0:ind_tmp+2]
                current_time[ii,1]=cells['time'][ii][-1]; 
                
###############################################################################
###############################################################################
###############################################################################

#if the time is > to observed_cells_time, we don't take into account new cells, only the fate of the cells that are already born.

            if current_time[ii,0]<run_time and current_time[ii,0]>=observed_cells_time and cells['recording'][ii]:#if there's time left and the cell has not died/or no fate after observed_cells_time
                #print('after')
                T = 41 - (current_time[ii,0]>(run_time-40))*(40 - run_time + current_time[ii,0])
                tmp_t = np.linspace(0, T, max(3,int(T*10)))
                
                
                if len(model_population.species)>0 :
                    Initial=[]
                    for j in cells['molecules'][ii] :
                        Initial=Initial+[cells['molecules'][ii][j][-1]]
                    simres = sim.run(initials=Initial, tspan=tmp_t); # stop_if = 'obsBH3>200'
                    #print(parameter_vector_original)
                    
                else :
                    simres = sim.run(verbose=False); # stop_if = 'obsBH3>200'
                    #print(parameter_vector_original)
                #print('stoch sim 2')
                
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
                
                death = np.array([i for i,x in enumerate(cells['molecules'][ii]['obsBAXm']) if x > 20]) if len([i for i,x in enumerate(cells['molecules'][ii]['obsBAXm']) if x > 25])>20 else np.array([])
                
                

###############################################################################
                #Cell apoptosis
                death_before_division = len(max_tmp_cdc)==0 or np.min(death) < max_tmp_cdc[min([cells['division'][ii], len(max_tmp_cdc)-1])] if len(death)>0 else False
                if death.size>0 and death_before_division:
                    cells['recording'][ii]=False; #the record stops
                    cells['living_cell'][ii]=False; #the cell dies
                    death=[i for i,x in enumerate(cells['molecules'][ii]['obsBAXm']) if x > 20];
                    current_time[ii,0]=cells['time'][ii][death[0]]; #death at time t
                    for obs in list(cells['molecules'][ii].keys()) :
                        cells['molecules'][ii][obs]=cells['molecules'][ii][obs][0:death[0]+1]
                    cells['time'][ii]=cells['time'][ii][0:death[0]+1];
                    current_time[ii,1]=0; #cell is death
                    cell_count_tmp = np.insert(cell_count_tmp, len(cell_count_tmp), [cells['time'][ii][death[0]], -1], axis=0); #one cell to be removed to the count
                    
###############################################################################

                elif len(max_tmp_cdc) > cells['division'][ii]:
                    #print('DIVISION 2')
                    divisions_tmp = np.count_nonzero(np.min(death) > np.array(max_tmp_cdc)) if  death.size != 0 else len(max_tmp_cdc)
                    
                    cells['recording'][ii]=False; #the record stops

                    for index in range(cells['division'][ii], divisions_tmp) :
                        
                        div_ind = (max_tmp_cdc[index]);
                        #print('index')
                        

                        
                        if index==0: #first cell, there is no need to check id CycE has cycled, it will be checked between two mitosis
                            #print('DIVISION 2 0')
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
                            #print('DIVISION 2 1')                               
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
                            else: 
                                current_time[ii,0]=cells['time'][ii][-1]; #no change because this is the time and cell_count_tmp for the cell
                    if death.size != 0: #apoptosis right after mitosis
                        #print('apoptosis after DIVISION 2')
                        cells['living_cell'][ii]=False; #the cell dies
                        cells['recording'][ii]=False; #the record stops
                        
                        current_time[ii,0] = cells['time'][ii][death[0]]; #death at time t
                        for obs in list(cells['molecules'][ii].keys()) :
                            cells['molecules'][ii][obs]=cells['molecules'][ii][obs][0:death[0]+1]
                        cells['time'][ii]=cells['time'][ii][0:death[0]+1];
                        current_time[ii,1]=0; #cell is dead
                        
                        cell_count_tmp = np.insert(cell_count_tmp, len(cell_count_tmp), [cells['time'][ii][death[0]], -1], axis=0); #one cell to be removed to the count
###############################################################################
                #Quiescent cell
                
                else : #size(max_tmp_cdc,1)==cells{ii,4} && length(find(cells{ii,2}(:,23)>25))<20 #quiescent cell$
                    #print('quiescence')
                    #print('quiescence 2')
                    current_time[ii,0]=cells['time'][ii][-1];
                    cell_count_tmp = np.insert(cell_count_tmp, len(cell_count_tmp), [cells['time'][ii][-1], 0], axis=0);
                    
###############################################################################

            if current_time[ii,0]>=run_time : 
                ##if the cell reached the end of the experiment the record stops
            
                #print('stop cell recording')
                cells['recording'][ii]=False; #Stop the recording
                ind_tmp = np.where(np.array(cells['time'][ii])<=run_time)[-1][-1]
                cells['time'][ii] = cells['time'][ii][0:ind_tmp+2]
                for obs in list(cells['molecules'][ii].keys()) : cells['molecules'][ii][obs]=cells['molecules'][ii][obs][0:ind_tmp+2]
                current_time[ii,1]=cells['time'][ii][-1]; 
                        
###############################################################################
###############################################################################

            else :
                cell_count_tmp = np.insert(cell_count_tmp, len(cell_count_tmp), [cells['time'][ii][-1], 0], axis=0);
                
            combo = [ii, cells['time'][ii][-1]];

###############################################################################
###############################################################################
###############################################################################
###############################################################################

        try :
            #print('updating cell count')
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
        except :
            print("Cell population simulation error")

    return (cells, cell_count)




if __name__=="__main__":
    # import model
    from cell_cycle_death_model import model
    import pandas as pd
    import os

    # import pysb sampled parameters
    sampled_parameter_names = np.loadtxt(os.path.dirname(__file__) + '/../model_data/sampled_parameter_names.txt',
                                         dtype=str)

    # Importing the parameters:
    new_param = pd.read_csv(os.path.dirname(__file__) + '/../model_data/parameters.csv')
    del new_param['Stress']
    for key in new_param.keys():
        try:
            if model.parameters[key].value != new_param[key][0]:
                model.parameters[key].value = new_param[key]
        except:
            print(key)
    model.parameters['omega'].value = 200

    # Setting initial concentrations '/../model_data/initial_concentrations.csv'
    initial_concentrations = pd.read_csv(os.path.dirname(__file__) + '/../model_data/initial_concentrations.csv')
    for key in initial_concentrations.keys():
        try:
            model.parameters[key].value = initial_concentrations[key][0] / model.parameters['omega'].value
        except:
            print(key)

    # # Import ECM_conditions.csv
    ECM_conditions = pd.read_csv(os.path.dirname(__file__) + '/../model_data/ECM_conditions_2.csv')
    condition = '4 kPa erlotinib'

    # Set the species AKT, ERK, FAK with species and value col
    # for input in ["AKT", "ERK", "FAK"]:
    #     model.parameters[input].value = ECM_conditions[
    #         (ECM_conditions["condition"] == condition) & (ECM_conditions["species"] == input)].value

    # Set inputs for runAKTERK function
    init = 1
    run_time = 40
    observed_cells_time = 40

    # Run the model
    # Simulate CheckModelConsistency
    # (variable_range, check, tmp_score) = CheckModelConsistency(model, thr=0.26, return_plot=True)

    (cells, cell_count) = runAKTERK(init, run_time, observed_cells_time, model, ECM_conditions, condition, sampled_parameter_names)
    # print(cells)
    # print(cell_count)
