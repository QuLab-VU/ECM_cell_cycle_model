# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 18:19:31 2020

@author: Grim.M
Title: runAKTERK
"""

# Import the basic libraries
import sys
import os
import numpy as np
import pandas as pd

# Import functions related to the model
sys.path.insert(0, os.path.dirname(__file__) + "/../cellular_cycle")
from Check_Model_Consistency import Check_Model_Consistency
from pysb.simulator import BngSimulator
from peakfinder import peakfinder

# Import functions related to the plots
sys.path.insert(0, os.path.dirname(__file__) + "/../figures")
from plots_cell_population_simulation import plot_cell_count
from plots_cell_population_simulation import break_down_cell_life
from plots_cell_population_simulation import plot_digraph_population_history
from plots_cell_population_simulation import plot_tree_from_population_df
from plots_cell_population_simulation import plot_event_population_history
from SMs_G1s_IMTs_extraction import SMs_G1s_IMTs_extraction
import fracprolif as fp

# Import functions related to the data processing
from retrieve_cellState_dw_2 import retrieve_cellState
from call_cell_events import call_events


#######################################################################

def Cell_Population_Simulation(init, run_time, observed_cells_time, model, ECM_conditions, condition, sampled_parameter_names, thr = [0.26, 0.18], return_plot=False, bypass=False, seed=None):
    """
    init: number of initial cells
    run_time: time of the experiment
    observed_cells_time: time of observation
    model: pysb model
    ECM_conditions: inputs for FAK, ERK, AKT signaling pathways related to the condition (stiffness and drug treatment)
    condition: condition (stiffness and drug treatment)
    cond_split_DMSO: condition for determinstic simulation, using DMSO treatment (initially, cells are untreated)
    """
    print(condition)

    initial_concentrations = pd.read_csv(os.getcwd() + '/model_data/initial_concentrations.csv')
    for key in initial_concentrations.keys():
        try:
            model.parameters[key].value = initial_concentrations[key][0] / model.parameters['omega'].value
        except:
            print(key)

    parameter_vector_original=[]
    for pname in sampled_parameter_names:
        #Change model parameter values to current location in parameter space
        parameter_vector_original+=[model.parameters[pname].value]

    cond_split_DMSO = condition.split(' ')[0:2]
    cond_split_DMSO.append('DMSO')
    condition_untreated = ' '.join(cond_split_DMSO)
    #Set initial condition to DMSO treatment
    for signal_tmp in ['ERK', 'FAK', 'AKT'] :
        # print(ECM_conditions)
        # ECM_conditions has col condition to be matched with the input splt DMSO, for each species from species col, and value col
        model.parameters[signal_tmp].value = ECM_conditions[(ECM_conditions["condition"] == condition_untreated) & (ECM_conditions["species"] == signal_tmp)]["value"]

    variable_range, check, tmp_score = Check_Model_Consistency(model, thr = thr, return_plot = return_plot)

    if not bypass:
        if check==False or np.any(np.array([variable_range[key] for key in variable_range.keys()])<0) or tmp_score>0.2: #the model doesn't oscillate or the model allows negativ concentrations (possible for the apoptosis module) or tmpscore>0.1
            print("check didn't pass");
            print("The model doesn't oscillate or the model allows too high concentrations")
            return (False, 0,0,0)

    #Set death module
    death_mod={'obsBAX': 186 ,'obsBAXmBCL2': 13, 'obsBH3': 0, 'obsBH3BCL2': 7, 'obsBAXm': 0, 'obsBCL2': 149}
    for prot in death_mod:
        variable_range[prot] = len(variable_range['obsCdc20'])*[death_mod[prot]]

    #Set the experiment to real condition
    for signal_tmp in ['ERK', 'FAK', 'AKT'] :
        model.parameters[signal_tmp].value = ECM_conditions[(ECM_conditions["condition"] == condition) & (ECM_conditions["species"] == signal_tmp)]["value"]

    cells={}
    cells['molecules']=[]; #Number of molecules taken randomly at one point during the cell cycle
    cells['time']=[]; #
    cells['division']=[]; #
    cells['living_cell']=[]; #living cell?
    cells['recording']=[];
    cells["ID"] = list(range(init))  # ID of the cell
    cells["Lineage"] = ["Initial cell"] * init  # Lineage of the cell
    cells["event"] = [None] * init  # event of the end of recording
    cells["Birth time"] = ["Initial cell"] * init  # time of the birth
    cells["divisions history"] = [[]] * init  # time of the divisions
    cells["division lineage history"] = [[]] * init  # daughter cells lineage
    cells["end of recording"] = [None] * init  # time of the end of recording
    cells["excluded Cdc20 peak"] = [0] * init  # Excluded Cdc20 peaks because of the check of the CycE oscillations

    for i in range(init) :
        cells['molecules']=cells['molecules']+[{}]
        a=np.random.randint(len(variable_range['obsSkp2']))
        for var in list(variable_range.keys()) :
            # print(variable_range[var][a])
            cells['molecules'][i][var] = [int(variable_range[var][a])]
        cells['time']=cells['time'] + [[0]]
        cells['division']=cells['division'] + [0]
        cells['living_cell']=cells['living_cell'] + [True] #living cell = True
        cells['recording']=cells['recording'] + [True] #recording the cell = True (if it isn't dead and if no event didn't occur after observed_cells_time)

    cell_count=np.array([[0., len(cells['molecules'])]]);   #start the cell counter at time zero and with number of initial cells
    current_time=np.transpose(np.concatenate([np.zeros((1,len(cells['molecules']))), np.ones((1,len(cells['molecules'])))])); #keeps track of the current time of each cell and whether it's dead or not
    sim = BngSimulator(model, verbose=False);

    while any(cells['recording']) : # Run the code only if there are cells that are recorded
        cell_count_tmp=np.zeros((1,2)); # Temporary cell count with current cell time and whether a cell was added, taken away, or not
        for ii in range(len(cells['molecules'])) : # For each cell
            if current_time[ii,0]<observed_cells_time and cells['living_cell'][ii] :#if there's time left and the cell has not died
                T = min(30, run_time+1-current_time[ii,0])  #80
                tmp_t = np.linspace(0, T, max(3,int(T*10)))

                Initial=[]
                for j in cells['molecules'][ii] :
                    Initial=Initial+[cells['molecules'][ii][j][-1]]
                simres = sim.run(initials=Initial, tspan = tmp_t); #seed removed

                tmp_t = simres.tout[0]
                yout = simres.observables

                for obs in list(cells['molecules'][ii].keys()) :
                    cells['molecules'][ii][obs]=cells['molecules'][ii][obs][:-1]+yout[obs].tolist(); #evolution of different species
                cells['time'][ii] = cells['time'][ii] + [x+cells['time'][ii][-1] for x in tmp_t[1:].tolist()]; #time step
                # current_ind_length=len(cells['molecules'][ii]['obsCdc20']);
                a, max_tmp_cdc = peakfinder(cells['molecules'][ii]['obsCdc20'], thr[0]*model.parameters['omega'].value, thr[0]*model.parameters['omega'].value, 1, includeEndpoints=1);      #cdc20
                a, max_tmp_cycE = peakfinder((cells['molecules'][ii]['obsCycE']), thr[1]*model.parameters['omega'].value, thr[1]*model.parameters['omega'].value, 1, includeEndpoints=1);  #cycE
                death = np.array([i for i,x in enumerate(cells['molecules'][ii]['obsBAXm']) if x > 20]) if len([i for i,x in enumerate(cells['molecules'][ii]['obsBAXm']) if x > 25])>20 else np.array([])
###############################################################################
                # Cell apoptosis
                death_before_division = len(max_tmp_cdc)==0 or np.min(death) < max_tmp_cdc[min([cells['division'][ii], len(max_tmp_cdc)-1])] if len(death)>0 else False

                if death.size>0 and death_before_division: # or np.any(np.array(cells['molecules'][ii]['obsBH3'])>=200) :#dying cell   #baxm >25 then apoptosis for 20 time steps or BAXm>=40 (it means it has started with more than 40 molecules of BAXm) # or np.any(np.array(cells['molecules'][ii]['obsBAXm'])>=40)

                    cells['living_cell'][ii]=False; #the cell dies
                    cells['recording'][ii]=False; #the record stops

                    current_time[ii,0] = cells['time'][ii][death[0]]; #death at time t
                    for obs in list(cells['molecules'][ii].keys()) :
                        cells['molecules'][ii][obs]=cells['molecules'][ii][obs][0:death[0]+1]
                    cells['time'][ii]=cells['time'][ii][0:death[0]+1];
                    current_time[ii,1]=0; #cell is dead

                    cell_count_tmp = np.insert(cell_count_tmp, len(cell_count_tmp), [cells['time'][ii][death[0]], -1], axis=0); #one cell to be removed to the count
                    cells["event"][ii] = "Apoptosis"
                    cells["end of recording"][ii] = cells['time'][ii][death[0]]  # time of the end of recording
                    #print('death 1')

###############################################################################
                #Dividing cell

                elif len(max_tmp_cdc) > (cells['division'][ii] + cells["excluded Cdc20 peak"][ii]):

                    divisions_tmp = np.count_nonzero(np.min(death) > np.array(max_tmp_cdc)) if  death.size != 0 else len(max_tmp_cdc) #divisions before apoptosis happens

                    #print('division 1')
                    for index in range(cells['division'][ii] + cells["excluded Cdc20 peak"][ii], divisions_tmp) :
                        div_ind = (max_tmp_cdc[index]);

                        if index==0 : # first cell, there is no need to check if CycE has cycled (it can have cycled before recording), it will be checked between two mitosis

                            cells['division'][ii] = cells['division'][ii]+1;
                            current_time[ii,0]=cells['time'][ii][-1]; #no change because this is the time and cell_count_tmp for the cell
                            cell_count_tmp = np.insert(cell_count_tmp, len(cell_count_tmp), [cells['time'][ii][div_ind], 1], axis=0);

                            if cells['time'][ii][div_ind] < observed_cells_time : # record the new cell only if it appears befor time of observation
                                initial_new_cell={}
                                for obs in list(cells['molecules'][ii].keys()) :
                                    initial_new_cell[obs]=[cells['molecules'][ii][obs][div_ind]]; #evolution of different species but at time=div_ind
                                cell_new_ID = cells["ID"][-1] + 1
                                cells['molecules'] = cells['molecules'] + [initial_new_cell]; #to change
                                cells['time'] = cells['time'] + [[cells['time'][ii][div_ind]]];
                                cells['division'] = cells['division'] + [1];
                                cells['living_cell'] = cells['living_cell'] + [True];
                                cells['recording'] = cells['recording'] + [True];
                                cells["ID"] = cells["ID"] + [cell_new_ID] # ID of the cell
                                cells["Lineage"] = cells["Lineage"] + [cells["ID"][ii]] # Lineage of the cell
                                cells["event"] = cells["event"] + [None] # event of the end of recording
                                cells["Birth time"] = cells["Birth time"] + [cells['time'][ii][div_ind]]  # time of the birth
                                cells["divisions history"] = cells["divisions history"] + [[]]  # time of the divisions
                                cells["division lineage history"] = cells["division lineage history"] + [[]]  # daughter cells lineage
                                cells["end of recording"] = cells["end of recording"] + [None]  # time of the end of recording
                                cells["excluded Cdc20 peak"] = cells["excluded Cdc20 peak"] + [0]
                                # Add this new cell to the lineage of cell ii
                                cells["divisions history"][ii] = cells["divisions history"][ii] + [cells['time'][ii][div_ind]]
                                cells["division lineage history"][ii] = cells["division lineage history"][ii] + [cell_new_ID]

                                current_time = np.concatenate((current_time, np.array([[cells['time'][ii][div_ind], 1]])), axis=0); #changed to the real time of the division

                        else :
                            #print('division 1 1')
                            if np.any(np.logical_and(max_tmp_cdc[index-1] < np.array(max_tmp_cycE) ,np.array(max_tmp_cycE) < max_tmp_cdc[index])) : #just checking if Cyc2 is actully oscillating between Cdc20

                                cells['division'][ii] = cells['division'][ii]+1; #it divides so let keep this
                                current_time[ii,0]=cells['time'][ii][-1]; #no change because this is the time and cell_count_tmp for the cell
                                cell_count_tmp = np.insert(cell_count_tmp, len(cell_count_tmp), [cells['time'][ii][div_ind], 1], axis=0);
                                if cells['time'][ii][div_ind] < observed_cells_time: #record the new cell only if the cdc20 appears after the cyc E peak
                                    cell_new_ID = cells["ID"][-1] + 1
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
                                    cells["ID"] = cells["ID"] + [cell_new_ID]  # ID of the cell
                                    cells["Lineage"] = cells["Lineage"] + [cells["ID"][ii]]  # Lineage of the cell
                                    cells["event"] = cells["event"] + [None]  # event of the end of recording
                                    cells["Birth time"] = cells["Birth time"] + [cells['time'][ii][div_ind]]  # time of the birth
                                    cells["divisions history"] = cells["divisions history"] + [[]]  # time of the divisions
                                    cells["division lineage history"] = cells["division lineage history"] + [[]]  # daughter cells lineage
                                    cells["end of recording"] = cells["end of recording"] + [None]  # time of the end of recording
                                    cells["excluded Cdc20 peak"] = cells["excluded Cdc20 peak"] + [0]

                                    # Add this new cell to the lineage of cell ii
                                    cells["divisions history"][ii] = cells["divisions history"][ii] + [cells['time'][ii][div_ind]]
                                    cells["division lineage history"][ii] = cells["division lineage history"][ii] + [cell_new_ID]

                                    current_time = np.concatenate((current_time, np.array([[cells['time'][ii][div_ind], 1]])), axis=0); #changed to the real time of the division
                            else:
                                current_time[ii,0]=cells['time'][ii][-1]; #no change because this is the time and cell_count_tmp for the cell
                                cells["excluded Cdc20 peak"][ii] += 1

                    if death.size != 0: #apoptosis right after mitosis

                        cells['living_cell'][ii]=False; #the cell dies
                        cells['recording'][ii]=False; #the record stops
                        cells["event"][ii] = "Apoptosis"
                        current_time[ii,0] = cells['time'][ii][death[0]]; #death at time t
                        for obs in list(cells['molecules'][ii].keys()) :
                            cells['molecules'][ii][obs]=cells['molecules'][ii][obs][0:death[0]+1]
                        cells['time'][ii]=cells['time'][ii][0:death[0]+1];
                        current_time[ii,1]=0; #cell is dead
                        cells["end of recording"][ii] = cells['time'][ii][death[0]]  # time of the end of recording

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
                cells['event'][ii] = "End of experiment"
                ind_tmp = np.where(np.array(cells['time'][ii])<=run_time)[-1][-1]
                cells['time'][ii] = cells['time'][ii][0:ind_tmp+2]
                for obs in list(cells['molecules'][ii].keys()) : cells['molecules'][ii][obs]=cells['molecules'][ii][obs][0:ind_tmp+2]
                current_time[ii,1]=cells['time'][ii][-1];
                cells["end of recording"][ii] = cells['time'][ii][-1]  # time of the end of recording

###############################################################################
###############################################################################
###############################################################################

#if the time is > to observed_cells_time, we don't take into account new cells, only the fate of the cells that are already born.

            if current_time[ii,0]<run_time and current_time[ii,0]>=observed_cells_time and cells['recording'][ii]:#if there's time left and the cell has not died/or no fate after observed_cells_time
                #print('after')
                T = 41 - (current_time[ii,0]>(run_time-40))*(40 - run_time + current_time[ii,0])
                tmp_t = np.linspace(0, T, max(3,int(T*10)))

                Initial=[]
                for j in cells['molecules'][ii] :
                    Initial=Initial+[cells['molecules'][ii][j][-1]]
                simres = sim.run(initials=Initial, tspan=tmp_t); # stop_if = 'obsBH3>200'

                tmp_t = simres.tout[0]
                yout = simres.observables

                for obs in list(cells['molecules'][ii].keys()) :
                    cells['molecules'][ii][obs]=cells['molecules'][ii][obs][:-1]+yout[obs].tolist(); #evolution of different species

                cells['time'][ii] = cells['time'][ii] + [x+cells['time'][ii][-1] for x in tmp_t[1:].tolist()]; #time step
                # current_ind_length=len(cells['molecules'][ii]['obsCdc20']);
                a, max_tmp_cdc = peakfinder(cells['molecules'][ii]['obsCdc20'],0.15*200,0.15*200);      #cdc20
                a, max_tmp_cycE = peakfinder(cells['molecules'][ii]['obsCycE'],0.15*200,0.15*200);  #cycE
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
                    cells["event"][ii] = "Apoptosis"

###############################################################################

                elif len(max_tmp_cdc) > (cells['division'][ii] + cells["excluded Cdc20 peak"][ii]):
                    divisions_tmp = np.count_nonzero(np.min(death) > np.array(max_tmp_cdc)) if  death.size != 0 else len(max_tmp_cdc)
                    cells['recording'][ii]=False; #the record stops
                    for index in range(cells['division'][ii] + cells["excluded Cdc20 peak"][ii], divisions_tmp) :
                        div_ind = (max_tmp_cdc[index]);
                        if index==0: #first cell, there is no need to check id CycE has cycled, it will be checked between two mitosis
                            #print('DIVISION 2 0')
                            cells['division'][ii] = cells['division'][ii]+1;
                            current_time[ii,0]=cells['time'][ii][-1]; #no change because this is the time and cell_count_tmp for the cell
                            cell_count_tmp = np.insert(cell_count_tmp, len(cell_count_tmp), [cells['time'][ii][div_ind], 1], axis=0);

                            if cells['time'][ii][div_ind] < observed_cells_time: #record the new cell only if it appears after the
                                initial_new_cell={}
                                cell_new_ID = cells["ID"][-1] + 1
                                for obs in list(cells['molecules'][ii].keys()) :
                                    initial_new_cell[obs]=[cells['molecules'][ii][obs][div_ind]]; #evolution of different species but at time=div_ind

                                cells['molecules']=cells['molecules'] + [initial_new_cell]; #to change
                                cells['time']=cells['time']+[[cells['time'][ii][div_ind]]];
                                cells['division']=cells['division']+[1];
                                cells['living_cell']=cells['living_cell']+[True];
                                cells['recording']=cells['recording']+[False];
                                cells["ID"] = cells["ID"] + [cell_new_ID]  # ID of the cell
                                cells["Lineage"] = cells["Lineage"] + [cells["ID"][ii]]  # Lineage of the cell
                                cells["event"] = cells["event"] + [None]  # event of the end of recording
                                cells["Birth time"] = cells["Birth time"] + [cells['time'][ii][div_ind]]  # time of the birth
                                cells["divisions history"] = cells["divisions history"] + [[]]  # time of the divisions
                                cells["division lineage history"] = cells["division lineage history"] + [[]]  # daughter cells lineage
                                cells["end of recording"] = cells["end of recording"] + [None]  # time of the end of recording
                                cells["excluded Cdc20 peak"] = cells["excluded Cdc20 peak"] + [0]
                                cells["divisions history"][ii] = cells["divisions history"][ii] + [cells['time'][ii][div_ind]]
                                cells["division lineage history"][ii] = cells["division lineage history"][ii] + [cell_new_ID]
                                current_time=np.concatenate((current_time, np.array([[cells['time'][ii][div_ind], 1]])), axis=0); #changed to the real time of the division

                        else :
                            if np.any(np.logical_and(max_tmp_cdc[index-1] < np.array(max_tmp_cycE) ,np.array(max_tmp_cycE) < max_tmp_cdc[index])) : #if cycE oscillates between two CDC20 peaks

                                #cells['time'][ii]=cells['time'][ii][0:div_ind];
                                cells['division'][ii] = cells['division'][ii]+1; #it divides so let keep this
                                current_time[ii,0]=cells['time'][ii][-1]; #no change because this is the time and cell_count_tmp for the cell
                                cell_count_tmp = np.insert(cell_count_tmp, len(cell_count_tmp), [cells['time'][ii][div_ind], 1], axis=0);
                                #add a new cell : Change to the time and concentrations of the previous cell but at time = div_ind
                                if cells['time'][ii][div_ind] < observed_cells_time: #record the new cell only if it appears after the
                                    cell_new_ID = cells["ID"][-1] + 1
                                    initial_new_cell={}
                                    for obs in list(cells['molecules'][ii].keys()) :
                                        initial_new_cell[obs]=[cells['molecules'][ii][obs][div_ind]]; #evolution of different species but at time=div_ind
                                    cells['molecules']=cells['molecules'] + [initial_new_cell]; #to change
                                    cells['time']=cells['time']+[[cells['time'][ii][div_ind]]];
                                    cells['division']=cells['division']+[1];
                                    cells['living_cell']=cells['living_cell']+[True];
                                    cells['recording']=cells['recording']+[False];
                                    cells["ID"] = cells["ID"] + [cell_new_ID]  # ID of the cell
                                    cells["Lineage"] = cells["Lineage"] + [cells["ID"][ii]]  # Lineage of the cell
                                    cells["event"] = cells["event"] + [None]  # event of the end of recording
                                    cells["Birth time"] = cells["Birth time"] + [cells['time'][ii][div_ind]]  # time of the birth
                                    cells["divisions history"] = cells["divisions history"] + [[]]  # time of the divisions
                                    cells["division lineage history"] = cells["division lineage history"] + [[]]  # daughter cells lineage
                                    cells["end of recording"] = cells["end of recording"] + [None]  # time of the end of recording
                                    cells["excluded Cdc20 peak"] = cells["excluded Cdc20 peak"] + [0]
                                    cells["divisions history"][ii] = cells["divisions history"][ii] + [cells['time'][ii][div_ind]]
                                    cells["division lineage history"][ii] = cells["division lineage history"][ii] + [cell_new_ID]


                                    current_time=np.concatenate((current_time, np.array([[cells['time'][ii][div_ind], 1]])), axis=0); #changed to the real time of the division
                            else:
                                current_time[ii,0]=cells['time'][ii][-1]; #no change because this is the time and cell_count_tmp for the cell
                                cells["excluded Cdc20 peak"][ii] += 1  # if the cell divides but CycE doesn't oscillate between two CDC20 peaks, we don't take into account this division, but it should be recorded to count the number of cdc20 peaks

                    if death.size != 0: #apoptosis right after mitosis
                        #print('apoptosis after DIVISION 2')
                        cells['living_cell'][ii]=False; #the cell dies
                        cells['recording'][ii]=False; #the record stops
                        current_time[ii,0] = cells['time'][ii][death[0]]; #death at time t
                        for obs in list(cells['molecules'][ii].keys()) :
                            cells['molecules'][ii][obs]=cells['molecules'][ii][obs][0:death[0]+1]
                        cells['time'][ii]=cells['time'][ii][0:death[0]+1];
                        current_time[ii,1]=0; #cell is dead
                        cells["event"][ii] = "Apoptosis"
                        cells["end of recording"][ii] = cells['time'][ii][death[0]]  # time of the end of recording
                        cell_count_tmp = np.insert(cell_count_tmp, len(cell_count_tmp), [cells['time'][ii][death[0]], -1], axis=0); #one cell to be removed to the count

###############################################################################
                #Quiescent cell

                else : #size(max_tmp_cdc,1)==cells{ii,4} && length(find(cells{ii,2}(:,23)>25))<20 #quiescent cell$
                    current_time[ii,0]=cells['time'][ii][-1];
                    cell_count_tmp = np.insert(cell_count_tmp, len(cell_count_tmp), [cells['time'][ii][-1], 0], axis=0);

###############################################################################

            if current_time[ii,0]>=run_time :
                ##if the cell reached the end of the experiment the record stops
                cells['event'][ii] = "End of experiment"
                cells['recording'][ii]=False; #Stop the recording
                ind_tmp = np.where(np.array(cells['time'][ii])<=run_time)[-1][-1]
                cells['time'][ii] = cells['time'][ii][0:ind_tmp+2]
                for obs in list(cells['molecules'][ii].keys()) : cells['molecules'][ii][obs]=cells['molecules'][ii][obs][0:ind_tmp+2]
                current_time[ii,1]=cells['time'][ii][-1];
                cells["end of recording"][ii] = cells['time'][ii][-1]  # time of the end of recording

###############################################################################
###############################################################################

            else :
                cell_count_tmp = np.insert(cell_count_tmp, len(cell_count_tmp), [cells['time'][ii][-1], 0], axis=0);


###############################################################################
###############################################################################
###############################################################################
###############################################################################

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

    #### Processing the cell lineage

    selected_columns = ['ID', 'Lineage', 'event', 'division', 'Birth time', 'divisions history',
                        'division lineage history', "end of recording", "living_cell", "recording"]
    subset_dict = {key: cells[key] for key in selected_columns}
    cells_df = pd.DataFrame(subset_dict)
    population_df = break_down_cell_life(cells_df)

    #### Plotting the results
    if return_plot :
        plot_cell_count(cell_count, run_time)
        plot_digraph_population_history(population_df)
        g = plot_event_population_history(population_df, run_time, condition)
        # save g in pdf with savefig
        # g.savefig("/Users/jeremie/Downloads/" + condition.replace(" ", "_") + '_population_events.pdf')



    return (True, population_df, cells, cell_count)




if __name__=="__main__":
    # import model
    from cell_cycle_death_model import model
    import pandas as pd
    import matplotlib.pyplot as plt

    # import pysb sampled parameters
    sampled_parameter_names = np.loadtxt(os.path.dirname(__file__) + '/../model_data/sampled_parameter_names.txt', dtype=str)

    # Importing the parameters:
    new_param = pd.read_csv(os.path.dirname(__file__) + '/../model_data/parameters.csv')
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

    # Import ECM_conditions.csv
    ECM_conditions = pd.read_csv(os.path.dirname(__file__) + '/../model_data/ECM_conditions_2.csv')
    condition = "4 kPa erlotinib"

    # Set the species AKT, ERK, FAK with species and value col
    for input in ["AKT", "ERK", "FAK"]:
        model.parameters[input].value = ECM_conditions[
            (ECM_conditions["condition"] == condition) & (ECM_conditions["species"] == input)].value

    # Set inputs for Cell_Population_Simulation function
    init = 2
    run_time = 120
    observed_cells_time = 120
    
    path_savefig = "/Users/jeremie/Downloads"

    (check, population_df, cells, cell_count) = Cell_Population_Simulation(init, run_time, observed_cells_time, model, ECM_conditions, condition, sampled_parameter_names, thr=[0.26, 0.14], return_plot=True)
    print(cell_count)
    plot_tree_from_population_df(population_df, path_savefig + "/tree-plot_" + condition.replace(" ", "_") + ".svg")

    condition = "0.6 kPa DMSO"
    (check, population_df, cells, cell_count) = Cell_Population_Simulation(init, run_time, observed_cells_time, model, ECM_conditions, condition, sampled_parameter_names, thr=[0.26, 0.14], return_plot=True)
    plot_tree_from_population_df(population_df, path_savefig + "/tree-plot_" + condition.replace(" ", "_") + ".svg")

    condition = "0.6 kPa erlotinib"
    (check, population_df, cells, cell_count) = Cell_Population_Simulation(init, run_time, observed_cells_time, model, ECM_conditions, condition, sampled_parameter_names, thr=[0.26, 0.14], return_plot=True)
    plot_tree_from_population_df(population_df, path_savefig + "/tree-plot_" + condition.replace(" ", "_") + ".svg")

    condition = "4 kPa DMSO"
    (check, population_df, cells, cell_count) = Cell_Population_Simulation(init, run_time, observed_cells_time, model, ECM_conditions, condition, sampled_parameter_names, thr=[0.26, 0.14], return_plot=True)
    plot_tree_from_population_df(population_df, path_savefig + "/tree-plot_" + condition.replace(" ", "_") + ".svg")

    condition = "55 kPa DMSO"
    (check, population_df, cells, cell_count) = Cell_Population_Simulation(init, run_time, observed_cells_time, model, ECM_conditions, condition, sampled_parameter_names, thr=[0.26, 0.14], return_plot=True)
    plot_tree_from_population_df(population_df, path_savefig + "/tree-plot_" + condition.replace(" ", "_") + ".svg")

    condition = "55 kPa erlotinib"
    (check, population_df, cells, cell_count) = Cell_Population_Simulation(init, run_time, observed_cells_time, model, ECM_conditions, condition, sampled_parameter_names, thr=[0.26, 0.14], return_plot=True)
    plot_tree_from_population_df(population_df, path_savefig + "/tree-plot_" + condition.replace(" ", "_") + ".svg")


    print(cells["ID"])
    print(cells["Lineage"])
    print(cells["event"])
    print(cells["division"])
    print(cells["Birth time"])
    print(cells["divisions history"])
    print(cells["division lineage history"])
    print(cells["end of recording"])

    IMT_eoe_counts = retrieve_cellState(run_time, cells, return_all_prot=False,
                                        cells_starting_at_mitosis=False)  # cell count
    all_events_simulated = call_events(run_time, IMT_eoe_counts);
    all_events_simulated = all_events_simulated[all_events_simulated['event'] != 'Starting cell']
    all_events_simulated = all_events_simulated[np.array(all_events_simulated['Birthtime'] < observed_cells_time)]
    mitotic_lifespans = all_events_simulated['Lifespan'][all_events_simulated['event'] == 'mitosis']
    censored_lifespans = all_events_simulated['Lifespan'][all_events_simulated['event'] == 'End_of_Expt']

    # # Remove Starting cells from population_df, ie population_df['Starting Cell'] should be False
    # population_df = population_df[population_df['Starting Cell'] == False]
    # mitotic_lifespans_new = population_df['Lifetime'][population_df['Event'] == 'Mitosis']
    # censored_lifespans_new = population_df['Lifetime'][population_df['Event'] == 'End of experiment']
    #
    # # print("Mitotic lifespans: ", mitotic_lifespans)
    # # print("Mitotic lifespans new method: ", mitotic_lifespans_new)
    # # print("Censored lifespans: ", censored_lifespans)
    # # print("Censored lifespans new method: ", censored_lifespans_new)
    #
    #
    # simulated_IMTs, simulated_SMs, simulated_G1s = SMs_G1s_IMTs_extraction(cells) #G1, SM, IMT
    #
    # parameters_sim = {}
    #
    # if len(population_df) > 0:
    #     number_of_deaths = len(population_df) - (population_df['Event']=='Apoptosis').value_counts()[0]
    #     cell_time_total = sum(population_df['Lifetime'])
    #     parameters_sim['a'] = number_of_deaths/cell_time_total
    # else:
    #     parameters_sim['a'] = np.nan
    #
    # # mitosis and quiescence rates
    # if len(mitotic_lifespans) > 0:
    #     m = fp.mle(mitotic_lifespans, censored_lifespans)
    #     parameters_sim['mu'] = m['x'][0]
    #     parameters_sim['sigma'] = m['x'][1]
    #     parameters_sim['lambda'] = m['x'][2]
    #     parameters_sim.update(fp.q_rates(m))
    # else:
    #     parameters_sim['d'] = np.nan
    #     parameters_sim['q'] = np.nan
    #     parameters_sim['mu'] = np.nan
    #     parameters_sim['sigma'] = np.nan
    #     parameters_sim['lambda'] = np.nan
    #
    plot_tree_from_population_df(population_df, path_savefig + "/tree-plot.svg")

