# -*- coding: utf-8 -*-
"""
Created on Mon Aug 10 11:07:54 2020

@author: Grim.M
"""
import pandas as pd

def all_events(run_time, IMT_eoe_counts):
    IMT_eoe_counts = IMT_eoe_counts.reset_index()
    all_events=pd.DataFrame(columns=['event' ,'Birthtime','Lifespan'])
    for ii in range(IMT_eoe_counts.shape[0]) :
        if IMT_eoe_counts['End_of_Expt'][ii] == False and IMT_eoe_counts['Death'][ii] == False:
            cell =  {'event' : 'mitosis', 'Birthtime' : IMT_eoe_counts['Birthtime'][ii], 'Lifespan' : IMT_eoe_counts['Lifespan'][ii]};
            all_events = all_events.append(cell, ignore_index=True)

        
        elif IMT_eoe_counts['End_of_Expt'][ii] == True and IMT_eoe_counts['Death'][ii] == False:
            if IMT_eoe_counts['Birthtime'][ii]==0:
                cell =  {'event' : 'quiescence', 'Birthtime' : IMT_eoe_counts['Birthtime'][ii], 'Lifespan' : IMT_eoe_counts['Lifespan'][ii]};
                all_events = all_events.append(cell, ignore_index=True)
            else:
                cell =  {'event' : 'End_of_Expt', 'Birthtime' : IMT_eoe_counts['Birthtime'][ii], 'Lifespan' : IMT_eoe_counts['Lifespan'][ii]};
                all_events = all_events.append(cell, ignore_index=True)
            
        elif IMT_eoe_counts['End_of_Expt'][ii] == False and IMT_eoe_counts['Death'][ii] == True:
           cell =  {'event' : 'apoptosis', 'Birthtime' : IMT_eoe_counts['Birthtime'][ii], 'Lifespan' : IMT_eoe_counts['Lifespan'][ii]};
           all_events = all_events.append(cell, ignore_index=True)
    
    all_events = all_events[all_events['Birthtime']+all_events['Lifespan']<=run_time] #mitosis after 120 hours in total doesn't count
    return all_events
