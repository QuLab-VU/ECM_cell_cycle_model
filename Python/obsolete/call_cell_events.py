# -*- coding: utf-8 -*-
"""
Modified from plot_frac_div_q_ap1.py
@author: Darren Tyson
Function designed to filter single-cell cell cycle position data for fracprolif
run_time = duration of experiment (e.g., end time)
IMT_eoe_counts = Pandas DataFrame with columns: End_of_Expt, Death Birthtime, Lifespan
new column (event) is cell fate decision call
save output with only 'event','Birthtime','Lifespan' columns
"""
import pandas as pd
import numpy as np

def call_events(run_time, IMT_eoe_counts):
    IMT_eoe_counts = IMT_eoe_counts.reset_index()
    IMT_eoe_counts['event'] = "undefined"
    
    IMT_eoe_counts.loc[IMT_eoe_counts.Death == True, 'event'] = "apoptosis"
    IMT_eoe_counts.loc[(IMT_eoe_counts.Death == False) & 
                       (IMT_eoe_counts.End_of_Expt == False) & 
                       (IMT_eoe_counts.Birthtime<=run_time), 'event'] = "mitosis"
    IMT_eoe_counts.loc[(IMT_eoe_counts.Death == False) & 
                       (IMT_eoe_counts.End_of_Expt == True) & 
                       (IMT_eoe_counts.Birthtime==0), 'event'] = "quiescence"
    IMT_eoe_counts.loc[(IMT_eoe_counts.Death == False) & 
                       (IMT_eoe_counts.End_of_Expt == True) & 
                       (IMT_eoe_counts.Birthtime!=0), 'event'] = "End_of_Expt"
    IMT_eoe_counts.loc[IMT_eoe_counts['Starting cell'] == True, 'event'] = "Starting cell"

    out = IMT_eoe_counts[['event','Birthtime','Lifespan']]
    out = out[out['Birthtime']+out['Lifespan'] <= run_time]
    return(out)