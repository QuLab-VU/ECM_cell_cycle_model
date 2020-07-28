#!/bin/sh
import subprocess
import os
import matlab.engine
import math 
# from StochCell import runCell

os.remove("tmp_state.xls")
os.remove("tmp_total.xls")

AKTs = [0.25]
ERKs = [0.25]
init = 2
FAK = 0.4
eng=matlab.engine.start_matlab()


for AKT in AKTs:
    for ERK in ERKs:
        eng.runStochCell(AKT,ERK,FAK,init,nargout=0)
        subprocess.check_call(['Rscript', 'getFits.R'], shell=True)
        akt_str=str(math.floor(AKT*100)%100).split('.')[0]
        erk_str=str(math.floor(ERK*100)%100).split('.')[0]
        os.rename('tmpplot.pdf',akt_str+erk_str+'.pdf')
        os.rename('tmp_r.Rdata',akt_str+erk_str+'.Rdata')