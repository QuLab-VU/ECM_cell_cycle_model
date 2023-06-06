#-*- coding: utf-8 -*-
"""
Created on Mon Feb 3 14:57:38 2020

@author: Jérémie Marlhens
Title: Austin_proliferation_model
"""

path="/home/jeremie"
#path=r"C:\Users\Grim.M\Documents\git_repos\ECM_cell_cycle_model\Python"
import sys
sys.path.insert(0, path+"/cellular_cycle") #add a pathway to the model
from Austin_proliferation_model_6 import model


from pysb import * 

###Initials
Parameter('Skp2_0', int(round(0.06 * 200)))
Parameter('p21_0', int(round(0.005 * 200)))
Parameter('p21p_0', int(round(0.08 * 200)))
Parameter('cdh1act_0', int(round(0.03 * 200)))
Parameter('cdh1inact_0', int(round((1-0.03) * 200))) #cdh1tot - cdh1act
Parameter('Myc_0', int(round(0.695 * 200)))
Parameter('CycD_0', int(round(0.235 * 200)))
Parameter('CycDp21_0', int(round(0.005 * 200)))
Parameter('pRb_0', int(round(0.025 * 200)))
Parameter('pRbp_0', int(round(1.21 * 200)))
Parameter('E2F_0', int(round(0.62 * 200)))
Parameter('pRbE2F_0', int(round(0.095 * 200)))
Parameter('CycE_0', int(round(0.815 * 200)))
Parameter('CycEp21_0', int(round(0.015 * 200)))
Parameter('CycA_0', int(round(0.51 * 200)))
Parameter('CycAp21_0', int(round(0.01  * 200)))
Parameter('CycB_0', int(round(0.085 * 200)))
Parameter('CycBp21_0', int(round(0.01 * 200)))
Parameter('Cdc20_0', int(round(0.03  * 200)))
Parameter('BAX_0', int(round(0.999000159839911  * 200)))
Parameter('BAXmBCL2_0', int(round(0.000997752220497874  * 200)))
Parameter('BH3_0', int(round(5.81398691338436e-05 * 200)))
Parameter('BH3BCL2_0', int(round(0.0332751934642001 * 200)))
Parameter('BAXm_0', int(round((1 - 0.999000159839911 - 0.000997752220497874) * 200))) #BAXt - BAX - BAXmBCL2
Parameter('BCL2_0', int(round((0.95 - 0.0332751934642001 - 0.000997752220497874)  * 200))) #BCL2t - BH3BCL2 -BAXmBCL2



Initial(model.monomers['Skp2'], model.parameters['Skp2_0'])
Initial(model.monomers['p21'](b=None, S='u'), model.parameters['p21_0'])
Initial(model.monomers['p21'](b=None, S='p'), model.parameters['p21p_0'])
Initial(model.monomers['cdh1'](A='yes'), model.parameters['cdh1act_0'])
Initial(model.monomers['cdh1'](A='no'), model.parameters['cdh1inact_0'])
Initial(model.monomers['Myc'](), model.parameters['Myc_0'])
Initial(model.monomers['CycD'](b=None), model.parameters['CycD_0'])
Initial(model.monomers['CycD'](b=1)%model.monomers['p21'](b=1, S='u'), model.parameters['CycDp21_0'])
Initial(model.monomers['pRb'](b=None, S='u'), model.parameters['pRb_0'])
Initial(model.monomers['pRb'](b=None, S='p'), model.parameters['pRbp_0'])
Initial(model.monomers['E2F'](b=None), model.parameters['E2F_0'])
Initial(model.monomers['pRb'](b=1, S='u')%model.monomers['E2F'](b=1), model.parameters['pRbE2F_0'])
Initial(model.monomers['CycE'](b=None), model.parameters['CycE_0'])
Initial(model.monomers['CycE'](b=1)%model.monomers['p21'](b=1, S='u'), model.parameters['CycEp21_0'])
Initial(model.monomers['CycA'](b=None), model.parameters['CycA_0'])
Initial(model.monomers['CycA'](b=1)%model.monomers['p21'](b=1, S='u'), model.parameters['CycAp21_0'])
Initial(model.monomers['CycB'](b=None), model.parameters['CycB_0'])
Initial(model.monomers['CycB'](b=1)%model.monomers['p21'](b=1, S='u'), model.parameters['CycBp21_0'])
Initial(model.monomers['Cdc20'](), model.parameters['Cdc20_0'])
Initial(model.monomers['BAX'](m='no', b=None), model.parameters['BAX_0'])
Initial(model.monomers['BAX'](m='yes', b=1)%model.monomers['BCL2'](b=1), model.parameters['BAXmBCL2_0'])
Initial(model.monomers['BH3'](b=None), model.parameters['BH3_0'])
Initial(model.monomers['BH3'](b=1)%model.monomers['BCL2'](b=1), model.parameters['BH3BCL2_0'])
Initial(model.monomers['BAX'](m='yes', b=None), model.parameters['BAXm_0'])
Initial(model.monomers['BCL2'](b=None), model.parameters['BCL2_0'])


