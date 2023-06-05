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
from Austin_proliferation_ODE_model_6 import model


from pysb import * 

###Initials
Parameter('Skp2_0', (0.0199387776457489 * 1))
Parameter('p21_0', (0.0212324326339528* 1))
Parameter('p21p_0', (0.104907517295920* 1))
Parameter('cdh1act_0', (0.970596980815461* 1))
Parameter('cdh1inact_0', ((1-0.970596980815461)* 1)) #cdh1tot - cdh1act
Parameter('Myc_0', (0.321067821067829* 1))
Parameter('CycD_0', (0.175953608509162* 1))
Parameter('CycDp21_0', (0.0138939800293163* 1))
Parameter('pRb_0', (0.301757475166483* 1))
Parameter('pRbp_0', (0.545620559818762* 1))
Parameter('E2F_0', (0.0403080965235495* 1))
Parameter('pRbE2F_0', (0.554487916752975* 1))
Parameter('CycE_0', (0.0640517424993558* 1))
Parameter('CycEp21_0', (0.00487528808175496* 1))
Parameter('CycA_0', (0.0182005140710040* 1))
Parameter('CycAp21_0', (0.00234453838504852 * 1))
Parameter('CycB_0', (0.0234757871244827* 1))
Parameter('CycBp21_0', (0.00435874152378834* 1))
Parameter('Cdc20_0', (0.0204771947463057 * 1))
Parameter('BAX_0', (0.999000159839911 * 1))
Parameter('BAXmBCL2_0', (0.000997752220497874 * 1))
Parameter('BH3_0', (5.81398691338436e-05* 1))
Parameter('BH3BCL2_0', (0.0332751934642001* 1))
Parameter('BAXm_0', ((1 - 0.999000159839911 - 0.000997752220497874)* 1)) #BAXt - BAX - BAXmBCL2
Parameter('BCL2_0', ((0.95 - 0.0332751934642001 - 0.000997752220497874) * 1)) #BCL2t - BH3BCL2 -BAXmBCL2


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

