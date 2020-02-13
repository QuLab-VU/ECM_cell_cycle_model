#-*- coding: utf-8 -*-
"""
Created on Mon Feb 3 14:57:38 2020

@author: Jérémie Marlhens
Title: Austin_proliferation_model
"""
from pysb import *

import sys
sys.path.insert(0, r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code\cellular_cycle") #add a pathway to the model
from Austin_proliferation_model import model


###Initials
Expression('Skp2_0', (0.0758759637154584 * model.parameters['omega']))
Expression('p21_0', (0.0179573231225730* model.parameters['omega']))
Expression('p21p_0', (0.0304907956309426* model.parameters['omega']))
Expression('cdh1act_0', (0.0868238125152800* model.parameters['omega']))
Expression('cdh1inact_0', ((1-0.0868238125152800)* model.parameters['omega'])) #cdh1tot - cdh1act
Expression('Myc_0', (0.358070500927645* model.parameters['omega']))
Expression('CycD_0', (0.246007623346968* model.parameters['omega']))
Expression('CycDp21_0', (0.0114732033609201* model.parameters['omega']))
Expression('pRb_0', (0.209618782280972* model.parameters['omega']))
Expression('pRbp_0', (0.596265732766763* model.parameters['omega']))
Expression('E2F_0', (0.0663628730831172* model.parameters['omega']))
Expression('pRbE2F_0', (0.659538609880404* model.parameters['omega']))
Expression('CycE_0', (0.0823229552545997* model.parameters['omega']))
Expression('CycEp21_0', (0.00601176670298263* model.parameters['omega']))
Expression('CycA_0', (0.0781546778431587* model.parameters['omega']))
Expression('CycAp21_0', (0.0110421142038035 * model.parameters['omega']))
Expression('CycB_0', (0.200858887607416* model.parameters['omega']))
Expression('CycBp21_0', (0.00209763607378452* model.parameters['omega']))
Expression('Cdc20_0', (0.367494461337630 * model.parameters['omega']))
Expression('BAX_0', (0.992676117362238 * model.parameters['omega']))
Expression('BAXmBCL2_0', (0.00730440429784349 * model.parameters['omega']))
Expression('BH3_0', (0.000428591174323906* model.parameters['omega']))
Expression('BH3BCL2_0', (0.199169986463653* model.parameters['omega']))
Expression('BAXm_0', ((1 - 0.992676117362238 - 0.00730440429784349)* model.parameters['omega'])) #BAXt - BAX - BAXmBCL2
Expression('BCL2_0', ((0.95 - 0.199169986463653 -0.00730440429784349) * model.parameters['omega'])) #BCL2t - BH3BCL2 -BAXmBCL2


Initial(model.monomers['Skp2'], model.expressions['Skp2_0'])
Initial(model.monomers['p21'](b=None, S='u'), model.expressions['p21_0'])
Initial(model.monomers['p21'](b=None, S='p'), model.expressions['p21p_0'])
Initial(model.monomers['cdh1'](A='yes'), model.expressions['cdh1act_0'])
Initial(model.monomers['cdh1'](A='no'), model.expressions['cdh1inact_0'])
Initial(model.monomers['Myc'](), model.expressions['Myc_0'])
Initial(model.monomers['Cdc20'](), model.expressions['Cdc20_0'])
Initial(model.monomers['CycA'](b=None), model.expressions['CycA_0'])
Initial(model.monomers['CycA'](b=1)%model.monomers['p21'](b=1, S='u'), model.expressions['CycAp21_0'])
Initial(model.monomers['CycB'](b=None), model.expressions['CycB_0'])
Initial(model.monomers['CycB'](b=1)%model.monomers['p21'](b=1, S='u'), model.expressions['CycBp21_0'])
Initial(model.monomers['CycD'](b=None), model.expressions['CycD_0'])
Initial(model.monomers['CycD'](b=1)%model.monomers['p21'](b=1, S='u'), model.expressions['CycDp21_0'])
Initial(model.monomers['CycE'](b=None), model.expressions['CycE_0'])
Initial(model.monomers['CycE'](b=1)%model.monomers['p21'](b=1, S='u'), model.expressions['CycEp21_0'])
Initial(model.monomers['pRb'](b=None, S='u'), model.expressions['pRb_0'])
Initial(model.monomers['pRb'](b=None, S='p'), model.expressions['pRbp_0'])
Initial(model.monomers['pRb'](b=1, S='u')%model.monomers['E2F'](b=1), model.expressions['pRbE2F_0'])
Initial(model.monomers['E2F'](b=None), model.expressions['E2F_0'])
Initial(model.monomers['BAX'](m='no', b=None), model.expressions['BAX_0'])
Initial(model.monomers['BAX'](m='yes', b=None), model.expressions['BAXm_0'])
Initial(model.monomers['BAX'](m='yes', b=1)%model.monomers['BCL2'](b=1), model.expressions['BAXmBCL2_0'])
Initial(model.monomers['BH3'](b=None), model.expressions['BH3_0'])
Initial(model.monomers['BH3'](b=1)%model.monomers['BCL2'](b=1), model.expressions['BH3BCL2_0'])
Initial(model.monomers['BCL2'](b=None), model.expressions['BCL2_0'])

