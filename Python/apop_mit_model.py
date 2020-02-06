#-*- coding: utf-8 -*-
"""
Created on Mon Feb 3 14:57:38 2020

@author: Jérémie Marlhens
Title: Austin_proliferation_model
"""

from pysb import*

Model()

###Monomers:

Monomer('Skp2', )
Monomer('p21', ['b', 'S'], {'S':['u', 'p']})
Monomer('cdh1', ['A'], {'A':['yes','no']}) #A=activated
Monomer('Myc', )
Monomer('CycD', ['b'])
Monomer('pRb', ['b', 'S'], {'S':['u', 'p']})
Monomer('E2F', ['b'])
Monomer('CycE', ['b'])
Monomer('CycA', ['b'])
Monomer('CycB',['b'] )
Monomer('Cdc20', )
Monomer('BAX', )
Monomer('BCL2', )
Monomer('BH3', )

'''
Monomer('YAP', )
Monomer('AKT', )
Monomer('ERK', )
'''

###Parameter

Parameter('AKT',1)
Parameter('ERK',1)
Parameter('FAK', 1)
Parameter('omega', 200)
Parameter('m', 3)
Parameter('BCL2t', 0.95) 
Parameter('BAXt', 1)
#Parameter(Stress, 2-AKT-ERK)
Parameter('kf1', 0.01)
Parameter('kf2', 2.4)
Parameter('kb', 0.14)
Parameter('kasbaxmbcl2', 75)
Parameter('kdsbaxmbcl2', 0.05)
Parameter('kasbh3bcl2',25)
Parameter('kdsbh3bcl2',0.01)
Parameter('k1sbh3',0.001)
Parameter('k2sbh3',0.009) #0.0095 #0.012 #0.025
Parameter('kdbh3',0.03)
Parameter('a', 1)  #basal term for activation of cyclin E/Cdk2 (without PF loop)) uM
Parameter('b', 1)  #basal term for cyclin B/Cdk1 (without PF loop)) uM
Parameter('s', 0.25) #0.4 #scaling parameter
Parameter('Hippo_tot', 1)
Parameter('Yap_tot', 1)
Parameter('AKT_tot', 1)
Parameter('ERK_tot',1)
Parameter('Cdc20_tot', 1) #total concentration of cdc20) uM
Parameter('FAK_tot', 1)
Parameter('skp2_tot', 1)
Parameter('p21_tot', 1)
Parameter('cdh1_tot', 1)
#Parameter(cmet, cmet) 
#Parameter(egfr, egfr)
#Parameter(stiff, stiff)
Parameter('CycD_tot', 1)
Parameter('CycE_tot', 1) #total concentration of cyclinE/Cdk2) uM
Parameter('CycA_tot', 1) #total concentration of cyclinA/Cdk2) uM
Parameter('CycB_tot', 1) #total concentration of cyclinB/Cdk1) uM
#K_dd, 0.1) #Michaelis constant for degradation of cyclinD/Cdk4-6 complex) uM
#K_gf, 0.1) #Michaelis constant for synthesis of cyclinD/Cdk4-6 induced by GF) uM
Parameter('K_1E2F', 0.5) #1  #Michaelis constant for E2F activation by cyclinD/Cdk4-6 and cyclinE/Cdk 2 complexes) uM
Parameter('K_2E2F', 0.5)  #0.5 #Michaelis constant for E2f inactivation by cyclinA/Cdk2 complex) uM
Parameter('K_1Me', 0.2) #0.1 #Michaelis constant for cyclinE/Cdk2 activation by E2F) uM
Parameter('K_2Me', 0.12) #0.1 #Michaelis constant for cyclinE/Cdk2 inactivation by cyclinA/Cdk2) uM
Parameter('K_1Ma', 0.2) #0.5 #Michaelis constant for cyclinA/Cdk2 activation by E2F) uM
Parameter('K_2Ma', 0.1) #0.1 #Michaelis constant for cyclinA/Cdk2 inactivation by Cdc20) uM
Parameter('K_1Mb', 0.3) #0.1 #Michaelis constant for cyclinB/Cdk1 activation by cyclinA/Cdk2) uM
Parameter('K_2Mb', 0.12) #0.1 #Michaelis constant for cyclinB/Cdk1 inactivation by Cdc20) uM
Parameter('K_1cdc20', 0.1)  #0.1 #Michaelis constant for Cdc20 activation through phosphorylation by cyclinB/Cdk1) uM
Parameter('K_2cdc20', 0.1)  #0.1 #MIchaelis constant for Cdc20 inactivation through dephosphorylation) uM
Parameter('K_ie', 0.2) #Michaelis constant for cyclinB/Cdk1 inactivation by cyclinE/Cdk2, representing the indeirect ecffecof the DNA relplication checkpoint) uM
#V_dd, 0.245)  #Maximum degradation rate of cyclinD/Cdk4-6 complex) uM/h
#V_sd, 0.175)  #Maximum rate for synthesis of cyclinD/Cdk4-6 induced by growth factor) uM/h
Parameter('V_1cdc20', 1.8)  #1.8  #Rate constant for activation of Cdc20 through phosphorylation by cyclinB/Cdk1) /h
Parameter('V_2cdc20', 0.6)  #0.6  #Rate for inactivation of Cdc20 through dephosphorylation) uM/h
Parameter('V_1E2F', 0.4) #0.4 #Rate constant for activation of E2F by cyclinD/Cdk4-6 and cyclinE/Cdk2 complexe) /h
Parameter('V_2E2F', 0.7) #0.7  #Rate constant for inactivation of E2F by cyclinA/Cdk2 complex) /h
Parameter('V_1Me', 1.1) #0.7 #1.1  #Rate for activation of cyclinE/Cdk2 by E2F) /uM*h
Parameter('V_2Me', 1.6) #0.5 #1.6  #Rate constant for inactivation fo clyclinE/Cdk2 by cyclinA/Cdk2) /h
Parameter('V_1Ma', 0.6) #0.6 #Rate constant for activation of cyclinA/Cdk2 by E2F) /h
Parameter('V_2Ma', 0.6) #0.6 #Rate constant for inactivation of cyclinA/Cdk2 by Cdc20) /h
Parameter('V_1Mb', 1)  #1 #Rate for activation of cyclinB/Cdk1 by cyclinA/Cdk2) /uM*h
Parameter('V_2Mb', 0.8) #0.8 #Rate constant for inactivation of cyclinB/Cdk1 by Cdc20) /h
Parameter('k_cd4', 0.08)
Parameter('k_myc1', 0.08)
Parameter('k_E', 0.4) #0.4
Parameter('k_M', 0.07) #0.27#0.1
Parameter('k_EMY', 0.3) #0.2 #0.12 #0.37#0.1
Parameter('k_ED', 0.22) #0.23#0.08
Parameter('k_CDF', 0.21) #0.18
Parameter('k_CD', 0.03) #0.03
Parameter('k_CDS', 0.2) #0.25#0.5
Parameter('k_R', 0.07) #.03 #0.18
Parameter('k_RE', 180) #180
Parameter('k_b', 0.003) #0.003
Parameter('K_S', 0.8) #0.5#0.1
Parameter('K_DF', 0.8)
Parameter('k_CE', 0.35)
Parameter('d_M', 0.7)
Parameter('d_E', 0.25)
Parameter('d_CD', 1.5) #1.5
Parameter('d_CE', 1.5) #1.5
Parameter('d_R', 0.06) #0.03
Parameter('d_RP', 0.06)
Parameter('d_RE', 0.03) #0.03
Parameter('k_P1', 18) #18
Parameter('k_P2', 18) #18
Parameter('k_DP', 3.6) #3.6
Parameter('K_M', 0.2) #0.5
Parameter('K_MY', 0.1) #0.15#0.1
Parameter('K_E', 0.15) #0.1
Parameter('K_CD', 0.92) #0.92
Parameter('K_CE', 0.92) #0.92
Parameter('K_RP', 0.01) #0.01
Parameter('K_i13', 0.8) #0.1
Parameter('K_i14', 5) #2
Parameter('K_dcescp2', 2)
Parameter('V_sac', 0.4) #0.4
Parameter('V_da', 0.35)
Parameter('V_sec', 0.7) #0.8
Parameter('V_de', 0.75)
Parameter('V_sf', 0.4) #0.4
Parameter('V_fs', 1.6)
Parameter('V_dfa', 0.35)
Parameter('V_sskp2', 0.01) #0.15
Parameter('V_s2skp', 0.01) #0.07
Parameter('V_dskp2', 1.1)
Parameter('V_1p21', 100)
Parameter('V_2p21', 0.1)
Parameter('V_dp21p', 5)
Parameter('K_dbcdc20', 0.3)
Parameter('K_dbcdh1', 0.3)
Parameter('K_db', 0.08) #0.07 #0.0005 #Michaelis constant for CycB deg
Parameter('V_1cdh1', 1.25) #1.25
Parameter('V_2cdh1', 8) #8
Parameter('V_db', 0.8) #0.08 #0.15 #0.06 #new value for cyc B deg
Parameter('v_scdh1a', 0.11)
Parameter('v_2fak', 0.1) #0.05
Parameter('v_s1p21', 0.5) #3#0.8
Parameter('v_s2p21', 0.1) #2#0.1
Parameter('K_1cdh1', 0.1)
Parameter('K_2cdh1', 0.1)
Parameter('K_dp21p', 0.1)
Parameter('K_dp21skp2', 0.1)
Parameter('K_cdh1', 0.4)
Parameter('K_dskp2', 0.5)
Parameter('K_AFAK', 0.1)
Parameter('K_dfa', 0.1)
Parameter('K_fs', 0.1)
Parameter('K_cf', 0.1)
Parameter('K_cac', 0.1)
Parameter('K_cec', 0.1)
Parameter('K_da', 0.1)
Parameter('K_de', 0.1)
Parameter('K_1p21', 0.5)
Parameter('K_2p21', 0.5)
Parameter('k_ddp21', 0.06)
Parameter('k_ddskp2', 0.005)
Parameter('k_dcdh1a', 0.1)
Parameter('k_c1', 3) #0.15
Parameter('k_c2', 1) #0.05
Parameter('k_c3', 8) #0.2
Parameter('k_c4', 2) #0.1
Parameter('k_c5', 5) #0.15
Parameter('k_c6', 1) #0.125
Parameter('k_c7', 0.12)
Parameter('k_c8', 0.2)
Parameter('V_Fp', 1.5)
Parameter('K_Fp', 0.1)

##Positive feedback loops on the G1/S and G2/M transitions
Parameter('b1', 0) #term representing the self-activation of cyclinE/Cdk2 through the mutual activation between cyclinE/Cdk2 and its phosphatase Cdc25(PF loop in the G1/S transition)) without=0 with=1
Parameter('b2', 0) #term representing the self-activation of cyclinB/Cdk1 through mutual activation between cyclin B/cdk1 and its phosphatase cdk25(PF loop in the G2/M transition)) without=0 with=1
Parameter('K_ib', 0.5) #0.5

###Observables
Observable('obscdh1tot', cdh1())
Observable('obsSkp2', Skp2())
Observable('obsp21', p21(b=None, S='u'))
Observable('obsp21p', p21(b=None, S='p'))
Observable('obscdh1act', cdh1(A='yes'))
Observable('obscdh1inact', cdh1(A='no'))
Observable('obsMyc', Myc())
Observable('obspRbp', pRb(b=None, S='p'))
Observable('obspRb', pRb(b=None, S='u'))
Observable('obspRbE2F', pRb(b=1, S='u')%E2F(b=1))
Observable('obsE2F', E2F(b=None))


###Expressions
##Skp2
Expression('ks_Skp2', s * V_sskp2 * omega)
Expression('kd_Skp2', s * k_ddskp2)
Expression('ks_Skp2_FAK', s * V_s2skp * FAK)
Expression('kd_Skp2_cdh1', s * V_dskp2 * omega / (( K_dskp2 * omega + obsSkp2 )*( K_cdh1 * omega + obscdh1act)))

##p21
Expression('ks_p21', s * v_s1p21 * omega)
Expression('kasso_p21_CycD', s * (k_c1 / omega))
Expression('kdiss_p21_CycD', s * k_c2)
Expression('kasso_p21_CycE', s * (k_c3 / omega))
Expression('kdiss_p21_CycE', s * k_c4)
Expression('kasso_p21_CycA', s * (k_c5 / omega))
Expression('kdiss_p21_CycA', s * k_c6)
Expression('kasso_p21_CycB', s * (k_c7 / omega))
Expression('kdiss_p21_CycB', s * k_c8)
Expression('kp_p21_CycE', s * V_1p21 / (K_1p21 * omega + obsp21))
Expression('kunph_p21p', s * ((V_2p21 * omega) / (K_2p21 * omega + obsp21p)))
Expression('kp_p21_FAK', s * ((V_Fp * FAK)/ ((K_Fp * omega) + obsp21)))
Expression('kd_p21', s * (k_ddp21))
Expression('kd_p21p_Skp2', s * (V_dp21p * omega * (1 / (K_dp21skp2 * omega + obsSkp2)) * (1 / (K_dp21p * omega + obsp21p))))

##cdh1
Expression('kautoact_cdh1', s * (V_1cdh1 * omega * ((1) / (K_1cdh1 * omega + (obscdh1inact)))))
Expression('k_inact_CycAorB', s * (V_2cdh1 * (1 / ( K_2cdh1 * omega + obscdh1act))))

##Myc
Expression('ks_Myc_AKT', s * ((k_M * AKT * omega) / ((K_S * omega) + AKT)))
Expression('ks_Myc_ERK', s * ((k_EMY * ERK * omega) / ((K_MY * omega) + ERK)))
Expression('ks_Myc_FAK', s * ((k_CDF * FAK * omega) / ((K_DF * omega) + FAK)))
Expression('kd_Myc', s * d_M)

##CycD
Expression('ks_CycD_AKT', s * ((k_CDS * AKT * omega) / ((K_S * omega) + AKT)))
Expression('ks_CycD_ERK', s * ((k_ED * AKT * omega) / ((K_MY * omega) + AKT)))
Expression('ks_CycD_Myc', s * ((k_CD * omega) / ((K_M * omega) + obsMyc)))
Expression('kd_CycD', s * d_CD)

##pRb and pRbp

Expression('ks_pRb', s * (k_R * omega))
Expression('kunph_pRbp', s * ((k_DP * omega) / ((K_RP * omega) + obspRbp)))
Expression('kf_RbE2F', s * (k_RE / omega))
Expression('kp_pRb_CycD', s * ((k_P1) / ((omega * K_CD) + obspRb)))
Expression('kp_pRb_CycE', s * ((k_P2) / ((omega * K_CD) + obspRb)))
Expression('kp_pRbE2F_CycD', s * ((k_P1) / ((omega * K_CD) + obspRbE2F)))
Expression('kp_pRbE2F_CycE', s * ((k_P2) / ((omega * K_CD) + obspRbE2F)))
Expression('kd_pRb', s * d_R)
Expression('kd_pRbp', s * d_RP)

##E2F 
Expression('ks_E2F_E2FMyc', s * (k_E * omega * (1 /((K_M * omega) + obsMyc)) * (1 / ((K_E * omega) + obsE2F))))
Expression('ks_E2F_Myc', s * ((k_b * omega) / ((K_M * omega) + obsMyc)))
Expression('kd_E2F', s * d_E)
Expression('kd_E2F_CycA', s * (((V_2E2F) / ((K_2E2F * omega) + obsE2F))))

#CycE
#Expression('ks_CycE', s * (V_1Me / omega * E2F * (a * omega + b1 * CycE) * ((CycE_tot * omega) - CycE - CycEp21) / (K_1Me * omega + ((CycE_tot * omega) - CycE - CycEp21))))


###Rules
##Skp2
Rule('Skp2_synthesis', None >> Skp2(), ks_Skp2) #Sytnhesis of Skp2
Rule('Skp2_degradation', Skp2() >> None, kd_Skp2) #Degradation of Skp2
Rule('Skp2_synthesis_FAK', None >> Skp2() , ks_Skp2_FAK) #Skp2 synthesis from FAK
Rule('Skp2_degradation_cdh1', Skp2() + cdh1() >> cdh1(), kd_Skp2_cdh1) #degradation of Skp2 by cdh1

##p21 and p21p
Rule('p21_synthesis', None >> p21(b=None, S='u'), ks_p21) #Sytnhesis of p21
Rule('p21_CycD_association', p21(b=None, S='u') + CycD(b=None) | p21(b=1, S='u')%CycD(b=1), kasso_p21_CycD, kdiss_p21_CycD) #p21 and CycD association and dissociation
Rule('p21_CycA_association', p21(b=None, S='u') + CycA(b=None) | p21(b=1, S='u')%CycA(b=1), kasso_p21_CycA, kdiss_p21_CycA) #p21 and CycA association and dissociation
Rule('p21_CycB_association', p21(b=None, S='u') + CycB(b=None) | p21(b=1, S='u')%CycB(b=1), kasso_p21_CycB, kdiss_p21_CycB) #p21 and CycB association and dissociation
Rule('p21_phosphorylation_CycE', p21(b=None, S='u') + CycE(b=None) >> p21(b=None, S='p') + CycE(b=None), kp_p21_CycE) #Phosphorylation of p21 by CycE
Rule('p21p_unphosphorylation', p21(b=None, S='p') >> p21(b=None, S='u'), kunph_p21p) #Unphosphorylation of p21
Rule('p21_phosphorylation_FAK', p21(b=None, S='u') >> p21(b=None, S='p'), kp_p21_FAK) #Phosphorylation of p21 by FAK
Rule('p21_degradation', p21(b=None, S='u') >> None, kd_p21) #p21 non specific degradation
Rule('p21p_degradation', p21(b=None, S='p') >> None, kd_p21) #p21p non specific degradation
Rule('p21p_degradation_Skp2',  p21(b=None, S='p') + Skp2() >> Skp2(), kd_p21p_Skp2) #p21p degradation from skp2

##cdh1
Rule('cdh1_autoactivation', cdh1(A='no') >> cdh1(A='yes'), kautoact_cdh1) #cdh1 autoactivation
Rule('cdh1_inhibiton_CycA', cdh1(A='yes') + CycA(b=None) >> cdh1(A='no') + CycA(b=None) ,k_inact_CycAorB) #cdh1 inactivation by Cyclin A
Rule('cdh1_inhibiton_CycB', cdh1(A='yes') + CycB(b=None) >> cdh1(A='no') + CycB(b=None) ,k_inact_CycAorB) #cdh1 inactivation by Cyclin B

##Myc
Rule('Myc_synthesis_AKT', None >> Myc(), ks_Myc_AKT) #Sytnhesis of Myc by AKT
Rule('Myc_synthesis_FAK', None >> Myc(), ks_Myc_FAK) #Sytnhesis of Myc by FAK
Rule('Myc_synthesis_ERK', None >> Myc(), ks_Myc_ERK) #Sytnhesis of Myc by ERK
Rule('Myc_degradation', Myc() >> None, kd_Myc) #Degradation of Myc

##CycD
Rule('Cyc_synthesis_AKT', None >> CycD(b=None), ks_CycD_AKT) #Sytnhesis of Cyclin D by AKT
Rule('Cyc_synthesis_ERK', None >> CycD(b=None), ks_CycD_ERK) #Sytnhesis of Cyclin D by ERK
Rule('Cyc_synthesis_Myc', Myc() >> CycD(b=None) + Myc(), ks_CycD_Myc) #Sytnhesis of Cyclin D by Myc
Rule('CycD_degradation', CycD(b=None) >> None, kd_CycD) #Degradation of CycD

##pRb and pRbp
Rule('pRb_synthesis', None >> pRb(b=None, S='u'), ks_pRb) #Sytnhesis of Rb
Rule('pRbp_unphosphorylation', pRb(b=None, S='p') >> pRb(b=None, S='p'), kunph_pRbp) #Unphosphorylation of pRb
Rule('pRbE2F_complex_formation', pRb(b=None, S='u') + E2F(b=None) >> pRb(b=1, S='u')%E2F(b=1), kf_RbE2F)
Rule('pRb_phosphorylation_CycE', pRb(b=None, S='u') + CycE(b=None) >> pRb(b=None, S='p') + CycE(b=None), kp_pRb_CycE) #Phosphorylation of pRb by Cyclin E
Rule('pRb_phosphorylation_CycD', pRb(b=None, S='u') + CycD(b=None) >> pRb(b=None, S='p') + CycD(b=None), kp_pRb_CycE) #Phosphorylation of pRb by Cyclin D
Rule('pRbE2F_phosphorylation_CycE', pRb(b=1, S='u')%E2F(b=1) + CycE(b=None) >> pRb(b=1, S='p')%E2F(b=1) + CycE(b=None), kp_pRbE2F_CycE) #Phosphorylation of pRbE2F by Cyclin E
Rule('pRbE2F_phosphorylation_CycD', pRb(b=1, S='u')%E2F(b=1) + CycD(b=None) >> pRb(b=1, S='p')%E2F(b=1) + CycD(b=None), kp_pRbE2F_CycE) #Phosphorylation of pRbE2F by Cyclin D
Rule('pRb_degradation', pRb(b=None, S='u') >> None, kd_pRb) #Degradation of pRb
Rule('pRbp_degradation', pRb(b=None, S='p') >> None, kd_pRbp) #Degradation of pRb

##E2F
Rule('E2F_synthesis_E2FMyc', Myc() + E2F(b=None) >> Myc() + E2F(b=None) + E2F(b=None), ks_E2F_E2FMyc) #e2f synthesis from myc e2f autoactivation
Rule('E2F_synthesis_Myc', Myc() >> Myc() + E2F(b=None), ks_E2F_Myc) #e2f synthesis from myc e2f autoactivation
Rule('E2F_degradation', E2F(b=None) >> None, kd_E2F) #Degradation of E2F
Rule('E2F_degradation_CycA', E2F(b=None) + CycA(b=None) >> CycA(b=None), kd_E2F_CycA) ##Degradation of E2F by Cyclin A

##CycE






###Initials

Parameter('Skp2_0', 0.0758759637154584)
Parameter('p21_0', 0.0179573231225730)
Parameter('p21p_0', 0.0304907956309426)
Parameter('cdh1act_0', 0.0868238125152800)
Parameter('cdh1inact_0', 0.0868238125152800)
Parameter('Myc_0', 0.358070500927645)
Parameter('CycD_0', 0.246007623346968)
Parameter('CycDp21_0', 0.0114732033609201)
Parameter('pRb_0', 0.209618782280972)
Parameter('pRbp_0', 0.596265732766763)
Parameter('E2F_0', 0.0663628730831172)
Parameter('pRbE2F_0', 0.659538609880404)
Parameter('CycE_0', 0.0823229552545997)
Parameter('CycEp21_0', 0.00601176670298263)
Parameter('CycA_0', 0.0781546778431587)
Parameter('CycAp21_0', 0.0110421142038035)
Parameter('CycB_0', 0.200858887607416)
Parameter('CycBp21_0', 0.00209763607378452)
Parameter('Ccdc20_0', 0.367494461337630)
Parameter('BAX_0', 0.992676117362238)
Parameter('BAXmBCL2_0', 0.00730440429784349)
Parameter('BH3_0', 0.000428591174323906)
Parameter('BH3BCL2_0', 0.199169986463653)


Initial(Skp2(), Skp2_0)
Initial(cdh1(A='yes'), cdh1act_0)

###Simulation

