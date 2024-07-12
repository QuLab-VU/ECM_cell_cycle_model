"""
Modified from Austin_proliferation_model_6.py
@author: Jérémie Marlhens
@author: Darren Tyson
@author: Leonard Harris
Title: cell_cycle_death_model
"""

from pysb import *
from pysb.simulator import ScipyOdeSimulator, BngSimulator
from pylab import linspace, plot, xlabel, ylabel, legend, figure, show
import sympy as sp

Model()

# MONOMERS
Monomer('Skp2')
Monomer('p21', ['b', 'S'], {'S': ['u', 'p']})
Monomer('cdh1', ['A'], {'A': ['yes', 'no']})  # A=activated
Monomer('Myc')
Monomer('CycD', ['b'])
Monomer('pRb', ['b', 'S'], {'S': ['u', 'p']})
Monomer('E2F', ['b'])
Monomer('CycE', ['b'])
Monomer('CycA', ['b'])
Monomer('CycB', ['b'])
Monomer('Cdc20')
Monomer('BAX', ['b', 'm'], {'m': ['yes', 'no']})
Monomer('BCL2', ['b'])
Monomer('BH3', ['b'])

Monomer('CycDi')
Monomer('CycEi')
Monomer('CycAi')
Monomer('CycBi')
Monomer('E2Fi')
Monomer('Cdc20i')
# PARAMETERS
# Stiffness, Drug: [AKT, ERK, FAK]
# Stiff, DMSO: 1, 1 ,1
# Interm, DMSO: 1, 1, 0.6
# Soft, DMSO: 0.6, 1, 0.4
# Stiff, ERL: 0.6, 0.2, 1
# Interm, ERL: 0.6, 0.2, 0.6
# Soft, ERL: 0.3, 0.2, 0.4
Parameter('AKT', 0.6)  # µM
Parameter('ERK', 0.2)  # µM
Parameter('FAK', 0.6)  # µM

Parameter('omega', 200)  # 1/µM
Parameter('BCL2t', 0.95)  # µM
Parameter('BAXt', 1)  # µM
# TODO: Is the scaling parameter necessary?
Parameter('s', 0.4)  # 0.4 # scaling parameter
# TODO: Is the 'm' parameter necessary? Are 's' and 'm' related to cytosolic and mitochondrial volumes?
Parameter('m', 3)  # 3 # unitless
Parameter('Cdc20_tot', 1)  # total concentration of cdc20 (µM)
Parameter('CycE_tot', 1)  # total concentration of cyclinE/Cdk2 (µM)
Parameter('CycA_tot', 1)  # total concentration of cyclinA/Cdk2 (µM)
Parameter('CycB_tot', 1)  # total concentration of cyclinB/Cdk1 (µM)

# INITIAL CONDITIONS (all µM units)
Parameter('Skp2_0', 0.06)
Parameter('p21_0', 0.005)
Parameter('p21p_0', 0.08)
Parameter('cdh1act_0', 0.03)
Parameter('cdh1inact_0', (1-0.03))  # cdh1tot - cdh1act
Parameter('Myc_0', 0.695)
Parameter('CycD_0', 0.235)
Parameter('CycDp21_0', 0.005)
Parameter('pRb_0', 0.025)
Parameter('pRbp_0', 1.21)
Parameter('E2F_0', 0.62)
Parameter('pRbE2F_0', 0.095)
Parameter('CycE_0', 0.815)
Parameter('CycEp21_0', 0.015)
Parameter('CycA_0', 0.51)
Parameter('CycAp21_0', 0.01)
Parameter('CycB_0', 0.085)
Parameter('CycBp21_0', 0.01)
Parameter('Cdc20_0', 0.03)
Parameter('BAX_0', 0.999000159839911)
Parameter('BAXmBCL2_0', 0.000997752220497874)
Parameter('BH3_0', 5.81398691338436e-05)
Parameter('BH3BCL2_0', 0.0332751934642001)
Parameter('CycDi_0', 0)
Parameter('CycEi_0', 0)
Parameter('CycAi_0', 0)
Parameter('CycBi_0', 0)
Parameter('E2Fi_0', 0)
Parameter('Cdc20i_0', 0)
# TODO: What are BAX, BAXmBCL2, and BH3BCL2? Steady state values? What parameter values?
Expression('BAXm_0_expr', (BAXt - 0.999000159839911 - 0.000997752220497874) * omega)  # BAXt - BAX - BAXmBCL2
Expression('BCL2_0_expr', (BCL2t - 0.0332751934642001 - 0.000997752220497874) * omega)  # BCL2t - BH3BCL2 -BAXmBCL2

Initial(Skp2(), Expression('Skp2_0_expr', Skp2_0 * omega))
Initial(p21(b=None, S='u'), Expression('p21_0_expr', p21_0 * omega))
Initial(p21(b=None, S='p'), Expression('p21p_0_expr', p21p_0 * omega))
Initial(cdh1(A='yes'), Expression('cdh1act_0_expr', cdh1act_0 * omega))
Initial(cdh1(A='no'), Expression('cdh1inact_0_expr', cdh1inact_0 * omega))
Initial(Myc(), Expression('Myc_0_expr', Myc_0 * omega))
Initial(CycD(b=None), Expression('CycD_0_expr', CycD_0 * omega))
Initial(CycD(b=1) % p21(b=1, S='u'), Expression('CycDp21_0_expr', CycDp21_0 * omega))
Initial(pRb(b=None, S='u'), Expression('pRb_0_expr', pRb_0 * omega))
Initial(pRb(b=None, S='p'), Expression('pRbp_0_expr', pRbp_0 * omega))
Initial(E2F(b=None), Expression('E2F_0_expr', E2F_0 * omega))
Initial(pRb(b=1, S='u') % E2F(b=1), Expression('pRbE2F_0_expr', pRbE2F_0 * omega))
Initial(CycE(b=None), Expression('CycE_0_expr', CycE_0 * omega))
Initial(CycE(b=1) % p21(b=1, S='u'), Expression('CycEp21_0_expr', CycEp21_0 * omega))
Initial(CycA(b=None), Expression('CycA_0_expr', CycA_0 * omega))
Initial(CycA(b=1) % p21(b=1, S='u'), Expression('CycAp21_0_expr', CycAp21_0 * omega))
Initial(CycB(b=None), Expression('CycB_0_expr', CycB_0 * omega))
Initial(CycB(b=1) % p21(b=1, S='u'), Expression('CycBp21_0_expr', CycBp21_0 * omega))
Initial(Cdc20(), Expression('Cdc20_0_expr', Cdc20_0 * omega))
Initial(BAX(m='no', b=None), Expression('BAX_0_expr', BAX_0 * omega))
Initial(BAX(m='yes', b=1) % BCL2(b=1), Expression('BAXmBCL2_0_expr', BAXmBCL2_0 * omega))
Initial(BH3(b=None), Expression('BH3_0_expr', BH3_0 * omega))
Initial(BH3(b=1) % BCL2(b=1), Expression('BH3BCL2_0_expr', BH3BCL2_0 * omega))
Initial(BAX(m='yes', b=None), BAXm_0_expr)
Initial(BCL2(b=None), BCL2_0_expr)
Initial(CycDi(), Expression('CycDi_0_expr', CycDi_0 * omega))
Initial(CycEi(), Expression('CycEi_0_expr', CycEi_0 * omega))
Initial(CycAi(), Expression('CycAi_0_expr', CycAi_0 * omega))
Initial(CycBi(), Expression('CycBi_0_expr', CycBi_0 * omega))
Initial(E2Fi(), Expression('E2Fi_0_expr', E2Fi_0 * omega))
Initial(Cdc20i(), Expression('Cdc20i_0_expr', Cdc20i_0 * omega))

# OBSERVABLES
Observable('obsSkp2', Skp2())
Observable('obsp21', p21(b=None, S='u'))
Observable('obsp21p', p21(b=None, S='p'))
Observable('obscdh1act', cdh1(A='yes'))
Observable('obscdh1inact', cdh1(A='no'))
Observable('obsMyc', Myc())
Observable('obsCycD', CycD(b=None))
Observable('obsCycDp21', p21(b=1, S='u') % CycD(b=1))
Observable('obspRb', pRb(b=None, S='u'))
Observable('obspRbp', pRb(b=None, S='p'))
Observable('obsE2F', E2F(b=None))
Observable('obspRbE2F', pRb(b=1, S='u') % E2F(b=1))
Observable('obsCycE', CycE(b=None))
Observable('obsCycEp21', p21(b=1, S='u') % CycE(b=1))
Observable('obsCycA', CycA(b=None))
Observable('obsCycAp21', p21(b=1, S='u') % CycA(b=1))
Observable('obsCycB', CycB(b=None))
Observable('obsCycBp21', p21(b=1, S='u') % CycB(b=1))
Observable('obsCdc20', Cdc20())
Observable('obsBAX', BAX(b=None, m='no'))
Observable('obsBAXmBCL2', BAX(b=1, m='yes') % BCL2(b=1))
Observable('obsBH3', BH3(b=None))
Observable('obsBH3BCL2', BH3(b=1) % BCL2(b=1))
Observable('obsBAXm', BAX(b=None, m='yes'))
Observable('obsBCL2', BCL2(b=None))
Observable('obsCycDi', CycDi())
Observable('obsCycEi', CycEi())
Observable('obsCycAi', CycAi())
Observable('obsCycBi', CycBi())
Observable('obsE2Fi', E2Fi())
Observable('obsCdc20i', Cdc20i())

# RULES

# ############
#    Skp2
# ############
Parameter('V_sskp2', 0.01)  # 0.15 (µM/h)
Expression('ks_Skp2', s * V_sskp2 * omega)
Rule('Skp2_synthesis', None >> Skp2(), ks_Skp2)  # µM/h

Parameter('k_ddskp2', 0.005)  # 1/h
Expression('kd_Skp2', s * k_ddskp2)
Rule('Skp2_degradation', Skp2() >> None, kd_Skp2)  # µM/h

Parameter('V_s2skp', 0.01)  # 0.07 (1/h)
Expression('ks_Skp2_FAK', s * V_s2skp * (FAK * omega))
Rule('Skp2_synthesis_FAK', None >> Skp2(), ks_Skp2_FAK)

Parameter('V_dskp2', 1.1)  # µM/h
Parameter('K_dskp2', 0.5)  # µM
Parameter('K_cdh1', 0.4)  # µM
Expression('kd_Skp2_cdh1', s * V_dskp2 * omega / ((K_dskp2 * omega + obsSkp2) * (K_cdh1 * omega + obscdh1act)))  # /µM-h
Rule('Skp2_degradation_cdh1', Skp2() + cdh1(A='yes') >> cdh1(A='yes'), kd_Skp2_cdh1)

# ############
# p21 and p21p
# ############
Parameter('v_s1p21', 0.7)  # µM^[(n+1)/n]/h # updated 20230808
Parameter('K_p21fak', 0.72)  # µM # updated 20230919
Parameter('F_n', 8)  # unitless # updated 20230808
# TODO: Ask about the commented out equation below. I changed it to the equation on the subsequent line.
# vv Austin's original equation vv
# Expression('ks_p21', s * ((v_s1p21 * omega) * (omega ** F_n) / (omega ** F_n + K_p21fak * (FAK * omega) ** F_n)))  # µM/h
# vv Leonard's modified equation vv
Expression('ks_p21', s * (v_s1p21 ** F_n * omega / (K_p21fak ** F_n + FAK ** F_n)))  # µM/h
Rule('p21_synthesis', None >> p21(b=None, S='u'), ks_p21)

Parameter('k_c1', 3)  # 0.15 # 1/µM-h
Parameter('k_c2', 1)  # 0.05 # /h
Expression('kasso_p21_CycD', s * (k_c1 / omega))
Expression('kdiss_p21_CycD', s * k_c2)
Rule('p21_CycD_association', p21(b=None, S='u') + CycD(b=None) | p21(b=1, S='u') % CycD(b=1),
     kasso_p21_CycD, kdiss_p21_CycD)

Parameter('k_c5', 5)  # 0.15 # 1/µM-h
Parameter('k_c6', 1)  # 0.125 # /h
Expression('kasso_p21_CycA', s * (k_c5 / omega))
Expression('kdiss_p21_CycA', s * k_c6)
Rule('p21_CycA_association', p21(b=None, S='u') + CycA(b=None) | p21(b=1, S='u') % CycA(b=1),
     kasso_p21_CycA, kdiss_p21_CycA)

Parameter('k_c7', 0.12)  # 1/µM-h
Parameter('k_c8', 0.2)  # /h
Expression('kasso_p21_CycB', s * (k_c7 / omega))
Expression('kdiss_p21_CycB', s * k_c8)
Rule('p21_CycB_association', p21(b=None, S='u') + CycB(b=None) | p21(b=1, S='u') % CycB(b=1),
     kasso_p21_CycB, kdiss_p21_CycB)

Parameter('k_c3', 8)  # 0.2 # 1/µM-h
Parameter('k_c4', 2)  # 0.1 # /h
Expression('kasso_p21_CycE', s * (k_c3 / omega))
Expression('kdiss_p21_CycE', s * k_c4)
Rule('p21_CycE_association', p21(b=None, S='u') + CycE(b=None) | p21(b=1, S='u') % CycE(b=1),
     kasso_p21_CycE, kdiss_p21_CycE)

Parameter('V_1p21', 100)  # /h
Parameter('K_1p21', 0.5)  # µM
Expression('kp_p21_CycE', s * V_1p21 / (K_1p21 * omega + obsp21))
Rule('p21_phosphorylation_CycE', p21(b=None, S='u') + CycE(b=None) >> p21(b=None, S='p') + CycE(b=None), kp_p21_CycE)

Parameter('V_2p21', 0.1)  # µM/h
Parameter('K_2p21', 0.5)  # µM
Expression('kunph_p21p', s * ((V_2p21 * omega) / (K_2p21 * omega + obsp21p)))
Rule('p21p_unphosphorylation', p21(b=None, S='p') >> p21(b=None, S='u'), kunph_p21p)

Parameter('k_ddp21', 0.06)  # /h
Expression('kd_p21', s * k_ddp21)
Rule('p21_degradation', p21(b=None) >> None, kd_p21)  # p21 and p21p non specific degradation

Parameter('V_dp21p', 5)  # µM/h
Parameter('K_dp21skp2', 0.1)  # µM
Parameter('K_dp21p', 0.1)  # µM
Expression('kd_p21p_Skp2',
           s * (V_dp21p * omega * (1 / (K_dp21skp2 * omega + obsSkp2)) * (1 / (K_dp21p * omega + obsp21p))))
Rule('p21p_degradation_Skp2',  p21(b=None, S='p') + Skp2() >> Skp2(), kd_p21p_Skp2)  # p21p degradation from skp2

# ############
#    cdh1
# ############
Parameter('V_1cdh1', 1.25)  # 1.25 # µM/h
Parameter('K_1cdh1', 0.1)  # µM
Expression('kautoact_cdh1', s * (V_1cdh1 * omega / (K_1cdh1 * omega + obscdh1inact)))
Rule('cdh1_autoactivation', cdh1(A='no') >> cdh1(A='yes'), kautoact_cdh1)

Parameter('V_2cdh1', 8)  # 8 # /h
Parameter('K_2cdh1', 0.1)  # µM
Expression('k_inact_CycAorB', s * (V_2cdh1 * (1 / (K_2cdh1 * omega + obscdh1act))))
Rule('cdh1_inhibiton_CycA', cdh1(A='yes') + CycA(b=None) >> cdh1(A='no') + CycA(b=None), k_inact_CycAorB)
Rule('cdh1_inhibiton_CycB', cdh1(A='yes') + CycB(b=None) >> cdh1(A='no') + CycB(b=None), k_inact_CycAorB)

# ############
#     Myc
# ############
Parameter('k_M', 0.27)  # updated 20230808 # µM/h
Parameter('K_S', 0.8)  # 0.5 # 0.1 # µM
Expression('ks_Myc_AKT', s * ((k_M * omega) * (AKT * omega) / ((K_S * omega) + (AKT * omega))))
Rule('Myc_synthesis_AKT', None >> Myc(), ks_Myc_AKT)

Parameter('k_EMY', 0.37)  # updated 20230808 # µM/h
Parameter('K_MY', 0.1)  # 0.15 # 0.1 # µM
Expression('ks_Myc_ERK', s * ((k_EMY * omega) * (ERK * omega) / ((K_MY * omega) + (ERK * omega))))
Rule('Myc_synthesis_ERK', None >> Myc(), ks_Myc_ERK)

Parameter('d_M', 0.7)  # /h
Expression('kd_Myc', s * d_M)
Rule('Myc_degradation', Myc() >> None, kd_Myc)

# ############
#    CycD
# ############
Parameter('k_CDS', 0.25)  # updated 20230808 # µM/h
Expression('ks_CycD_AKT', s * ((k_CDS * omega) * (AKT * omega) / ((K_S * omega) + (AKT * omega))))
Rule('CycD_synthesis_AKT', None >> CycD(b=None), ks_CycD_AKT)

Parameter('k_ED', 0.2)  # 0.23 # 0.08 # µM/h
Expression('ks_CycD_ERK', s * ((k_ED * omega) * (ERK * omega) / ((K_MY * omega) + (ERK * omega))))
Rule('CycD_synthesis_ERK', None >> CycD(b=None), ks_CycD_ERK)

Parameter('k_CD', 0.03)  # 0.03 # µM/h
Parameter('K_M', 0.2)  # updated 20240711 to match parameters.csv # µM # NOTE that this is distinct from k_M in Myc above
Expression('ks_CycD_Myc', s * ((k_CD * omega) / ((K_M * omega) + obsMyc)))
Rule('CycD_synthesis_Myc', Myc() >> CycD(b=None) + Myc(), ks_CycD_Myc)

Parameter('d_CD', 1.5)  # 1.5 # /h
Expression('kd_CycD', s * d_CD)
Rule('CycD_degradation', CycD(b=None) >> None, kd_CycD)

Parameter('V_Di', 0.6)
Parameter('K_Di', 0.2)
Expression('kd_CycD_i', s * (V_Di * (1 / (K_Di * omega + obsCycD))))
Rule('CycD_inhibition', CycD(b=None) + CycDi() >> CycDi(), kd_CycD_i)

# ############
# pRb and pRbp
# ############
Parameter('k_R', 0.07)  # 0.03 # 0.18 # µM/h
Expression('ks_pRb', s * (k_R * omega))
Rule('pRb_synthesis', None >> pRb(b=None, S='u'), ks_pRb)

Parameter('k_DP', 3.6)  # 3.6 # µM/h
Parameter('K_RP', 0.01)  # 0.01 # µM
Expression('kunph_pRbp', s * ((k_DP * omega) / ((K_RP * omega) + obspRbp)))
Rule('pRbp_unphosphorylation', pRb(b=None, S='p') >> pRb(b=None, S='u'), kunph_pRbp)

Parameter('k_RE', 180)  # 180 # 1/µM-h
Expression('kf_RbE2F', s * (k_RE / omega))
Rule('pRbE2F_complex_formation', pRb(b=None, S='u') + E2F(b=None) >> pRb(b=1, S='u') % E2F(b=1), kf_RbE2F)

Parameter('d_RE', 0.03)  # 0.03 # 1/h
Expression('kd_pRbE2F', s * d_RE)
Rule('pRbE2F_degradation', pRb(b=1, S='u') % E2F(b=1) >> None, kd_pRbE2F)

Parameter('k_P2', 18)  # 18 # 1/h
Parameter('K_CE', 0.55)  # 0.92 # µM
Expression('kp_pRb_CycE', s * (k_P2 / ((K_CE * omega) + obspRb)))
Expression('kp_pRbE2F_CycE', s * (k_P2 / ((K_CE * omega) + obspRbE2F)))
Rule('pRb_phosphorylation_CycE', pRb(b=None, S='u') + CycE(b=None) >> pRb(b=None, S='p') + CycE(b=None), kp_pRb_CycE)
Rule('pRbE2F_phosphorylation_CycE', pRb(b=1, S='u') % E2F(b=1) + CycE(b=None) >>
     pRb(b=None, S='p') + E2F(b=None) + CycE(b=None), kp_pRbE2F_CycE)

Parameter('k_P1', 18)  # 18 # 1/h
Parameter('K_CD', 0.55)  # 0.92 # µM
Expression('kp_pRb_CycD', s * (k_P1 / ((K_CD * omega) + obspRb)))
Expression('kp_pRbE2F_CycD', s * (k_P1 / ((K_CD * omega) + obspRbE2F)))
Rule('pRb_phosphorylation_CycD', pRb(b=None, S='u') + CycD(b=None) >> pRb(b=None, S='p') + CycD(b=None), kp_pRb_CycD)
Rule('pRbE2F_phosphorylation_CycD', pRb(b=1, S='u') % E2F(b=1) + CycD(b=None) >>
     pRb(b=None, S='p') + E2F(b=None) + CycD(b=None), kp_pRbE2F_CycD)

# TODO: Do the two rules below need different rate constants, or can we combine them into one rule?
Parameter('d_R', 0.06)  # 0.03 # 1/h
Expression('kd_pRb', s * d_R)
Rule('pRb_degradation', pRb(b=None, S='u') >> None, kd_pRb)

Parameter('d_RP', 0.06) # 1/h
Expression('kd_pRbp', s * d_RP)
Rule('pRbp_degradation', pRb(b=None, S='p') >> None, kd_pRbp)

# ############
#     E2F
# ############
# E2F synthesis from Myc E2F autoactivation
Parameter('k_E', 0.4)  # 0.4 # µM/h
Parameter('K_E', 0.15)  # 0.1 # µM
Expression('ks_E2F_E2FMyc', s * (k_E * omega * (1 / ((K_M * omega) + obsMyc)) * (1 / ((K_E * omega) + obsE2F))))
Rule('E2F_synthesis_E2FMyc', Myc() + E2F(b=None) >> Myc() + E2F(b=None) + E2F(b=None), ks_E2F_E2FMyc)

# E2F synthesis from Myc E2F autoactivation
Parameter('k_b', 0.003)  # 0.003 # µM/h
Expression('ks_E2F_Myc', s * (k_b * omega / ((K_M * omega) + obsMyc)))
Rule('E2F_synthesis_Myc', Myc() >> Myc() + E2F(b=None), ks_E2F_Myc)

Parameter('d_E', 0.25)  # 1/h
Expression('kd_E2F', s * d_E)
Rule('E2F_degradation', E2F(b=None) >> None, kd_E2F)

Parameter('V_2E2F', 0.7)  # 0.7  # Rate constant for inactivation of E2F by cyclinA/Cdk2 complex (1/h)
Parameter('K_2E2F', 0.5)   # 0.5 # Michaelis constant for E2f inactivation by cyclinA/Cdk2 complex (µM)
Expression('kd_E2F_CycA', s * (V_2E2F / ((K_2E2F * omega) + obsE2F)))
Rule('E2F_degradation_CycA', E2F(b=None) + CycA(b=None) >> CycA(b=None), kd_E2F_CycA)

Parameter('V_E2Fi', 0.6)
Parameter('K_E2Fi', 0.2)
Expression('kd_E2F_i', s * (V_E2Fi * (1 / (K_E2Fi * omega + obsE2Fi))))
Rule('E2F_inhibition', E2F(b=None) + E2Fi() >> E2Fi(), kd_E2F_i)

# ############
#    CycE
# ############
Parameter('V_1Me', 1.1)  # 0.7 # 1.1 # Rate for activation of cyclinE/Cdk2 by E2F (1/µM-h)
Parameter('a', 1)  # basal term for activation of cyclin E/Cdk2 (without PF loop)) (µM)
# Term representing the self-activation of cyclinE/Cdk2 through the mutual activation between cyclinE/Cdk2 and its
# phosphatase Cdc25 (PF loop in the G1/S transition)) without=0 with=1
Parameter('b1', 1)  # 1/µM
Parameter('K_1Me', 0.2)  # 0.1 # Michaelis constant for cyclinE/Cdk2 activation by E2F (µM)
Expression('kact_CycE_E2F', s * ((V_1Me / omega) * ((a * omega) + b1 * obsCycE) *
                                 ((CycE_tot * omega) - obsCycE - obsCycEp21) /
                                 ((K_1Me * omega) + (CycE_tot * omega) - obsCycE - obsCycEp21)))
Rule('CycE_activation_E2F', E2F(b=None) >> E2F(b=None) + CycE(b=None), kact_CycE_E2F)  # 1/h

Parameter('V_2Me', 1.6)  # 0.5 # 1.6 # Rate constant for inactivation fo clyclinE/Cdk2 by cyclinA/Cdk2 (1/h)
Parameter('K_2Me', 0.12)  # 0.1 # Michaelis constant for cyclinE/Cdk2 inactivation by cyclinA/Cdk2 (µM)
Expression('kinact_CycE_CycA', s * (V_2Me * (1 / (K_2Me * omega + obsCycE))))
Rule('CycE_inactivation_CycA', CycA(b=None) + CycE(b=None) >> CycA(b=None), kinact_CycE_CycA)

Parameter('V_Ei', 0.6)
Parameter('K_Ei', 0.2)
Expression('kd_CycE_i', s * (V_Ei * (1 / (K_Ei * omega + obsCycE))))
Rule('CycE_inhibition', CycE(b=None) + CycEi() >> CycEi(), kd_CycE_i)

# ############
#    CycA
# ############

Parameter('V_1Ma', 0.6)  # 0.6 # Rate constant for activation of cyclinA/Cdk2 by E2F (1/h)
Parameter('K_1Ma', 0.2)  # 0.5 # Michaelis constant for cyclinA/Cdk2 activation by E2F (µM)
Expression('kact_CycA_E2F', s * (V_1Ma * ((CycA_tot * omega) - obsCycA - obsCycAp21) /
                                 ((K_1Ma * omega) + ((CycA_tot * omega) - obsCycA - obsCycAp21))))
Rule('CycA_activation_E2F', E2F(b=None) >> E2F(b=None) + CycA(b=None), kact_CycA_E2F)

Parameter('V_2Ma', 0.6)  # 0.6 # Rate constant for inactivation of cyclinA/Cdk2 by Cdc20 (1/h)
Parameter('K_2Ma', 0.1)  # 0.1 # Michaelis constant for cyclinA/Cdk2 inactivation by Cdc20 (µM)
Expression('kinact_CycA_Cdc20', s * (V_2Ma * (1 / (K_2Ma * omega + obsCycA))))
Rule('CycA_inactivation_Cdc20', CycA(b=None) + Cdc20() >> Cdc20(), kinact_CycA_Cdc20)

Parameter('V_Ai', 0.6)
Parameter('K_Ai', 0.2)
Expression('kd_CycA_i', s * (V_Ai * (1 / (K_Ai * omega + obsCycA))))
Rule('CycA_inhibition', CycA(b=None) + CycAi() >> CycAi(), kd_CycA_i)

# ############
#    CycB
# ############

Parameter('V_1Mb', 1)  # 1 # Rate for activation of cyclinB/Cdk1 by cyclinA/Cdk2 (1/µM-h)
# Term representing the self-activation of cyclinB/Cdk1 through mutual activation between cyclin B/cdk1 and its
# phosphatase cdk25 (PF loop in the G2/M transition)) without=0 with=1
Parameter('b2', 1)  # 1/µM
Parameter('K_ie', 0.2)  # µM
Parameter('K_1Mb', 0.3)  # 0.1 # Michaelis constant for cyclinB/Cdk1 activation by cyclinA/Cdk2 (µM)
Expression('kact_CycB_CycA', s * ((V_1Mb / omega) * ((a * omega) + b2 * obsCycB) *
                                  ((K_ie * omega) / ((K_ie * omega) + obsCycE)) *
                                  ((CycB_tot * omega) - obsCycBp21 - obsCycB) /
                                  ((K_1Mb * omega) + (CycB_tot * omega) - obsCycBp21 - obsCycB)))
Rule('CycB_activation_CycA', CycA(b=None) >> CycA(b=None) + CycB(b=None), kact_CycB_CycA)  # 1/h

Parameter('V_2Mb', 0.8)  # 0.8 # Rate constant for inactivation of cyclinB/Cdk1 by Cdc20 (1/h)
Parameter('K_ib', 0.5)  # 0.5 # (µM)
Parameter('K_2Mb', 0.12)  # 0.1 # Michaelis constant for cyclinB/Cdk1 inactivation by Cdc20 (µM)
Expression('kinact_CycB_Cdc20', s * (V_2Mb * ((K_ib * omega) / ((K_ib * omega) + obsCycB)) *
                                     (1 / ((K_2Mb * omega) + obsCycB))))
Rule('CycB_inactivation_Cdc20', CycB(b=None) + Cdc20() >> Cdc20(), kinact_CycB_Cdc20)  # 1/µM-h

Parameter('V_Bi', 0.6)
Parameter('K_Bi', 0.2)
Expression('kd_CycB_i', s * (V_Bi * (1 / (K_Bi * omega + obsCycB))))
Rule('CycB_inhibition', CycB(b=None) + CycBi() >> CycBi(), kd_CycB_i)

# ############
#    Cdc20
# ############

Parameter('V_1cdc20', 1.8)  # 1.8  # Rate constant for activation of Cdc20 through phosphorylation by cyclinB/Cdk1 (1/h)
Parameter('K_1cdc20', 0.1)  # µM
Expression('kact_Cdc20_CycB', s * (V_1cdc20 * (((Cdc20_tot * omega) - obsCdc20) /
                                               ((K_1cdc20 * omega) + ((Cdc20_tot * omega) - obsCdc20)))))
Rule('Cdc20_activation_CycB', CycB(b=None) >> CycB(b=None) + Cdc20(), kact_Cdc20_CycB)

Parameter('V_2cdc20', 0.6)  # 0.6  # Rate for inactivation of Cdc20 through dephosphorylation (µM/h)
Parameter('K_2cdc20', 0.1)  # 0.1 # Michaelis constant for Cdc20 inactivation through dephosphorylation (µM)
Expression('kinact_Cdc20', s * ((V_2cdc20 * omega) * (1 / ((K_2cdc20 * omega) + obsCdc20))))
Rule('Cdc20_inactivation', Cdc20() >> None, kinact_Cdc20)

Parameter('V_cdc20i', 0.6)
Parameter('K_cdc20i', 0.2)
Expression('kd_Cdc20_i', s * (V_cdc20i * (1 / (K_cdc20i * omega + obsCdc20i))))
Rule('Cdc20_inhibition', Cdc20() + Cdc20i() >> Cdc20i(), kd_Cdc20_i)

# ################
# Apoptotic switch
# ################

# BAX activation from BH3 and BAXm inactivation
Parameter('kf1', 0.01)  # 1/h
Parameter('kf2', 2.4)  # 1/µM-h
Parameter('kb', 0.14)  # 1/h
Expression('kact_BAXm', m * (kf1 + (kf2 / omega) * obsBH3))
Expression('kinact_BAXm', m * kb)
Rule('BAX_activation_inactivation', BAX(m='no', b=None) | BAX(m='yes', b=None), kact_BAXm, kinact_BAXm)
Rule('BAXmBCL2_inactivation', BAX(m='yes', b=1) % BCL2(b=1) >> BAX(m='no', b=None) + BCL2(b=None), kinact_BAXm)

# BAXm and BCL2 association and dissociation
Parameter('kasbaxmbcl2', 75)  # 1/µM-h
Parameter('kdsbaxmbcl2', 0.05)  # 1/h
Expression('kasso_BAXm_BCL2', m * (kasbaxmbcl2 / omega))
Expression('kdiss_BAXmBCL2', m * kdsbaxmbcl2)
Rule('BAXmBCL2_association_dissociation', BAX(m='yes', b=None) + BCL2(b=None) | BAX(m='yes', b=1) % BCL2(b=1),
     kasso_BAXm_BCL2, kdiss_BAXmBCL2)

Parameter('k1sbh3', 0.001)  # µM/h
Parameter('k2sbh3', 0.025)  # 0.0095 # 0.012 # 0.025 # 1/h
Parameter('a_stress', 0.6515)  # µM updated 20230926
Parameter('b_stress_AKT', 0.112)  # unitless
Parameter('b_stress_ERK', 0.0111)  # unitless
Expression('ks_BH3', m * ((k1sbh3 * omega) +
                          k2sbh3 * ((a_stress * omega) - b_stress_AKT * (AKT * omega) - b_stress_ERK * (ERK * omega))))
Rule('BH3_synthesis', None >> BH3(b=None), ks_BH3)  # µM/h

Parameter('kdbh3', 0.03)  # 1/h
Expression('kd_BH3', m * kdbh3)
Rule('BH3_degradation', BH3(b=None) >> None, kd_BH3)
Rule('BH3BCL2_degradation', BH3(b=1) % BCL2(b=1) >> BCL2(b=None), kd_BH3)

Parameter('kasbh3bcl2', 25)  # 1/µM-h
Parameter('kdsbh3bcl2', 0.01)  # 1/h
Expression('kasso_BH3_BCL2', m * (kasbh3bcl2 / omega))
Expression('kdiss_BH3BCL2', m * kdsbh3bcl2)
Rule('BH3mBCL2_association_dissociation', BH3(b=None) + BCL2(b=None) | BH3(b=1) % BCL2(b=1),
     kasso_BH3_BCL2, kdiss_BH3BCL2)

if __name__ == '__main__':
    # Plot trajectory of CDC20
    def plot_results(sim_result, label):
        # check if more than one sim_result
        try:
            iter(sim_result)
        except TypeError:
            sim_result = [sim_result]
            label = [label]
        # plot results
        figure()
        for i in range(len(sim_result)):
            plot(time, sim_result[i].observables['obsCdc20'], label=label[i])
        xlabel('Time (h)')
        ylabel('Amount of CDC20 (# of molecules)')
        legend(loc=0)

    # Simulate the model for 240 hours
    time = linspace(0, 240, 2401)
    print("Simulating...")
    sim_result_ODE = ScipyOdeSimulator(model, time, verbose=True).run()
    sim_result_SSA = BngSimulator(model, time, verbose=True).run(seed=33)
    plot_results([sim_result_ODE, sim_result_SSA], ['ODE', 'SSA'])

    # Import most up-to-date model and compare time courses
    from Austin_proliferation_model_6 import model as model_austin
    param_values = {}
    for p in model.parameters:
        if p.name not in ['AKT', 'ERK', 'FAK']:
            param_values[p.name] = p.value
            # print('%s, %g, %g' % (p.name, p.value, model_austin.parameters[p.name].value))
        else:
            param_values[p.name] = p.value * omega.value
    # quit()
    # Use same parameter values and initial species amounts to compare models
    sim_result_Austin = ScipyOdeSimulator(model_austin, time, verbose=True).\
        run(param_values=param_values, initials=sim_result_ODE.species[0])
    plot_results([sim_result_ODE, sim_result_Austin], ['ODE', 'Austin'])

    show()
