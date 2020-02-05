function [dydt]=apop_mit_UP_ODE_AKTERK_10_7(t,y,CI,p)
 
dydt = zeros(size(y));
FAK=p.FAK;

Yap=0;  %changed from the Goldbeter value of 0.7588.  We are not considering input from the YAP pathway
AKT=p.AKT;
ERK=p.ERK;
%FAK=0; %.4;  %0.975;
Skp2=y(1);
p21=y(2);
p21p=y(3);  
Cdh1=y(4);  %constrained by a total
Myc=y(5);
CycD=y(6);  %
CycDp21=y(7);
Rb=y(8);
pRb=y(9);
E2F=y(10);
RbE2F=y(11);
CycE=y(12);
CycEp21=y(13);
CycA=y(14);
CycAp21=y(15);
CycB=y(16);
CycBp21=y(17);
Cdc20=y(18);
BAX = y(19);
BAXmBCL2 = y(20);
BH3=y(21);
BH3BCL2 = y(22);
BAXm=p.BAXt-BAX-BAXmBCL2;
BCL2=p.BCL2t-BH3BCL2-BAXmBCL2;
Stress=2-AKT-ERK;
 
 
%RUN THIS IN THE PRESENCE AND ABSENCE OF FAK on P1  so dialing p.V_Fp from
%0 to a reasonabley mix.  Then Look at P21 LEVELS.   
 
 
 
%checked same as Austin's Version
%Skp2 has a natural synthesis and the rate is increased by FAK.   Has a
%nature degradation and one that is stimulated as part of a complex with
%Cdh1
dydt(1) = p.s*((p.V_sskp2) + (p.V_s2skp*FAK) - (p.V_dskp2*(Skp2/(p.K_dskp2+Skp2))*(Cdh1/(p.K_cdh1+Cdh1)))...
    -(p.k_ddskp2*Skp2));
 
%Austin's follow Goldbeter and includes a FAK dependent term
%p21 has natural synthesis.  It forms a reversible complex with each of the
%cyclins.  Cyclin E phosphorylates p21 and targets it for degradation Austin mimiced this term with  does FAK
dydt(2) = p.s*((p.v_s1p21)+(p.v_s2p21*E2F*(p.K_i13/(p.K_i13+Rb))*(p.K_i14/(p.K_i14+pRb)))...
    -(p.k_c1*CycD*p21)+(p.k_c2*CycDp21)-(p.k_c3*CycE*p21)+(p.k_c4*CycEp21)...
    -(p.k_c5*CycA*p21)+(p.k_c6*CycAp21)-(p.k_c7*CycB*p21)...
    +(p.k_c8*CycBp21)-(p.V_1p21*CycE*p21/(p.K_1p21+p21))...
    -(p.V_Fp*p21*FAK/(p.K_Fp+p21))+p.V_2p21*p21p/(p.K_2p21+p21p)-(p.k_ddp21*p21));
 
% phosphorylate p21 is increased due to CycE and FAK, both complexes break up.  Skp degrades p21p specifically. 
dydt(3) = p.s*((p.V_1p21*CycE*p21/(p.K_1p21+p21))+(p.V_Fp*p21*FAK/(p.K_Fp+p21))...
    -(p.V_2p21*p21p/(p.K_2p21+p21p))-(p.V_dp21p*Skp2/(p.K_dp21p+Skp2)*p21p/(p.K_dp21skp2+p21p))-(p.k_ddp21*p21p));
 
% needed for the inhibition of s2skp ....
dydt(4) = p.s * ((p.V_1cdh1*((p.cdh1_tot-Cdh1)/(p.K_1cdh1+(p.cdh1_tot-Cdh1))))...
    -(p.V_2cdh1*(Cdh1/(p.K_2cdh1+Cdh1))*(CycA+CycB)));
 
%reversible E2F Rb switch driven by growth factors ..... 
dydt(5) = p.s*((p.k_M*AKT)/((p.K_S+AKT)) + (p.k_EMY*ERK)/(p.K_MY+ERK) -p.d_M*Myc);
 
% Cyclin D is synthesized when AKT and ERK or Myc are High Are high / forms a reversable complex with
% p21 and is degraded
 
dydt(6) = p.s*(p.k_CDS*AKT/(p.K_S+AKT)+ p.k_CDF*FAK/(p.K_DF+FAK)...
    -(p.k_c1*CycD*p21)+(p.k_c2*CycDp21)...
    + p.k_ED*ERK/(p.K_MY+ERK)+ (p.k_CD*Myc)/(p.K_M+Myc)-(p.d_CD*CycD));
dydt(7) = p.s*((p.k_c1*CycD*p21)-(p.k_c2*CycDp21));
 
%RB some natural synthesis, return from phosphorylated state, -
%binds to E2F, CycD and CYCE phosphorylates RB, natural degradation 
dydt(8) = p.s*(((p.k_R)+((p.k_DP*pRb)/((p.K_RP)+pRb)))...
    -(1*(((p.k_RE)*Rb*E2F)+((p.k_P1*CycD*Rb)/((p.K_CD)+Rb))...
    + ((p.k_P2*CycE*Rb)/((p.K_CE)+Rb))+(p.d_R*Rb))));
 
%RBp inceased by CycD and CycE, some returned -> also increased through RB
%E2F complex
dydt(9) = p.s*(((p.k_P1*CycD*Rb)/((p.K_CD)+Rb)) + ((p.k_P2*CycE*Rb)/((p.K_CE)+Rb))...
    +((p.k_P1*CycD*RbE2F)/((p.K_CD)+RbE2F))+((p.k_P2*CycE*RbE2F)/((p.K_CE)+RbE2F))...
    -((p.k_DP*pRb)/((p.K_RP)+pRb))-(p.d_RP*pRb));
 
%From reversible switch paper until the last term which is the Goldbeter
%term linking E2F and cyclin A   -  Note E2F is not constrained in this
%model because it is synthesized by growth factors...   COMMENT ON THIS
%FORM ON THE LINK TO GROWTH FACTORS ....
 
%ON THE TO DO LIST IS TO EXPLORE WHAT HAPPENS WITH THIS THAT E2F is made by
%growth factors....
dydt(10) = p.s*(((p.k_E*(Myc/((p.K_M)+Myc))*(E2F/((p.K_E)+E2F)))+((p.k_b*Myc)/((p.K_M)+Myc)))+((p.k_P1*CycD*RbE2F)/((p.K_CD)+RbE2F))+((p.k_P2*CycE*RbE2F)/((p.K_CE)+RbE2F))-((p.k_RE)*Rb*E2F)-(p.d_E*E2F)-(((p.V_2E2F*E2F*CycA)/((p.K_2E2F)+E2F))));
 
%RB E2F bind complex is degraded , removed from phosphorylation by cyc D
%and cycE 
dydt(11) = p.s*(((p.k_RE)*Rb*E2F)-(1*(((p.k_P1*CycD*RbE2F)/((p.K_CD)+RbE2F))+((p.k_P2*CycE*RbE2F)/((p.K_CE)+RbE2F))+(p.d_RE*RbE2F))));
 
 
%You want this constrained by the total form - matching the Skeleton
dydt(12) = p.s*(p.V_1Me *E2F*(p.a+p.b1*CycE)*(p.CycE_tot-CycE-CycEp21)/(p.K_1Me+p.CycE_tot-CycE-CycEp21)  -...
    (p.k_c3*CycE*p21)+(p.k_c4*CycEp21)-(p.V_2Me*CycA*(CycE/(p.K_2Me+CycE)))-(p.V_de*(Skp2/(p.K_dcescp2+Skp2))*(CycE/(p.K_de+CycE))));
%Inactive Form is CycE_tot-CycE-CycEp21
dydt(13) = p.s*((p.k_c3*CycE*p21)-(p.k_c4*CycEp21));
 
%Cyclin A  - You want this constrained by the form - matching the Skeleton
%Model
dydt(14) = p.s*(p.V_1Ma * E2F*(p.CycA_tot-CycAp21-CycA)/(p.K_1Ma+p.CycA_tot-CycAp21-CycA) -(p.k_c5*CycA*p21)+(p.k_c6*CycAp21)-(p.V_2Ma*Cdc20*(CycA/(p.K_2Ma+CycA))));
dydt(15) = p.s*((p.k_c5*CycA*p21)-(p.k_c6*CycAp21));
 
%Cyclin B should be constrained by the total form - matching Skeleton
dydt(16) = p.s*(p.V_1Mb *CycA*(p.CycB_tot-CycBp21-CycB)/(p.K_1Mb+p.CycB_tot-CycBp21-CycB) -(p.k_c7*CycB*p21)...
    +(p.k_c8*CycBp21)-(p.V_db*(Cdc20)*(CycB/(p.K_db+CycB))));   %%%changed this
%(p.V_db*((Cdc20/(p.K_dbcdc20+Cdc20))+(Cdh1/(p.K_dbcdh1+Cdh1)))*(CycB/(p.K_db+CycB))) 
%(p.V_2Mb*Cdc20*(CycB/(p.K_2Mb+CycB))))  .1351
%I think the plan was to keep track of 2 cyclin B's a nuclear and a
%cytoplasmic portion, and to make their interaction slightly different -
%do this later ...
 
dydt(17) = p.s*((p.k_c7*CycB*p21)-(p.k_c8*CycBp21));
dydt(18)= p.s*((p.V_1cdc20 * CycB *(p.Cdc20_tot - Cdc20)/(p.K_1cdc20+(p.Cdc20_tot-Cdc20)))-(p.V_2cdc20*(Cdc20/(p.K_2cdc20+Cdc20))));

dydt(19)= p.m*(-(p.kf1+p.kf2*BAX)*BH3+p.kb*BAXm+p.kb*BAXmBCL2);
dydt(20)= p.m*(p.kasbaxmbcl2*BAXm*BCL2-p.kdsbaxmbcl2*BAXmBCL2-p.kb*BAXmBCL2);
dydt(21)= p.m*(p.k1sbh3+p.k2sbh3*Stress-p.kdbh3*BH3-p.kasbh3bcl2*BH3*BCL2+p.kdsbh3bcl2*BH3BCL2);
dydt(22)= p.m*(p.kasbh3bcl2*BH3*BCL2-p.kdsbh3bcl2*BH3BCL2-p.kdbh3*BH3BCL2);
 
 
%does anything change in the stability analysis if there is a constraint on
%Cyclin A total / Cylin B, E, and D.   Where is the data for the totals 
 
 
%Austin's going to check to see if the system oscillates as AKT and ERK
%increase, but does not oscillate when AKT and ERK are below a reasonable
%threshold.
 
%then the next step would be to see what happens as FAK goes from 0-1  / we
%would expect to see also have an oscillation transition
 
%Examine CycD and E2F levels to make sure they are well constrained in
%reasonable values
 
%We should look at p21 as well, to see how it adjusts as a function of FAK.