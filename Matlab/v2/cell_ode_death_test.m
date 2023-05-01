function [dydt]=cell_ode_death_test(t,y,p)
 
dydt = zeros(size(y));
FAK=p.FAK;
AKT=p.AKT;
ERK=p.ERK;
Skp2=y(1);
p21=y(2);
p21p=y(3);  
Cdh1=y(4);  
Myc=y(5);
CycD=y(6);  
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
BH3 = y(21);
BH3BCL2 = y(22);
BAXm=p.BAXt-BAX-BAXmBCL2;
BCL2=p.BCL2t-BH3BCL2-BAXmBCL2;
Stress=0.6815-0.112*p.AKT-0.0111*p.ERK;
%Stress=1-0.8*AKT-0.2*ERK;
 
%Skp2 
dydt(1) = p.s*((p.V_sskp2)... %Basal skp2 synthesis
    + (p.V_s2skp*FAK)... %Skp2 synthesis from FAK
    - (p.V_dskp2*(Skp2/(p.K_dskp2+Skp2))*(Cdh1/(p.K_cdh1+Cdh1)))... %skp2 degradation from cdh1
    - (p.k_ddskp2*Skp2)); %skp2 degradation
 
%p21
dydt(2) = p.s*((p.v_s1p21/(1+p.K_p21fak*FAK^p.F_n))... %p21 synthesis
    -(p.k_c1*CycD*p21)... %p21 cycD association
    +(p.k_c2*CycDp21)... %p21 cycD dis
    -(p.k_c3*CycE*p21)... %p21 cycE association
    +(p.k_c4*CycEp21)... %p21 cycE dis
    -(p.k_c5*CycA*p21)... %p21 cycA association
    +(p.k_c6*CycAp21)... %p21 cycA dis
    -(p.k_c7*CycB*p21)... %p21 cycB association
    +(p.k_c8*CycBp21)... %p21 cycB dis
    -(p.V_1p21*CycE*p21/(p.K_1p21+p21))...%p21 inact from cycE
    +(p.V_2p21*p21p/(p.K_2p21+p21p))... %p21 autoact from p21p 
    -(p.k_ddp21*p21)); %p21 non specific inact  %these prob shouldn't be the same
 
%p21p
dydt(3) = p.s*((p.V_1p21*CycE*p21/(p.K_1p21+p21))... %p21 inact from cycE
    +(p.V_Fp*p21*FAK/(p.K_Fp+p21))... %p21 autoact from p21p 
    -(p.V_2p21*p21p/(p.K_2p21+p21p))... %p21 autoact from p21p 
    -(p.V_dp21p*Skp2/(p.K_dp21p+Skp2)*p21p/(p.K_dp21skp2+p21p))... %p21p deg from skp2
    -(p.k_ddp21*p21p)); %p21p non specific inact  %these prob shouldn't be the same
 
%Cdh1
dydt(4) = p.s*((p.V_1cdh1*((p.cdh1_tot-Cdh1)/(p.K_1cdh1+(p.cdh1_tot-Cdh1))))... %cdh1 autoact
    -(p.V_2cdh1*(Cdh1/(p.K_2cdh1+Cdh1))*(CycA+CycB)));%cdh1 inact from cycA and cycB
 
%Myc
dydt(5) = p.s*((p.k_M*AKT)/((p.K_S+AKT))...%myc synthesis from AKT
    +(p.k_EMY*ERK)/(p.K_MY+ERK) ...%myc synthesis from ERK
    -p.d_M*Myc);%myc degradation
 
%CycD 
dydt(6) = p.s*(p.k_CDS*AKT/(p.K_S+AKT)...%cyclin D synthesis from AKT
    -(p.k_c1*CycD*p21)...%p21 cycD association
    +(p.k_c2*CycDp21)...%p21 cycD dis
    +p.k_ED*ERK/(p.K_MY+ERK)...%cyclin D synthesis from ERK
    +(p.k_CD*Myc)/(p.K_M+Myc)...%cycD synth from myc
    -(p.d_CD*CycD));%cycD degradation
dydt(7) = p.s*((p.k_c1*CycD*p21)...%p21 cycD association
    -(p.k_c2*CycDp21));%p21 cycD dis
 
%Rb
dydt(8) = p.s*(((p.k_R)...%Rb synthesis 
    +((p.k_DP*pRb)/((p.K_RP)+pRb)))...%pRb dephosphorylation
    -((p.k_RE)*Rb*E2F)...%Rb-E2F complex formation
    -((p.k_P1*CycD*Rb)/((p.K_CD)+Rb))...%Rb phospho from cycD
    -((p.k_P2*CycE*Rb)/((p.K_CE)+Rb))...%Rb phosphorylation from cyc E
    -(p.d_R*Rb));%Rb degradation
 
%pRb
dydt(9) = p.s*(((p.k_P1*CycD*Rb)/((p.K_CD)+Rb))...%Rb phospho from cycD
    +((p.k_P2*CycE*Rb)/((p.K_CE)+Rb))...%Rb phosphorylation from cyc E
    +((p.k_P1*CycD*RbE2F)/((p.K_CD)+RbE2F))...%rb e2f phospho from cycD
    +((p.k_P2*CycE*RbE2F)/((p.K_CE)+RbE2F))...%Rb-e2f phosphorylation from cyc E
    -((p.k_DP*pRb)/((p.K_RP)+pRb))...%pRb dephosphorylation
    -(p.d_RP*pRb));%pRb degradation
 
%E2F
dydt(10) = p.s*((p.k_E*(Myc/((p.K_M)+Myc))*(E2F/((p.K_E)+E2F)))...%e2f synthesis from myc e2f autoactivation
    +((p.k_b*Myc)/((p.K_M)+Myc))...%e2f synthesis from myc
    +((p.k_P1*CycD*RbE2F)/((p.K_CD)+RbE2F))...%rb e2f phospho from cycD
    +((p.k_P2*CycE*RbE2F)/((p.K_CE)+RbE2F))...%Rb-e2f phosphorylation from cyc E
    -((p.k_RE)*Rb*E2F)...%Rb-E2F complex formation
    -(p.d_E*E2F)...%e2f degradation
    -(((p.V_2E2F*E2F*CycA)/((p.K_2E2F)+E2F))));%e2f degradation from cycA
 
%Rb-E2F
dydt(11) = p.s*(((p.k_RE)*Rb*E2F)...%Rb-E2F complex formation
    -((p.k_P1*CycD*RbE2F)/((p.K_CD)+RbE2F))...%rb e2f phospho from cycD
    -((p.k_P2*CycE*RbE2F)/((p.K_CE)+RbE2F))...%Rb-e2f phosphorylation from cyc E
    -(p.d_RE*RbE2F));%Rb-e2f degradation
 
 
%CycE
dydt(12) = p.s*(p.V_1Me *E2F*(p.a+p.b1*CycE)*(p.CycE_tot-CycE-CycEp21)/(p.K_1Me+p.CycE_tot-CycE-CycEp21)...%Cyc E activation
    -(p.k_c3*CycE*p21)...%p21 cycE association
    +(p.k_c4*CycEp21)...%p21 cycE dis
    -(p.V_2Me*CycA*(CycE/(p.K_2Me+CycE))));%Cyc E inactivation
dydt(13) = p.s*((p.k_c3*CycE*p21)...%p21 cycE association
    -(p.k_c4*CycEp21));%p21 cycE dis
 
%CycA
dydt(14) = p.s*(p.V_1Ma * E2F*(p.CycA_tot-CycAp21-CycA)/(p.K_1Ma+p.CycA_tot-CycAp21-CycA)...%Cyc A activation
    -(p.k_c5*CycA*p21)...%p21 cycA association
    +(p.k_c6*CycAp21)...%p21 cycA dis
    -(p.V_2Ma*Cdc20*(CycA/(p.K_2Ma+CycA))));%Cyc A inactivation
dydt(15) = p.s*((p.k_c5*CycA*p21)...%p21 cycA association
    -(p.k_c6*CycAp21));%p21 cycA dis
 
%CycB
dydt(16) = p.s*((p.V_1Mb *CycA * (p.b+p.b2*CycB)*(p.K_ie/(p.K_ie+CycE))*(p.CycB_tot-CycBp21-CycB)/(p.K_1Mb+(p.CycB_tot-CycBp21-CycB)))...%Cyc B activation
    -(p.k_c7*CycB*p21)...%p21 cycB association
    +(p.k_c8*CycBp21)...%p21 cycB dis
    -(p.V_2Mb*(p.K_ib/(p.K_ib+CycB))*Cdc20*(CycB/(p.K_2Mb+CycB))));%Cyc B inactivation
dydt(17) = p.s*((p.k_c7*CycB*p21)...%p21 cycB association
    -(p.k_c8*CycBp21));%p21 cycB dis

%Cdc20
dydt(18)= p.s*((p.V_1cdc20 * CycB *(p.Cdc20_tot - Cdc20)/(p.K_1cdc20+(p.Cdc20_tot-Cdc20)))...%cdc20 activation
    -(p.V_2cdc20*(Cdc20/(p.K_2cdc20+Cdc20))));%cdc20 inactivation

%BAX
dydt(19)= p.m*(-(p.kf1+p.kf2*BH3)*BAX...%bax activation from bh3
    +p.kb*BAXm...%baxm inactivation to bax
    +p.kb*BAXmBCL2);%baxmblc2 inactivation to bax


%BAXmBCL2
dydt(20)= p.m*(p.kasbaxmbcl2*BAXm*BCL2...%baxm-bcl2 association 
    -p.kdsbaxmbcl2*BAXmBCL2...%baxmbcl2 dissociation
    -p.kb*BAXmBCL2);%baxmblc2 inactivation to bax

%BH3
dydt(21)= p.m*(p.k1sbh3...%bh3 synthesis
    +p.k2sbh3*Stress...%bh3 synthesis
    -p.kdbh3*BH3...%bh3 degradation
    -p.kasbh3bcl2*BH3*BCL2...%bh3 bcl2 association
    +p.kdsbh3bcl2*BH3BCL2);%bh3bcl2 dissociation

%BH3BCL2
dydt(22)= p.m*(p.kasbh3bcl2*BH3*BCL2...%bh3 bcl2 association
    -p.kdsbh3bcl2*BH3BCL2...%bh3bcl2 dissociation
    -p.kdbh3*BH3BCL2);%bh3 degradation