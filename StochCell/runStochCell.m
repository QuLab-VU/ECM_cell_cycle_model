function [ret] = runStochCell(AKT,ERK,FAK,init)
clear cells
init_cell=init;
v = 0.05;     
sys_size=200;
run_time=120;
p.AKT=1;p.ERK=1;

CI=0;
p.FAK=FAK;

% run ode model to get initial starting protein levels ie. every cell is
% going to be at a different point of the cell cycle

p.omega=sys_size;
p.m = 3;
p.BCL2t = 0.95;
p.BAXt = 1;
p.Stress = 2-p.AKT-p.ERK;
p.kf1 = 0.01;
p.kf2 = 2.4;
p.kb = 0.14;
p.kasbaxmbcl2 = 75;
p.kdsbaxmbcl2 = 0.05;
p.kasbh3bcl2=25;
p.kdsbh3bcl2=0.01;
p.k1sbh3=0.001;
p.k2sbh3=0.025;  %0.009%0.025
p.kdbh3=0.03;

%%%%%%%%%%%
p.K_i13 = 0.8;  %0.1
p.K_i14 = 5;  %2
p.K_dcescp2 = 2;
p.V_sac = 0.4;  %0.4
p.V_da = 0.35;
p.V_sec = 0.7;  %0.8
p.V_de = 0.75;
p.V_sf = 0.4;  %0.4
p.V_fs = 1.6;
p.V_dfa = 0.35;
p.V_sskp2 = 0.01; %0.15
p.V_s2skp = 0.01;%5.20103250e-02; %0.01; %0.07
p.V_dskp2 = 1.1;
p.V_1p21 = 100;
p.V_2p21 = 0.1;
p.V_dp21p = 5;
p.K_dbcdc20 = 0.3;    %%%%%
p.K_dbcdh1 = 0.3;     %%%%%%
p.K_db = 0.08; %0.07  %0.0005  %%%%%% Michaelis constant for cycB deg
p.V_1cdh1 = 1.25; %1.25
p.V_2cdh1 = 8;  %8
p.V_db = 0.8;  %0.08  %0.15 %0.06   %%%%%new value for cyc B deg
p.v_scdh1a = 0.11;
p.v_2fak = 0.1; %0.05
p.v_s1p21 = 0.5;  %3%0.8
p.v_s2p21 = 0.1;  %2%0.1
p.K_1cdh1 = 0.1;
p.K_2cdh1 = 0.1;
p.K_dp21p = 0.1;
p.K_dp21skp2 = 0.1;
p.K_cdh1 = 0.4;
p.K_dskp2 = 0.5;
p.K_AFAK = 0.1;
p.K_dfa = 0.1;
p.K_fs = 0.1;
p.K_cf = 0.1;
p.K_cac = 0.1;
p.K_cec = 0.1;
p.K_da = 0.1;
p.K_de = 0.1;
p.K_1p21 = 0.5;
p.K_2p21 = 0.5;
p.k_ddp21 = 0.06;
p.k_ddskp2 = 0.005;
p.k_dcdh1a = 0.1;
p.k_c1 = 3;  %0.15
p.k_c2 = 1;  %0.05
p.k_c3 = 8;  %0.2
p.k_c4 = 2;  %0.1
p.k_c5 = 5;  %0.15
p.k_c6 = 1; %0.125
p.k_c7 = 0.12;
p.k_c8 = 0.2;
p.V_Fp = 1.5;%4.42837548e-02;%1.5;
p.K_Fp = 0.1;
p.k_CDF = 0.21;%1.71177745e-04;%0.21;
p.K_DF = 0.8;%4.27401052e-01;%0.8;
%%%%%%%%%%%%%%%%%%%%%

%p.omega = sys_size;
p.a = 0.8;                        % basal term for activation of cyclin E/Cdk2 and cyclin B/Cdk1 (without PF loop); uM
p.s = 0.37; %0.4;     %0.43
% scaling parameter
p.Cdc20_tot = 1;        % total concentration of cdc20; uM
p.FAK_tot = 1;
p.skp2_tot = 1;
p.p21_tot = 1;
p.cdh1_tot = 1;
%e2f_tot = 1;            % total concentration of e2f; uM
p.GF = p.AKT;                      % growth factor; uM
p.CycE_tot = 1;              % total concentration of cyclinE/Cdk2; uM
p.CycA_tot = 1;              % total concentration of cyclinA/Cdk2; uM
p.CycB_tot = 1;              % total concentration of cyclinB/Cdk1; uM
%K_dd = 0.1;                % Michaelis constant for degradation of cyclinD/Cdk4-6 complex; uM
%K_gf = 0.1;                % Michaelis constant for synthesis of cyclinD/Cdk4-6 induced by GF; uM
p.K_1E2F = 1;   %1         % Michaelis constant for E2F activation by cyclinD/Cdk4-6 and cyclinE/Cdk 2 complexes; uM
p.K_2E2F = 0.5;    %0.5        % Michaelis constant for E2f inactivation by cyclinA/Cdk2 complex; uM
p.K_1Me = 0.1;     %0.1         % Michaelis constant for cyclinE/Cdk2 activation by E2F; uM
p.K_2Me = 0.1;     %0.1         % Michaelis constant for cyclinE/Cdk2 inactivation by cyclinA/Cdk2; uM
p.K_1Ma = 0.5;     %0.5         % Michaelis constant for cyclinA/Cdk2 activation by E2F; uM
p.K_2Ma = 0.1;     %0.1         % Michaelis constant for cyclinA/Cdk2 inactivation by Cdc20; uM
p.K_1Mb = 0.1;     %0.1         % Michaelis constant for cyclinB/Cdk1 activation by cyclinA/Cdk2; uM
p.K_2Mb = 0.1;     %0.1         % Michaelis constant for cyclinB/Cdk1 inactivation by Cdc20; uM
p.K_1cdc20 = 0.1;  %0.1      % Michaelis constant for Cdc20 activation through phosphorylation by cyclinB/Cdk1; uM
p.K_2cdc20 = 0.1;  %0.1      % MIchaelis constant for Cdc20 inactivation through dephosphorylation; uM
p.K_ie = 0.2;                % Michaelis constant for cyclinB/Cdk1 inactivation by cyclinE/Cdk2, representing the indeirect ecffecof the DNA relplication checkpoint; uM
%V_dd = 0.245;              % Maximum degradation rate of cyclinD/Cdk4-6 complex; uM/h
%V_sd = 0.175;              % Maximum rate for synthesis of cyclinD/Cdk4-6 induced by growth factor; uM/h
p.V_1cdc20 = 1.8;        % Rate constant for activation of Cdc20 through phosphorylation by cyclinB/Cdk1; /h
p.V_2cdc20 = 0.6;        % Rate for inactivation of Cdc20 through dephosphorylation; uM/h
p.V_1E2F = 0.4;            % Rate constant for activation of E2F by cyclinD/Cdk4-6 and cyclinE/Cdk2 complexe; /h
p.V_2E2F = 0.7;            % Rate constant for inactivation of E2F by cyclinA/Cdk2 complex; /h
p.V_1Me = 0.9; %1.1;        %1.1      % Rate for activation of cyclinE/Cdk2 by E2F; /uM*h
p.V_2Me = 1.6;              % Rate constant for inactivation fo clyclinE/Cdk2 by cyclinA/Cdk2; /h
p.V_1Ma = 0.6; %0.6;              % Rate constant for activation of cyclinA/Cdk2 by E2F; /h
p.V_2Ma = 0.6;              % Rate constant for inactivation of cyclinA/Cdk2 by Cdc20; /h
p.V_1Mb = 1;      %1          % Rate for activation of cyclinB/Cdk1 by cyclinA/Cdk2; /uM*h
p.V_2Mb = 0.8;              % Rate constant for inactivation of cyclinB/Cdk1 by Cdc20; /h

p.k_E = 0.4;
p.k_M = 0.13;%2.39078059e-01;%0.13;
p.k_CD = 0.03;
p.k_CDS = 0.1;%1.23105202e-04;%0.1; %0.3;
p.k_R = 0.18;
p.k_RE = 180;
p.k_b = 0.003;
p.K_S = 0.5;%1.12304753e+00;%0.5;
p.k_CE = 0.35;
p.d_M = 0.6;
p.d_E = 0.25;
p.d_CD = 1.5;
p.d_CE = 1.5;
p.d_R = 0.06;
p.d_RP = 0.06;
p.d_RE = 0.03;
p.k_P1 = 18;
p.k_P2 = 18;
p.k_DP = 3.6;
p.K_M = 0.13;
p.K_E = 0.1;
p.K_CD = 0.92;
p.K_CE = 0.92;
p.K_RP = 0.01;
p.K_MY=0.1;%1.59868274e-02;%0.1;
p.k_EMY = 0.15;%3.91191419e-02;% 0.15;  %0.2
p.k_ED = 0.17;%4.61787386e-01;%0.17;  %0.22


p.b=0.8;
% Positive feedback loops on the G1/S and G2/M transitions

p.b1 = 1;                      % term representing the self-activation of cyclinE/Cdk2 through the mutual activation between cyclinE/Cdk2 and its phosphatase Cdc25(PF loop in the G1/S transition); without=0 with=1
p.b2 = 1;                      % term representing the self-activation of cyclinB/Cdk1 through mutual activation between cyclin B/cdk1 and its phosphatase cdk25(PF loop in the G2/M transition); without=0 with=1
p.K_ib = 0.5; 


%this is the starting for a normal oscillation

Initial=[0.0758759637154584,0.0179573231225730,0.0304907956309426,0.0868238125152800...
    ,0.358070500927645,0.246007623346968,0.0114732033609201,0.209618782280972,...
    0.596265732766763,0.0663628730831172,0.659538609880404,0.0823229552545997,...
    0.00601176670298263,0.0781546778431587,0.0110421142038035,0.200858887607416,...
    0.00209763607378452,0.367494461337630,0.992676117362238,0.00730440429784349,...
    0.000428591174323906,0.199169986463653];

BAX = Initial(19);
BAXmBCL2 = Initial(20);
BH3 = Initial(21);
BH3BCL2 = Initial(22);

BAXm = p.BAXt - BAX - BAXmBCL2;
BCL2 = p.BCL2t - BH3BCL2 -BAXmBCL2;
Initial = horzcat(Initial,BAXm,BCL2);
TIMES = (0:1/3:1000);

[t,y] = ode23(@(t,y) apop_mit_UP_ODE_AKTERK_10_7(t,y,CI,p), TIMES, Initial); 


%set up initial conditions
maxima=peakfinder(y(:,18),[],0.2);   %find the peak of cdc20 - last protein to peak
initial_divided=round(y(maxima(6),:).*sys_size);
variable_range=y(maxima(6):maxima(7),:);
%disp('Finished ODE')
p.AKT=AKT;
p.ERK=ERK;

runAKTERK
filename1 = strcat('cells_10_16_',num2str(10*rem(AKT,1)),num2str(10*rem(ERK,1)));
save(filename1,'cells','cell_count');
retrieve_cellState

end



