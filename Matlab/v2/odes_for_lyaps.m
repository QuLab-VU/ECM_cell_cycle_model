[p] = set_params(200,1,1);

Initial=[0.0199387776457489,0.0212324326339528,0.104907517295920,0.970596980815461,...
    0.321067821067829,0.175953608509162,0.0138939800293163,0.301757475166483,...
    0.545620559818762,0.0403080965235495,0.554487916752975,0.0640517424993558,...
    0.00487528808175496,0.0182005140710040,0.00234453838504852,0.0234757871244827,...
    0.00435874152378834,0.0204771947463057,0.999000159839911,0.000997752220497874,...
    5.81398691338436e-05,0.0332751934642001];

BAX = Initial(19);
BAXmBCL2 = Initial(20);
BH3 = Initial(21);
BH3BCL2 = Initial(22);

BAXm = p.BAXt - BAX - BAXmBCL2;
BCL2 = p.BCL2t - BH3BCL2 -BAXmBCL2;
Initial = horzcat(Initial,BAXm,BCL2);
TIMES = (0:0.1:1000);

AKT = 1;
ERK = 1;
FAK = 1;
[t,y] = ode23(@(t,y) cell_ode_death_test(t,y,p,AKT,ERK,FAK), TIMES, Initial); 
maxima=peakfinder(y(:,18),[],0.2); 
variable_range=y(maxima(6):maxima(7),:);
r = randi([1 length(variable_range)],1,1);
for_init = [1 50 100 150 200 250];

Initial = variable_range(1,:);

options=odeset('RelTol',1.0e-10);
T=1000;

%erl
FAKs = linspace(0.4,0.6,500);
tmp = linspace(0.6,1,500);
FAKs = [FAKs,tmp];
AKTs = linspace(0.3,0.6,500);
tmp = ones(1,500).*0.6;
AKTs = [AKTs, tmp];
ERKs = [0.2];

M=11;
win_len = 40;
erl_IMT = zeros(1000,1);

parfor (ii = 1:1000, M)
    AKT=AKTs(ii);
    ERK=ERKs(1);
    FAK=FAKs(ii);

    [t,y] = ode23(@(t,y) cell_ode_death_test(t,y,p,AKT,ERK,FAK), TIMES, Initial);
    max_cdc=peakfinder(y(:,18),[],0.2);

    variable_range=y(max_cdc(9)-(win_len/2):max_cdc(11)+(win_len/2),1:18);
    variable_range=variable_range.*200;
    
    erl_IMT(ii) = t(max_cdc(11)-max_cdc(9))/2;

    dlmwrite(strcat('C:\Users\olesk\ECM_mod_output\lyap_odes\erl_range_9-23_', sprintf('%04d',ii),'.csv'),variable_range)
end
dlmwrite(strcat('C:\Users\olesk\ECM_mod_output\lyap_odes\erl_IMT.csv'),erl_IMT)

%dmso
FAKs = linspace(0.4,0.6,500);
tmp = linspace(0.6,1,500);
FAKs = [FAKs,tmp];
AKTs = linspace(0.6,1,500);
tmp = ones(1,500).*1;
AKTs = [AKTs, tmp];
ERKs = [1];
dmso_IMT = zeros(1000,1);

parfor (ii = 1:1000, M)
    AKT=AKTs(ii);
    ERK=ERKs(1);
    FAK=FAKs(ii);

    [t,y] = ode23(@(t,y) cell_ode_death_test(t,y,p,AKT,ERK,FAK), TIMES, Initial);
    max_cdc=peakfinder(y(:,18),[],0.2);

    variable_range=y(max_cdc(9)-(win_len/2):max_cdc(11)+(win_len/2),1:18);
    variable_range=variable_range.*200;

    dmso_IMT(ii) = t(max_cdc(11)-max_cdc(9))/2;

    dlmwrite(strcat('C:\Users\olesk\ECM_mod_output\lyap_odes\dmso_range_9-23_', sprintf('%04d',ii),'.csv'),variable_range)
end
dlmwrite(strcat('C:\Users\olesk\ECM_mod_output\lyap_odes\dmso_IMT.csv'),dmso_IMT)

