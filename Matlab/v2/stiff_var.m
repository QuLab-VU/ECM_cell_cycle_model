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
TIMES = (0:0.01:2000);

p.AKT = 1;
p.ERK = 1;
p.FAK = 1;

[t,y] = ode23(@(t,y) cell_ode(t,y,p), TIMES, Initial);

sys_size=200;
maxima=peakfinder(y(:,18),[],0.2);   %find the peak of cdc20 - last protein to peak
initial_divided=round(y(maxima(6),:).*sys_size);
variable_range=y(maxima(6):maxima(7),:);

Initial_s = round(datasample(variable_range,1,1).*sys_size);
p.run_time=120;

FAKs = linspace(0.4,0.6,100);
tmp = linspace(0.602,1,100);
FAKs = [FAKs,tmp];
AKTs = linspace(0.3,0.6,100);
tmp = ones(1,100).*0.6;
AKTs = [AKTs, tmp];
ERKs = [0.2];
s1p21s = 0.3:0.02:0.8;
IMTs = [];

for ii = 1:length(s1p21s)
    p.v_s1p21 = s1p21s(ii);
    Initial=[0.0199387776457489,0.0212324326339528,0.104907517295920,0.970596980815461,...
        0.321067821067829,0.175953608509162,0.0138939800293163,0.301757475166483,...
        0.545620559818762,0.0403080965235495,0.554487916752975,0.0640517424993558,...
        0.00487528808175496,0.0182005140710040,0.00234453838504852,0.0234757871244827,...
        0.00435874152378834,0.0204771947463057,0.999000159839911,0.000997752220497874,...
        5.81398691338436e-05,0.0332751934642001];
    IMTs_tmp = [];
    for jj = 1:length(FAKs)
    
        p.AKT=AKTs(jj);
        p.ERK=ERKs(1);
        p.FAK=FAKs(jj);
        [t,y] = ode23(@(t,y) cell_ode_death_test(t,y,p), TIMES, Initial);
        maxima=peakfinder(y(:,18),[],0.2);
        Initial = y(end,:);
        IMTs_tmp(jj)=t(maxima(10))-t(maxima(9));
        ii
        jj
    end
    IMTs(ii,:)=IMTs_tmp;
end

figure;hold on;
for ii = 1:26
    plot(IMTs(ii,:))
end

%%
CDS_vals = 0.11:0.02:0.25;
IMTs = [];

for ii = 1:length(CDS_vals)
    p.k_CDS = CDS_vals(ii);
    Initial=[0.0199387776457489,0.0212324326339528,0.104907517295920,0.970596980815461,...
        0.321067821067829,0.175953608509162,0.0138939800293163,0.301757475166483,...
        0.545620559818762,0.0403080965235495,0.554487916752975,0.0640517424993558,...
        0.00487528808175496,0.0182005140710040,0.00234453838504852,0.0234757871244827,...
        0.00435874152378834,0.0204771947463057,0.999000159839911,0.000997752220497874,...
        5.81398691338436e-05,0.0332751934642001];
    IMTs_tmp = [];
    for jj = 1:length(FAKs)
    
        p.AKT=AKTs(jj);
        p.ERK=ERKs(1);
        p.FAK=FAKs(jj);
        [t,y] = ode23(@(t,y) cell_ode_death_test(t,y,p), TIMES, Initial);
        maxima=peakfinder(y(:,18),[],0.2);
        Initial = y(end,:);
        IMTs_tmp(jj)=t(maxima(7))-t(maxima(6));
        ii
        jj
    end
    IMTs(ii,:)=IMTs_tmp;
end

figure;hold on;
for ii = 1:8
    plot(IMTs(ii,:))
end
%%

FAKs = linspace(0.4,0.6,100);
tmp = linspace(0.6,1,100);
FAKs = [FAKs,tmp];
AKTs = linspace(0.6,1,100);
tmp = ones(1,100).*1;
AKTs = [AKTs, tmp];
ERKs = [1];

for jj = 1:length(FAKs)

    p.AKT=AKTs(jj);
    p.ERK=ERKs(1);
    p.FAK=FAKs(jj);
    [t,y] = ode23(@(t,y) cell_ode_death_test(t,y,p), TIMES, Initial);
    Initial = y(end,:);
    dlmwrite(strcat('dmso_range_4-29_', sprintf('%04d',jj),'.txt'),y)
   
    jj
end