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
TIMES = (0:0.01:1000);

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

p.AKT = 0.6;
p.ERK = 0.2;
p.FAK = 0.4;
dt = 0.01;

FAKs = 0:0.1:1;

figure;hold on
for ii = 1:length(FAKs)
    p.FAK = FAKs(ii);
    [t,y] = ode23(@(t,y) cell_ode_death_test(t,y,p), TIMES, Initial);
    plot(t,y(:,12))
end

[p] = set_params(200,1,1);
p.AKT = 1;
p.ERK = 1;
p.FAK = 0.6;
figure;hold on
[t,y] = ode23(@(t,y) cell_ode_death_test(t,y,p), TIMES, Initial);
plot(t,y(:,12))
xlim([500,1000])
ylim([0, 1])
p.AKT = 0.6;
p.ERK = 0.2;
p.FAK = 0.6;
[t,y] = ode23(@(t,y) cell_ode_death_test(t,y,p), TIMES, Initial);
plot(t,y(:,12))

[p] = set_params(200,1,1);
p.AKT = 0.6;
p.ERK = 0.2;
p.FAK = 0.6;
figure;hold on
[t,y] = ode23(@(t,y) cell_ode_death_test(t,y,p), TIMES, Initial);
plot(t,y(:,12))
xlim([500,1000])
ylim([0, 1])
p.V_sskp2 = 0;
p.FAK=0;
p.AKT = 0.6;
p.ERK = 0.2;
[t,y] = ode23(@(t,y) cell_ode_death_test(t,y,p), TIMES, Initial);
plot(t,y(:,12))

[p] = set_params(200,1,1);
p.AKT = 0.6;
p.ERK = 0.2;
p.FAK = 0.6;
figure;hold on
[t,y] = ode23(@(t,y) cell_ode_death_test(t,y,p), TIMES, Initial);
plot(t,y(:,12))
xlim([500,1000])
ylim([0, 1])
p.V_dp21p = 0;
[t,y] = ode23(@(t,y) cell_ode_death_test(t,y,p), TIMES, Initial);
plot(t,y(:,12))

[p] = set_params(200,1,1);
p.AKT = 0.3;
p.ERK = 0.2;
p.FAK = 0.4;
figure;hold on
[t,y] = ode23(@(t,y) cell_ode_death_test(t,y,p), TIMES, Initial);
plot(t,y(:,12))
xlim([500,1000])
ylim([0, 1])
p.V_1p21 = 0.9;
[t,y] = ode23(@(t,y) cell_ode_death_test(t,y,p), TIMES, Initial);
plot(t,y(:,12))



[p] = set_params(200,1,1);
p.AKT = 0.6;
p.ERK = 0.2;
p.FAK = 0.6;
figure;hold on
[t,y] = ode23(@(t,y) cell_ode_death_test(t,y,p), TIMES, Initial);
plot(t,y(:,12))
%plot(t,y(:,2))
xlim([500,1000])
ylim([0, 1])
p.v_s1p21 = 1.6;
%figure;hold on
[t,y] = ode23(@(t,y) cell_ode_death_test(t,y,p), TIMES, Initial);
plot(t,y(:,12))
%plot(t,y(:,2))
xlim([500,1000])
ylim([0, 1])



v_s1p21s = 1:0.1:2;
[p] = set_params(200,1,1);
p.AKT = 0.6;
p.ERK = 0.2;
p.FAK = 0.6;
figure;hold on
[t,y] = ode23(@(t,y) cell_ode_fakp21_test(t,y,p), TIMES, Initial);
plot(t,y(:,12),"LineWidth",1.5,"Color",'k')
%plot(t,y(:,2))
xlim([500,1000])
ylim([0, 1])
for ii = 1:length(v_s1p21s)
    p.v_s1p21 = v_s1p21s(ii);
    %figure;hold on
    [t,y] = ode23(@(t,y) cell_ode_fakp21_test(t,y,p), TIMES, Initial);
    plot(t,y(:,12))
    %plot(t,y(:,2))
end



[p] = set_params(200,1,1);
p.AKT = 0.6;
p.ERK = 0.2;
p.FAK = 0.6;
figure;hold on
[t,y] = ode23(@(t,y) cell_ode_death_test(t,y,p), TIMES, Initial);
plot(t,y(:,12))
plot(t,y(:,2))


p.v_s1p21 = 2.5;
figure;hold on
[t,y] = ode23(@(t,y) cell_ode_death_test(t,y,p), TIMES, Initial);
plot(t,y(:,12))
plot(t,y(:,2))


[p] = set_params(200,1,1);
p.AKT = 0.3;
p.ERK = 0.2;
p.FAK = 0.4;
figure;hold on
[t,y] = ode23(@(t,y) cell_ode_p21_test(t,y,p), TIMES, Initial);
plot(t,y(:,12))
plot(t,y(:,2))

p.v_s1p21 = 5;
figure;hold on
[t,y] = ode23(@(t,y) cell_ode_p21_test(t,y,p), TIMES, Initial);
plot(t,y(:,12))
plot(t,y(:,2))

p.v_s1p21=0.8;
p.F_n=8;
p.K_p21fak=5;
FAKs = 0:0.01:1;
p21r = [];
p21r2=[];
ns = 1:10;

figure; hold on
for jj = 1:length(ns)
    
    for ii=1:length(FAKs)
        FAK = FAKs(ii);
        p21r(ii)=(p.v_s1p21/(1+p.K_p21fak*FAK^p.F_n));
        %p21r2(ii )=(p.v_s1p21-((p.v_s1p21*FAK^p.F_n)/(p.K_p21fak^p.F_n+FAK^p.F_n)));
    end
    plot(FAKs,p21r);
end
figure; hold on
plot(FAKs,p21r);
plot(FAKs,p21r2)

%p.K_p21fak = 0.1;
%p.v_s1p21 = 0.5; %1.5
K_p21faks = [0:0.2:1];
[p] = set_params(200,1,1);
p.AKT = 0.6;
p.ERK = 0.2;
p.FAK = 0.6;
figure;hold on
[t,y] = ode23(@(t,y) cell_ode_fakp21_test(t,y,p), TIMES, Initial);
plot(t,y(:,12),"LineWidth",1.5,"Color",'k')
%plot(t,y(:,2))
xlim([500,1000])
ylim([0, 1])
for ii = 1:length(K_p21faks)
    p.K_p21fak = K_p21faks(ii);
    %figure;hold on
    [t,y] = ode23(@(t,y) cell_ode_fakp21_test(t,y,p), TIMES, Initial);
    plot(t,y(:,12))
    %plot(t,y(:,2))
end


p.v_s1p21 = 0.5;
p.F_n=2;
p.AKT = 0.3;
p.ERK = 0.2;
p.FAK = 0.4;
figure;hold on
[t,y] = ode23(@(t,y) cell_ode_fakp21_test(t,y,p), TIMES, Initial);
plot(t,y(:,12))
%plot(t,y(:,2))
xlim([500,1000])
ylim([0, 1])
p.AKT = 1;
p.ERK = 1;
p.FAK = 0.6;
[t,y] = ode23(@(t,y) cell_ode_fakp21_test(t,y,p), TIMES, Initial);
plot(t,y(:,12))

p.v_s1p21 = 1;
p.AKT = 0.6;
p.ERK = 0.2;
p.FAK = 0.6;
FAKs = 0:0.2:1;
figure;hold on
for ii = 1:length(FAKs)
    p.FAK = FAKs(ii);
    [t,y] = ode23(@(t,y) cell_ode_fakp21_test(t,y,p), TIMES, Initial);
    plot(t,y(:,12))
end

%%
AKT_DMSO = [0.6, 1, 1];
ERK_DMSO = [1, 1, 1];
FAK_DMSO = [0.4, 0.6, 1];

AKT_ERL = [0.3, 0.6, 0.6];
ERK_ERL = [0.2, 0.2, 0.2];
FAK_ERL = [0.4, 0.6, 1];

AKTs = [0.6, 1, 1, 0.3, 0.6, 0.6];
ERKs = [1, 1, 1, 0.2, 0.2, 0.2];
stress = [0.645, 0.62, 0.62, 0.67, 0.65, 0.65];

DMSO_stress = [0.66, 0.62, 0.62];
ERL_stress = [0.67, 0.65, 0.65];

p.v_s1p21 = 1;
p.F_n=2;
p.K_p21fak=2;
ii=2;
figure;hold on
p.AKT = AKT_DMSO(ii);
p.ERK = ERK_DMSO(ii);
p.FAK = FAK_DMSO(ii);
[t,y] = ode23(@(t,y) cell_ode_fakp21_test(t,y,p), TIMES, Initial);
plot(t,y(:,12))
p.AKT = AKT_ERL(ii);
p.ERK = ERK_ERL(ii);
p.FAK = FAK_ERL(ii);
[t,y] = ode23(@(t,y) cell_ode_fakp21_test(t,y,p), TIMES, Initial);
plot(t,y(:,12))

Kp21s = [0.1:0.1:2];
FAKs = 0:0.01:1;
p21r = [];
Fs = 1:10;


figure;hold on
for ii =1:length(Fs)
%    p.K_p21fak=Kp21s(ii);
    p.F_n = Fs(ii);
    for jj = 1:length(FAKs)
        FAK=FAKs(jj);
        p21r(jj) = p.s*((p.v_s1p21/(1+p.K_p21fak*FAK^p.F_n)));
    end
    plot(FAKs,p21r)
end

p.K_p21fak = 5;
p.F_n = 7;
p.v_s1p21 = 0.9;
FAKs = 0:0.01:1;
p21r = [];
figure;hold on
for jj = 1:length(FAKs)
    FAK=FAKs(jj);
    p21r(jj) = ((p.v_s1p21/(1+p.K_p21fak*FAK^p.F_n)));
end
plot(FAKs,p21r)

[t,x] = cell_stoch_fakp21_test(sys_size,p,Initial_s);
plot(t,x(:,23))


[cells, cell_count] = runStochCell_p21FAK(1,1,0.6,5,0.01)
figure; hold on
plot(cell_count(:,1), log(cell_count(:,2)./cell_count(1,2)))
[cells, cell_count] = runStochCell_p21FAK(0.6,0.2,0.6,5,0.01)
plot(cell_count(:,1), log(cell_count(:,2)./cell_count(1,2)))

figure; hold on
set(gcf,'color','w')
axis([0 120 0 5])
death_rates_low=[];
for ii = 1:3
    [cells, cell_count] = runStochCell_death(0.6,1,0.4,5,0.01)
    plot(cell_count(:,1), log(cell_count(:,2)./cell_count(1,2))/log(2),"Color",[0 0.4470 0.7410])
    death_rates_low(ii,1) = calc_death(cells);
    [cells, cell_count] = runStochCell_death(0.3,0.2,0.4,5,0.01)
    plot(cell_count(:,1), log(cell_count(:,2)./cell_count(1,2))/log(2),"Color",[0.6350 0.0780 0.1840])
    death_rates_low(ii,2) = calc_death(cells);
end

figure; hold on
set(gcf,'color','w')
axis([0 120 0 5])
death_rates_int=[];
for ii = 1:3
    [cells, cell_count] = runStochCell_death(1,1,0.6,5,0.01)
    plot(cell_count(:,1), log(cell_count(:,2)./cell_count(1,2))/log(2),"Color",[0 0.4470 0.7410])
    death_rates_int(ii,1) = calc_death(cells);
    [cells, cell_count] = runStochCell_death(0.6,0.2,0.6,5,0.01)
    plot(cell_count(:,1), log(cell_count(:,2)./cell_count(1,2))/log(2),"Color",[0.6350 0.0780 0.1840])
    death_rates_int(ii,2) = calc_death(cells);
end


figure; hold on
set(gcf,'color','w')
axis([0 120 0 5])
death_rates_hi=[];
for ii = 1:3
    [cells, cell_count] = runStochCell_death(1,1,1,5,0.01)
    plot(cell_count(:,1), log(cell_count(:,2)./cell_count(1,2))/log(2),"Color",[0 0.4470 0.7410])
    death_rates_hi(ii,1) = calc_death(cells);
    [cells, cell_count] = runStochCell_death(0.6,0.2,1,5,0.01)
    plot(cell_count(:,1), log(cell_count(:,2)./cell_count(1,2))/log(2),"Color",[0.6350 0.0780 0.1840])
    death_rates_hi(ii,2) = calc_death(cells);
end


%%
p.v_s1p21 = 1;
p.F_n=7;
p.K_p21fak=10;

figure; hold on
set(gcf,'color','w')
axis([0 120 0 5])
death_rates_si=[];
for ii = 1:3
    [cells, cell_count] = runStochCell_death(0.6,0.2,0.6,5,0.01)
    plot(cell_count(:,1), log(cell_count(:,2)./cell_count(1,2))/log(2),"Color",[0 0.4470 0.7410])
    death_rates_int(ii,2) = calc_death(cells);
    [cells, cell_count] = runStochCell_death(0.6,0.2,1,5,0.01)
    plot(cell_count(:,1), log(cell_count(:,2)./cell_count(1,2))/log(2),"Color",[0.6350 0.0780 0.1840])
    death_rates_si(ii,2) = calc_death(cells);
end

%%
%p.v_s1p21 = 0.15;
%p.v_s1p21 = 0.15;
figure; hold on
set(gcf,'color','w')
axis([0 120 0 5])
death_rates_hi=[];
for ii = 1
    [cells, cell_count] = runStochCell_death(1,1,1,0.57,5,0.01) %0.56
    plot(cell_count(:,1), log(cell_count(:,2)./cell_count(1,2))/log(2),"Color",[0 0.4470 0.7410])
    death_rates_hi(ii,1) = calc_death(cells);
    [cells, cell_count] = runStochCell_death(0.6,0.2,1,0.615,5,0.01) %0.61
    plot(cell_count(:,1), log(cell_count(:,2)./cell_count(1,2))/log(2),"Color",[0.6350 0.0780 0.1840])
    death_rates_hi(ii,2) = calc_death(cells);
end

%p.v_s1p21 = 0.75;


figure; hold on
set(gcf,'color','w')
axis([0 120 0 5])
death_rates_int=[];
for ii = 1
    [cells, cell_count] = runStochCell_death(1,1,0.6,0.61,5,0.01) %0.56
    plot(cell_count(:,1), log(cell_count(:,2)./cell_count(1,2))/log(2),"Color",[0 0.4470 0.7410])
    death_rates_int(ii,1) = calc_death(cells);
    [cells, cell_count] = runStochCell_death(0.6,0.2,0.6,0.615,5,0.01) %0.61
    plot(cell_count(:,1), log(cell_count(:,2)./cell_count(1,2))/log(2),"Color",[0.6350 0.0780 0.1840])
    death_rates_int(ii,2) = calc_death(cells);
end


%p.v_s1p21 = 0.8;

figure; hold on
set(gcf,'color','w')
axis([0 120 0 5])
death_rates_low=[];
for ii = 1
    [cells, cell_count] = runStochCell_death(0.6,1,0.4,0.61,5,0.01) %0.6
    plot(cell_count(:,1), log(cell_count(:,2)./cell_count(1,2))/log(2),"Color",[0 0.4470 0.7410])
    death_rates_low(ii,1) = calc_death(cells);
    [cells, cell_count] = runStochCell_death(0.3,0.2,0.4,0.62,10,0.01) %0.65
    plot(cell_count(:,1), log(cell_count(:,2)./cell_count(1,2))/log(2),"Color",[0.6350 0.0780 0.1840])
    death_rates_low(ii,2) = calc_death(cells);
end

%%



figure; hold on
set(gcf,'color','w')
axis([0 120 0 5])

for ii = 1
    [cells, cell_count] = runStochCell_death(0.6,0.2,1,0.65,5,0.01) %0.61
    plot(cell_count(:,1), log(cell_count(:,2)./cell_count(1,2))/log(2))
    [cells, cell_count] = runStochCell_death(0.6,0.2,0.6,0.65,5,0.01) %0.61
    plot(cell_count(:,1), log(cell_count(:,2)./cell_count(1,2))/log(2))
end

%%
p.v_s1p21 = 0.8;
p.F_n=2;
p.K_p21fak=6;

FAKs = 0:0.2:1;
plot_c = parula(length(FAKs));

figure; hold on
set(gcf,'color','w')
axis([0 120 0 5])

for ii = 1:length(FAKs)
    [cells, cell_count] = runStochCell_death(0.6,0.2,FAKs(ii),3,0.01)
    plot(cell_count(:,1), log(cell_count(:,2)./cell_count(1,2)),Color=plot_c(ii,:))
end
%%



p.K_p21fak
p.F_n

figure; hold on
set(gcf,'color','w')
axis([0 120 0 5])
death_rates_hi=[];
[cells, cell_count] = runStochCell_death(1,1,0.6,3,0.01)
plot(cell_count(:,1), log(cell_count(:,2)./cell_count(1,2)),"Color",[0 0.4470 0.7410])
[cells, cell_count] = runStochCell_death(0.6,0.2,0.6,3,0.01)
plot(cell_count(:,1), log(cell_count(:,2)./cell_count(1,2)),"Color",[0.6350 0.0780 0.1840])

[cells, cell_count] = runStochCell_death(1,1,1,3,0.01)
plot(cell_count(:,1), log(cell_count(:,2)./cell_count(1,2)),"Color",[0 0.4470 0.7410])
[cells, cell_count] = runStochCell_death(0.6,0.2,1,3,0.01)
plot(cell_count(:,1), log(cell_count(:,2)./cell_count(1,2)),"Color",[0.6350 0.0780 0.1840])





figure; hold on
set(gcf,'color','w')
axis([0 120 0 5])
death_rates_low=[];
for ii = 1:3
    [cells, cell_count] = runStochCell_p21FAK(0.6,1,0.4,5,0.01)
    plot(cell_count(:,1), log(cell_count(:,2)./cell_count(1,2)),"Color",[0 0.4470 0.7410])
    death_rates_low(ii,1) = calc_death(cells);
    [cells, cell_count] = runStochCell_p21FAK(0.3,0.2,0.4,5,0.01)
    plot(cell_count(:,1), log(cell_count(:,2)./cell_count(1,2)),"Color",[0.6350 0.0780 0.1840])
    death_rates_low(ii,2) = calc_death(cells);
end


%%
