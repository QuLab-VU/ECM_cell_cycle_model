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

p.AKT = 0.3;
p.ERK = 0.2;
p.FAK = 0.4;
dt = 0.01;

% p.AKT = 0.01;
% p.ERK = 0.01;
% p.FAK = 0.01;

p.Stress=0.75;
figure; hold on;
title('Low Input')

for ii = 1:10
    rng(ii)
    
    [t,x] = cell_stoch_fakp21_test(sys_size,p,Initial_s);
    plot(t,x(:,23))
        
end

p.AKT = 0.3;
p.ERK = 0.2;
p.FAK = 0.4;

[t,y] = ode23(@(t,y) cell_ode(t,y,p), TIMES, Initial); 
figure;hold on
plot(t,p.BAXt-y(:,19)-y(:,20),'DisplayName','BAXm')
plot(t,y(:,19),'DisplayName','BAX')
plot(t,y(:,20),'DisplayName','BAXmBCL2')
plot(t,y(:,21),'DisplayName','BH3')
plot(t,y(:,22),'DisplayName','BH3BCL2')
plot(t,p.BCL2t-y(:,22)-y(:,20),'DisplayName','BCL2')
legend()

p.kf2 = 2.4;
p.k2sbh3=0.025;
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
p.AKT = 0.3;
p.ERK = 0.2;
p.FAK = 0.4;

[t,x] = cell_stoch(sys_size,p,Initial_s);
figure;hold on
plot(t,x(:,23),'DisplayName','BAXm')
plot(t,x(:,19),'DisplayName','BAX')
plot(t,x(:,20),'DisplayName','BAXmBCL2')
plot(t,x(:,21),'DisplayName','BH3')
plot(t,x(:,22),'DisplayName','BH3BCL2')
plot(t,x(:,24),'DisplayName','BCL2')
legend()

[cells, cell_count, cell_table] = runStochCell(0.3,0.2,0.4,15,0.01,false)



[cells, cell_count, cell_table] = runStochCell_death(0.3,0.2,0.4,10,0.01,0.56,true)
[cells, cell_count, cell_table] = runStochCell_death(0.3,0.2,0.4,10,0.01,0.6,true)
[cells, cell_count, cell_table] = runStochCell_death(0.3,0.2,0.4,10,0.01,0.7,true)

[cells, cell_count, cell_table] = runStochCell_death(0.6,0.2,0.6,10,0.01,0.56,true)
[cells, cell_count, cell_table] = runStochCell_death(0.6,0.2,0.6,10,0.01,0.6,true)
[cells, cell_count, cell_table] = runStochCell_death(0.6,0.2,0.6,10,0.01,0.7,true)


stress = 0.55:0.01:0.66;
plot_c = parula(length(stress));
death_rate_lo_dmso = zeros(length(stress),3);
figure; hold on
set(gcf,'color','w')
for ii = 1:11
    for jj = 1:3
        [cells, cell_count] = runStochCell_death(0.6,1,0.4,10,0.01,stress(ii)) %low dmso
        plot(cell_count(:,1), log(cell_count(:,2)./cell_count(1,2)),Color=plot_c(ii,:))
        death_rate_lo_dmso(ii,jj) = calc_death(cells);
        ii
        jj
    end
end
axis([0 120 0 5])


death_rate_lo_erl = zeros(length(stress),3);
figure; hold on
set(gcf,'color','w')
axis([0 120 0 5])
for ii = 1:11
    for jj = 1:3
        [cells, cell_count] = runStochCell_death(0.3,0.2,0.4,10,0.01,stress(ii)) %low erl
        plot(cell_count(:,1), log(cell_count(:,2)./cell_count(1,2)),Color=plot_c(ii,:))
        death_rate_lo_erl(ii,jj) = calc_death(cells);
        ii
        jj
    end
end




death_rate_int_dmso = zeros(length(stress),3);
figure; hold on
set(gcf,'color','w')
axis([0 120 0 5])
for ii = 1:11
    for jj = 1:3
        [cells, cell_count] = runStochCell_death(1,1,0.6,5,0.01,stress(ii)) %int dmso
        plot(cell_count(:,1), log(cell_count(:,2)./cell_count(1,2)),Color=plot_c(ii,:))
        death_rate_int_dmso(ii,jj) = calc_death(cells);
        ii
        jj
    end
end


death_rate_int_erl = zeros(length(stress),3);
figure; hold on
set(gcf,'color','w')
axis([0 120 0 5])
for ii = 1:11
    for jj = 1:3
        [cells, cell_count] = runStochCell_death(0.6,0.2,0.6,5,0.01,stress(ii)) %int erl
        plot(cell_count(:,1), log(cell_count(:,2)./cell_count(1,2)),Color=plot_c(ii,:))
        death_rate_int_erl(ii,jj) = calc_death(cells);
        ii
        jj
    end
end

death_rate_hi_dmso = zeros(length(stress),3);
figure; hold on
set(gcf,'color','w')
axis([0 120 0 5])
for ii = 1:11
    for jj = 1:3
        [cells, cell_count] = runStochCell_death(1,1,1,5,0.01,stress(ii)) %hi dmso
        plot(cell_count(:,1), log(cell_count(:,2)./cell_count(1,2)),Color=plot_c(ii,:))
        death_rate_hi_dmso(ii,jj) = calc_death(cells);
        ii
        jj
    end
end

death_rate_hi_erl = zeros(length(stress),3);
figure; hold on
set(gcf,'color','w')
axis([0 120 0 5])
for ii = 1:11
    for jj = 1:3
        [cells, cell_count] = runStochCell_death(0.6,0.2,1,5,0.01,stress(ii)) %hi erl
        plot(cell_count(:,1), log(cell_count(:,2)./cell_count(1,2)),Color=plot_c(ii,:))
        death_rate_hi_erl(ii,jj) = calc_death(cells);
        ii
        jj
    end
end

figure; hold on
set(gcf,'color','w')
for ii = 1:11
    scatter(stress(ii).*[1,1,1],death_rate_lo_dmso(ii,:),10,'blue')
    scatter(stress(ii).*[1,1,1],death_rate_int_dmso(ii,:),10,'green')
    scatter(stress(ii).*[1,1,1],death_rate_hi_dmso(ii,:),10,'cyan')
    scatter(stress(ii).*[1,1,1],death_rate_lo_erl(ii,:),10,'red','filled')
    scatter(stress(ii).*[1,1,1],death_rate_int_erl(ii,:),10,'magenta','filled')
    scatter(stress(ii).*[1,1,1],death_rate_hi_erl(ii,:),10,'yellow','filled')
end
legend('Low, DMSO', 'Int, DMSO', 'High, DMSO', 'Low, ERL', 'Int, ERL', 'High, ERL')


figure; hold on
set(gcf,'color','w')
for ii = 1:11
    scatter(stress(ii).*[1,1,1],death_rate_lo_dmso(ii,:),10,'blue')
    scatter(stress(ii).*[1,1,1],death_rate_int_dmso(ii,:),10,'green')
    scatter(stress(ii).*[1,1,1],death_rate_hi_dmso(ii,:),10,'cyan')
    scatter(stress(ii).*[1,1,1],death_rate_lo_erl(ii,:),10,'red','filled')
    scatter(stress(ii).*[1,1,1],death_rate_int_erl(ii,:),10,'magenta','filled')
    scatter(stress(ii).*[1,1,1],death_rate_hi_erl(ii,:),10,'yellow','filled')
end


ttest_lo = zeros(11,2);
for ii = 1:11
    [ttest_lo(ii,1),ttest_lo(ii,2)] = ttest2(death_rate_lo_dmso(ii,:),death_rate_lo_erl(ii,:));
end

ttest_int = zeros(11,2);
for ii = 1:11
    [ttest_int(ii,1),ttest_int(ii,2)] = ttest2(death_rate_int_dmso(ii,:),death_rate_int_erl(ii,:));
end

ttest_hi = zeros(11,2);
for ii = 1:11
    [ttest_hi(ii,1),ttest_hi(ii,2)] = ttest2(death_rate_hi_dmso(ii,:),death_rate_hi_erl(ii,:));
end

stress_fit=[1-0.64;1-0.6;1-0.54;1-0.65;1-0.63;1-0.55];
preds = [0.6, 1, 0.4;
    1, 1, 0.6;
    1, 1, 1;
    0.3, 0.2, 0.4;
    0.6, 0.2, 0.6;
    0.6, 0.2, 0.55];
b = regress(stress_fit, preds);
b_1 = regress(stress_fit, cat(2,x,z));

x = [0.6; 1; 1; 0.3; 0.6; 0.6];
y = [1; 1; 1; 0.2; 0.2; 0.2];
z = [0.4; 0.6; 1; 0.4; 0.6; 1];

out = [0.64, 0.6, 0.54, 0.65, 0.63, 0.55];

p=polyfitn([x,y,z],stress_fit,'constant x, y, z');

polyvaln(p,[x(1),y(1),z(1)]);

1-polyvaln(p,[1,1]);

stress_out=zeros(6,1);
for ii=1:6
    stress_out(ii) = 1-polyvaln(p,[x(ii),y(ii),z(ii)]);
end

stress_fit=[1-0.565;1-0.56;1-0.56;1-0.59;1-0.57;1-0.57];
out = [0.6, 0.56, 0.56, 0.65, 0.61, 0.61];
x = [0.6; 1; 1; 0.3; 0.6; 0.6];
y = [1; 1; 1; 0.2; 0.2; 0.2];
p=polyfitn([x,y],out,'constant x, y');
coeffs = [0.6815, -0.112, -0.0111];


cell_count_all = {};
for ii = 1:6
    figure; hold on
    set(gcf,'color','w')
    axis([0 120 0 5])
    for jj = 1:5
        [cells, cell_count] = runStochCell_death(x(ii),y(ii),z(ii),10,0.01,out(ii),false) %hi erl
        plot(cell_count(:,1), log(cell_count(:,2)./cell_count(1,2)),Color=plot_c(ii,:))
        cell_count_all{ii,jj}=cell_count;
        ii
        jj
    end
end

figure; hold on
set(gcf,'color','w')
axis([0 120 0 5])
for ii = 1:5
    plot(cell_count_all{3,ii}(:,1), log(cell_count_all{3,ii}(:,2)./cell_count_all{3,ii}(1,2)),Color=plot_c(3,:))
end

stress_out=zeros(6,1);
for ii=1:6
    stress_out(ii) = 1-polyvaln(p,[x(ii),y(ii),z(ii)]);
end

figure;plot(cell_count(:,1), log(cell_count(:,2)./cell_count(1,2)))

sys_size=200;
p.run_time=120;
p.AKT = 0.3;
p.ERK = 0.2;
p.FAK = 0.4;
Stress = [0.55:0.01:0.65];
all_dead = [];
num = 10;
t_dead_all = [];

for ii = 1:length(Stress)
    p.Stress = Stress(ii);
    dead_cells = 0;
    t_dead = [];
    figure; hold on;
    for n = 1:num
        [t,x] = cell_stoch_death_test(sys_size,p,Initial_s);
        dead = find((x(:,23)>25));
        if length(dead) > 0
            dead_cells = dead_cells + 1;
            t_dead = [t_dead, t(dead(1))];
        end
        plot(t,x(:,23),"Color",[.7 .7 .7 .5])
    end
    t_dead_all = [t_dead_all, mean(t_dead)];
    all_dead(ii) = dead_cells/num;
    ii
    dead_cells
end

figure;scatter(Stress,all_dead)
figure;scatter(Stress,t_dead_all)
   
