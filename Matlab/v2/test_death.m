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

figure; hold on;
title('Low Input')

for ii = 1:10
    rng(ii)
    
    [t,x] = cell_stoch(sys_size,p,initial_divided);
    plot(t,x(:,23))
        
end

p.AKT = 0.6;
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

figure;plot(cell_count(:,1), log(cell_count(:,2)./cell_count(1,2)))
