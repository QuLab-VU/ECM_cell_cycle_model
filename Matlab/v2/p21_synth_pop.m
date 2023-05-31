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

p21s = [0.5, 0.6, 0.7, 0.76];


for jj = 1:length(p21s)
    p.v_s1p21 = p21s(jj);


    figure; hold on
    set(gcf,'color','w')
    axis([0 120 0 5])
    title(strcat('High stiff', num2str(p21s(jj))))
    death_rates_hi=[];
    for ii = 1:3
        [cells, cell_count] = runStochCell_death(1,1,1,p21s(jj),5,0.01) %0.56
        plot(cell_count(:,1), log(cell_count(:,2)./cell_count(1,2))/log(2),"Color",[0 0.4470 0.7410])
        death_rates_hi(ii,1) = calc_death(cells);
        [cells, cell_count] = runStochCell_death(0.6,0.2,1,p21s(jj),5,0.01) %0.61
        plot(cell_count(:,1), log(cell_count(:,2)./cell_count(1,2))/log(2),"Color",[0.6350 0.0780 0.1840])
        death_rates_hi(ii,2) = calc_death(cells);
    end
    
    %p.v_s1p21 = 0.75;
    
    
    figure; hold on
    set(gcf,'color','w')
    axis([0 120 0 5])
    title(strcat('Int. stiff', num2str(p21s(jj))))
    death_rates_int=[];
    for ii = 1:3
        [cells, cell_count] = runStochCell_death(1,1,0.6,p21s(jj),5,0.01) %0.56
        plot(cell_count(:,1), log(cell_count(:,2)./cell_count(1,2))/log(2),"Color",[0 0.4470 0.7410])
        death_rates_int(ii,1) = calc_death(cells);
        [cells, cell_count] = runStochCell_death(0.6,0.2,0.6,p21s(jj),5,0.01) %0.61
        plot(cell_count(:,1), log(cell_count(:,2)./cell_count(1,2))/log(2),"Color",[0.6350 0.0780 0.1840])
        death_rates_int(ii,2) = calc_death(cells);
    end
    
    
    %p.v_s1p21 = 0.8;
    
    figure; hold on
    set(gcf,'color','w')
    axis([0 120 0 5])
    title(strcat('Low stiff', num2str(p21s(jj))))
    death_rates_low=[];
    for ii = 1:3
        [cells, cell_count] = runStochCell_death(0.6,1,0.4,p21s(jj),5,0.01) %0.6
        plot(cell_count(:,1), log(cell_count(:,2)./cell_count(1,2))/log(2),"Color",[0 0.4470 0.7410])
        death_rates_low(ii,1) = calc_death(cells);
        [cells, cell_count] = runStochCell_death(0.3,0.2,0.4,p21s(jj),10,0.01) %0.65
        plot(cell_count(:,1), log(cell_count(:,2)./cell_count(1,2))/log(2),"Color",[0.6350 0.0780 0.1840])
        death_rates_low(ii,2) = calc_death(cells);
    end
end