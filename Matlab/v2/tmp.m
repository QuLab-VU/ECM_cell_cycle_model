figure; hold on
for ii = 1:3
    [cells, cell_count] = runStochCell_death(0.6,0.2,1,5,0.01,0.5)
    plot(cell_count(:,1), log(cell_count(:,2)./cell_count(1,2))/log(2),Color='red')
    [cells, cell_count] = runStochCell_death(1,1,1,3,0.01,0.5) %hi dmso
    plot(cell_count(:,1), log(cell_count(:,2)./cell_count(1,2))/log(2),Color='blue')
end

figure; hold on
for ii = 1:3
    [cells, cell_count] = runStochCell_death(0.6,0.2,0.6,5,0.01,0.5)
    plot(cell_count(:,1), log(cell_count(:,2)./cell_count(1,2))/log(2),Color='red')
    [cells, cell_count] = runStochCell_death(1,1,0.6,5,0.01,0.5) %int dmso
    plot(cell_count(:,1), log(cell_count(:,2)./cell_count(1,2))/log(2),Color='blue')
end

figure; hold on
for ii = 1:3
    [cells, cell_count] = runStochCell_death(0.3,0.2,0.4,10,0.01,0.5) %low erl
    plot(cell_count(:,1), log(cell_count(:,2)./cell_count(1,2))/log(2),Color='red')
    [cells, cell_count] = runStochCell_death(0.6,1,0.4,10,0.01,0.5) %low dmso
    plot(cell_count(:,1), log(cell_count(:,2)./cell_count(1,2))/log(2),Color='blue')
end