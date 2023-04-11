function death_rate = calc_death(cells)
tot_time = 0;
dead_cells = 0;

for ii = 1:length(cells)
    tot_time = tot_time + (cells{ii,3}(end)-cells{ii,3}(1));
    if cells{ii,5} == 0
        dead_cells = dead_cells + 1;
    end
end

death_rate = dead_cells/tot_time;


