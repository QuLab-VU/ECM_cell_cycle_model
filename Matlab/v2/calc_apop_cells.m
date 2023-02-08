init_cell=init;
v = 0.05;     
sys_size=200;
run_time=120;

% run ode model to get initial starting protein levels ie. every cell is
% going to be at a different point of the cell cycle

[p] = set_params(sys_size,AKT,ERK);
p.AKT=1;
p.ERK=1;
p.FAK=FAK;

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

[t,y] = ode23(@(t,y) cell_ode(t,y,p), TIMES, Initial); 


%set up initial conditions
maxima=peakfinder(y(:,18),[],0.2);   %find the peak of cdc20 - last protein to peak
initial_divided=round(y(maxima(6),:).*sys_size);
variable_range=y(maxima(6):maxima(7),:);
%disp('Finished ODE')
p.AKT=AKT;
p.ERK=ERK;
p.Stress=Stress;
p.run_time=120;

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
        [t,x] = cell_stoch_death_test(sys_size,p,round(datasample(variable_range,1,1).*sys_size));
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