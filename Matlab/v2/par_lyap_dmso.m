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
TIMES = (0:0.001:2000);

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

%dmso
FAKs = linspace(0.4,0.6,50);
tmp = linspace(0.6,1,50);
FAKs = [FAKs,tmp];
AKTs = linspace(0.6,1,50);
tmp = ones(1,50).*1;
AKTs = [AKTs, tmp];
ERKs = [1];
lyap = [];
ep = 0.001;
ep_all = [-1e-2 -1e-3 1e-3 1e-2];

dt = 0.001;

c=1;
d={};
d1={};

M=10;

parfor (ii = 1:5, M)
    AKT=AKTs(ii);
    ERK=ERKs(1);
    FAK=FAKs(ii);

    
    for kk = 1:6

        [t,y] = ode23(@(t,y) cell_ode_death_test(t,y,p,AKT,ERK,FAK), TIMES, Initial);
        max_cdc=peakfinder(y(:,18),[],0.2);

        new_y = y(max_cdc(end-11):max_cdc(end-1),:);
        new_t = t(max_cdc(end-11):max_cdc(end-1),:);
        
        variable_range=y(max_cdc(end-11):max_cdc(end-10),:);        
        tmp_loc = 1+(round(length(variable_range)/kk)*(kk-1));

        Initial_1 = new_y(tmp_loc,:)+ep;
        new_y_ab = new_y(tmp_loc:end,:);

        [t,x]=ode23(@(t,y) (cell_ode_death_test(t,y,p,AKT,ERK,FAK)),[0:dt:T],Initial_1,options);

        z=[];
        z1=[];
        n=1;

        for i=1:2000:length(x)
            dists = [];
            %L1_norm = [];
            for j = 1:length(new_y_ab)
                dists(j) = norm(new_y_ab(j,:)-x(i,:));
                %L1_norm(j) = norm(new_y_ab(j,:)-x(i,:),1);
            end
            z(n) = min(dists);
            %z1(n) = min(L1_norm);
            n=n+1;
            % n
            % kk
            % ii
        end

        d{ii,kk} = log(z);
        %d1{ii,kk} = log(z1);
    end
end

save lyap_erl_160923.mat d
