AKT_DMSO = [0.6, 1, 1];
ERK_DMSO = [1, 1, 1];
FAK_DMSO = [0.4, 0.6, 1];

AKT_ERL = [0.3, 0.6, 0.6];
ERK_ERL = [0.2, 0.2, 0.2];
FAK_ERL = [0.4, 0.6, 1];

init = 5;  %50 for ERL + low stiffness
dt = 0.01;
save_t = false;

SAVEDIR = strcat(matlabroot, "\ECM_UMAP_output");

v_names = {'id', 'eoe', 'death', 'time', 'Skp2', 'p21', 'p21p', 'Cdh1', 'Myc', 'CycD', 'CycDp21',...
    'Rb', 'pRb', 'E2F', 'RbE2F', 'CycE', 'CycEp21', 'CycA', 'CycAp21', 'CycB', 'CycBp21', 'Cdc20',...
    'BAX', 'BAXmBCL2', 'BH3', 'BH3BCL2', 'BAXm', 'BCL2'};




for ii = 1:3
    divs = 0;
    count = 1;
    all_divs = {};
    AKT = AKT_DMSO(ii);
    ERK = ERK_DMSO(ii);
    FAK = FAK_DMSO(ii);
    while divs < 200
        

        [cells, cell_count, cell_table] = runStochCell(AKT,ERK,FAK,init,dt,save_t);
        %[cell_table] = restructure_for_UMAP(cells,dt,p,save_t);

        tmp_t = cell_table(cell_table.eoe == 0, :);
        tmp_t1 = tmp_t(tmp_t.death == 0, :);

        [C,~,~] = unique(tmp_t1.id);
        
        if ~isempty(C)
            for jj = 1:length(C)
                tmp_div = tmp_t1(tmp_t1.id == C(jj), :);
                tmp_cell = [];
                tmp_cell(:,1) = ones(height(tmp_div),1).*(jj+divs);
                tmp_cell(:,2:28) = table2array(tmp_div(:,2:28));
                all_divs = [all_divs; num2cell(tmp_cell)];    
            end
        end

        divs = divs + length(C);
        divs

    end
    cell_table = cell2table(all_divs, "VariableNames",v_names);
    writetable(cell_table,strcat(SAVEDIR, '\cells_divs_',num2str(AKT*100),'_',num2str(ERK*100),...
        '_',num2str(FAK*100),'.csv'));
end

for ii = 1:5
    divs = 0;
    count = 1;
    all_divs = {};
    AKT = AKT_ERL(1);
    ERK = ERK_ERL(1);
    FAK = FAK_ERL(1);
    
    [cells, cell_count, cell_table] = runStochCell(AKT,ERK,FAK,25,dt,save_t);

    
    writetable(cell_table,strcat(SAVEDIR, '\all_cells_0904_',num2str(ii),'_',num2str(AKT*100),'_',num2str(ERK*100),...
        '_',num2str(FAK*100),'.csv'));
    csvwrite(strcat(SAVEDIR, '\cells_count_0904_',num2str(ii),'_',num2str(AKT*100),'_',num2str(ERK*100),...
    '_',num2str(FAK*100),'.csv'),cell_count)
end

for ii = 1:5
    divs = 0;
    count = 1;
    all_divs = {};
    AKT = AKT_ERL(2);
    ERK = ERK_ERL(2);
    FAK = FAK_ERL(2);
    
    [cells, cell_count, cell_table] = runStochCell(AKT,ERK,FAK,5,dt,save_t);

    
    writetable(cell_table,strcat(SAVEDIR, '\all_cells_0904_',num2str(ii),'_',num2str(AKT*100),'_',num2str(ERK*100),...
        '_',num2str(FAK*100),'.csv'));
    csvwrite(strcat(SAVEDIR, '\cells_count_0904_',num2str(ii),'_',num2str(AKT*100),'_',num2str(ERK*100),...
    '_',num2str(FAK*100),'.csv'),cell_count)
end

for ii = 1:5
    divs = 0;
    count = 1;
    all_divs = {};
    AKT = AKT_ERL(3);
    ERK = ERK_ERL(3);
    FAK = FAK_ERL(3);
    
    [cells, cell_count, cell_table] = runStochCell(AKT,ERK,FAK,5,dt,save_t);

    
    writetable(cell_table,strcat(SAVEDIR, '\all_cells_0904_',num2str(ii),'_',num2str(AKT*100),'_',num2str(ERK*100),...
        '_',num2str(FAK*100),'.csv'));
    csvwrite(strcat(SAVEDIR, '\cells_count_0904_',num2str(ii),'_',num2str(AKT*100),'_',num2str(ERK*100),...
    '_',num2str(FAK*100),'.csv'),cell_count)
end



%%
p21_vals = [0.05, 0.25, 0.5, 1, 2];
V_2Me_vals = [0.16, 0.8, 1.6, 2, 2.5];

AKT = AKT_ERL(2);
ERK = ERK_ERL(2);
FAK = FAK_ERL(2);
sys_size = 200;

for ii = 1:5
    [p] = set_params(sys_size,AKT,ERK);
    p.FAK=FAK;
    p.AKT=AKT;
    p.ERK=ERK;
    p.v_s1p21 = p21_vals(ii);

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
    TIMES = (0:1/3:4000);
    
    [t,y] = ode23(@(t,y) cell_ode(t,y,p), TIMES, Initial); 
    maxima=peakfinder(y(:,18),[],0.2);

    if  ~isempty(maxima)
        initial_divided=round(y(maxima(5),:).*sys_size);
        variable_range=y(maxima(4):maxima(5),:);

        initial_s = round(datasample(variable_range,1,1).*sys_size);
        p.run_time = 5000;

        [t,x] = cell_stoch(sys_size,p,initial_divided);
        [new_t, new_x] = rest_stoch(t,x,dt);
        new_x_n = new_x./200;
        
        csvwrite(strcat(SAVEDIR, '\prot_p21_',num2str(ii),'csv'), new_x_n);
    end
end


V_2Me_vals = [0.8, 1, 1.6, 2, 2.5];
for ii = 1:5
    [p] = set_params(sys_size,AKT,ERK);
    p.FAK=FAK;
    p.AKT=AKT;
    p.ERK=ERK;
    p.V_2Me = V_2Me_vals(ii);

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
    TIMES = (0:1/3:4000);
    
    [t,y] = ode23(@(t,y) cell_ode(t,y,p), TIMES, Initial); 
    maxima=peakfinder(y(:,18),[],0.2);

    if  ~isempty(maxima)
        initial_divided=round(y(maxima(5),:).*sys_size);
        variable_range=y(maxima(4):maxima(5),:);

        initial_s = round(datasample(variable_range,1,1).*sys_size);
        p.run_time = 5000;

        [t,x] = cell_stoch(sys_size,p,initial_divided);
        [new_t, new_x] = rest_stoch(t,x,dt);
        new_x_n = new_x./200;
        
        csvwrite(strcat(SAVEDIR, '\prot_V2Me_',num2str(ii),'csv'), new_x_n);
    end
end
    
for ii = 1:3
    [p] = set_params(sys_size,1,1);
    p.FAK=FAK_ERL(ii);
    p.AKT=AKT_ERL(ii);
    p.ERK=ERK_ERL(ii);

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
    TIMES = (0:0.01:1000);
    
    [t,y] = ode23(@(t,y) cell_ode(t,y,p), TIMES, Initial);

    csvwrite(strcat(SAVEDIR, '\ode_cond_',num2str(ii),'_',num2str(p.AKT*100),'_',num2str(p.ERK*100),...
    '_',num2str(p.FAK*100),'.csv'),y)
end

for ii = 1:3
    [p] = set_params(sys_size,1,1);
    p.FAK=FAK_DMSO(ii);
    p.AKT=AKT_DMSO(ii);
    p.ERK=ERK_DMSO(ii);

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
    TIMES = (0:0.01:1000);
    
    [t,y] = ode23(@(t,y) cell_ode(t,y,p), TIMES, Initial);

    csvwrite(strcat(SAVEDIR, '\ode_cond_',num2str(ii),'_',num2str(p.AKT*100),'_',num2str(p.ERK*100),...
    '_',num2str(p.FAK*100),'.csv'),y)
end