function [cell_table] = restructure_for_UMAP(cells,dt,p)
SAVEDIR = strcat(matlabroot, "\ECM_UMAP_output");
if ~isdir(SAVEDIR)
    mkdir(SAVEDIR)
end


all_cells = {};
v_names = {'id', 'eoe', 'death', 'time', 'Skp2', 'p21', 'p21p', 'Cdh1', 'Myc', 'CycD', 'CycDp21',...
    'Rb', 'pRb', 'E2F', 'RbE2F', 'CycE', 'CycEp21', 'CycA', 'CycAp21', 'CycB', 'CycBp21', 'Cdc20',...
    'BAX', 'BAXmBCL2', 'BH3', 'BH3BCL2', 'BAXm', 'BCL2'};

count = 1;
for ii = 1:size(cells,1) 
    new_t = cells{ii,3}(1):dt:cells{ii,3}(end);
    new_t = new_t';
    new_x = zeros(length(new_t),size(cells{ii,3},2));

    approx_step = length(find(cells{ii,3}<dt))+100;
    
    nind = zeros(length(new_t),1);
    oind=1;
    for jj = 1:length(new_t)
        if oind+approx_step<length(cells{ii,3})
            tmp = find(cells{ii,3}(oind:oind+approx_step)<=new_t(jj));
        else
            tmp = find(cells{ii,3}(oind:end)<=new_t(jj));
        end

        if jj ==1 
            nind(jj) = tmp(end);
        else
            nind(jj) = tmp(end)+(oind-1);
        end
        oind=nind(jj);

    end
    new_x=cells{ii,2}(nind,:);
    
    max_tmp_cdc=peakfinder(new_x(:,18),0.2*p.omega,0.2*p.omega,'false');
    new_x_n = new_x./200;

    
    if length(max_tmp_cdc) > 1 && new_x(end,23) < 20
        for jj = 1:length(max_tmp_cdc)-1
            tmp_cell = [];
            tmp_cell(:,1) = ones(length(max_tmp_cdc(jj):max_tmp_cdc(jj+1)),1).*count;
            tmp_cell(:,2) = zeros(length(max_tmp_cdc(jj):max_tmp_cdc(jj+1)),1);
            tmp_cell(:,3) = zeros(length(max_tmp_cdc(jj):max_tmp_cdc(jj+1)),1);
            tmp_cell(:,4) = new_t(max_tmp_cdc(jj):max_tmp_cdc(jj+1));
            tmp_cell(:,5:28) = new_x_n(max_tmp_cdc(jj):max_tmp_cdc(jj+1),:);
            all_cells = [all_cells; num2cell(tmp_cell)];
            count = count+1;
        end
    elseif length(max_tmp_cdc) == 1 && new_x(end,23) < 20
        tmp_cell = [];
        tmp_cell(:,1) = ones(length(new_t(max_tmp_cdc(end):end)),1).*count;
        tmp_cell(:,2) = ones(length(new_t(max_tmp_cdc(end):end)),1);
        tmp_cell(:,3) = zeros(length(new_t(max_tmp_cdc(end):end)),1);
        tmp_cell(:,4) = new_t(max_tmp_cdc(end):end);
        tmp_cell(:,5:28) = new_x_n(max_tmp_cdc(end):end,:);
        all_cells = [all_cells; num2cell(tmp_cell)];
        count = count+1;
    elseif isempty(max_tmp_cdc) && new_x(end,23) < 20
        tmp_cell = [];
        tmp_cell(:,1) = ones(length(new_t(:)),1).*count;
        tmp_cell(:,2) = ones(length(new_t(:)),1);
        tmp_cell(:,3) = zeros(length(new_t(:)),1);
        tmp_cell(:,4) = new_t(:);
        tmp_cell(:,5:28) = new_x_n(:,:);
        all_cells = [all_cells; num2cell(tmp_cell)];
        count = count+1;

    elseif length(max_tmp_cdc) > 1 && new_x(end,23) > 20
        for jj = 1:length(max_tmp_cdc)-1
            if jj < length(max_tmp_cdc)-1
                tmp_cell = [];
                tmp_cell(:,1) = ones(length(max_tmp_cdc(jj):max_tmp_cdc(jj+1)),1).*count;
                tmp_cell(:,2) = zeros(length(max_tmp_cdc(jj):max_tmp_cdc(jj+1)),1);
                tmp_cell(:,3) = zeros(length(max_tmp_cdc(jj):max_tmp_cdc(jj+1)),1);
                tmp_cell(:,4) = new_t(max_tmp_cdc(jj):max_tmp_cdc(jj+1));
                tmp_cell(:,5:28) = new_x_n(max_tmp_cdc(jj):max_tmp_cdc(jj+1),:);
                all_cells = [all_cells; num2cell(tmp_cell)];
                count = count+1;
            else
                tmp_cell = [];
                tmp_cell(:,1) = ones(length(max_tmp_cdc(jj):max_tmp_cdc(jj+1)),1).*count;
                tmp_cell(:,2) = zeros(length(max_tmp_cdc(jj):max_tmp_cdc(jj+1)),1);
                tmp_cell(:,3) = ones(length(max_tmp_cdc(jj):max_tmp_cdc(jj+1)),1);
                tmp_cell(:,4) = new_t(max_tmp_cdc(jj):max_tmp_cdc(jj+1));
                tmp_cell(:,5:28) = new_x_n(max_tmp_cdc(jj):max_tmp_cdc(jj+1),:);
                all_cells = [all_cells; num2cell(tmp_cell)];
                count = count+1;
            end
        end

    elseif length(max_tmp_cdc) == 1 && new_x(end,23) > 20
        tmp_cell = [];
        tmp_cell(:,1) = ones(length(new_t(max_tmp_cdc(end):end)),1).*count;
        tmp_cell(:,2) = zeros(length(new_t(max_tmp_cdc(end):end)),1);
        tmp_cell(:,3) = ones(length(new_t(max_tmp_cdc(end):end)),1);
        tmp_cell(:,4) = new_t(max_tmp_cdc(end):end);
        tmp_cell(:,5:28) = new_x_n(max_tmp_cdc(end):end,:);
        all_cells = [all_cells; num2cell(tmp_cell)];
        count = count+1;
    elseif isempty(max_tmp_cdc) && new_x(end,23) > 20
        tmp_cell = [];
        tmp_cell(:,1) = ones(length(new_t(:)),1).*count;
        tmp_cell(:,2) = zeros(length(new_t(:)),1);
        tmp_cell(:,3) = ones(length(new_t(:)),1);
        tmp_cell(:,4) = new_t(:);
        tmp_cell(:,5:28) = new_x_n(:,:);
        all_cells = [all_cells; num2cell(tmp_cell)];
        count = count+1;
    end   
    
end



cell_table = cell2table(all_cells, "VariableNames",v_names);
writetable(cell_table,strcat(SAVEDIR, '\cells_',num2str(p.AKT*100),'_',num2str(p.ERK*100),...
    '_',num2str(p.FAK*100),'.csv'));
    