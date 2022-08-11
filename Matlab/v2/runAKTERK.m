clear cells

v = 0.05;     
sys_size=200;
run_time=120;
cells=cell(init_cell,5);
 
%%
for ii = 1:init_cell
    cells{ii,1}=p;
    cells{ii,2}=round(datasample(variable_range,1,1).*sys_size);
    cells{ii,3}=0;
    cells{ii,4}=0;
    cells{ii,5}=true;
end
 
%%
cell_count=[0 size(cells,1)];   %start the cell counter at time zero and with number of initial cells
current_time=ones(size(cells,1),2); %keeps track of the current time of each cell and whether it's dead or not
while isempty(find(current_time(:,2)==1 & current_time(:,1)<run_time,1))==false %run the code only if there are live cells under the run time
    cell_count_tmp=zeros(size(current_time,1),2); %temporary cell count with current cell time and whether a cell was added, taken away, or not
    for ii = 1:size(cells,1)
       if current_time(ii,1)<120 && cells{ii,5}==true
            if run_time-cells{ii,3}(end)<30
                p.run_time=run_time-cells{ii,3}(end);
            else
                p.run_time=30;
            end
            [tmp_t,tmp_x] = cell_stoch_07_06_21(sys_size,p,cells{ii,2}(end,:));
            cells{ii,2}=[cells{ii,2}; tmp_x(2:end,:)];
            cells{ii,3}=[cells{ii,3}; tmp_t(2:end)+current_time(ii,1)];
            current_ind_length=size(cells{ii,2},1);
 
            max_tmp_cdc=peakfinder(cells{ii,2}(:,18),0.15*p.omega,0.15*p.omega,'false');      %cdc20
            max_tmp_cycE=peakfinder(cells{ii,2}(:,12),0.15*p.omega,0.15*p.omega,'false');  %cycE
            %min_tmp=peakfinder(cells{ii,2}(:,14),[],30,-1);   %cdc20
 
            if length(find(cells{ii,2}(:,23)>25))>20 %dying cell   %baxm >20
                cells{ii,5}=false;
                death=find(cells{ii,2}(:,23)>25);
                current_time(ii,1)=cells{ii,3}(death(1));
                cells{ii,2}=cells{ii,2}(1:death(1),:);
                cells{ii,3}=cells{ii,3}(1:death(1));
                current_time(ii,2)=0;
                cell_count_tmp(ii,:)=[cells{ii,3}(death(1)), -1];                
            elseif size(max_tmp_cdc,1)>cells{ii,4} && length(find(cells{ii,2}(:,23)>25))<20 %dividing cell
%                 if isempty(max_tmp_cycE)
%                      current_time(ii)=cells{ii,3}(end);
%                      cell_count_tmp(ii,:)=[cells{ii,3}(end), 0];
                if size(max_tmp_cdc,1)>cells{ii,4}% && isempty(max_tmp_cycE)==false%make sure cycE has cycled
                    div_ind = (max_tmp_cdc(cells{ii,4}+1));
                    if isempty(div_ind)
                        current_time(ii,1)=cells{ii,3}(end);
                        cell_count_tmp(ii,:)=[cells{ii,3}(end), 0];
                    elseif cells{ii,4} == length(max_tmp_cycE)
                        current_time(ii)=cells{ii,3}(end);
                        cell_count_tmp(ii,:)=[cells{ii,3}(end), 0];
                    else
                        cells{ii,2}=cells{ii,2}(1:div_ind(1),:);
                        cells{ii,3}=cells{ii,3}(1:div_ind(1));
                        cells{ii,4}=cells{ii,4}+1;
                        current_time(ii,1)=cells{ii,3}(end);
                        cell_count_tmp(ii,:)=[cells{ii,3}(end), 1];
                        %add a new cell
                        tmp_num=size(cells,1);
                        cells{tmp_num+1,1}=p;
                        cells{tmp_num+1,2}=initial_divided;
                        cells{tmp_num+1,3}=current_time(ii,1);
                        cells{tmp_num+1,4}=0;
                        cells{tmp_num+1,5}=true;
                        current_time(tmp_num+1,:)=[current_time(ii,1), 1];
                    end  
                else               
                    current_time(ii)=cells{ii,3}(end);
                    cell_count_tmp(ii,:)=[cells{ii,3}(end), 0];
                end

            else %size(max_tmp_cdc,1)==cells{ii,4} && length(find(cells{ii,2}(:,23)>25))<20 %quiescent cell
                current_time(ii)=cells{ii,3}(end);
                cell_count_tmp(ii,:)=[cells{ii,3}(end), 0];
            end
       else
           cell_count_tmp(ii,:)=[cells{ii,3}(end), 0];
       end
       combo = [ii, cells{ii,3}(end)];
       disp(combo)
    end
    for ii=1:size(cell_count_tmp,1)
        if cell_count_tmp(ii,2)==1
            ins_ind=find(cell_count_tmp(ii,1)>cell_count(:,1));
            current_cell_num=cell_count(ins_ind(end),2);
            cell_count=insertrows(cell_count,cell_count_tmp(ii,:),ins_ind(end));
            cell_count(ins_ind(end)+1,2)=current_cell_num;
            cell_count(ins_ind(end)+1:end,2)=cell_count(ins_ind(end)+1:end,2)+1;
        elseif cell_count_tmp(ii,2)==-1
            ins_ind=find(cell_count_tmp(ii,1)>cell_count(:,1));
            current_cell_num=cell_count(ins_ind(end),2);
            cell_count=insertrows(cell_count,cell_count_tmp(ii,:),ins_ind(end));
            cell_count(ins_ind(end)+1,2)=current_cell_num;
            cell_count(ins_ind(end)+1:end,2)=cell_count(ins_ind(end)+1:end,2)-1;
        elseif cell_count_tmp(ii,2)==0
            ins_ind=find(cell_count_tmp(ii,1)>cell_count(:,1));
            current_cell_num=cell_count(ins_ind(end),2);
            cell_count=insertrows(cell_count,cell_count_tmp(ii,:),ins_ind(end));
            cell_count(ins_ind(end)+1,2)=current_cell_num;
        end
         
    end 
    [~,ia]=unique(cell_count(:,1));
    cell_count=cell_count(ia,:);
    %disp(cell_count(end,:))    
     
end