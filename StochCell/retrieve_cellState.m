IMT_eoe_counts = cell(1,4);
IMT_eoe_counts{1,1}='End_of_Expt';
IMT_eoe_counts{1,2}='Death';
IMT_eoe_counts{1,3}='Birthtime';
IMT_eoe_counts{1,4}='Lifespan';

divSum=0;
for ii = 1:size(cells,1)
    divSum=divSum+cells{ii,4};
end


for ii = 1:size(cells,1)
    end_index=find(cells{ii,3}<run_time);
    tmp_counts=cell(1,4);
    if isempty(end_index)==false
        max_tmp_cdc=peakfinder(cells{ii,2}(:,18),[],0.15*p.omega,1);  
        max_tmp_cycE=peakfinder(cells{ii,2}(:,12),[],0.15*p.omega,1);

        div_ind = max_tmp_cdc;
        if ii<=init_cell   %initial cells don't include first point
            if length(div_ind)>=2 && cells{ii,5}==1
                for jj = 1:length(div_ind)
                    
                    if jj<length(div_ind)
                        tmp_counts{1,1}=false;
                        tmp_counts{1,2}=false;
                        tmp_counts{1,3}=cells{ii,3}(div_ind(jj));
                        tmp_counts{1,4}=cells{ii,3}(div_ind(jj+1))-cells{ii,3}(div_ind(jj));
                    else
                        tmp_counts{1,1}=true;
                        tmp_counts{1,2}=false;
                        tmp_counts{1,3}=cells{ii,3}(div_ind(jj));
                        tmp_counts{1,4}=run_time-cells{ii,3}(div_ind(jj));
                    end
                    IMT_eoe_counts=[IMT_eoe_counts; tmp_counts];                    
                end
            elseif length(div_ind)==1 && cells{ii,5}==1
                tmp_counts=cell(1,4);
                tmp_counts{1,1}=true;
                tmp_counts{1,2}=false;
                tmp_counts{1,3}=cells{ii,3}(div_ind(1));
                tmp_counts{1,4}=run_time-cells{ii,3}(div_ind(1));
                IMT_eoe_counts=[IMT_eoe_counts; tmp_counts]; 
            elseif isempty(div_ind) && cells{ii,5}==1
                tmp_counts{ii+1,1}=true;
                tmp_counts{ii+1,2}=false;
                tmp_counts{ii+1,3}=cells{ii,3}(1);
                tmp_counts{ii+1,4}=run_time-cells{ii,3}(1);
                IMT_eoe_counts=[IMT_eoe_counts; tmp_counts]; 
            elseif length(div_ind)>=2 && cells{ii,5}==0
                for jj = 1:length(div_ind)
                    tmp_counts=cell(1,4);
                    if jj<length(div_ind)
                        tmp_counts{1,1}=false;
                        tmp_counts{1,2}=false;
                        tmp_counts{1,3}=cells{ii,3}(div_ind(jj));
                        tmp_counts{1,4}=cells{ii,3}(div_ind(jj+1))-cells{ii,3}(div_ind(jj));
                    else
                        tmp_counts{1,1}=false;
                        tmp_counts{1,2}=true;
                        tmp_counts{1,3}=cells{ii,3}(div_ind(jj));
                        tmp_counts{1,4}=cells{ii,3}(end)-cells{ii,3}(div_ind(jj));
                    end
                    IMT_eoe_counts=[IMT_eoe_counts; tmp_counts];                    
                end
            elseif length(div_ind)==1 && cells{ii,5}==0
                tmp_counts=cell(1,4);
                tmp_counts{1,1}=false;
                tmp_counts{1,2}=true;
                tmp_counts{1,3}=cells{ii,3}(div_ind(1));
                tmp_counts{1,4}=cells{ii,3}(end)-cells{ii,3}(div_ind(1));
                IMT_eoe_counts=[IMT_eoe_counts; tmp_counts]; 
%            else %isempty(div_ind) && cells{ii,5}==0
%                 tmp_counts{ii+1,1}=false;
%                 tmp_counts{ii+1,2}=true;
%                 tmp_counts{ii+1,3}=cells{ii,3}(1);
%                 tmp_counts{ii+1,4}=cells{ii,3}(end)-cells{ii,3}(1);
%                 IMT_eoe_counts=[IMT_eoe_counts; tmp_counts];             
            end
        else
            if length(div_ind)>=1 && cells{ii,5}==1
                for jj = 1:length(div_ind)
                    tmp_counts=cell(1,4);
                    if jj==1
                        tmp_counts{1,1}=false;
                        tmp_counts{1,2}=false;
                        tmp_counts{1,3}=cells{ii,3}(1);
                        tmp_counts{1,4}=cells{ii,3}(div_ind(jj))-cells{ii,3}(1);
                    elseif jj==length(div_ind)
                        tmp_counts{1,1}=true;
                        tmp_counts{1,2}=false;
                        tmp_counts{1,3}=cells{ii,3}(div_ind(jj));
                        tmp_counts{1,4}=run_time-cells{ii,3}(div_ind(jj));                        
                    else
                        tmp_counts{1,1}=false;
                        tmp_counts{1,2}=false;
                        tmp_counts{1,3}=cells{ii,3}(div_ind(jj));
                        tmp_counts{1,4}=cells{ii,3}(div_ind(jj+1))-cells{ii,3}(div_ind(jj));
                    end
                    IMT_eoe_counts=[IMT_eoe_counts; tmp_counts];                    
                end
            elseif length(div_ind)==1 && cells{ii,5}==1                
                tmp_counts{1,1}=false;
                tmp_counts{1,2}=false;
                tmp_counts{1,3}=cells{ii,3}(1);
                tmp_counts{1,4}=cells{ii,3}(div_ind(1))-cells{ii,3}(1);
                IMT_eoe_counts=[IMT_eoe_counts; tmp_counts]; 
                tmp_counts{1,1}=true;
                tmp_counts{1,2}=false;
                tmp_counts{1,3}=cells{ii,3}(div_ind(1));
                tmp_counts{1,4}=run_time-cells{ii,3}(div_ind(1));
                IMT_eoe_counts=[IMT_eoe_counts; tmp_counts]; 
            elseif isempty(div_ind) && cells{ii,5}==1
                tmp_counts{ii+1,1}=true;
                tmp_counts{ii+1,2}=false;
                tmp_counts{ii+1,3}=cells{ii,3}(1);
                tmp_counts{ii+1,4}=run_time-cells{ii,3}(1);
                IMT_eoe_counts=[IMT_eoe_counts; tmp_counts]; 
            elseif length(div_ind)>=2 && cells{ii,5}==0
                for jj = 1:length(div_ind)
                    tmp_counts=cell(1,4);
                    if jj==1
                        tmp_counts{1,1}=false;
                        tmp_counts{1,2}=false;
                        tmp_counts{1,3}=cells{ii,3}(1);
                        tmp_counts{1,4}=cells{ii,3}(div_ind(jj))-cells{ii,3}(1);
                    elseif jj==length(div_ind)
                        tmp_counts{1,1}=false;
                        tmp_counts{1,2}=true;
                        tmp_counts{1,3}=cells{ii,3}(div_ind(jj));
                        tmp_counts{1,4}=cells{ii,3}(end)-cells{ii,3}(div_ind(jj));                        
                    else
                        tmp_counts{1,1}=false;
                        tmp_counts{1,2}=false;
                        tmp_counts{1,3}=cells{ii,3}(div_ind(jj));
                        tmp_counts{1,4}=cells{ii,3}(div_ind(jj+1))-cells{ii,3}(div_ind(jj));
                    end
                    IMT_eoe_counts=[IMT_eoe_counts; tmp_counts];                    
                end
            elseif length(div_ind)==1 && cells{ii,5}==0
                tmp_counts=cell(1,4);
                tmp_counts{1,1}=false;
                tmp_counts{1,2}=false;
                tmp_counts{1,3}=cells{ii,3}(1);
                tmp_counts{1,4}=cells{ii,3}(div_ind(1))-cells{ii,3}(1);
                IMT_eoe_counts=[IMT_eoe_counts; tmp_counts]; 
                tmp_counts{1,1}=false;
                tmp_counts{1,2}=true;
                tmp_counts{1,3}=cells{ii,3}(div_ind(1));
                tmp_counts{1,4}=cells{ii,3}(end)-cells{ii,3}(div_ind(1));
                IMT_eoe_counts=[IMT_eoe_counts; tmp_counts]; 
            else %isempty(div_ind) && cells{ii,5}==0
                tmp_counts{ii+1,1}=false;
                tmp_counts{ii+1,2}=true;
                tmp_counts{ii+1,3}=cells{ii,3}(1);
                tmp_counts{ii+1,4}=cells{ii,3}(end)-cells{ii,3}(1);
                IMT_eoe_counts=[IMT_eoe_counts; tmp_counts];             
            end
        end            
    end
end

IMT_eoe_counts(all(cellfun('isempty',IMT_eoe_counts),2),:)=[];

cell_count_new = cell(122,2);
cell_count_new{1,1}='Time_h';
cell_count_new{1,3}='AKTERK_max';
tmp_AKT=100*rem(AKT,1);
tmp_ERK=100*rem(ERK,1);
cell_count_new{1,2}=sprintf('AKTERK_treat');
cell_count_new{2,2} = init_cell;
for ii = 0:120
    cell_count_new{ii+2,1}=ii;
end

for ii = 3:122
    tmp = cell_count(find(cell_count(:,1)>ii-1 & cell_count(:,1)<ii),2);
    if isempty(tmp)
        cell_count_new{ii,2}=cell_count_new{ii-1,2};
    else
        cell_count_new{ii,2}=tmp(1);
    end
end
for ii = 2:122
    cell_count_new{ii,2}=log2(double(cell_count_new{ii,2}));
end
x1=(0:1:120)';
a=init_cell;
b=0.03913;
y1=a*exp(b*x1);

for ii = 2:122
    cell_count_new{ii,3}=log2(double(y1(ii-1)));
end

tots=sprintf('AKTERK_%02d_%02d',tmp_AKT,tmp_ERK);

xlswrite(strcat('C:\Users\Austin\Documents\pc9_cell_runs\',tots,'_total'),cell_count_new);
xlswrite(strcat('C:\Users\Austin\Documents\pc9_cell_runs\',tots,'_state'),IMT_eoe_counts);
xlswrite(strcat('C:\Users\Austin\Documents\pc9_cell_runs\tmp_total'),cell_count_new);
xlswrite(strcat('C:\Users\Austin\Documents\pc9_cell_runs\tmp_state'),IMT_eoe_counts);
