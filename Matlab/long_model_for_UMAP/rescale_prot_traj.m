st_tx_erl = cell(3);
st_tx_dmso = cell(3);
max_all = [];
dt=0.01;
for ii = 1:3
    for jj = 1
        load(strcat('erl_long_4_29_',num2str(ii),'_',num2str(jj),'.mat'));
        [new_t, new_x] = rest_stoch(t,x,dt);
        new_x=new_x./200;
        
        dlmwrite(strcat('erl_long_4_29_',num2str(ii),'_',num2str(jj),'.txt'),new_x(:,1:18))
        ii
        jj
    end
    for jj =1
        load(strcat('dmso_long_4_29_',num2str(ii),'_',num2str(jj),'.mat'));
        [new_t, new_x] = rest_stoch(t,x,dt);
        new_x=new_x./200;

        dlmwrite(strcat('dmso_long_4_29_',num2str(ii),'_',num2str(jj),'.txt'),new_x(:,1:18))
    end
end 

dlmwrite(strcat('t_4-9_',num2str(ii),'_',num2str(jj),'.txt'),new_t)