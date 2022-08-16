function [IMTs, g1, s_m, G1_S, mit] = getIMTs_long(t,y,p)
IMTs=[];
IMT_tmp=[];
g1 = [];
g1_tmp =[];
s_m = [];
s_m_tmp = [];


max_tmp_cycE=peakfinder(y(:,12),[],0.2*p.omega,'false');
max_tmp_cdc=peakfinder(y(:,18),[],0.15*p.omega,'false');
G1_S = max_tmp_cycE;
mit = max_tmp_cdc;
if ~isempty(max_tmp_cdc) && ~isempty(max_tmp_cycE)            
    div_ind = max_tmp_cdc;
    if length(div_ind)>=2
        for jj =2:length(div_ind)
            IMT_tmp(jj-1,1)=t(div_ind(jj))-t(div_ind(jj-1));
        end
    end

    for jj=1:length(max_tmp_cycE)
        tmp = max_tmp_cdc(max_tmp_cdc > max_tmp_cycE(jj));
        if ~isempty(tmp)
            s_m_tmp(jj) = t(tmp(1)) - t(max_tmp_cycE(jj));
        end
    end

    for jj=1:length(max_tmp_cdc)
        tmp = max_tmp_cycE(max_tmp_cycE > max_tmp_cdc(jj));
        if ~isempty(tmp)
            g1_tmp(jj) = t(tmp(1)) - t(max_tmp_cdc(jj));
        end
    end
end
IMTs=[IMTs; IMT_tmp];
g1 = [g1, g1_tmp];
s_m = [s_m, s_m_tmp];
IMT_tmp=[];
g1_tmp = [];
s_m_tmp = [];

end