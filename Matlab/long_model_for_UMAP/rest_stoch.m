function [new_t, new_x] = rest_stoch(t,x,dt)
new_t = t(1):dt:t(end);
new_x = zeros(length(new_t),size(x,2));

approx_step = length(find(t<dt))+100;

nind = zeros(length(new_t),1);
oind=1;
for ii = 1:length(new_t)
    if oind+approx_step<length(t)
        tmp = find(t(oind:oind+approx_step)<=new_t(ii));
    else
        tmp = find(t(oind:end)<=new_t(ii));
    end
    
    if ii ==1 
        nind(ii) = tmp(end);
    else
        nind(ii) = tmp(end)+(oind-1);
    end
    oind=nind(ii);
       
end
    
new_x=x(nind,:);