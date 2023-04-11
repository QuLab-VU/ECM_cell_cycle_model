Initial = [5,0];
TIMES = (0:0.01:100);
S=5;

[t,y] = ode23(@(t,y) test(t,y,S), TIMES, Initial); 

sig=0:0.01:2;

Initial = [5,0];
resp=[];
for ii=1:length(sig)
    S=sig(ii);
    [t,y] = ode23(@(t,y) test(t,y,S), TIMES, Initial); 
    resp(ii)=y(end,1);
end
figure; hold on
plot(resp)
Initial = [0,0];
resp=[];
for ii=1:length(sig)
    S=sig(ii);
    [t,y] = ode23(@(t,y) test(t,y,S), TIMES, Initial); 
    resp(ii)=y(end,1);
end
plot(resp)