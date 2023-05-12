[p] = set_params(200,1,1);
p.v_s1p21=0.3;
p.AKT = 0.3;
p.ERK = 0.2;
p.FAK = 0.4;
figure;hold on
[t,y] = ode23(@(t,y) cell_ode_death_test(t,y,p), TIMES, Initial);
plot(t,y(:,12))
p.AKT = 0.6;
p.ERK = 0.2;
p.FAK = 0.6;
[t,y] = ode23(@(t,y) cell_ode_death_test(t,y,p), TIMES, Initial);
plot(t,y(:,12))
p.AKT = 0.6;
p.ERK = 0.2;
p.FAK = 1;
[t,y] = ode23(@(t,y) cell_ode_death_test(t,y,p), TIMES, Initial);
plot(t,y(:,12))

[p] = set_params(200,1,1);
p.AKT = 0.3;
p.ERK = 0.2;
p.FAK = 0.4;
figure;hold on
[t,y] = ode23(@(t,y) cell_ode(t,y,p), TIMES, Initial);
plot(t,y(:,12))
p.AKT = 0.6;
p.ERK = 0.2;
p.FAK = 0.6;
[t,y] = ode23(@(t,y) cell_ode(t,y,p), TIMES, Initial);
plot(t,y(:,12))
p.AKT = 0.6;
p.ERK = 0.2;
p.FAK = 1;
[t,y] = ode23(@(t,y) cell_ode(t,y,p), TIMES, Initial);
plot(t,y(:,12))