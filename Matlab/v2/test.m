function [dydt]=test(t,y,S)
dydt = zeros(size(y));

k_1r = 0.5;
k_2r = 0.4;
k_1e = 0.15;
k_2e = 0.09;
K_m1e = 0.01;
K_m2e = 0.01;
E_tot = 1;


R=y(1);
E=y(2);

dydt(1) = k_1r*(E_tot - E) + k_1r*S - k_2r*R;
dydt(2) = k_1e*(E_tot - E)/(K_m1e+(E_tot - E)) - (k_2e*E*R)/(K_m2e+E);