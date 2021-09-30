clear
clc
close all


wmax = 2; %rad/s
m = 2.6; %kg
b = 10/100; %cm/100 = meters
c = 20/100; %cm/100 = meters
Imax = m/12*(b^2+c^2);


dt = 0.000001;
L = dt:dt:0.1; %N-m

tstar = wmax*Imax./(L);

plottool(1,'Ttumble',18,'Time (min)','Max Moment (N-m)')
plot(tstar/60,L)


