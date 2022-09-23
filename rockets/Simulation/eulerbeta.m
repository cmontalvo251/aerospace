function [veuler, teuler] = eulerbeta(time,V,Tk)
close all
clc,clear
data = dlmread('udata.csv');

x = data(:,1);
y = data(:,2);

xfit = x(x>0.01);
yfit = y(x>0.01);

xfit = xfit(xfit < 25.8);
yfit = yfit(xfit < 25.8);

time = xfit;
V = yfit;


g = 9.81; %m/s^2
m = 2.316; %kg
d = 3.125*.0254; %m
rho = 1.226601593; %kg/m^3
theta = 45*(pi/180); %launch angle
mfuel = 0.640699223; %kg
totaltburn = 2.56;%s
mdot = mfuel/totaltburn;



Cdguess = 0.23;
t=0.00;
tfinal=25.8;
dt = 0.05;
nsteps = round((tfinal-t)/dt);
veuler(1) = V(1,1);
teuler = t:dt:(tfinal-dt);

tburn = [0, .176, .399, .744, 1.214, 1.601, 2.563,5,7.5,10,15,25];
Tvec = [0, 1980.928, 1160.584, 1260.273,1143.485, 1093.96,200,114,94,84,78,0];


for k = 1:nsteps -1
    if teuler(k) <= 25.75
        T = interp1(tburn,Tvec,teuler(k));
	 m = m-mdot*dt;
    else 
	T = 0;
    end

    vprime = T/m -veuler(k)^2*((pi*rho*d^2*Cdguess)/(8*m))-g*sin(theta);
    veuler(k+1) = veuler(k) + vprime*dt;

end
  


