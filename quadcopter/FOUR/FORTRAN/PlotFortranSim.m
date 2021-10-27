clear
clc 
close all

system('del Run.exe');
system('del SimulationResults.txt');
system('gfortran.exe -o Run.exe QuadSim[2].f95 -w');
system('Run.exe');

data=dlmread('SimulationResults.txt');

t = data(:,1);
xpos = data(:,2);
ypos = data(:,3);
zpos = data(:,4);
phi  = data(:,5);
theta= data(:,6);
psi=   data(:,7);
u=     data(:,8);
v=     data(:,9);
w=     data(:,10);
p=     data(:,11);
q=     data(:,12);
r=     data(:,13);
omega1= data(:,14);
omega2= data(:,15);
omega3= data(:,16);
omega4= data(:,17);

xc = [0 10 10 0];
yc = [10 10 0 0];
zc = [-10 -20 -30 0];

%%% Plot xyz
fig = figure();
set(fig,'color','white');
plot(t,xpos,'k-')
hold on
plot(t,ypos,'r--')
plot(t,zpos,'b-')
%axis([0 10 0 12])
set(gca,'FontSize',18)
xlabel('Time, sec'); ylabel('Distance, m'); title('x,y,z')
legend('x','y','z')
grid on

%%% Plot y versus x 
fig = figure();
set(fig,'color','white');
plot(xpos,ypos,'k-')
hold on
for idx=1:4
plot(xc(idx), yc(idx),'*r','Linewidth',1)
end
plot(0,0,'bs','Linewidth',2')
set(gca,'FontSize',18)
xlabel('X Distance, m'); ylabel('Y Distance, m');

grid on

% Position of the Quad, x y and z 
fig = figure();
set(fig,'color','white');
plot3(xpos, ypos, zpos, 'k-','Linewidth',2)
hold on
for idx=1:4
plot3(xc(idx), yc(idx), zc(idx),'*r','Linewidth',2)
end
plot3(0,0,0,'bs','Linewidth',2)
set(gca,'Fontsize',18)
xlabel( 'X  (m)')
ylabel( 'Y  (m)')
title('Position')
grid on

%%% Plot phi,theta,psi
fig = figure();
set(fig,'color','white');
plot(t,phi*180/(pi),'k-')
hold on
plot(t,theta*180/(pi),'r--')
plot(t,psi*180/(pi),'b-')
set(gca,'FontSize',18)
xlabel('Time, sec'); ylabel('Angle, deg'); title('\phi,\theta,\psi')
legend('\phi','\theta','\psi')
grid on

%%% Plot uvw
fig = figure();
set(fig,'color','white');
plot(t,u,'k-')
hold on
plot(t,v,'r--')
plot(t,w,'b-')
set(gca,'FontSize',18)
xlabel('Time, sec'); ylabel('Velocity, m/s'); title('u,v,w')
legend('u','v','w')
grid on

%%% Plot pqr
fig = figure();
set(fig,'color','white');
plot(t,p*180/(pi),'k-')
hold on
plot(t,q*180/(pi),'r--')
plot(t,r*180/(pi),'b-')
set(gca,'FontSize',18)
xlabel('Time, sec'); ylabel('Rate, deg/s'); title('p,q,r')
legend('p','q','r')
grid on

%%% Plot omega
fig = figure();
set(fig,'color','white');
plot(t,omega1,'g-')
hold on
plot(t,omega2,'r-')
plot(t,omega3,'b-')
plot(t,omega4,'k-')
!axis([0 10 0 1000])
set(gca,'FontSize',18)
xlabel('Time, sec'); ylabel('Angular Velocity, rad/s');
legend('\Omega_1','\Omega_2','\Omega_3','\Omega_4')
grid on