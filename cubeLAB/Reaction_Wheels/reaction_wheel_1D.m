%%%Reaction Wheel 1D Code
clear
clc
close all

%%%TIME STUFF
dt = 0.1;
tfinal = 100;
tout = 0:dt:tfinal;

%%%INITIAL CONDITIONS
p_sat0 = 0;
phi_sat0 = 0;
phi_rw0 = 0;
p_rw0 = 0;
xinitial = [phi_sat0;p_sat0;phi_rw0;p_rw0];

%%%PRE-ALLOCATE
xout = zeros(length(xinitial),length(tout));
xout(:,1) = xinitial;

%%%RK4
for idx = 1:length(tout)-1
  xk = xout(:,idx);
  tk = tout(idx);
  k1 = Derivatives(tk,xk);
  k2 = Derivatives(tk+dt/2,xk+k1*dt/2); 
  k3 = Derivatives(tk+dt/2,xk+k2*dt/2); 
  k4 = Derivatives(tk+dt,xk+k3*dt); 
  phi = (1/6)*(k1 + 2*k2 + 2*k3 + k4);
  xk1 = xk + phi*dt;
  xout(:,idx+1) = xk1;
end

%%%PLOT STUFF
phi_sat = xout(1,:);
p_sat = xout(2,:);
phi_rw = xout(3,:);
p_rw = xout(4,:);

plottool(1,'Rotation Angles',12,'Time (sec)','Angles (deg)')
plot(tout,phi_sat*180/pi,'b-')
plot(tout,phi_rw*180/pi,'r-')
legend('Satellite','Reaction Wheel')

plottool(1,'Anglular Rates',12,'Time (sec)','Anglular Rates (rad/s)')
plot(tout,p_sat*180/pi,'b-')
plot(tout,p_rw*180/pi,'r-')
legend('Satellite','Reaction Wheel')


%%%Plot total momentum of system
%%%Masses and Inertias (kg-m^2)
Irw_cg = 0.5; %%%Moment of inertia of the reaction wheel about the center of mass of the reaction wheel
mass_rw = 0.1; %%Mass of reaction wheel
d = 0.2; %%%Distance from Center of Mass
Irw_sat = Irw_cg + mass_rw*d^2;
I_sat = 1; %%%Moment of inertia of the cubesat with out reaction wheels
I_sys = I_sat + Irw_sat; %%%Total moment of inertia of system

H = I_sys*p_sat + Irw_cg*p_rw;

plottool(1,'Total Momentum',12,'Time (sec)','Momentum (N-m-s)')
plot(tout,H)




