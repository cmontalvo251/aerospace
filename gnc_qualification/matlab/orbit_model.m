function [x,y,z,t,r,T_orbit,vx,vy,vz,v] = orbit_model(apogee,perigee,N)
constants

%%%COnvert
apogee = REarth + apogee*1000;
perigee = REarth + perigee*1000;

%%Integrate the equations of motion
semi_major = (apogee+perigee)/2.0;
T_orbit = 2*pi*semi_major^(3./2)/sqrt(G*MEarth);
t = linspace(0,T_orbit,N);
dt = t(2)-t(1);

%%%Assume we start at perigee
x0 = perigee;
y0 = 0;
z0 = 0;
velx = 0;
ecc = (apogee-perigee)/(apogee+perigee);
par = apogee*(1-ecc);
hmom = sqrt(par*MEarth*G);
vely = hmom/perigee;
velz = 0;
stateinitial = [x0;y0;z0;velx;vely;velz];

%functionHandle,tspan,xinitial,timestep,extraparameters,next,quat
[t,stateout] = odeK4(@Derivatives,[t(1) t(end)],stateinitial,dt,[],10,0);
x = stateout(1,:);
y = stateout(2,:);
z = stateout(3,:);
r = sqrt(x.^2+y.^2+z.^2);
vx = stateout(4,:);
vy = stateout(5,:);
vz = stateout(6,:);
v = sqrt(vx.^2+vy.^2+vz.^2);

%%%Plot the orbit
plot3(x,y,z,'LineWidth',2)
set(gcf,'color','white')
xlabel('X (m)')
ylabel('Y (m)')
zlabel('Z (m)')
grid on
title('Orbit')
axis equal
%%%Plot the velocity
figure()
set(gcf,'color','white')
plot(t,vx,'LineWidth',2)
hold on
plot(t,vy,'LineWidth',2)
plot(t,vz,'LineWidth',2)
plot(t,sqrt(vx.^2+vy.^2+vz.^2),'LineWidth',2)
xlabel('Time (sec)')
ylabel('Velocity (m/s)')
legend('Vx','Vy','Vz','V')
grid on
title('Velocity of Orbit')
%%%Plot Altitude
figure()
set(gcf,'color','white')
plot(t,(r-REarth)/1000,'LineWidth',2)
xlabel('Time (sec)')
ylabel('Altitude (km)')
grid on

end

%#Equations of Motion
function dstatedt = Derivatives(t,state)
constants

%%Extract State
x = state(1);
y = state(2);
z = state(3);
velx = state(4);
vely = state(5);
velz = state(6);

%#Kinematics
xdot = velx;
ydot = vely;
zdot = velz;

%#Dynamics
%#Gravitational Acceleration
rSat = sqrt(x^2 + y^2 + z^2);
if rSat < REarth
  xdbldot = 0;
  ydbldot = 0;
  zdbldot = 0;
else
  xdbldot = (-G*MEarth*(x/rSat))/(rSat^2);
  ydbldot = (-G*MEarth*(y/rSat))/(rSat^2);
  zdbldot = (-G*MEarth*(z/rSat))/(rSat^2);
end

dstatedt = [xdot;ydot;zdot;xdbldot;ydbldot;zdbldot];

end