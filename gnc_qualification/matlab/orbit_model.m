function [x,y,z,t,T_orbit] = orbit_model(apogee,perigee,N,flag)
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

if flag
  plot3(x,y,z)
  axis equal
end

end

%#Equations of Motion
function dstatedt = Derivatives(t,state)
MEarth = 5.9736e24;% #earth mass
REarth = 6357000.0;% #earth radius
G = 6.6742e-11; %%Gravitational constant

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