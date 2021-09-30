function Controller_Tuning()
purge

[tout,xout] = ode45(@Derivatives,[0 30],zeros(6,1));
yaw = xout(:,3);
xw = 10;
yw = 10;

x = xout(:,1);
y = xout(:,2);

delW = sqrt((xw-x).^2+(yw-y).^2);

plot(tout,delW)
xlabel('Time (sec)')
ylabel('Distance to Waypoint (m)')

figure
plot(tout,yaw*180/pi)
xlabel('Time (sec)')
ylabel('Yaw Angle (deg)')

figure
plot(x,y)
xlabel('X (m)')
ylabel('Y (m)')

function statedot = Derivatives(t,state)
x = state(1);
y = state(2);
yaw = state(3); 
yawdot = state(4);
u = state(5);
v = state(6);

b = (38/12)/3.28; %%%diameter of blimp
Vol = (4/3)*pi*(b/2)^3; %%%volume
mstructure = 32/1000; %%%%mass of structure
mgondola = 68/1000; %%%mass of gondola and batteries
density_helium = 0.18; %%%%kg/m^3 of helium
m = density_helium*Vol + mstructure + mgondola; %%%total mass
I =  (2/5)*m*(b/2)^2; %%%Inertia of a sphere
d = (5/12)/3.28; %%%Moment arm
k = (8/3)*(b/2)^3;
k1 = k;
k2 = k;
k3 = 0;
rho = 1.225;

sig = .045;
Cnr = .023;

Na = -0.5*rho*(Vol)^(5/3)*abs(yawdot)*Cnr*yawdot;

xw = 10;
yw = 10;

u_command = 1;
kp = 0.008;

kpN = -0.009;
kdN = -0.014;

psic = atan2(yw - y,xw-x);
%psic = 45*pi/180;
spsi = sin(yaw);
cpsi = cos(yaw);
cpsic = cos(psic);
spsic = sin(psic);
delpsi = atan2(spsi*cpsic-cpsi*spsic,cpsi*cpsic+spsi*spsic);
TN = kpN*delpsi + kdN*yawdot;
T1 = (pi-abs(delpsi))/pi*(kp*(u_command-u)) - TN/d;
T2 = (pi-abs(delpsi))/pi*(kp*(u_command-u)) + TN/d;

if T1 > sig
    T1 = sig;
end
if T2 > sig
    T2 = sig;
end

N = (T2-T1)*d + Na;    % The total moment

yawdbldot = N/(k3 + I);    % EOMs

C_D = 0.125;
rho = 1.225;
V = sqrt(u^2+v^2);
X = -0.5*rho*abs(V)*Vol^(2/3)*C_D*u + T1 + T2;
Y = -0.5*rho*abs(V)*Vol^(2/3)*C_D*v;

udot = X/(m+k1) + yawdot*v;
vdot = Y/(m+k2) - yawdot*u;

xdot = cos(yaw)*u - sin(yaw)*v;
ydot = sin(yaw)*u + cos(yaw)*v;

statedot = [xdot;ydot;yawdot;yawdbldot;udot;vdot];
