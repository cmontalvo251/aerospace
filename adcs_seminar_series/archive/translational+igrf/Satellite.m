function dstatedt = Satellite(t,state)
%%%stateinitial = [x0;y0;z0;xdot0;ydot0;zdot0];
global BxI ByI BzI

x = state(1);
y = state(2);
z = state(3);
%xdot = state(4);
%ydot = state(5);
%zdot = state(6);

%%%inertia and mass
m = 2.6; %%kilograms

%%%Kinematics
vel = state(4:6);

%%%Gravity Model
planet
r = state(1:3); %% r = [x;y;z]
rho = norm(r);
rhat = r/rho;
Fgrav = -(G*M*m/rho^2)*rhat;

%%%Call the magnetic field model
%%%Convert Cartesian x,y,z into Lat,Lon, Alt
phiE = 0;
thetaE = acos(z/rho);
psiE = atan2(y,x);
latitude = 90-thetaE*180/pi;
longitude = psiE*180/pi;
rhokm = (rho)/1000;
[BN,BE,BD] = igrf('01-Jan-2020',latitude,longitude,rhokm,'geocentric');
%%%Convert NED (North East Down to X,Y,Z in ECI frame)
%%%First we need to create a rotation matrix from the NED frame to the 
%%%inertial frame
BNED = [BN;BE;BD]; 
BI = TIB(phiE,thetaE+pi,psiE)*BNED;
%BI = eye(3)*BNED;
BxI = BI(1);
ByI = BI(2);
BzI = BI(3);

%%%Dynamics
F = Fgrav;
accel = F/m;

%%%Return derivatives vector
dstatedt = [vel;accel];

