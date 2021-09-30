function dstatedt = Satellite(t,state)
%%%stateinitial = [x0;y0;z0;xdot0;ydot0;zdot0];
global BB invI I m nextMagUpdate lastMagUpdate lastSensorUpdate 
global nextSensorUpdate BfieldMeasured pqrMeasured
x = state(1);
y = state(2);
z = state(3);
%xdot = state(4);
%ydot = state(5);
%zdot = state(6);
q0123 = state(7:10);
p = state(11);
q = state(12);
r = state(13);
pqr = state(11:13);

%%%Translational Kinematics
vel = state(4:6);
%%%Rotational Kinematics

PQRMAT = [0 -p -q -r;p 0 r -q;q -r 0 p;r q -p 0];
q0123dot = 0.5*PQRMAT*q0123;

%%%Gravity Model
planet
r = state(1:3); %% r = [x;y;z]
rho = norm(r);
rhat = r/rho;
Fgrav = -(G*M*m/rho^2)*rhat;

%%%Call the magnetic field model
if t >= lastMagUpdate
    lastMagUpdate = lastMagUpdate + nextMagUpdate;
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
    BB = TIBquat(q0123)'*BI;
    %%%Convert to Tesla
    BB = BB*1e-9;
end

if t >= lastSensorUpdate
    lastSensorUpdate = lastSensorUpdate + nextSensorUpdate;
    [BfieldMeasured,pqrMeasured] = Sensor(BB,pqr); 
end

%%%Translational Dynamics
F = Fgrav;
accel = F/m;

%%%Magtorquer Model
LMN_magtorquers = [0;0;0];

%%%Rotational Dynamics
H = I*pqr;
pqrdot = invI*(LMN_magtorquers - cross(pqr,H));

%%%Return derivatives vector
dstatedt = [vel;accel;q0123dot;pqrdot];

