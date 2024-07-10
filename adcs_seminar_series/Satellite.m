function dstatedt = Satellite(t,state)
%%%stateinitial = [x0;y0;z0;xdot0;ydot0;zdot0];
global BB invI I m mu nextMagUpdate lastMagUpdate lastSensorUpdate maxSpeed
global nextSensorUpdate BfieldMeasured pqrMeasured BfieldNav pqrNav
global BfieldNavPrev pqrNavPrev current Is Ir1Bcg Ir2Bcg Ir3Bcg n1 n2 n3
global maxAlpha Ir1B Ir2B Ir3B ptpMeasured ptpNavPrev ptpNav rwalphas
global fsensor MagFieldBias AngFieldBias EulerBias
global MagFieldNoise AngFieldNoise EulerNoise R Amax lmax CD

x = state(1);
y = state(2);
z = state(3);
%xdot = state(4);
%ydot = state(5);
%zdot = state(6);
q0123 = state(7:10);
ptp = Quaternions2EulerAngles(q0123')';
p = state(11);
q = state(12);
r = state(13);
pqr = state(11:13);
w123 = state(14:16);

%%%Translational Kinematics
vel = state(4:6);
%%%Rotational Kinematics

PQRMAT = [0 -p -q -r;p 0 r -q;q -r 0 p;r q -p 0];
q0123dot = 0.5*PQRMAT*q0123;

%%%Gravity Model
%planet
r = state(1:3); %% r = [x;y;z]
rho = norm(r);
rhat = r/rho;
Fgrav = -(mu*m/rho^2)*rhat;

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
    BNED = [BN;BE;-BD]; %%PCI has Down as Up
    BI = TIB(phiE,thetaE+pi,psiE)*BNED;
    %BI = eye(3)*BNED;    
    BB = TIBquat(q0123)'*BI;
    %%%Convert to Tesla
    BB = BB*1e-9;
end


if t >= lastSensorUpdate
    %%%%SENSOR BLOCK
    lastSensorUpdate = lastSensorUpdate + nextSensorUpdate;
    [BfieldMeasured,pqrMeasured,ptpMeasured] = Sensor(BB,pqr,ptp); 
    
    %%%NAVIGATION BLOCK
    [BfieldNav,pqrNav,ptpNav] = Navigation(BfieldMeasured,pqrMeasured,ptpMeasured);   
end


%%%CONTROL BLOCK
[current,rwalphas] = Control(BfieldNav,pqrNav,ptpNav);

%%%Magtorquer Model
magtorquer_params
%%%Add in saturation filter
if sum(abs(current)) > maxCurrent/1000
   current = (current/sum(abs(current)))*maxCurrent/1000; 
end
muB = current*n*A;
LMN_magtorquers = cross(muB,BB);

%%%Reaction Wheels
w123dot = [0;0;0];
for idx = 1:3
    if abs(w123(idx)) > maxSpeed
        w123dot(idx) = 0;
    else
        if abs(rwalphas(idx)) > maxAlpha
            rwalphas(idx) = sign(rwalphas(idx))*maxAlpha;
        end
        w123dot(idx) = rwalphas(idx);
    end
end
LMN_RWs = Ir1B*w123dot(1)*n1 + Ir2B*w123dot(2)*n2 + Ir3B*w123dot(3)*n3;

%%%Compute Disturbance Forces and Moments
[XYZD,LMND] = Disturbance(rho-R,Amax,lmax,vel,CD,BB);

%%%TOTAL MOMENTS
LMN = LMN_magtorquers - LMN_RWs + LMND;

%%%Translational Dynamics
F = Fgrav + XYZD;
accel = F/m;

%%%Compute the total angular momentum
w1 = w123(1);
w2 = w123(2);
w3 = w123(3);
H = Is*pqr + Ir1B*w1*n1 + Ir2B*w2*n2 + Ir3B*w3*n3;

%%%Rotational Dynamics
pqrdot = invI*(LMN - cross(pqr,H));

%%%Return derivatives vector
dstatedt = [vel;accel;q0123dot;pqrdot;w123dot];

