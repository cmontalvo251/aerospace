function [LMNControl,xyzcurrent] = ComputeControl(state,time)

%%%Extract State Vector
q0123 = state(1:4);
pqr = state(5:7);
p = pqr(1);
q = pqr(2);
r = pqr(3);

%%%%Magnetorquer parameters
n = 84; %%%number of turns
A = 0.02; %%%area(m^2)
max_current = 0.04; %%%Amps

%%%First thing we need is the magnetic field strength
EarthB_Inertial = [-1.001e4;-2053;2.066e4]*(1e-9); %%Teslas

%%%Then we need to rotate the magnetic field into the body frame
EarthB_Body = RQUAT(q0123)'*EarthB_Inertial;

%%%%%%%%%%%%%%%%CONTROLLER DERIVATION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%Then we need to run our controllers and determine the current in x,y,z
CONTROLLER_TYPE = 3;
if CONTROLLER_TYPE == 1
    xcurrent = max_current/3;
    ycurrent = max_current/3;
    zcurrent = max_current/3;
    %disp('Max Current Applied')
elseif CONTROLLER_TYPE == 2
    %disp('RREF Controller')
    k = 100;
    Mx = -k*p;
    Mz = -k*r;
    Mxbar = Mx/(n*A);
    Mzbar = Mz/(n*A);
    bx = EarthB_Body(1);by = EarthB_Body(2);bz = EarthB_Body(3);
    Mybar = -Mxbar*bx/by - Mzbar*bz/by;
    zcurrent = 0;
    ycurrent = (Mxbar+by*zcurrent)/bz;
    xcurrent = (Mzbar+bx*ycurrent)/by;
elseif CONTROLLER_TYPE == 3
    %disp('Magnetic Field Cross Omega')
    k = 1e5;
    pqr = [p;q;r];
    muIDEAL = -k*(cross(EarthB_Body,pqr));
    xyzcurrentIDEAL = muIDEAL/(n*A);
    xcurrent = xyzcurrentIDEAL(1);
    ycurrent = xyzcurrentIDEAL(2);
    zcurrent = xyzcurrentIDEAL(3);
else
    xcurrent = 0;
    ycurrent = 0;
    zcurrent = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%Saturation Block
mag = abs(xcurrent)+abs(ycurrent)+abs(zcurrent);
if mag > max_current
    xcurrent = xcurrent*max_current/mag;
    ycurrent = ycurrent*max_current/mag;
    zcurrent = zcurrent*max_current/mag;
end

%%%Once we have the current we can create magnetic moment vector in the
%%%body frame
xyzcurrent = [xcurrent;ycurrent;zcurrent];
muXYZ = n*A*xyzcurrent;
LMNControl = cross(muXYZ,EarthB_Body);
