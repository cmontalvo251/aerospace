function [XYZD,LMND] = Disturbance(altitude,Amax,lmax,vel,CD,BI_Tesla)
%%%Compute the disturbance forces and moments on a Satellite

%%%Aerodynamic Drag
V = norm(vel);
d = density(altitude);
Fdrag = 0.5*d*V^2*Amax*CD;
vhat = vel / V;
XYZAERO = -Fdrag*vhat;
Maero = Fdrag*lmax/2;
LMNAERO = [Maero;Maero;Maero];

%%%Solar Radiation Pressure
solar_pressure = 4.5e-6; %%Pa
Fpressure = solar_pressure*Amax;
shat = vhat; %%%We need to code the direction of the sun (In a later video)
XYZSRP = -Fpressure*shat;
Mpressure = Fpressure*lmax/2;
LMNSRP = [Mpressure;Mpressure;Mpressure];

%%%Magnetic Dipole Moment
dconstant = 2.64e-3;
LMNMDM = dconstant*BI_Tesla;

%%%Add it all up
XYZD = XYZSRP + XYZAERO;
LMND = LMNSRP + LMNAERO + LMNMDM;