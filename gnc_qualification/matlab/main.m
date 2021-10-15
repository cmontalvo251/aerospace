clear
clc
close all

disp('GNC QUALIFICATION')

%%USer defined parameters
apogee = 600; %%km above surface
perigee = 600; %km above surface
N = 1000; %%%precision of orbit
magnetic_moment = 0.6;  %%number of turns times current time area of magnetorquers
L = 36.6; %%%length of satellite in centimeters
W  = 22.6; %%%width of satellite in centirmeter
D = 22.6; %%depth of satellite in centimeters
CD = 1.0; %%drag coefficient of satellite
m = 25.0; %%%mass of satellite in kg
w0 = 10.0; %%%initial angular velocity of sat in deg/s
SF = 2.0; %%%safety factor of detumbling
HxRW = 0.1; %%%Size of reaction wheel in Nms
HyRW = 0.1; %%%Size of reaction wheel in Nms
HzRW = 0.1; %%%Size of reaction wheel in Nms
Power = 6.0; %%%Total wattage of reaction wheel during max torque
MaxT = 0.025; %%%Maximum Torque of RW in N-m
mission_duration = 12; %%%mission duration in months
maneuver_angle = 180; %%%manuever_angle in degrees
f = 0.5; %%%a factor from 0(non-inclusive) to 0.5 which dictate the speed of the manuever (0 is not moving and 0.5 is as fast as possible)
constants %%Other constants

%%%% Inertia Calculator
[Ixx,Iyy,Izz,max_moment_arm,max_area] = inertia(L,W,D,m);

%%%Run the orbit model
[x,y,z,t,r,T_orbit,vx,vy,vz,v] = orbit_model(apogee,perigee,N);
disp(['Orbit Time = ',num2str(T_orbit)])

%%%Call the magnetic field model
addpath('../../igrf') %%%Hey you need to make sure you download igrf from mathworks
%%https://www.mathworks.com/matlabcentral/fileexchange/34388-international-geomagnetic-reference-field-igrf-model
[BxI,ByI,BzI,B500,Btotal] = magnetic_field(x,y,z,r,t);
disp(['Magnetic Field Strength @ 500 km (nT) = ',num2str(B500)])

%%%%Initial Tumbling
[Hx0,Hy0,Hz0] = initial_angular_momentum(Ixx,Iyy,Izz,w0,SF);
disp(['Initial Angular Momentum (N-m-s) ',num2str(Hx0),' ',num2str(Hy0),' ',num2str(Hz0)])

%%%DISTURBANCE TORQUES
Maero = aerodynamics(t,r,v,max_area,max_moment_arm,rhosl,CD);
Mgrav = gravity(r,max_moment_arm,m);
Mrad = solar_radiation(r,max_area,max_moment_arm);
Mdipole = dipole(Mrad,Btotal,B500);

%%%%Compute Total Disturbance Per Orbit
Hdist = disturbance(t,Maero,Mgrav,Mrad,Mdipole);
disp(['Total Disturbance Momentum Per Orbit (N-m-s) = ',num2str(Hdist)])

%%%%Compute Magnetorquer Effectiveness 
magnetorquers(t,magnetic_moment,Btotal,Hdist,Hx0,Hy0,Hz0,HxRW,HyRW,HzRW);

%%%Compute Reaction Wheel Effectiveness
reaction_wheels(t,Hx0,Hy0,Hz0,HxRW,HyRW,HzRW,Hdist,mission_duration);

%%%Maneuver time with reaction wheels
maneuver_time(Power,MaxT,maneuver_angle,f,Ixx,Iyy,Izz);