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
HxRW = 200; %%%Size of reaction wheel in mNms
HyRW = 100; %%%Size of reaction wheel in mNms
HzRW = 100; %%%Size of reaction wheel in mNms
Power = 6.0; %%%Total wattage of reaction wheel during max torque
MaxT = 0.025; %%%Maximum Torque of RW in N-m
mission_duration = 12; %%%mission duration in months
maneuver_angle = 180; %%%manuever_angle in degrees
f = 0.5; %%%a factor from 0(non-inclusive) to 0.5 which dictate the speed of the manuever (0 is not moving and 0.5 is as fast as possible)

%%%% Inertia Calculator
[Ixx,Iyy,Izz,max_moment_arm,max_area] = inertia(L,W,D,m);

%%%Run the orbit model
[x,y,z,t,r,T_orbit,vx,vy,vz,v] = orbit_model(apogee,perigee,N,0);
disp(['Orbit Time = ',num2str(T_orbit)])

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
constants
plot(t,(r-REarth)/1000,'LineWidth',2)
xlabel('Time (sec)')
ylabel('Altitude (km)')
grid on

%%%Call the magnetic field model
addpath('../../igrf') %%%Hey you need to make sure you download igrf from mathworks
%%https://www.mathworks.com/matlabcentral/fileexchange/34388-international-geomagnetic-reference-field-igrf-model
[BxI,ByI,BzI,B500,B] = magnetic_field(x,y,z,r);
disp(['Magnetic Field Strength @ 500 km (nT) = ',num2str(B500)])

%%%Plot the magnetic field
figure()
set(gcf,'color','white')
plot(t,BxI,'LineWidth',2)
hold on
plot(t,ByI,'LineWidth',2)
plot(t,BzI,'LineWidth',2)
legend('X','Y','Z')
xlabel('Time (sec)')
ylabel('Magnetic Field (nano Tesla)')
grid on

%%%%Initial Tumbling
[Hx0,Hy0,Hz0] = initial_angular_momentum(Ixx,Iyy,Izz,w0,SF);
disp(['Initial Angular Momentum (N-m-s) ',num2str(Hx0),' ',num2str(Hy0),' ',num2str(Hz0)])

%%%DISTURBANCE TORQUES
Maero = aerodynamics(t,r,v,max_area,max_moment_arm,rhosl,CD);
Mgrav = gravity(r,max_moment_arm,m);
Mrad = solar_radiation(r,max_area,max_moment_arm);
Mdipole = dipole(Mrad,B,B500);

%%%%Compute Total Disturbance Per Orbit
Hdist = disturbance(t,Maero,Mgrav,Mrad,Mdipole);
disp(['Total Disturbance Momentum Per Orbit (N-m-s) = ',num2str(Hdist)])
