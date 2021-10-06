clear
clc
close all

disp('GNC QUALIFICATION')

%%USer defined parameters
apogee = 600; %%km above surface
perigee = 600; %km above surface
N = 1000; %%%precision of orbit
magnetic_moment = 1.176;  %%number of turns times current time area of magnetorquers
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
[Ixx,Iyy,Izz,max_moment_arm] = inertia(L,W,D,m);

%%%Run the orbit model
[x,y,z,t,T_orbit] = orbit_model(apogee,perigee,N,0);
disp(['Orbit Time = ',num2str(T_orbit)])

%%%Plot the orbit
plot3(x,y,z)
set(gcf,'color','white')
xlabel('X (m)')
ylabel('Y (m)')
zlabel('Z (m)')
grid on
title('Orbit')
axis equal
%%%Plot Altitude
altitude = sqrt(x.^2 + y.^2 + z.^2);
figure()
set(gcf,'color','white')
plot(t,altitude)
xlabel('Time (sec)')
ylabel('Altitude (m)')
grid on

%%%Call the magnetic field model
addpath('../../igrf') %%%Hey you need to make sure you download igrf from mathworks
%%https://www.mathworks.com/matlabcentral/fileexchange/34388-international-geomagnetic-reference-field-igrf-model
[BxI,ByI,BzI] = magnetic_field(x,y,z);

%%%Plot the magnetic field
figure()
set(gcf,'color','white')
plot(t,BxI,'b-')
hold on
plot(t,ByI,'r-')
plot(t,BzI,'g-')
legend('X','Y','Z')
xlabel('Time (sec)')
ylabel('Magnetic Field (nano Tesla)')
grid on