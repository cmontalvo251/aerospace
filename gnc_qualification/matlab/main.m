clear
clc
close all

disp('GNC QUALIFICATION')

%%USer defined parameters
apogee = 600; %%km above surface
perigee = 600; %km above surface
N = 1000;

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
addpath('../igrf') %%%Hey you need to make sure you download igrf from mathworks
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