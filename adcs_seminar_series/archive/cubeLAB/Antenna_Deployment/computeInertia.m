function inertia = computeInertia(theta)
%%%One antenna
a = .2; %h
b = .1; %l
c = .1; %w
%%% Dimensions of the CubeSat

M = 1;

Ixx = (b^2 + c^2);
Iyy = (a^2 + c^2);
Izz = (a^2 + b^2);

I_matrix = [Ixx 0 0; 0 Iyy 0; 0 0 Izz];

I_body = (M/12) * I_matrix; 
%%%%%Inertia of the satellite body

l = .5; %%%%dummy number = 5
m = .1; %%% 0.5
r = .01; %%% 0.5
%%%length, estimated mass and radius of the antenna. Radius of the "disc"
%%%end of the "rod" antenna. The antenna is treated as a rod, for which Ixx
%%%and Iyy are known (mass*length^2/3). However the Izz of the rod is
%%%treated as a disc. (mass*radius^2/2)

I_antEA = [(m*(l^2)/3) 0 0; 0 (m*(l^2)/3) 0;0 0 (m*(r^2)/2)];
%%%%%Inertia of the antenna about the end, in the antenna frame

T0 = [1 0 0;0 cos(theta) -sin(theta);0 sin(theta) cos(theta)];
%%%%%Transform matrix 

I_antES = T0*I_antEA*T0';
%%%Inertia of the antenna about the end in the satellite frame

x = 0;
y = .05;
z = .1;

Re = [0 -z y; z 0 -x; -y x 0];
%%%% radius vector from satellite body origin to end of antenna

I_antOS = I_antES + m*Re*(Re');
%%% parallel axis theorem -- Inertia of antenna about the satellite
%%% origin, within the frame of the antenna


%%%Total inertia: satellitel body plus antenna
inertia = I_body + I_antOS;