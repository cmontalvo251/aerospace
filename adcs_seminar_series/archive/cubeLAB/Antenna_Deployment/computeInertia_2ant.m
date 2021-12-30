function inertia = computeInertia_2ant(theta)
a = .2;
b = .1;
c = .1;
%%% Dimensions of the CubeSat

M = 1;

Ixx = (b^2 + c^2);
Iyy = (a^2 + c^2);
Izz = (a^2 + b^2);

I_matrix = [Ixx 0 0; 0 Iyy 0; 0 0 Izz];

I_body = (M/12) * I_matrix; 
%%%%%Inertia of the satellite body

l = .5;
m = .1;
r = .01;
%%%length, estimated mass and radius of the antenna. Radius of the "disc"
%%%end of the "rod" antenna. The antenna is treated as a rod, for which Ixx
%%%and Iyy are known (mass*length^2/3). However the Izz of the rod is
%%%treated as a disc. (mass*radius^2/2)

%%%%%%%%%%%%%First antenna 

I_ant1EA = [(m*(l^2)/3) 0 0; 0 (m*(l^2)/3) 0;0 0 (m*(r^2)/2)];
%%%%%Inertia of the first antenna about the end, in the antenna frame

T01 = [1 0 0;0 cos(theta) -sin(theta);0 sin(theta) cos(theta)];
%%%%%Transform matrix 

I_ant1ES = T01*I_ant1EA*T01';
%%%Inertia of the antenna about the end in the satellite frame

x = 0;
y = .05;
z = .1;

Re1 = [0 -z y; z 0 -x; -y x 0];
%%%% radius vector from satellite body origin to end of first antenna

I_ant1OS = I_ant1ES + m*Re1*(Re1');
%%% parallel axis theorem -- Inertia of antenna about the satellite
%%% origin, within the frame of the antenna

%%%%%%%%%%%%%Second antenna 

I_ant2EA = [(m*(l^2)/3) 0 0; 0 (m*(l^2)/3) 0;0 0 (m*(r^2)/2)];
%%% Inertia of the second antenna about its end, in the antenna frame

T02 = [1 0 0;0 cos(theta) -sin(theta);0 sin(theta) cos(theta)];
%%%Transform matrix

I_ant2ES = T02*I_ant2EA*T02';
%%%Inertia of the second antenna about the end in the satellite frame

x2 = 0;
y2 = -.05;
z2 = .1;

Re2 = [0 -z2 y2;z2 0 -x2; -y2 x2 0];
%%%Radius vector from sat body origin to end of antenna
I_ant2OS = I_ant2ES + m*Re2*Re2';

%%%Total inertia: satellite body plus antennas
inertia = I_body + I_ant1OS + I_ant2OS;