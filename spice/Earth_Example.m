%%% Construct a meta kernel, "standard.tm”, which will be used to load the needed
%%% generic kernels: "naif0009.tls," "de421.bsp,” and "pck0009.tpc.”
%%% Load the generic kernels using the meta kernel, and a Cassini spk.
clear
clc
close all
%%%Add the path of the mice toolkit
addpath('mice/src/mice')
addpath('mice/lib/')
%%LOCATION OF KERNELS https://naif.jpl.nasa.gov/pub/naif/
%%%DE Planets
%%%https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440_tech-comments.txt
%%%https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440_and_de441.pdf
%%%https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/
%%%https://naif.jpl.nasa.gov/pub/naif/generic_kernels/fk/planets/earth_assoc_itrf93.tf
root = '/home/files/Docs/Git_Repos/Github/aerospace/spice/kernels/';

cspice_furnsh({[root,'gen/lsk/naif0012.tls'],[root,'earth/de440.bsp']})%,[root,'earth/earth_000101_220118_211025.bpc']})

%%%% Define the number of divisions of the time interval and the time interval.
STEP = 10000;
et = cspice_str2et( {'Sep 2, 1987', 'Oct 27, 2021'} );
times = (0:STEP-1) * ( et(2) - et(1) )/STEP + et(1);
AU = 1.496e+8; %%Km

%%%GET EARTH
[pos, ltime]= cspice_spkpos( 'EARTH BARYCENTER', times, 'J2000', 'NONE', 'SUN' );
%Plot the resulting trajectory.
xE = pos(1,:)/AU;
yE = pos(2,:)/AU;
zE = pos(3,:)/AU;
plot3(xE,yE,zE)
hold on
plot3(xE(1),yE(1),zE(1),'b*')
%%%GET SATURN
[pos, ltime]= cspice_spkpos( 'SATURN BARYCENTER', times, 'J2000', 'NONE', 'SUN' );
%Plot the resulting trajectory.
xS = pos(1,:)/AU;
yS = pos(2,:)/AU;
zS = pos(3,:)/AU;
plot3(xS,yS,zS)
plot3(xS(1),yS(1),zS(1),'y*')

xlabel('X')
ylabel('Y')
zlabel('Z')



grid on
cspice_kclear