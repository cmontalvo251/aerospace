%%% Construct a meta kernel, "standard.tm”, which will be used to load the needed
%%% generic kernels: "naif0009.tls," "de421.bsp,” and "pck0009.tpc.”
%%% Load the generic kernels using the meta kernel, and a Cassini spk.

%%%Add the path of the mice toolkit
addpath('mice/src/mice')
addpath('mice/lib/')
%%LOCATION OF KERNELS https://naif.jpl.nasa.gov/pub/naif/
root = '/home/files/Docs/Git_Repos/Github/aerospace/spice/kernels/';

cspice_furnsh({[root,'gen/lsk/naif0012.tls'],[root,'cassini/spk/030201AP_SK_SM546_T45.bsp']})

%%%% Define the number of divisions of the time interval and the time interval.
STEP = 1000;
et = cspice_str2et( {'Jun 20, 2004', 'Dec 1, 2005'} );
times = (0:STEP-1) * ( et(2) - et(1) )/STEP + et(1);

[pos, ltime]= cspice_spkpos( 'Cassini', times, 'J2000', 'NONE', 'SATURN BARYCENTER' );

%Plot the resulting trajectory.

AU = 1.496e+11;

x = pos(1,:)/AU;
y = pos(2,:)/AU;
z = pos(3,:)/AU;
plot3(x,y,z)

cspice_kclear