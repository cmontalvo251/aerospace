%%%%BREAK SYSTEM INTO TWO SEPARATE SIMS%%%%%
SIMXYZ = 0;
SIMPTP = 1;

%%%%%%%
CPPDATA = 0; %%%1 = Plot c++ simulation results. Use this to debug C++ code

%%%%%%%%Hardware in the Loop%%%%%%
%%%0 = Simulate everything in MATLAB
%%%1 = Only simulate satellite and sensors in MATLAB (Don't
%simulate Kalman Filter)
%%%2 = Only simulate satellite in MATLAB (Don't simulate Kalman
%filter or sensors)
HIL = 2; 

%%%%%%%Connection Type%%%%%%%%%%
%%%0 = Communicate with ADCS on Computer
%%%1 = Communicate with ADCS over serial
SERIALTYPE = 0; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%Time step of simulation
%%%%FOR ROTATIONAL DYNAMICS
if (SIMPTP)
  timestep = 0.01; %%seconds (1000 Hz = 0.001, 100 Hz = 0.01)
  tskip = 1;
  tfinal = 10; %%%seconds
else
  %%%%FOR ORBIT DYNAMICS
  timestep = 0.01*3600;
  tskip = 1;
  tfinal = 1.5*3600; %1.5 hrs
end

%%%%%%%%%%SENSOR INFORMATION%%%%%%%%%%%%%%%%

gpsUpdateRate = 10; %%Hz
controlUpdateRate = 100; %%Hz
noiseTypeGPS = 2;

%%%%Mass and Inertia of Satellite
mass = 2.6; %kg
a_x = 10/100; %meters
b_y = 10/100; %meters
c_z = 20/100; %meters

Ixx = (1/12)*mass*(b_y^2+c_z^2);
Iyy = (1/12)*mass*(a_x^2+c_z^2);
Izz = (1/12)*mass*(a_x^2+b_y^2);

I_vec = [Ixx,Iyy,Izz];
I = diag(I_vec);
Iinv = diag(1./I_vec);

%%%NON-DIMENSIONAL PARAMETERS FOR ORBIT
L = 42164*(1000); %%%meters Radius of earth in meters
T = 86400; %%%The length of 1 LEO orbit - seconds (86400 seconds = 1 day)
M = 5.972*10^24 ; %%%Mass of the earth in kg
G = 6.6742*10^-11; %%Gravitational constant - units?
muSat = -G*(mass+M);

%%%Rotational Control Parameters
Kp = diag([1.0,1.0,1.0]);
Kd = diag([0.5,0.5,0.5]);
Km = 0.01*diag([1.0,1.0,0.4]);
%Km = 0;

