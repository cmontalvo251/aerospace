function [Hx0,Hy0,Hz0] = initial_angular_momentum(Inertia,w0,SF)
Ixx = Inertia(1,1);
Iyy = Inertia(2,2);
Izz = Inertia(3,3);
%%Convert w0 to rad/s
wrad = w0*pi/180.0;
%%%Compute Angular Momentum
Hx0 = Ixx*wrad*SF;
Hy0 = Iyy*wrad*SF;
Hz0 = Izz*wrad*SF;