function [Hx0,Hy0,Hz0] = initial_angular_momentum(Ixx,Iyy,Izz,w0,SF)
%%Convert w0 to rad/s
wrad = w0*pi/180.0;
%%%Compute Angular Momentum
Hx0 = Ixx*wrad*SF;
Hy0 = Iyy*wrad*SF;
Hz0 = Izz*wrad*SF;