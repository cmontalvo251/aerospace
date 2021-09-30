function [dwdt,d] = derivs(t,w,td)
%td = 15;
theta = thetafunct(t,td);
I = computeInertia_2ant(theta);

thetadot = thetadotfunct(t,td);
I_dot = computeIdot(theta)*thetadot;

dwdt = -inv(I)*(I_dot*w + cross(w,I*w));

d = td;