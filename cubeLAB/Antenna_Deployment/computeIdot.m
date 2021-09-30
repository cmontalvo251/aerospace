function Idot = computeIdot(theta)

l = .5;
m = .1;
r = .01;

A = (m*l^2)/3;
B = (m*r^2)/2;

Idot = [0 0 0;
    0 (2*B-2*A)*sin(theta)*cos(theta) (A-B)*((cos(theta)^2)-(sin(theta)^2));
    0 (A-B)*((cos(theta)^2)-(sin(theta)^2)) ((2*A-2*B)*sin(theta)*cos(theta))];

