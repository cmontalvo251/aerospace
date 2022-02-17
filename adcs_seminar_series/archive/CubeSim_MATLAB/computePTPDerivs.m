function statedot = computePTPDerivs(state,LMNcontrol)
%%%DERIVATIVE OF SATELLITE
parameters

q0 = state(1);
q1 = state(2);
q2 = state(3);
q3 = state(4);
pqr = state(5:7);
p = pqr(1);
q = pqr(2);
r = pqr(3);

eta = q0;
epsilon = [q1;q2;q3];
etaskew = [0 -q3 q2;q3 0 -q1;-q2 q1 0];
%%%%Pulled from Attitude Stabilization
%quatdot = 0.5*[-epsilon';(eta*eye(3) + etaskew)]*pqr; 
%%%Pulled from Boom 2010 - Both of these are the same.
quatdot = [-p*q1-q*q2-r*q3;p*q0+r*q2-q*q3;q*q0-r*q1+p*q3;r*q0+q*q1-p*q2]/2;

pqrskew = [0 -r q;r 0 -p;-q p 0];

LMNdrag = [0;0;0];

LMN = LMNdrag + LMNcontrol;

pqrdot = Iinv*(LMN-pqrskew*I*pqr);

statedot = [quatdot;pqrdot];
