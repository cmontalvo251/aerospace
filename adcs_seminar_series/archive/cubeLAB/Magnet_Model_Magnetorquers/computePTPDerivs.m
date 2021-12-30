function statedot = computePTPDerivs(state,LMN)

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

q0 = state(1);
q1 = state(2);
q2 = state(3);
q3 = state(4);
pqr = state(5:7);
p = pqr(1);
q = pqr(2);
r = pqr(3);

%%%%Pulled from Attitude Stabilization
%quatdot = 0.5*[-epsilon';(eta*eye(3) + etaskew)]*pqr; 
%%%Pulled from Boom 2010 - Both of these are the same.
quatdot = [-p*q1-q*q2-r*q3;p*q0+r*q2-q*q3;q*q0-r*q1+p*q3;r*q0+q*q1-p*q2]/2;

pqrskew = [0 -r q;r 0 -p;-q p 0];

pqrdot = Iinv*(LMN-pqrskew*I*pqr);

statedot = [quatdot;pqrdot];
