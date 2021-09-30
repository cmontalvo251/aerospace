function R = Rscrew(nhat)
%%%Assume that nhat is a 3x1 vector in R3 space. 
%%compute R such that v(inertial) = R v(body)
%%%There is an equation in Bachau's book on screw rotation but I'm at NASA
%%%right now and I can't find it. Instead I will compute phi theta and psi
%%%based on the vector nhat and run it through R123

x = nhat(1);
y = nhat(2);
z = nhat(3);

%%%yaw 
ps = atan2(y,x);
%%%Theta
theta = atan2(z,sqrt(x^2+y^2));
%%%Phi is always zero because it is a line
phi = 0;

R = R123(phi,theta,ps);