function LMN = thrusters(state,time)
global I q0
%%%This routine will fire thrusters to stabilize a satellite

%%%%Extract states
roll_rate = state(4);
pitch_rate = state(5);
yaw_rate = state(6);

%%%If we had 3 axis control we could just use -k*w and it would stabilize
kp = 2;
L = -kp*roll_rate;
M = -kp*pitch_rate;
N = -kp*yaw_rate;

%%%If we only had two axis control we would not be stable
L = 0;
M = -kp*pitch_rate;
N = -kp*yaw_rate;

%%%According to my Lyapunov Derivation though if we use the following
%%%We are stable
L = 0;
Mbar = -3*roll_rate*yaw_rate;
Nbar = -kp*yaw_rate;
Ixx = I(1,1);
Iyy = I(2,2);
Izz = I(3,3);
I1 = (Izz-Ixx)/Iyy;
I2 = (Ixx-Iyy)/Izz;
Mhat = Mbar*I1;
Nhat = Nbar*I2;
M = Mhat*Iyy;
N = Nhat*Izz;
%%%But if you investigate the zero dynamics this system is not
%%%stable so we need to phase this controller

I0 = (Iyy-Izz)/Ixx;

%%%Phase 1
pdotc = 0;
qdotc = 0;
rdotc = 0;
pc = 0;
qc = q0;
gam = pdotc - kp*(roll_rate-pc);
rc = gam/(I0*pitch_rate);

%%%Just program Phase 2 with dead reckoning
if (time > 10)
    qc = 0;
    rc = 0;
end

lam = qdotc - kp*(pitch_rate-qc);
bet = rdotc - kp*(yaw_rate-rc);
Mhat = lam - I1*roll_rate*yaw_rate;
Nhat = bet - I2*roll_rate*pitch_rate;
M = Mhat*Iyy;
N = Nhat*Izz;



LMN = [L;M;N];