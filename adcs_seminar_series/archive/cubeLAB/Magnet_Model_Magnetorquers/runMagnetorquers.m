purge

%%%Ok so here we are going to test all of the magnetorquer controllers

%%%Initial State of the satellite
phi = 0;
theta = 45*pi/180;
psi = 0;
p = 0.1;
q = 0;
r = 0;

q0123 = euler2quat([phi,theta,psi])';

state = [q0123;p;q;r];
timestep = 0.1;
tout = 1:timestep:3000;
stateout = zeros(length(tout),length(state));
current_out = zeros(length(tout),3);
for idx = 1:length(tout)
    [nextstate,xyzcurrent] = computeRK4Step(@computePTPDerivs,state,timestep);
    %%%Normalize quaternions
    q0 = nextstate(1);
    q1 = nextstate(2);
    q2 = nextstate(3);
    q3 = nextstate(4);
    normq = sqrt(q0^2+q1^2+q2^2+q3^2);
    nexstate(1) = nextstate(1)/normq;
    nexstate(2) = nextstate(2)/normq;
    nexstate(3) = nextstate(3)/normq;
    nexstate(4) = nextstate(4)/normq;
    %%%Propagate step
    stateout(idx,:) = state;
    state = nextstate;
    %%%Save currents
    current_out(idx,:) = xyzcurrent';
end    

q0123_vec = stateout(:,1:4);
ptp_vec = quat2euler(q0123_vec);
pqr_vec = stateout(:,5:7);

plottool(1,'Euler Angles',18,'Time(sec)','Angle (deg)');
plot(tout,ptp_vec*180/pi)
legend('\phi','\theta','\psi')

plottool(1,'Angular Velocity',18,'Time (sec)','Angular Velocity (rad/s)');
plot(tout,pqr_vec)
legend('P','Q','R')

plottool(1,'Currents',18,'Time (sec)','Current (Amps)');
plot(tout,current_out)
legend('X','Y','Z')


