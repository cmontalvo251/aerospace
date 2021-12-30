function LMNcontrol = computeControl(state_PTP_tilde)
parameters

q0123 = state_PTP_tilde(1:4);

pqr = state_PTP_tilde(5:7); %%rad/s

ptp = quat2euler(q0123')';

Kpptp = Kp*-ptp;
Kdpqr = Kd*-pqr;
pqrc =  Kpptp + Kdpqr;

pqrc_pqr = (pqrc-pqr);

LMNcontrol = Km*pqrc_pqr;






