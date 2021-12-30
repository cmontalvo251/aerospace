function dxdt = Derivatives(t,x)

phi_sat = x(1);
p_sat = x(2);
phi_rw = x(3);
p_rw = x(4);

%%%Masses and Inertias (kg-m^2)
Irw_cg = 0.5; %%%Moment of inertia of the reaction wheel about the center of mass of the reaction wheel
mass_rw = 0.1; %%Mass of reaction wheel
d = 0.2; %%%Distance from center of mass of satellite to center of mass or reaction wheel
Irw_sat = Irw_cg + mass_rw*d^2;
I_sat = 1; %%%Moment of inertia of the cubesat with out reaction wheels
I_sys = I_sat + Irw_sat; %%%Total moment of inertia of system

%%%Properties of Reaction Wheel
TMAX = 5; %%Maximum torque on the motor at 1 Volt? (N-m)
p_rw_max = 5; %%rad/s %%%Maximum angular speed of rw at 1 Volt?
sig1 = -TMAX/p_rw_max;

%%%Voltage Computation
phi_sat_command = 0.7; %%%Commands in radians
p_sat_command = 0.0;  %%%Command in rad/s
kp = 0.5; %%%Proportional Gains?
kd = 0.2; %%%Derivative Gains?
Voltage = kp*(phi_sat-phi_sat_command) + kd*(p_sat-p_sat_command); %%Volt applied

%%%Torque of Motor
Torque_Motor = Voltage*(sig1*p_rw + TMAX);

%%%Damping Motor
c_rw = 0.0; %(N-m*(s/rad)) %%%Damping constant of reaction wheel?
Damping_Motor = -c_rw*p_rw;

%%%Angular Acceleration of Reaction Wheel
pdot_rw = (Torque_Motor-Damping_Motor)/Irw_cg;

%%%External Moments
cb = 0.0; %%%Bearing friction of platform (nd)
L = -cb*I_sys*abs(p_sat)*p_sat; %%%Bearing Friction N-m -> kg * m^2/s^2 -> kg-m^2 

%%%Angular Acceleration of Satellite
pdot_sat = (L-Irw_cg*pdot_rw)/(I_sys);

dxdt = [p_sat;pdot_sat;p_rw;pdot_rw];
