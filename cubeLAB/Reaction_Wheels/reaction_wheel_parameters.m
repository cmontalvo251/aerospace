purge
%%%%What parameters do we need?

format long g

%%%%%%%%%%%%%INPUTS%%%%%%%%%%%%%%%%%%

%%%Position of reaction wheel from cg
rwx = 0.0;
rwy = 1.0;
rwz = 0.0;

%%%Maximum RPM
rpm_max = 7800;

%%%%Total Amount of inertia stored
inertia = 100; %%%mN-s

%%%Mass of reaction wheel
mass = 0.35;

%%%Unit vector in body frame
nhat = [1;0;0];

%%%Max torque
maxT = 0.007; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%COMPUTATIONS%%%%%%%%%%%%%%

rw_max = rpm_max*2*pi/60; %%Convert to rad/s
inertia_store = inertia/1000; %%%Convert to N-s
rotational_inertia = inertia_store/rw_max; %%%Compute rotational inertia
max_accel = maxT/rotational_inertia;
radius = sqrt(rotational_inertia*2/mass); %%Compute radius of reaction wheel
height = radius/1.2;  %%%Estimate height of disk

%%%Compute inertia of reaction wheel at cg of reaction wheel
Ixxcg = mass*radius^2/2;
Iyycg = (mass/12)*(3*radius^2 + height^2);
Izzcg = Iyycg;
Icg = [Ixxcg 0 0;0 Iyycg 0;0 0 Izzcg];

%%%Rotate this to the body frame
nhat = nhat./norm(nhat);
R = Rscrew(nhat);
I1 = R*Icg*R';

%%%Compute inertia of reaction wheel at cg of satellite
rRW = skew([rwx,rwy,rwz]);
IRW = I1 + mass*(rRW)*(rRW')

%%%At this point you have IRW of each reaction wheel
%%%you can then add IRW to the satellite to get Itotal.

