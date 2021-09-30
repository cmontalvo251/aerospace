function quad_animation(xk,uk,tk)
%%%In order to visualize properly
%%%We wait until the simulation
%%%Catches up to the real time
%%%Clock
real_time = toc;
%%%If more time has passed than the simulation
%%%Our timestep is too small. Our computer
%%%Can't keep up.
if real_time > tk
    disp('Timestep too small')
end
%%%If our computer is too fast we wait until
%%%The real world has caught up with our sim
while real_time < tk
    real_time = toc;
end

%%%%Extract state
x = xk(1);
z = xk(2);
theta = -xk(3); %%%Flipping axis
%%%Velocity states don't matter
%%%Control states do matter
dt = uk(1); %%0-100
de = uk(2); %%

%Debugging
%x = 10;
%y = 5;
%psi = 45*pi/180;
%dt = 100;
%dr = 45*pi/180;

%%%%First clear figure
cla;

%%%%Now plot a rectangle for fuselage
%%%%must use patch because of rotation
L = 2.04;
Lv = 0.7;
xfuse_I = [-3*L/4 L/4 L/4   0     -L/4 -L/2 -3*L/4];
zfuse_I = -[0       0  Lv/4  Lv/2  Lv/4  Lv/4   Lv];

%%%%Rotate to body frame
%Save variables for speed
ctheta  = cos(theta);
stheta = sin(theta);
xfuse_B = ctheta*xfuse_I - stheta*zfuse_I;
zfuse_B = stheta*xfuse_I + ctheta*zfuse_I;

%%%Plot and offset by x and y
patch(xfuse_B+x,zfuse_B+z,[20 100 42]./255)

aircraft_cg_I = [x;z];
r_cg_elevator_B = [-3*L/4;0];
r_cg_prop_B = [L/4;0];
r_cg_prop_LENGTH_B = [L/4+dt/100;0];
r_elevator_E = [-L/2;0]; %%exxagerate for effect

TIB = [ctheta -stheta;stheta ctheta];
TBE = [cos(de) sin(de);-sin(de) cos(de)];

r_prop_START_I = aircraft_cg_I + TIB*r_cg_prop_B;
r_prop_END_I = aircraft_cg_I + TIB*r_cg_prop_LENGTH_B;
r_elevator_START_I = aircraft_cg_I + TIB*r_cg_elevator_B;
r_elevator_END_I = r_elevator_START_I + TIB*TBE*r_elevator_E;
%%%Now plot thrust and elevator
hold on
plot([r_prop_START_I(1) r_prop_END_I(1)],[r_prop_START_I(2) r_prop_END_I(2)],'b-','LineWidth',2)
plot([r_elevator_START_I(1) r_elevator_END_I(1)],[r_elevator_START_I(2) r_elevator_END_I(2)],'r-','LineWidth',2)

%%%Make a title
title(tk)

%%%Make axis square
world_lim = 3*L;
xlim([-world_lim+x world_lim+x])
ylim([-world_lim+z world_lim+z])
axis equal

%%%%Plot everything now
drawnow
