function boat_animation(xk,uk,tk)
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
y = xk(2);
psi = xk(3);
%%%Velocity states don't matter
%%%Control states do matter
dt = uk(1); %%0-100 percent
dr = uk(2); %%angle in radians

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
xfuse_I = (([0 18 18 0]-9)./12)./3.28;
yfuse_I = (([0 0 7.75 7.75]-7.75/2)./12)./3.28;

%%%%Rotate to body frame
%Save variables for speed
cpsi = cos(psi);
spsi = sin(psi);
xfuse_B = cpsi*xfuse_I - spsi*yfuse_I;
yfuse_B = spsi*xfuse_I + cpsi*yfuse_I;

%%%Plot and offset by x and y
patch(xfuse_B+x,yfuse_B+y,[165 42 42]./255)

%%%Now plot the rudder
%%%First compute the front and back
%%%of the rudder. Remember that the rudder can
%%%Rotate by dr
lr = (10/12)/3.28;
c = (2.25/12)/3.28;

boat_cg_I = [x;y];
r_cg_rudder_B = [-lr;0];
rFRONT_rudder_R = [c/2;0];
rREAR_rudder_R = [-2*c;0]; %%Exagerate for effect

TIB = [cpsi -spsi;spsi cpsi];
TBR = [cos(dr) sin(dr);-sin(dr) cos(dr)];

rFRONT_I = boat_cg_I + TIB*r_cg_rudder_B + TIB*TBR*rFRONT_rudder_R;
rREAR_I = boat_cg_I + TIB*r_cg_rudder_B + TIB*TBR*rREAR_rudder_R;
%%%Now using the rotation of the rudder and the cg
hold on
plot([rFRONT_I(1) rREAR_I(1)],[rFRONT_I(2) rREAR_I(2)],'b-','LineWidth',3)

%%%Plot Thrust vector so we know how much thrust
%%%we are getting
thrust_vec = boat_cg_I + TIB*[((18/12)/3.28)*dt/100 0]';
plot([boat_cg_I(1) thrust_vec(1)],[boat_cg_I(2) thrust_vec(2)],'r-')

%%%Make a title
title(tk)

%%%Make axis square
world_lim = 2;
xlim([-world_lim+x world_lim+x])
ylim([-world_lim+y world_lim+y])
axis equal

%%%%Plot everything now
drawnow
