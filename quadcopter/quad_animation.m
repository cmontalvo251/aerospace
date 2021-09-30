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
phi = xk(3);
%%%Velocity states don't matter
%%%Control states do matter
r1 = uk(1); %%0-100
r2 = uk(2); %%

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
L = 1/3.28;
H = (2/12)/3.28;
xfuse_I = [-L/2 L/2 L/2 -L/2];
zfuse_I = [-H/2 -H/2 H/2 H/2];

%%%%Rotate to body frame
%Save variables for speed
cphi = cos(phi);
sphi = sin(phi);
xfuse_B = cphi*xfuse_I - sphi*zfuse_I;
zfuse_B = sphi*xfuse_I + cphi*zfuse_I;

%%%Plot and offset by x and y
patch(xfuse_B+x,zfuse_B+z,[165 42 42]./255)

quad_cg_I = [x;z];
r_cg_left_B = [-L/2;0];
r_cg_right_B = [L/2;0];
r_cg_left_LENGTH_B = [0;-r1]/100; %%Normalize by 100
r_cg_right_LENGTH_B = [0;-r2]/100;

TIB = [cphi -sphi;sphi cphi];

rLEFT_I = quad_cg_I + TIB*r_cg_left_B;
rRIGHT_I = quad_cg_I + TIB*r_cg_right_B;
rLEFT_TOP_I = rLEFT_I + TIB*r_cg_left_LENGTH_B;
rRIGHT_TOP_I = rRIGHT_I + TIB*r_cg_right_LENGTH_B;
%%%Now plot both thrust vectors
hold on
plot([rLEFT_I(1) rLEFT_TOP_I(1)],[rLEFT_I(2) rLEFT_TOP_I(2)],'b-','LineWidth',2)
plot([rRIGHT_I(1) rRIGHT_TOP_I(1)],[rRIGHT_I(2) rRIGHT_TOP_I(2)],'r-','LineWidth',2)

%%%Make a title
title(tk)

%%%Make axis square
world_lim = 3*L;
xlim([-world_lim+x world_lim+x])
ylim([-world_lim+z world_lim+z])
axis equal

%%%%Plot everything now
drawnow
