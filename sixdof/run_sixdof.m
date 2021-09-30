%%This function uses a 6dof fixed timestep RK4 integrator to
%propagate the state forward. The only thing the user must do is
%to include the script to add in forces and moments given to the
%body
%%The states are standard AE format x = [x y z phi theta psi u v w
%p q r];
%The form of the forces and moment function is as follows
%[XYZ,LMN] = Fext_Mext(time,state,inertia)
%state is a 12x1, and inertia is a cell array with inertia{1} =
%mass and inertia{2} = I(moment of inertia matrix)
%LMN and XYZ are column vectors of body moments and forces respectively
%%%An example script is in run_sixdof.m and an example force and moment model for a ball is in 
%%%aero_ball.m

clear
clc
close all

Fext_Mext = 'aero_ball';
tspan = [0 10];
xinitial = [0 0 0 0 0 0 0 0 -10 0 0 0]';
timestep = 0.01;
next = 10;
mass = 1;
I = eye(3);
inertia = {mass,I};

[tout,xout] = sixdof(Fext_Mext,tspan,xinitial,timestep,inertia,next);

ylabels = {'x','y','z','\phi','\theta','\psi','u','v','w','p','q','r'};

for idx = 1:12
    figure()
    plot(tout,xout(:,idx))
    xlabel('Time (sec)')
    grid on
    ylabel(ylabels{idx})
end