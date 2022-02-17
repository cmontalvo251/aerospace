%%%Simulate rotational Dynamics of a CubeSAT
%%%3/27/2020 - Commit 1 - Inertia diagonal - no controller - assuming perfect
%%%state feedback

clear
clc
close all

%%%SIMULATE FOR 0 to TMAX seconds
TMAX = 20;
tspan = [0 TMAX];
timestep = 0.001;

%%%Inertia
global I invI q0
%%%Principal axes of inertia
Ixx = 1;
Iyy = 2;
Izz = 3;
%%%Off Axes
Ixy = 0;
Iyz = 0;
Ixz = 0;
I = [Ixx Ixy Ixz;Ixy Iyy Iyz;Ixz Iyz Izz]; %%This doesn't have the cross products of inertia
invI = inv(I);

%%%INITIAL CONDITIONS
%phi,theta,psi are the angles
%p,q,r are the angular rates
%xinitial = [phi0,theta0,psi0,p0,q0,r0]
p0 = 1;
q0 = 1;
r0 = 1;
xinitial = [0;0;0;p0;q0;r0];

%%%%SETUP RK4
next = 10;
[N,flag] = size(xinitial);
tout = tspan(1):timestep:tspan(end);
integrationsteps = length(tout);
xout = zeros(N,integrationsteps);
x = xinitial;
threshold = 0;

%%%RK4 SIMULATION
for ii = 1:length(tout)
    
  %%%NOTIFY USER OF PROGRESS
  time = tout(ii);
  if next
    percent = (floor(100*time/tout(end)));
    if percent >= threshold
      disp(['Simulation ',num2str(percent),'% Complete'])
      threshold = threshold + next;
    end
  end
  
  %%Save States
  xout(:,ii) = x;
  
  %%Integrate 4 times like RK4
  xdot1 = dynamics(time, x);
  xdot2 = dynamics(time + (.5*timestep), x + (xdot1*.5*timestep));
  xdot3 = dynamics(time + (.5*timestep), x + (xdot2*.5*timestep));
  xdot4 = dynamics(time + timestep, x + (xdot3*timestep));
  xdotRK4 = (1/6) * (xdot1 + (2*xdot2) + (2*xdot3) + xdot4);
  nextstate = x + (timestep * xdotRK4);
  x = nextstate;
  
end

%%%GET READY FOR PLOTTING
xout = xout'; %%This changes everything to a more file friendly format
tout = tout';

%%%PLOT EVERYTHING
BodyName = 'Satellite';
dim1 = 'deg';
Names_XYZ = {['\phi ',BodyName],['\theta ',BodyName],['\psi ',BodyName],['P ',BodyName],['Q ',BodyName],['R ',BodyName]};
ylabels_XYZ = {['\phi (',dim1,') ',BodyName],['\theta (',dim1,') ',BodyName],['\psi (',dim1,') ',BodyName],['p (',dim1,'/s) ',BodyName],['q (',dim1,'/s) ',BodyName],['r (',dim1,'/s) ',BodyName]};

for idx = 1:6
  fancy_plotting(1,Names_XYZ{idx},18,'Time (sec)',ylabels_XYZ{idx});
  plot(tout,180/pi*xout(:,idx),'b-')
  title(['Final Value = ',num2str(xout(end,idx)*180/pi)])
  hold on
end