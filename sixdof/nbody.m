function [tout,xout] = nbody(Fext_Mext,tspan,xinitial,timestep,mass,MOI,next)
%%This function uses a 6dof fixed timestep RK4 integrator to
%propagate the state forward. The only thing the user must due is
%to include the script to add in forces given to the bodies
%%The states are standard AE format x = [x y z phi theta psi u v w
%p q r];
%The form of the forces and moment function is as follows
%[XYZ,LMN] = Fext_Mext(time,state,inertia)
%state is a 6x1, and inertia is a cell array with inertia{1} =
%mass and inertia{2} = I(moment of inertia matrix)
%LMN and XYZ are column vectors of body moments and forces respectively

%%Get size of vectors
[Numstates,numbodies] = size(xinitial);

%%Setup time
tout = tspan(1):timestep:tspan(end);
integrationsteps = length(tout);

%%Set up state vector
xout = zeros(Numstates,integrationsteps,numbodies);
x = xinitial;
threshold = 0;
for ii = 1:length(tout)
  time = tout(ii);
  if next
    percent = (floor(100*time/tout(end)));
    if percent >= threshold
      disp(['Simulation ',num2str(percent),'% Complete'])
      threshold = threshold + next;
    end
  end
  %%Save States
  xout(:,ii,:) = x;
  %Get mass properties
  %%Integrate
  xdot1 = dynamics( Fext_Mext,time, x,mass,MOI);
  xdot2 = dynamics( Fext_Mext,time + (.5*timestep), x + (xdot1*.5*timestep), mass,MOI);
  xdot3 = dynamics( Fext_Mext,time + (.5*timestep), x + (xdot2*.5*timestep), mass,MOI);
  xdot4 = dynamics( Fext_Mext,time + timestep, x + (xdot3*timestep), mass,MOI);
  xdotRK4 = (1/6) * (xdot1 + (2*xdot2) + (2*xdot3) + xdot4);
  nextstate = x + (timestep * xdotRK4);
  x = nextstate;
end
  
xout = xout; %%This changes everything to a more file friendly format
tout = tout;

function dstatedt = dynamics(Fext_Mext,time,xin,mass,MOI)

%initialize
dstatedt = 0.*xin;

%%Get Forces and Moments on all particles
[XYZ,LMN] = feval(Fext_Mext,time,xin);
Numparticles = length(mass);

for ii = 1:Numparticles

  %%%Unwrap states
  state = xin(:,ii);
  x = state(1);
  y = state(2);
  z = state(3);
  phi = state(4);
  theta = state(5);
  psi = state(6);
  u = state(7);
  v = state(8);
  w = state(9);
  p = state(10);
  q = state(11);
  r = state(12);

  %Compute sines and cosines
  ctheta = cos(theta);
  stheta = sin(theta);
  ttheta = tan(theta);
  sphi = sin(phi);
  cphi = cos(phi);
  spsi = sin(psi);
  cpsi = cos(psi);


  %%Unwrap Inertia
  m = mass(ii);
  I = MOI(:,:,ii);

  %%Derivatives

  %Kinematics
  Kxyz = [ctheta*cpsi sphi*stheta*cpsi-cphi*spsi cphi*stheta*cpsi+sphi*spsi;ctheta*spsi sphi*stheta*spsi+cphi*cpsi cphi*stheta*spsi-sphi*cpsi;-stheta sphi*ctheta cphi*ctheta];
  xyzdot = Kxyz*[u v w]';

  Kptp = [1 sphi*ttheta cphi*ttheta;0 cphi -sphi;0 sphi/ctheta cphi/ctheta];
  pqr = [p q r]';
  ptpdot = Kptp*pqr;

  %%Dynamics
  Kuvw = skew(pqr);
  uvwdot = -Kuvw*[u v w]'+XYZ(:,ii)./m;

  invI = inv(I);
  pqrdot =invI*LMN(:,ii)+invI*(Kuvw*(I*[p q r]'));

  dstatedt(:,ii) = [xyzdot;ptpdot;uvwdot;pqrdot];
end


% Copyright - Carlos Montalvo 2015
% You may freely distribute this file but please keep my name in here
% as the original owner
