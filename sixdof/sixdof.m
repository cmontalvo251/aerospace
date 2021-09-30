function [tout,xout] = sixdof(Fext_Mext,tspan,xinitial,timestep,inertia,next)
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

[N,flag] = size(xinitial);
if flag > 1
  disp('xinitial must be a column vector')
  tout = 0;
  xout = 0;
  return
end
if N ~= 12
   disp('Model must have 12 initial conditions (6 states and 6 derivatives)')
   tout = 0;
   xout = 0;
   return
end
tout = tspan(1):timestep:tspan(end);
integrationsteps = length(tout);
xout = zeros(N,integrationsteps);
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
  xout(:,ii) = x;
  %%Integrate
  xdot1 = dynamics( Fext_Mext,time, x, inertia);
  xdot2 = dynamics( Fext_Mext,time + (.5*timestep), x + (xdot1*.5*timestep), inertia);
  xdot3 = dynamics( Fext_Mext,time + (.5*timestep), x + (xdot2*.5*timestep), inertia);
  xdot4 = dynamics( Fext_Mext,time + timestep, x + (xdot3*timestep), inertia);
  xdotRK4 = (1/6) * (xdot1 + (2*xdot2) + (2*xdot3) + xdot4);
  nextstate = x + (timestep * xdotRK4);
  x = nextstate;
end

xout = xout'; %%This changes everything to a more file friendly format
tout = tout';

function dstatedt = dynamics(Fext_Mext,time,state,inertia)

%%Get Forces and Moments
[XYZ,LMN] = feval(Fext_Mext,time,state,inertia);

%%%Unwrap states
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
m = inertia{1};
I = inertia{2};

%%Derivatives

%Kinematics
TIB = [ctheta*cpsi sphi*stheta*cpsi-cphi*spsi cphi*stheta*cpsi+sphi*spsi;ctheta*spsi sphi*stheta*spsi+cphi*cpsi cphi*stheta*spsi-sphi*cpsi;-stheta sphi*ctheta cphi*ctheta];
xyzdot = TIB*[u v w]';

H = [1 sphi*ttheta cphi*ttheta;0 cphi -sphi;0 sphi/ctheta cphi/ctheta];
pqr = [p q r]';
ptpdot = H*pqr;

%%Dynamics
Kuvw = skew(pqr);
uvwdot = -Kuvw*[u v w]'+XYZ./m;

invI = inv(I);
pqrdot =invI*LMN-invI*(Kuvw*(I*[p q r]'));

dstatedt = [xyzdot;ptpdot;uvwdot;pqrdot];

% Copyright - Carlos Montalvo 2015
% You may freely distribute this file but please keep my name in here
% as the original owner
