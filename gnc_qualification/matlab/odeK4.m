function [tout,xout,other] = odeK4(functionHandle,tspan,xinitial,timestep,extraparameters,next,quat)
%function [tout,xout] = odeK4(functionHandle,tspan,xinitial,timestep,extraparameters,next)
%%This function operates much like ode45 only it uses a fixed step
%RK4 integration
[N,flag] = size(xinitial);
if flag > 1
  disp('xinitial must be a column vector')
  tout = 0;
  xout = 0;
  return
end
tout = tspan(1):timestep:tspan(end);
integrationsteps = length(tout);
xout = zeros(N,integrationsteps);
x = xinitial;
if ~exist('extraparameters','var')
  extraparameters = 0;
  y = 0;
elseif isempty(extraparameters)
  extraparameters = 0;
  y = 0;
else
  y = extraparameters;
  extraparameters = 1;
  other = zeros(length(y),integrationsteps);
end
if ~exist('next','var')
  next = 0;
end
if strcmp(next,'off')
    next = 0;
end
threshold = 0;
if ~exist('quat','var')
  quat = 'none';
end
for ii = 1:length(tout)
  time = tout(ii);
  if next
    percent = (floor(100*(time-tout(1))/(tout(end)-tout(1))));
    if percent >= threshold
      disp(['Simulation ',num2str(percent),'% Complete'])
      threshold = threshold + next;
    end
  end
  %%Save States
  xout(:,ii) = x;
  %xdot(:,ii) = xdotRK4;
  %%Integrate
  if extraparameters
    if find(y == 42)
      [xdot1,outs] = feval(functionHandle, time, x, y);
      other(:,ii) = outs;
      [xdot2,outs] = feval(functionHandle, time + (.5*timestep), x + (xdot1*.5*timestep), y);
      [xdot3,outs] = feval(functionHandle, time + (.5*timestep), x + (xdot2*.5*timestep), y);
      [xdot4,outs] = feval(functionHandle, time + timestep, x + (xdot3*timestep), y);
      xdotRK4 = (1/6) * (xdot1 + (2*xdot2) + (2*xdot3) + xdot4);
      nextstate = x + (timestep * xdotRK4);
    else
      xdot1 = feval(functionHandle, time, x, y);
      xdot2 = feval(functionHandle, time + (.5*timestep), x + (xdot1*.5*timestep), y);
      xdot3 = feval(functionHandle, time + (.5*timestep), x + (xdot2*.5*timestep), y);
      xdot4 = feval(functionHandle, time + timestep, x + (xdot3*timestep), y);
      xdotRK4 = (1/6) * (xdot1 + (2*xdot2) + (2*xdot3) + xdot4);
      nextstate = x + (timestep * xdotRK4);
    end
  else
    xdot1 = feval(functionHandle, time, x);
    xdot2 = feval(functionHandle, time + (.5*timestep), x + (xdot1*.5*timestep));
    xdot3 = feval(functionHandle, time + (.5*timestep), x + (xdot2*.5*timestep));
    xdot4 = feval(functionHandle, time + timestep, x + (xdot3*timestep));
    xdotRK4 = (1/6) * (xdot1 + (2*xdot2) + (2*xdot3) + xdot4);
    nextstate = x + (timestep * xdotRK4);
  end
  if strcmp(quat,'Quat')
    %%Normalize the quaternions
    q0 = nextstate(4);
    q1 = nextstate(5);
    q2 = nextstate(6);
    q3 = nextstate(7);
    normq = sqrt(q0^2+q1^2+q2^2+q3^2);
    nexstate(4) = nextstate(4)/normq;
    nexstate(5) = nextstate(5)/normq;
    nexstate(6) = nextstate(6)/normq;
    nexstate(7) = nextstate(7)/normq;
  end
  x = nextstate;
end

% Copyright - Carlos Montalvo 2015
% You may freely distribute this file but please keep my name in here
% as the original owner
