function [timestep, nextstate] = RK4(time, derivative, x, y)
% RK4 integration
% given time, derivative function handle, x{state vector} and y{vector of
% auxilary values}
%
% returns timestep and x{state vector}

timestep = y(1);
%timestep = 0.00005 ; % fixed timestep, to be set on case by case basis

xdot1 = feval(derivative, time, x, y);
xdot2 = feval(derivative, time + (.5*timestep), x + (xdot1*.5*timestep), y);
xdot3 = feval(derivative, time + (.5*timestep), x + (xdot2*.5*timestep), y);
xdot4 = feval(derivative, time + timestep, x + (xdot3*timestep), y);

xdotRK4 = (1/6) * (xdot1 + (2*xdot2) + (2*xdot3) + xdot4);

nextstate = x + (timestep * xdotRK4) ;

