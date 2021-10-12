function out = integrate_curve(x,y)

%%%The Orbit model right now is a constant step integrator which means we can just
%%Use a summing algorithm
dx = x(2)-x(1);
out = sum(dx*y);