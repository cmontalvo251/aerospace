function tout = throttle_curve_alt(tin,tmid,texp)

%%%%Assume that throttle_in varies from 1100 to 1900

%our polynomial is tout = p1 + p2*(tin - 1100) + p3*(tin-1100)^3

% First boundary condition
% tout(tin = 1100) = 1100 = p1
p1 = 1100;

% Second Boundary Condition
% toutprime(tin = tmid) = texp
% toutprime = p2 + 3*p3*(tin-1100)^2
% toutprime(tin = 1100) = p2 + 3*p3*(tmid-1100)^2 = texp
% texp = [1 3*(tmid-1100)^2]*[p2;p3];

% Third Boundary Condition 
% tout(tin = 1900) = 1900 = p1 + p2*800 + p3*800^3
% 1900 - p1 = [800 800^3]*[p2;p3]
% (1900-p1)/800 = [1 800^2]*[p2;p3];

% Setup a linear system of equations
matrix = [1 3*(tmid-1100)^2;1 800^2];
sol = [texp;(1900-p1)/800];
p2p3 = inv(matrix)*sol;
p2 = p2p3(1);
p3 = p2p3(2);

tout = p1 + p2*(tin - 1100) + p3*(tin-1100).^3;

