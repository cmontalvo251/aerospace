function tout = throttle_curve(tin,tmid,texp)

%%%%Assume that throttle_in varies from 1100 to 1900

%our polynomial is tout = p1 + p2*(tin - 1100) + p3*(tin-1100)^2 + p4*(tin-1100)^3

% First boundary condition
% tout(tin = 1100) = 1100 = p1
p1 = 1100;

% Second Boundary Condition
% toutprime(tin = tmid) = texp
% toutprime = p2 + 2*p3*(tin-1100) + 3*p4*(tin-1100)^2
% toutprime(tin = 1100) = p2 = texp
p2 = texp;

% Third Boundary Condition 
% tout(tin = 1900) = 1900 = p1 + p2*800 + p3*800^2 + p4*800^3;
% 1900 - p1 - p2*800 = [800^2 800^3]*[p3;p4]
% (1900-p1-p2*800)/800^2 = [1 800]*[p3;p4];

% Fourth Boundary Condition
% toutprime(tin = tmid) = 0
% toutprime(tin = tmid) = p2 + 2*p3*(tmid-1100) + 3*p4*(tmid - 1100)^2 = 0
% -p2 = [2*(tmid-1100) 3*(tmid-1100)^2]*[p3;p4]

% Setup a linear system of equations
matrix = [1 800;2*(tmid-1100) 3*(tmid-1100)^2];
sol = [(1900-p1-p2*800)/800^2;-p2];
p3p4 = inv(matrix)*sol;
p3 = p3p4(1);
p4 = p3p4(2);

tout = p1 + p2*(tin - 1100) + p3*(tin-1100).^2 + p4*(tin-1100).^3;

