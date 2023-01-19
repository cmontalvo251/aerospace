function tout = throttle_curve_FINAL(tin,tmid,texp)

%%%%Assume that throttle_in varies from 1100 to 1900

% Assume texp must be 1 or greater
if texp < 1
    texp = 1;
end

%our polynomial is 5th order with 6 coefficients

% First boundary condition
% tout(tin = 1100) = 1100 = p1
p1 = 1100;

% Second Boundary Condition
% toutprime(tin = tmid) = texp
% toutprime = p2 + 2*p3*(tin-1100) + 3*p4*(tin-1100)^2
% toutprime(tin = 1100) = p2 = texp
p2 = texp;

% Third Boundary Condition 
% tout(tin = 1900) = 1900

% Fourth Boundary Condition
% tout(tin = tmid) = tmid

% Fifth Boundary Condition
% toutprime(tin = 1900) = texp

% Sixth Boundary Condition
% toutprime(tin = tmid) = 2-texp --- if texp > 2, 0

texp_end = 2-texp;
if texp_end < 0
    texp_end = 0;
end

% Setup a linear system of equations
dt = (tmid-1100);
matrix = [800^2 800^3 800^4 800^5;dt^2 dt^3 dt^4 dt^5;2*800 3*800^2 4*800^3 5*800^4;2*dt 3*dt^2 4*dt^3 5*dt^4];
sol = [1900-p1-800*p2;tmid-p1-dt*p2;texp-p2;texp_end-p2];
coeff = inv(matrix)*sol;
p3 = coeff(1);
p4 = coeff(2);
p5 = coeff(3);
p6 = coeff(4);

% Compute Throttle Out
dt = tin-1100;
tout = p1 + p2*dt + p3*dt.^2 + p4*dt.^3 + p5*dt.^4 + p6*dt.^5;

% Detect Negative Slope
toutprime = p2 + 2*p3*dt + 3*p4*dt.^2 + 4*p5*dt.^3 + 5*p6*dt.^4;

if sum(toutprime < 0)
    disp('Negative Slope Detected - Reduce Value of texp')
    texp = 0;
    tout = throttle_curve_FINAL(tin,tmid,texp);
end

