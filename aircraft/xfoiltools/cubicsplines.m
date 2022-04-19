function cubicsplines(X,Y,p0)

%%%Linear Splines
for idx = 1:length(X)-1
   m = (Y(idx+1)-Y(idx) )/(X(idx+1)-X(idx));
   xspline = linspace(X(idx),X(idx+1),10);
   yspline = m*(xspline-X(idx)) + Y(idx);
   p1 = plot(xspline,yspline,'g-','LineWidth',2);
end

%%%Quadratic Splines
N = length(X);
n = N-1;

%%%%Create H0 = [2*(n-1) x 3*n]
H0 = zeros(2*n-2,3*n);
K0 = zeros(2*n-2,1);
for idx = 1:(n-1)
    col = idx;
    row = 2*(idx-1)+1;
    H0(row,col) = X(idx+1)^2;
    H0(row+1,col+1) = X(idx+1)^2;
    H0(row,n+col) = X(idx+1);
    H0(row+1,n+col+1) = X(idx+1);
    H0(row,2*n+col) = 1;
    H0(row+1,2*n+col+1) = 1;
    K0(row) = Y(idx+1);
    K0(row+1) = Y(idx+1);
end

%%%Create H1 = [n-1 x 3*n];
H1 = zeros(n-1,3*n);
for idx = 1:n-1
    H1(idx,idx) = 2*X(idx+1);
    H1(idx,idx+1) = -2*X(idx+1);
    H1(idx,n+idx) = 1;
    H1(idx,n+idx+1) = -1;
end
%%%Create HE
HE = zeros(2,3*n);
HE(1,1) = X(1)^2;
HE(1,n+1) = X(1);
HE(1,2*n+1) = 1;
HE(2,n) = X(end)^2;
HE(2,2*n) = X(end);
HE(2,end) = 1;

%%%Compute H
H = [H0;H1;HE];
[r,c] = size(H);

%%%Compute K
K = zeros(r,1);
K(1:(2*(n-1))) = K0;
K(end-1) = Y(1);
K(end) = Y(end);

%%%Assume now that a1 = 0
H = H(:,2:end);

coeffs = inv(H)*K;

%%%Get all your A's,B's, and C's
A = [0;coeffs(1:n-1)];
B = coeffs(n:2*n-1);
C = coeffs(2*n:end);

%%%%Plot the results
for idx = 1:n
    xspline = linspace(X(idx),X(idx+1),10);
    yspline = A(idx)*xspline.^2 + B(idx)*xspline + C(idx);
    p2 = plot(xspline,yspline,'r-','LineWidth',2);
end

%%%%%Cubic Splines
N = length(Y);
Num_eqns = N-1;
Unknowns = 4*Num_eqns;

%%%Create the Y vector
B = zeros(2*(N-2)+2 + (N-2) + (N-2) + 2,1);

%%%Creating H0 - This ensures that the spline passed through all data points
Num_row = 2*(N-2) + 2;
H0 = zeros(Num_row,Unknowns);
for idx = 1:(Num_row/2)
  col = ((idx-1)*4+1);
  row = (idx-1)*2 + 1;
  for jdx = col:(col+3)
    H0(row,jdx) = X(idx)^(3-jdx+col);
    H0(row+1,jdx) = X(idx+1)^(3-jdx+col);
  end
  B(row,1) = Y(idx);
  B(row+1,1) = Y(idx+1);
end

%%%Creating H1 - this ensures that all first derivatives are the same
H1 = zeros(N-2,Unknowns);
for idx = 2:(N-1)
  col = 1 + (idx-2)*4;
  H1(idx-1,col) = 3*X(idx)^2;
  H1(idx-1,col+1) = 2*X(idx);
  H1(idx-1,col+2) = 1;
  H1(idx-1,col+3) = 0;
  H1(idx-1,col+4) = -3*X(idx)^2;
  H1(idx-1,col+5) = -2*X(idx);
  H1(idx-1,col+6) = -1;
  H1(idx-1,col+7) = 0;
end

%%%Create H2 - This ensures that all second derivatives are the same
H2 = 0*H1;
for idx = 2:(N-1)
  col = 1 + (idx-2)*4;
  H2(idx-1,col) = 6*X(idx);
  H2(idx-1,col+1) = 2;
  H2(idx-1,col+2) = 0;
  H2(idx-1,col+3) = 0;
  H2(idx-1,col+4) = -6*X(idx);
  H2(idx-1,col+5) = -2;
  H2(idx-1,col+6) = 0;
  H2(idx-1,col+7) = 0;
end

%%%This ensures that the second derivative at the endpoints must be zero
Hepts = zeros(2,Unknowns);
%%%1st end point
Hepts(1,1) = 6*X(1);
Hepts(1,2) = 2*X(1);
Hepts(1,3) = 1;
Hepts(2,end) = 1;
Hepts(2,end-1) = 2*X(end);
Hepts(2,end-2) = 6*X(end);

%%%stack everything together
H = [H0;H1;H2;Hepts];

%%%Solve for coefficients
ABCs = inv(H)*B;

for idx = 1:Num_eqns
  row = 1 + (idx-1)*4;
  a = ABCs(row);
  b = ABCs(row+1);
  c = ABCs(row+2);
  d = ABCs(row+3);
  xspline = linspace(X(idx),X(idx+1),100);
  yspline = a*xspline.^3 + b*xspline.^2 + c*xspline + d;
  p3 = plot(xspline,yspline,'b-');
end
legend([p0 p1 p2 p3],'Input Data','Linear Spline','Quadratic Spline','Cubic Spline')