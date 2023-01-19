purge


x = 0:0.01:200;
u = 1*cos(x);
N = 300;
w = linspace(0.0,2*N/(N-1),N);

%Compute S
sigu = 1;
Lu = 1;
S1 = 0.*w;
S3 = S1;
dx = x(2)-x(1);
for ii = 1:length(w)
  for jj = 1:length(u)
    S1(ii) = S1(ii) + (1/(pi))*u(jj)*cos(w(ii)*x(jj))*dx;
  end
  S1(ii) = norm(S1(ii));
  % if abs(w(ii)-1) > 1e-2
  %   S3(ii) = 1;
  % else
  %   S3(ii) = 0;
  % end
  S3(ii) = 0.1*w(ii);
end

wlower = w(1);
wupper = w(end);
P = length(u);
N = length(S1);
dw = (wupper-wlower)./N;
dwp = 0.01.*dw;

%Regenerate u using S
%%Generate PHI
PHI = zeros(N,1);
dwp_rand = zeros(N,1);
for ii = 1:N(1)
  PHI(ii) = 2*pi*rand;
  dwp_rand(ii) = rand;
end
%%Generate f
f = 0.*u;
for idu = 1:P
  for k1 = 1:N
    w1k1 = wlower(1) + (k1-1/2)*dw;
    delw1 = -dwp(1)/2 + dwp(1)*dwp_rand(k1,1);
    wp = [w1k1+0*delw1];
    % if abs(wp-1) < 1e-3
    %   S = 1;
    % else
    %   S = 0;
    % end
    f(idu) = f(idu) + sqrt(S3(k1))*cos(wp*x(idu)+0*PHI(k1));
  end
end
plottool(1,'Name',12,'x','f');
plot(x,f)
plot(x,u,'r-')
drawnow

%%Regenerate S using f
S2 = 0.*w;
for ii = 1:length(w)
  for jj = 1:length(f)
    S2(ii) = S2(ii) + (1/(pi))*f(jj)*cos(w(ii)*x(jj))*dx;
  end
  S2(ii) = norm(S2(ii))./(N/4);
end
plottool(1,'Name',12,'w','S');
plot(w,S1)
plot(w,S2,'r-')
plot(w,S3,'g-')
drawnow