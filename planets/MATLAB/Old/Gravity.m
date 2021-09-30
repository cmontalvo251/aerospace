function dxdt = Gravity(t,x,y)
global m

[Numparticles,c] = size(x);
Numparticles = Numparticles/2;

xpos_vec = x(1:Numparticles,1);
ypos_vec = x(Numparticles+1:2*Numparticles,1);
xvel_vec = x(1:Numparticles,2);
yvel_vec = x(Numparticles+1:2*Numparticles,2);


%%Calculate Position Derivatives
xposdot = xvel_vec;
yposdot = yvel_vec;
xveldot = zeros(Numparticles,1);
yveldot = zeros(Numparticles,1);
%%Calculate Velocity Derivatives
for i = 1:length(xposdot)
  for j = 1:length(xposdot)
    if i ~= j
      %%Calculate Gravitational pull
      G = m(i)*m(j)*100000;
      r = (ypos_vec(i)-ypos_vec(j))^2 + (xpos_vec(i)-xpos_vec(j))^2;
      Fmag = G/r;
      direction = -[xpos_vec(i) - xpos_vec(j);ypos_vec(i) - ypos_vec(j)];
      u = direction/norm(direction);
      if abs(Fmag) > 10000000/2
	Fmag = sign(Fmag)*10000000/2;
      end
      F = Fmag*u;
      xveldot(i) = xveldot(i) + F(1)/m(i);
      yveldot(i) = yveldot(i) + F(2)/m(i);
    end
  end
end
xveldot;
yveldot;
%pause
dxdt = [[xposdot;yposdot],[xveldot;yveldot]];
