%%%Simulate Particles orbiting and combining
close all
clear
clc

global m

Numparticles = 100;

%%Generate Vectors
xpos_vec = [0];
ypos_vec = [0];
xvel_vec = [0];
yvel_vec = [0];
m = [0];
r = [0];
xposition = [-50 50 0];
yposition = [0 0 50];
%%Generate initial conditions
pos = 100;
vel = 100;
for i = 1:Numparticles
  xpos_vec(i) = -pos + 2*pos*rand(1,1);
  ypos_vec(i) = -pos + 2*pos*rand(1,1);
  xvel_vec(i) = -vel + 2*vel*rand(1,1);
  yvel_vec(i) = -vel + 2*vel*rand(1,1);
  %xpos_vec(i) = xposition(i);
  %ypos_vec(i) = yposition(i);
  %xvel_vec(i) = 0;
  %yvel_vec(i) = 0;
  m(i) = 1; %%kg
  r(i) = 1; %%m
end
scale = 2;
xpos_vec = xpos_vec';
ypos_vec = ypos_vec';
xvel_vec = xvel_vec';
yvel_vec = yvel_vec';
limit = 100;
%%Generate Plot
h1 = figure();
set(h1,'color','white')
%%Start Integrating
go = 1;
counter = 0;
dt = 0.001;
SKIP = 0;
step = 0.01;
frames = 1;
while go
  hold off
  if counter >= SKIP
    %%Plot States
    for i = 1:length(xpos_vec)
      plot([xpos_vec(i) xpos_vec(i)],[ypos_vec(i) ypos_vec(i)],'bo','LineWidth',r(i)*scale)
      hold on
    end
    if Numparticles < 10
      %step = dt*10;
    end
    if Numparticles < 5
      if (max(abs(xpos_vec))>limit) || (max(abs(ypos_vec)) > limit)
	limit = limit*1.2;
      end
    end
    axis([-limit limit -limit limit])
    title(num2str(counter));
    drawnow
    M(frames) = getframe;
    frames = frames + 1;
    SKIP = SKIP + step;
  end
  counter = counter + dt;
  %%Collision Detection
  loc = [];
  for i = 1:Numparticles-1
    for j = i+1:Numparticles
      dist = sqrt((xpos_vec(i)-xpos_vec(j))^2 + (ypos_vec(i)-ypos_vec(j))^2);
      if dist <= (r(i) + r(j))
	%disp('Collision')
	loc = [i j];
      end
    end
  end
  %%Combine Collided Ones
  if ~isempty(loc)
    deletes = [];
    combines = [0];
    for i = 1:2:length(loc)-1
      loc1 = loc(i);
      loc2 = loc(i+1);
      if isempty(find(combines == loc1,1))
	combines = [combines loc1];
	deletes = [deletes loc2];
	%%Use momentum conservation to combine particles
	mtotal = m(loc1) + m(loc2);
	xvel_vec(loc1) = (m(loc1)*xvel_vec(loc1)+m(loc2)*xvel_vec(loc2))/mtotal;
	yvel_vec(loc1) = (m(loc1)*yvel_vec(loc1)+m(loc2)*yvel_vec(loc2))/mtotal;
	%%Pick the position with the larger mass(if same then average)
	if m(loc1) < m(loc2)
	  xpos_vec(loc1) = xpos_vec(loc2);
	  ypos_vec(loc1) = ypos_vec(loc2);
	elseif m(loc1) == m(loc2)
	  xpos_vec(loc1) = (xpos_vec(loc1)+xpos_vec(loc2))/2;
	  ypos_vec(loc1) = (ypos_vec(loc1)+ypos_vec(loc2))/2;
	end
	Numparticles = Numparticles - 1
	m(loc1) = mtotal;
	%%For Radius assume spheres of constant volume
	Vtotal = (4/3)*pi*r(loc1)^3 + (4/3)*pi*r(loc2)^3;
	r(loc1) = ((Vtotal*(3/4))/pi)^(1/3);
	m(deletes) = [];
	r(deletes)= [];
	xpos_vec(deletes) = [];
	ypos_vec(deletes) = [];
	xvel_vec(deletes) = [];
	yvel_vec(deletes) = [];
      end
    end
    %%Delete Extras
  end
  %%Integrate
  xstate = [[xpos_vec;ypos_vec],[xvel_vec;yvel_vec]];
  [dt,xnext] = RK4(counter,@Gravity,xstate,dt);
  xpos_vec = xnext(1:Numparticles,1);
  ypos_vec = xnext(Numparticles+1:2*Numparticles,1);
  xvel_vec = xnext(1:Numparticles,2);
  yvel_vec = xnext(Numparticles+1:2*Numparticles,2);
  if counter > 200
   go = 0;
  end
end

  

