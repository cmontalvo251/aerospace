%%%BOIDS Algorithm
purge

global size_of_grid tnext xcommand_vec ycommand_vec tskip counter

num_boids = 3;
Velocity = 10;
size_of_grid = 500;
radius = 10;
beak_length = 2*radius;
dyn = 0.9;
dt = 0.5;
t = 0;
tskip = 200;
tnext = tskip;
xcommand_vec = [size_of_grid/2 size_of_grid/2 -size_of_grid/2 -size_of_grid/2];
ycommand_vec = [-size_of_grid/2 size_of_grid/2 size_of_grid/2 -size_of_grid/2];
counter = 1;

%%Create Boids (Assume Dubins Model)
xboids = struct('Position',[0;0],'Velocity',10,'Orientation',0);
for ii = 1:num_boids
  x = -size_of_grid + 2*size_of_grid*rand;
  y = -size_of_grid + 2*size_of_grid*rand;
  xboids(ii).Position = [x;y];
  xboids(ii).Velocity = 10;
  xboids(ii).Orientation = 2*pi*rand;
  xboids(ii).Radius = radius;
end

plottool(1,'Boids',12);
while 1
  %%%Plotting
  cla;
  for ii = 1:num_boids
    x = xboids(ii).Position(1);
    y = xboids(ii).Position(2);
    %%Plot Head
    if ii == 1
      col = 'r';
    else
      col = 'b';
    end
    circle(x,y,radius,col)
    %%Plot Beak
    theta = xboids(ii).Orientation;
    %plot([x x+beak_length*cos(theta)],[y y+beak_length*sin(theta)],'LineWidth',3)
  end
  axis([-size_of_grid size_of_grid -size_of_grid size_of_grid])
  axis square
  title(['T = ',num2str(t),' Tnext = ',num2str(tnext)])

  %%%Propagate
  for ii = 1:num_boids
    theta = xboids(ii).Orientation;
    %%Position
    xboids(ii).Position = xboids(ii).Position + dt*xboids(ii).Velocity*[cos(theta);sin(theta)];
    %%Control
    control = compute_control(t,xboids,ii);
    %%Orientation
    xboids(ii).Orientation = xboids(ii).Orientation*dyn + control*(1-dyn);
  end
  %%%Time
  t = t + dt;
  
  drawnow
end

    


