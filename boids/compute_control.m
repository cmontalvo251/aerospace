function control = compute_control(t,xboids,ii)
%%%Flocking Algorithm
global size_of_grid tnext xcommand_vec ycommand_vec tskip counter

x = xboids(ii).Position(1);
y = xboids(ii).Position(2);
psi = xboids(ii).Orientation;

if ii == 1
  %%%Main Boid so just follow waypoints
  xc = xcommand_vec(counter);
  yc = ycommand_vec(counter);
  dx = (x-xc);
  dy = (y-yc);
  d = sqrt(dx^2+dy^2);
  %%Waypoint tracking
  if d <= 2*xboids(ii).Radius
    counter = counter + 1;
    if counter > 4
      counter = 1;
    end
  end
else
  %%%Other boids so participate in flocking algorithm
  %%Main goal is to follow main boid
  xc = xboids(1).Position(1);
  yc = xboids(1).Position(2);
end

psicoord = atan2(yc-y,xc-x);
spsi = sin(psi);
cpsi = cos(psi);
spsic = sin(psicoord);
cpsic = cos(psicoord);

delpsi = atan2(spsi*cpsic-cpsi*spsic,cpsi*cpsic+spsi*spsic);

control = psi-delpsi;

if ii ~= 1
  %%%Avoidance routine
  %%First we must find our nearest_neighbors
  num_boids = length(xboids);
  dmin = 1e20;
  for jj = 1:num_boids
    if jj ~= ii
      dx = (x-xboids(jj).Position(1));
      dy = (y-xboids(jj).Position(2));
      d = sqrt(dx^2 + dy^2);
      if d < dmin
	dmin = d;
      end
    end
  end
  %%%Now we do a waiting sum to compute control
  %%%if dmin = 0 then control = control + pi
  %%%if dmin = inf then control = control + 0
  %%%This could work -> 2*(pi/2-atan(dmin))
  %%%This isn't quite working
  del = 2*(pi/2-atan((dmin-2*xboids(ii).Radius)^2));
  control = control + del;
  if dmin <= 2*xboids(ii).Radius
    disp(['Parameters(dmin,del) = ',num2str(dmin),' ',num2str(del)])
  end
end


