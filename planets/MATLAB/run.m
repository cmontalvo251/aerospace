purge

global m

Numparticles = 4;

dt = 0.001;
tout = 0:dt:5;

%%Generate Vectors
xpos_vec = zeros(Numparticles,length(tout));
ypos_vec = xpos_vec;
xvel_vec = xpos_vec;
yvel_vec = xpos_vec;
m = [0];
r = [0];
xposition = [-50 50 0];
yposition = [0 0 50];
%%Generate initial conditions
pos = 100;
vel = 10;
for i = 1:Numparticles
  xpos_vec(i,1) = -pos + 2*pos*rand(1,1);
  ypos_vec(i,1) = -pos + 2*pos*rand(1,1);
  xvel_vec(i,1) = -vel + 10*vel*rand(1,1);
  yvel_vec(i,1) = -vel + 10*vel*rand(1,1);
 end
m = [100 100 100 100 100 100];
r = [1 1 1 1 1];
limit = 100;
threshold = 0;
next = 1;
%%Start Integrating
for ii = 1:length(tout)
  if next
    percent = (floor(100*(tout(ii)-tout(1))/(tout(end)-tout(1))));
    if percent >= threshold
      disp(['Simulation ',num2str(percent),'% Complete'])
      threshold = threshold + next;
    end
  end
  %%Integrate
  xstate = [[xpos_vec(:,ii);ypos_vec(:,ii)],[xvel_vec(:,ii);yvel_vec(:,ii)]];
  [dt,xnext] = RK4(tout(ii),@Gravity,xstate,dt);
  if ii ~= length(tout)
     xpos_vec(:,ii+1) = xnext(1:Numparticles,1);
     ypos_vec(:,ii+1) = xnext(Numparticles+1:2*Numparticles,1);
     xvel_vec(:,ii+1) = xnext(1:Numparticles,2);
     yvel_vec(:,ii+1) = xnext(Numparticles+1:2*Numparticles,2);
  end
end
tscale = 10;

%Make out file
xpos_vec = xpos_vec';
ypos_vec = ypos_vec';
[r,c] = size(xpos_vec);
step = 100;
l = 1:step:r;
outfile = [];
for ii = 1:Numparticles
    outfile = [outfile,xpos_vec(:,ii),ypos_vec(:,ii)];   
end

dlmwrite('/Users/carlos/Dropbox/Projects/Planets/Processing/planets/position.txt',outfile,'delimiter',',');

%save('planets.mat','xpos_vec','ypos_vec','xvel_vec','yvel_vec','tout','Numparticles','tscale','r')
%make_sound
%buildsound

%play_movie;