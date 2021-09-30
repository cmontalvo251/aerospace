% Copyright - Carlos Montalvo 2015
% You may freely distribute this file but please keep my name in here
% as the original owner
clear
clc
close all

ANIMATE = 1;
ANIMATEORBIT = 1;
SIMULATE = 1;
skip = 10;

if SIMULATE
  system('make');
  system('./Orbit.exe');
end

%%%Now that these parameters are in input files we can pull these straight
%%%from the input files
fid = fopen('Input_Files/Satellite.SIM');
ln = fgetl(fid);
SIMXYZ = str2num(ln(1:(find(ln=='!')-1)));
ln = fgetl(fid);
SIMPTP = str2num(ln(1:(find(ln=='!')-1)));
ln = fgetl(fid);
SIMGRAV = str2num(ln(1:(find(ln=='!')-1)));
ln = fgetl(fid); %%%year
ln = fgetl(fid); 
SIMMAGN = str2num(ln(1:(find(ln=='!')-1)));
ln = fgetl(fid); %%update rate of mag field
ln = fgetl(fid); %%timestep
ln = fgetl(fid); %%tfinal
ln = fgetl(fid); %%print rate
ln = fgetl(fid); 
NUMSATS = str2num(ln(1:(find(ln=='!')-1)));
fclose(fid);

if ANIMATE
  animate_cubes();
  return
end

BodyName = '';
dim1 = 'km';
maincolor = 'k-';
colors = {'b-','r-','g-'};

%%%Generate plots
if SIMXYZ
  [forbit,axorbit] = plottool(1,'xy',18,'X Orbit (km)','Y Orbit (km)','Z Orbit (km)');

  for idx = 1:6
    Names_XYZ = {['X ',BodyName],['Y ',BodyName],['Z ',BodyName],['U ',BodyName],['V ',BodyName],['W ',BodyName]};
    ylabels_XYZ = {['x (',dim1,') ',BodyName],['y (',dim1,') ',BodyName],['z (',dim1,') ',BodyName],['u (',dim1,'/s) ',BodyName],['v (',dim1,'/s) ',BodyName],['w (',dim1,'/s) ',BodyName]};
    [fstates{idx},axstates{idx}] = plottool(1,Names_XYZ{idx},18,'Time (min)',ylabels_XYZ{idx});
  end
end

if SIMPTP
  ylabels_PTP = {['q0 (rad) ',BodyName],['q1 (rad) ',BodyName],['q2 (rad) ',BodyName],['q3 (rad) ',BodyName],['p (rad/s) ',BodyName],['q (rad/s) ',BodyName],['r (rad/s) ',BodyName]};
  Names_PTP = {['q0 ',BodyName],['q1 ',BodyName],['q2 ',BodyName],['q3 ',BodyName],['P ',BodyName],['Q ',BodyName],['R ',BodyName]};
  for idx = 1:7
    [fptp{idx},axptp{idx}] = plottool(1,Names_PTP{idx},18,'Time (min)',ylabels_PTP{idx});
  end

  eangle = {['\phi (deg) ',BodyName],['\theta (deg) ',BodyName],['\psi (deg) ',BodyName]};
  for idx = 1:3
    [feuler{idx},axeuler{idx}] = plottool(1,eangle{idx},18,'Time (min)',eangle{idx});
  end
end

if SIMGRAV == 1
  ylabels_GRAV = {['Gx (m/s^2) ',BodyName],['Gy (m/s^2) ',BodyName],['Gz (m/s^2) ',BodyName]};
  Names_GRAV = {['Gx ',BodyName],['Gy ',BodyName],['Gz ',BodyName]};
  for idx = 1:3
    [fgrav{idx},axgrav{idx}] = plottool(1,Names_GRAV{idx},18,'Time (min)',ylabels_GRAV{idx});
  end
  [fgravmag,axgravmag] = plottool(1,'Gravity Magnitude',12,'Time (min)','Magnitude of Gravity (m/s^2)');
end

if SIMMAGN == 1
  ylabels_MAGN = {['Bx (nT) ',BodyName],['By (nT) ',BodyName],['Bz (nT) ',BodyName]};
  Names_MAGN = {['Bx (Body Frame)',BodyName],['By (Body Frame)',BodyName],['Bz (Body Frame)',BodyName]};
  for idx = 1:3
    [fmagn{idx},axmagn{idx}] = plottool(1,Names_MAGN{idx},18,'Time (min)',ylabels_MAGN{idx});
  end
end


for nsat = 1:NUMSATS
  data = dlmread(['Output_Files/State',num2str(nsat-1),'.OUT']);  
  state_XYZ_vec = data(:,2:7)';
  state_PTP_vec = data(:,8:14)';
  t_vec = data(:,1);
  GVEC = data(:,15:17)';
  BECSPH = data(:,18:20)';
  BVECINE = data(:,21:23)';
  BVEC = data(:,24:26)'; %%%These are body frame components
  CURRENT = data(:,27:29)'; %%THese are the currents in amps

  %%%Display last XYZ coordinate
  state_XYZ_vec_last = state_XYZ_vec(:,end)'

  %%%%Convert to lat and lon
  x = state_XYZ_vec(1,:);
  y = state_XYZ_vec(2,:);
  z = state_XYZ_vec(3,:);
  rho = sqrt(x.^2 + y.^2 + z.^2);
  phi = (acos(z ./ rho));
  the = atan2(y , x);
  lat = 90 - phi*(180 / pi);
  lon = the*(180/pi);

  if (SIMXYZ)
    
    plottool(1,'GPS Coordinates',18,'Time (min)','Latitude/Longitude');
    plot(t_vec./60,lat,'b-')
    plot(t_vec./60,lon,'r-')
    legend('Latitude','Longitude')
    
    %%%Plot orbit
    plot3(axorbit,state_XYZ_vec(1,:)./1000,state_XYZ_vec(2,:)./1000,state_XYZ_vec(3,:)./1000,colors{nsat})

    %x,y,z,q0,q1,q2,q3,u,v,w ,p ,q , r
    %1,2,3,4 ,5 ,6 ,7 ,8,9,10,11,12,13
    for idx = 1:6
      plot(axstates{idx},t_vec./60,state_XYZ_vec(idx,:)./1000,colors{nsat})
    end
    
  end

  if (SIMPTP)
       
    for idx = 1:7
      plot(axptp{idx},t_vec./60,state_PTP_vec(idx,:),colors{nsat})
    end

    q0123_vec = state_PTP_vec(1:4,:);
    ptp_vec = quat2euler(q0123_vec'); %'


    for idx = 1:3
      plot(axeuler{idx},t_vec./60,ptp_vec(:,idx)*180/pi,colors{nsat})
    end
    
  end

  if (SIMGRAV==1)

    GVEC_MAG = 0*GVEC(1,:);
    for idx = 1:3
      plot(axgrav{idx},t_vec./60,GVEC(idx,:),colors{nsat})
      GVEC_MAG = GVEC_MAG + GVEC(idx,:).^2;
    end

    %%%Plot the Magnitude of GVEC
    GVEC_MAG = sqrt(GVEC_MAG);
    plot(axgravmag,t_vec./60,GVEC_MAG,colors{nsat})
  end

  if (SIMMAGN)

    for idx = 1:3
      plot(axmagn{idx},t_vec./60,BVEC(idx,:),colors{nsat})
    end

    plottool(1,'All Magnetic',18,'Time (min)','Magnetic Field (Body Frame) (nT)');
    plot(t_vec./60,BVEC,'LineWidth',2)
    legend('X_B','Y_B','Z_B')
    
    fig = figure();
    set(fig,'color','white')
    set(axes,'FontSize',18)
    plot(BVECINE','LineWidth',2)
    xlabel('Longitude (deg)')
    ylabel('Magnetic Field (nT)')
    grid on
    legend('X_I','Y_I','Z_I')

    fig = figure();
    set(fig,'color','white')
    set(axes,'FontSize',18)
    plot(t_vec./60,CURRENT,'LineWidth',2)
    xlabel('Time (min)')
    ylabel('Magnetorquer Current (Amps)')
    grid on
    legend('x','y','z')
  
  end %%%End MAGN
  
end

return

figure()
hold on
view([-27,50])
b_ine = BVEC;
[r,c] = size(BVEC);
% xlim([-1.01 1.01])
% ylim([-1.01 1.01])
% zlim([-2 2])
for idx = 1:c
  %%%Convert to cartesian coordinates
  phi = pi/180*(90 - lat(idx));
  theta = pi/180*lon(idx);
  x = sin(phi)*cos(theta);
  y = sin(phi)*sin(theta);
  z = cos(phi);
  % %%%Speherical
  % bx = BVEC(1,idx);
  % by = BVEC(2,idx);
  % bz = BVEC(3,idx);
  % %%%Convert to Inertial
  % BSPH = [bx;by;bz];
  % %%%Rotation Matrix
  % %%%This assumes a 3-2-1 sequence using standard Euler angle convention
  % psiE = -theta;
  % thetaE = -phi;
  % phiE = 0;
  % %Compute sines and cosines
  % ctheta = cos(thetaE);
  % stheta = sin(thetaE);
  % sphi = sin(phiE);
  % cphi = cos(phiE);
  % spsi = sin(psiE);
  % cpsi = cos(psiE);
  % %Kinematics
  % R = [ctheta*cpsi sphi*stheta*cpsi-cphi*spsi cphi*stheta*cpsi+sphi*spsi;ctheta*spsi sphi*stheta*spsi+cphi*cpsi cphi*stheta*spsi-sphi*cpsi;-stheta sphi*ctheta cphi*ctheta];
  % BINE = R*BSPH;
  BINE = BVECINE(:,idx);
  b_ine(:,idx) = BINE;
  bnorm = norm(BINE);
  BINE_hat = BINE/bnorm;
  f = 1;
  plot3([x x+BINE_hat(1)*f],[y y+BINE_hat(2)*f],[z z+BINE_hat(3)*f],'b-*','MarkerSize',2)
  axis equal
  %drawnow
end


