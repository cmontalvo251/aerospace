%%%%First let's read in all the datas
for nsat = 1:NUMSATS
  data = dlmread(['Output_Files/State',num2str(nsat-1),'.OUT']);  
  t_vec = data(:,1); %%%This should be the same for all
  %%%Now let's get xyz and ptp
  state_XYZ_vec{nsat} = data(:,2:7)';
  state_QUAT_vec = data(:,8:14)';
  %%%Finally let's convert quat to euler
  q0123_vec = state_QUAT_vec(1:4,:);
  ptp_vec{nsat} = quat2euler(q0123_vec'); %'
end

%%%Ok so now we have state_XYZ_vec{} and ptp_vec{} so animating
%should be fairly easy

%%%Make a world
screenwidth = 1600;
screenheight = 980;
if SIMXYZ && ANIMATEORBIT
  L = 42164*(1000);
  T = 86400;
  M = 5.972*10^24 ;
  G = 6.6742*10^-11;
  mu = G*M;
  REARTH = 6371*1000;
  REARTHND = REARTH/L;
  A = imread('objects/textures/earth_large.png');
  sky = imread('objects/textures/world-map.jpg');
  [XE,YE,ZE] = ellipsoid(0,0,0,REARTH/1000,REARTH/1000,REARTH/1000);
  [forbit,axorbit] = plottool(1,'xy',18,'X Orbit (km)','Y Orbit (km)','Z Orbit (km)');
  orbit_pos = get(gcf,'Position'); %%left bottom width height
  %left = orbit_pos(1) %%1041 - for some reason this returs 1320 sometimes which doesn't make sense.
  %bottom = orbit_pos(2); %%554
  width = orbit_pos(3)*0.8; %%560
  height = orbit_pos(4)*0.8; %%420
  bottom = screenheight - height;
  left = screenwidth - width;
  left0 = left;
  bottom0 = bottom;
  %get(0,'sreensize') returns 3200 and 1080 which is also wrong. Hmmmm.
  set(forbit,'Position',[left bottom width height])
else
  fig = figure();
  pos = get(gcf,'Position'); %%left bottom width height
  screenwidth = 1600;
  screenheight = 980;
  width = pos(3)*0.8; %%560
  height = pos(4)*0.8; %%420
  bottom = screenheight - height;
  left = screenwidth - width;
  left0 = left;
  bottom0 = bottom;  
  close(fig);
end
colors = {'r-','g-','b-'};
facecolors = {[1 0 0],[0 1 0],[0 0 1]};

if SIMPTP
  for nsat = 1:NUMSATS
    [fzoom{nsat},axzoom{nsat}] = plottool(1,['Satellite',num2str(nsat-1)],18,'X Orbit (km)','Y Orbit (km)','Z Orbit (km)');
    left = left - width;
    if left < 0
      left = left0;
      bottom = bottom - height*1.3;
    end
    set(fzoom{nsat},'Position',[left bottom width height])
  end
end

pause

for nframe = 1:skip:length(t_vec)
  %%Clear figures
  if SIMPTP
    for nsat = 1:NUMSATS
      cla(axzoom{nsat});
      hold(axzoom{nsat},'on');
    end
  end
  if SIMXYZ && ANIMATEORBIT
    cla(axorbit);
    hold(axorbit,'on');

    %%PLot the earth
    surf(axorbit,XE,YE,ZE);
    %set(0,'CurrentFigure',forbit)
    h = findobj('Type','surface');
    set(h,'CData',flipdim(sky,1),'FaceColor','texturemap','edgecolor','none','FaceLighting','Gouraud','Clipping','off')
  end

  %%%Plot individual satellites
  for nsat = 1:NUMSATS
    %%%Extract State
    state_XYZ = state_XYZ_vec{nsat}./1000;
    state_PTP = ptp_vec{nsat};
    %%%Plot on orbit
    if SIMXYZ && ANIMATEORBIT
      %set(0,'CurrentFigure',forbit)
      plot3(axorbit,state_XYZ(1,1:nframe),state_XYZ(2,1:nframe),state_XYZ(3,1:nframe),colors{nsat})
    end

    if SIMPTP
      %%Ok now plot a zoomed photo
      set(0,'CurrentFigure',fzoom{nsat})
      x = state_XYZ(1,nframe);
      y = state_XYZ(2,nframe);
      z = state_XYZ(3,nframe);
      phi = state_PTP(nframe,1);
      theta = state_PTP(nframe,2);
      psi = state_PTP(nframe,3);
      CubeDraw((10/100)/1000,(10/100)/1000,(20/100)/1000,x,y,z,phi,theta,psi,facecolors{nsat})
      view(axzoom{nsat},37,18)
      axis equal
      title(axzoom{nsat},num2str(t_vec(nframe)))
    end
  end
  
  %%%Set parameters for orbit plot
  if SIMXYZ && ANIMATEORBIT
    %set(0,'CurrentFigure',forbit)
    view(axorbit,37,18)
    axis equal
    title(axorbit,num2str(t_vec(nframe)))
  end
  
  %%%Draw everything now
  drawnow
  %pause(.01)
  %f = getfilename(n,5);
  %saveas(gcf,['Frames/Frame_',f,'.jpg'])
end 


