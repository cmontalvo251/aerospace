clear
clc
close all
%% File Path Setup
%dataloc = '/home/carlos/Files/GitLab_Repos/Research/MaaWM/WRF_Model/';
dataloc = '/home/carlos/Files/GitLab_Repos/Research/MaaWM/Auto_Generated_WindField/Wave1/';
UVWfrontend
%% U

[f1,ax1] = plottool(1,'Winds',18,'X Coordinate (m)','Y coordinate (m)','Velocity (m/s)');
[f2,ax2] = plottool(1,'Winds',18,'Frequency (1/m)','Y coordinate (m)','Cn (m/s)');

%%Grab a slice of data
cnu_mesh = [];
cnv_mesh = [];
cnw_mesh = [];
for y = ycoord
  %y
  %%%Plot the raw data
  ui = 0*xcoord;
  vi = 0*xcoord;
  wi = 0*xcoord;
  ii = 1;
  L = xcoord(end)-xcoord(1);
  w0 = 2*pi/L;
  ii = 1;
  for x = xcoord
    uvw = uvwout(x,y,200,0,dataloc,0);
    %uvw = [sin(w0*x);cos(w0*x);sin(2*w0*x)];
    ui(ii) = uvw(1);
    vi(ii) = uvw(2);
    wi(ii) = uvw(3);
    ii = ii + 1;
  end
  %%Run the Fourier Series Fit
  NMAX = 20;
  [d,anu,bnu,iters,frequencies,ufft] = myfft(ui,xcoord,NMAX,0,0);
  cnu = sqrt(anu.^2 + bnu.^2);
  [d,anv,bnv,iters,frequencies,vfft] = myfft(vi,xcoord,NMAX,0,0);
  cnv = sqrt(anv.^2 + bnv.^2);
  [d,anw,bnw,iters,frequencies,wfft] = myfft(wi,xcoord,NMAX,0,0);
  cnw = sqrt(anw.^2 + bnw.^2);
  
  cnu_mesh = [cnu_mesh;cnu];
  cnv_mesh = [cnv_mesh;cnv];
  cnw_mesh = [cnw_mesh;cnw];

  wavelength = 2*pi./(iters*w0);
  freqs = 1./wavelength;
  
  %%%Plot to make sure the fit is correct
  yslice = y*ones(1,length(xcoord));
  %plot3(ax1,xcoord,yslice,ui,'b*')
  %plot3(ax1,xcoord,yslice,ufft,'b-')
  %plot3(ax1,xcoord,yslice,vi,'r*')
  %plot3(ax1,xcoord,yslice,vfft,'r-')
  %plot3(ax1,xcoord,yslice,wi,'g*')
  %plot3(ax1,xcoord,yslice,wfft,'g-')
  
  %%Plot the Frequencies
  %yslice = y*ones(1,length(wavelength));
  %plot3(ax2,wavelength,yslice,cnu,'b-')
  %plot3(ax2,wavelength,yslice,cnv,'r-')
  %plot3(ax2,wavelength,yslice,cnw,'g-')
end

cnu_mean = mean(cnu_mesh);
cnv_mean = mean(cnv_mesh);
cnw_mean = mean(cnw_mesh);

[f8,ax8] = plottool(1,'Winds',18,'Frequency (1/m)','Cn (m/s)');
%freqs = wavelength;
plot(ax8,freqs,cnu_mean,'b-')
plot(ax8,freqs,cnv_mean,'r-')
plot(ax8,freqs,cnw_mean,'g-')
legend('U','V','W')

[ww,yy] = meshgrid(freqs,ycoord);
contourf(ax2,ww,yy,cnw_mesh)

[f3,ax3] = plottool(1,'Winds',18,'Frequency (1/m)','Y coordinate (m)','Cn (m/s)');
contourf(ww,yy,cnv_mesh)

[f4,ax4] = plottool(1,'Winds',18,'Frequency (1/m)','Y coordinate (m)','Cn (m/s)');
contourf(ww,yy,cnu_mesh)