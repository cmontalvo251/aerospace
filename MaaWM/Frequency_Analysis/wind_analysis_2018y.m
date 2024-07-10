clear
clc
close all
%% File Path Setup
dataloc = '/home/carlos/Files/GitLab_Repos/Research/MaaWM/WRF_Model/';
UVWfrontend
%% U


[f1,ax1] = plottool(1,'Winds',18,'Y Coordinate (m)','X coordinate (m)','Velocity (m/s)');
[f2,ax2] = plottool(1,'Winds',18,'Wavelength (m)','X coordinate (m)','Cn (m/s)');

%%Grab a slice of data
cnu_mesh = [];
cnv_mesh = [];
cnw_mesh = [];
for x = xcoord
  %y
  %%%Plot the raw data
  ui = 0*ycoord;
  vi = 0*ycoord;
  wi = 0*ycoord;
  ii = 1;
  L = ycoord(end)-ycoord(1);
  w0 = 2*pi/L;
  ii = 1;
  for y = ycoord
    uvw = uvwout(x,y,200,0,dataloc,0);
    %uvw = [sin(w0*x);cos(w0*x);sin(2*w0*x)];
    ui(ii) = uvw(1);
    vi(ii) = uvw(2);
    wi(ii) = uvw(3);
    ii = ii + 1;
  end
  %%Run the Fourier Series Fit
  NMAX = 10;
  [d,anu,bnu,iters,frequencies,ufft] = myfft(ui,ycoord,NMAX,0,0);
  cnu = sqrt(anu.^2 + bnu.^2);
  [d,anv,bnv,iters,frequencies,vfft] = myfft(vi,ycoord,NMAX,0,0);
  cnv = sqrt(anv.^2 + bnv.^2);
  [d,anw,bnw,iters,frequencies,wfft] = myfft(wi,ycoord,NMAX,0,0);
  cnw = sqrt(anw.^2 + bnw.^2);
  
  cnu_mesh = [cnu_mesh;cnu];
  cnv_mesh = [cnv_mesh;cnv];
  cnw_mesh = [cnw_mesh;cnw];

  wavelength = 2*pi./(iters*w0);
  
  %%%Plot to make sure the fit is correct
  xslice = x*ones(1,length(ycoord));
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
  %legend('U','V','W')
end

[ww,xx] = meshgrid(wavelength,xcoord);

contourf(ax2,ww,xx,cnw_mesh)

[f3,ax3] = plottool(1,'Winds',18,'Wavelength (m)','X coordinate (m)','Cn (m/s)');
contourf(ww,xx,cnv_mesh)

[f4,ax4] = plottool(1,'Winds',18,'Wavelength (m)','X coordinate (m)','Cn (m/s)');
contourf(ww,xx,cnu_mesh)