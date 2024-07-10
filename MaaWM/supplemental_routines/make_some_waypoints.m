purge

num_Aircraft = 1;
num_quads = 1;

[DYNSKIP,BOUNDARY,IKIND,K35] = read_boundary('Input_Files/MaaWM.SIM');

outfile = 'Input_Files/Waypoints.WAY';
iplot = 1;
freq = 20;
pathtype = 'spiral';
mult = 0;
ALT = 0;
inset = (BOUNDARY-0.5*BOUNDARY);

Create_Waypoints_AC_QUAD(num_Aircraft,num_quads,BOUNDARY,outfile,iplot,freq,pathtype,mult,ALT,inset);

