function freqx = grid_space(datax,boundary,GRIDTYPE,K35,IKIND,dataloc)

global markX markY markZ markT parameters zcoord tcoord
global U0 Udt V0 Vdt W0 Wdt xcoord ycoord tlength tstr
global bounds boundflag tcoord zcoord terrain
global mark1 mark2

%datax = ceil(sqrt(num_data_pts));

datay = datax;

total_area = (2*boundary)^2;

totalpts = datax*datay;

density = 1000*totalpts/total_area

X = linspace(-boundary,boundary,datax);
Y = linspace(-boundary,boundary,datay);
Z = -200*3.28;
freqx = 0;
%%%Get Ready for UVW stuff
switch GRIDTYPE
  case 'WRF'
   dataloc = getdataloc('Input_Files/MaaWM.ifiles') %%%This routine was written to just grab the dataloc from the .ATM file
   %dataloc = '/home/carlos/Documents/GitLab_Repos/MAAWM/Auto_Generated_WindField/Wave1/';
   UVWfrontend;
 case 'RBF'
  root_directory = 'Input_Files/Example_Fit/';
  [dataflag,xknot,alfu_uno,alfv_dos,alfw_tres,uveeze,vveeze,wveeze] = read_some_data(root_directory,K35);
 case 'SINE'
  freqx = 2;
  freqy = freqx;
end

%%%Open file
fid = fopen('Output_Files/Winds.OUT','wb');

%plottool(1,'Grid Space',18,'X (ft)','Y (ft)');
for idx = 1:length(X)
  for jdx = 1:length(Y)
    switch GRIDTYPE
     case 'WRF'
         %Position is in feet and needs to be in meters and output is in m/s and needs to be in ft/s
         uvwwrf = 3.28*uvwout(X(idx)/3.28,Y(jdx)/3.28,-Z/3.28,0,dataloc,1); 
     case 'RBF'
      state = [X(idx) Y(jdx) Z];
      ux = RBFhack2(state,xknot,alfu_uno,IKIND);
      ux = ux + Vbuild(state,uveeze,K35);
      vx = RBFhack2(state,xknot,alfv_dos,IKIND);
      vx = vx + Vbuild(state,vveeze,K35);
      wx = RBFhack2(state,xknot,alfw_tres,IKIND);
      wx = wx + Vbuild(state,wveeze,K35);
      uvwwrf = [ux,vx,wx];
     case 'SINE'
      lx = X(idx)/boundary;
      ly = Y(jdx)/boundary;
      ux = 5*sin(freqx*lx)*cos(freqy*ly);
      vx = 5*sin(freqx*lx)*cos(freqy*ly);
      wx = 5*sin(freqx*lx)*cos(freqy*ly);
      uvwwrf = [ux,vx,wx];
    end
    u = uvwwrf(1);
    v = uvwwrf(2);
    w = uvwwrf(3);
    %plot3(X(idx),Y(jdx),u,'b*')
    fprintf(fid,'%f %f %f %f %f %f %f %d \n',X(idx),Y(jdx),Z,0,u,v,w,1);
  end
end
%view(-27,30)
fclose(fid);
