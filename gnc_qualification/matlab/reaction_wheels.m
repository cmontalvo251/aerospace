function reaction_wheels(t,Hx0,Hy0,Hz0,HxRW,HyRW,HzRW,Hdist,mission_duration)

%%%Compute Amount of Reaction Wheels Needed to detumble
px = (Hx0/HxRW)*100;
py = (Hy0/HyRW)*100;
pz = (Hz0/HzRW)*100;

disp(['Percentage of RW Needed to Detumble (%) = ',num2str(px),[' ' ...
		    ' ',num2str(py),' ',num2str(pz)]])

%%%Compute Number of Orbits before saturation
Nx0 = HxRW/Hdist;
Ny0 = HyRW/Hdist;
Nz0 = HzRW/Hdist;

disp(['Number of Orbits to Saturate (nd) = ',num2str(Nx0),[' ' ...
		    ' ',num2str(Ny0),' ',num2str(Nz0)]])

%%Time to Saturate
tx0 = Nx0*t(end);
ty0 = Ny0*t(end);
tz0 = Nz0*t(end);

disp(['Time to saturate (sec) = ',num2str(tx0),[' ' ...
		    ' ',num2str(ty0),' ',num2str(tz0)]])

%%%Compute number of desats per mission
tmin = min([tx0,ty0,tz0]);
ndesats = mission_duration*30.416*24*60*60/tmin;
disp(['Number of Desaturizations Per Mission (sec) = ',num2str(ndesats)])


