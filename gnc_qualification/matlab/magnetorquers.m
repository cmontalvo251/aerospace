function magnetorquers(t,magnetic_moment,Btotal,Hdist,Hx0,Hy0,Hz0,HxRW,HyRW,HzRW)

%%%%First compute torque across an entire orbit
Mmag = magnetic_moment.*Btotal*1e-9;

figure()
set(gcf,'color','white')
plot(t,Mmag,'LineWidth',2)
xlabel('Time (sec)')
ylabel('Magnetorquer Torque (N-m)')
grid on

%%%THen compute magnetorquer momentum
Hmag = integrate_curve(t,Mmag);
disp(['Magnetorquer Momentum Per Orbit (N-m-s) = ',num2str(Hmag)])

%%THen total magnetorquer effectiveness
Htotal = Hmag - Hdist;
disp(['Magnetorquer Momentum with Disturbance Per Orbit (N-m-s) = ',num2str(Htotal)])

%%%Number of ORbits to detumble
Nx0 = (Hx0/Htotal);
Ny0 = (Hy0/Htotal);
Nz0 = (Hz0/Htotal);
disp(['Number of Orbits to detumble (nd) = ',num2str(Nx0),[' ' ...
		    ' ',num2str(Ny0),' ',num2str(Nz0)]])

%%Time to Detumble
tx0 = Nx0*t(end);
ty0 = Ny0*t(end);
tz0 = Nz0*t(end);

disp(['Time to detumble (sec) = ',num2str(tx0),[' ' ...
		    ' ',num2str(ty0),' ',num2str(tz0)]])

%%%Number of Orbits to desaturate Reaction Wheels
NxRW = (HxRW/Htotal);
NyRW = (HyRW/Htotal);
NzRW = (HzRW/Htotal);
disp(['Number of Orbits to desaturate (nd) = ',num2str(NxRW),[' ' ...
		    ' ',num2str(NyRW),' ',num2str(NzRW)]])

%%Time to Detumble
txRW = NxRW*t(end);
tyRW = NyRW*t(end);
tzRW = NzRW*t(end);

disp(['Time to desaturate (sec) = ',num2str(txRW),[' ' ...
		    ' ',num2str(tyRW),' ',num2str(tzRW)]])

