function Maero = aerodynamics(t,r,v,max_area,max_moment_arm,rhosl,CD)
constants

%%%Compute Density
rho = rhosl*exp(-scale_height*(r-REarth)/1000);
%%Plot Density
figure()
plot(t,rho,'LineWidth',2)
xlabel('Time (sec)')
ylabel('Density of Atmosphere (kg/m^3)')
set(gcf,'color','white')
grid on

%%%Compute Aero Force
Faero = 0.5*rho.*(v.^2)*max_area*CD;
%%%Aero Moment
Maero = Faero*max_moment_arm;	