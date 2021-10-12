function Hdist = disturbance(t,Maero,Mgrav,Mrad,Mdipole);

Mtotal = Maero + Mgrav + Mrad + Mdipole;

%%%Plot all disturbance torques
figure()
set(gcf,'color','white')
hold on
plot(t,log10(Maero),'LineWidth',2)
plot(t,log10(Mgrav),'LineWidth',2)
plot(t,log10(Mrad),'LineWidth',2)
plot(t,log10(Mdipole),'LineWidth',2)
plot(t,log10(Mtotal),'LineWidth',2)
legend('Aerodynamics','Gravity','Radiation Pressure','Dipole Moment','Total')
xlabel('Time (sec)')
ylabel('Logarithm (Base10) of Disturbance Torques (N-m)')
grid on


%%%%Integrate Disturbance Momentum
Hdist = integrate_curve(t,Mtotal);