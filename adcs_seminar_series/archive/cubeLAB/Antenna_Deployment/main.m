clear
clc
close all
%global td
td = 0.5;
w_init=[0; 0; 1];
t_interval=[0 1];
%timestep = 0.0001;

%%%Compute tout, wout and others is just td for now
%[tout, wout, others] = odeRK4(@derivs, t_interval, w_init, timestep, td,10);
[tout, wout] = ode45(@derivs, t_interval, w_init,[], td);

disp('Integration Complete')

dwdtout = zeros(3,length(tout));
H = zeros(length(tout),1);
Idotterm = dwdtout;
crossterm = dwdtout;

%%%%Plot the derivative as well and check for angular momentum
%%%Also plot effect of Idot*w and cross(w,I*w);
for idx = 1:length(tout)
    [dwdt,d] = derivs(tout(idx),wout(idx,:)',td);
    dwdtout(:,idx) = dwdt;
    
    theta = thetafunct(tout(idx),td);
    thetadot = thetadotfunct(tout(idx),td);
    
    I = computeInertia(theta);
    Idot = computeIdot(theta)*thetadot;
    
    H(idx) = norm(I*wout(idx,:)');
    
    Idotterm(:,idx) = -inv(I)*Idot*wout(idx,:)';
    crossterm(:,idx) = -inv(I)*cross(wout(idx,:)',I*wout(idx,:)');
end

disp('Plotting')

figure();
p1 = plot(tout,wout*180/pi);
grid on
xlabel('Time (s)')
ylabel('Angular velocity (deg/s)')
legend('Roll Rate','Pitch Rate','Yaw Rate')

figure();
y1 = thetafunct(tout,td)';
y2 = thetadotfunct(tout,td)';
plot(tout,[y1;y2], 'LineWidth',4)
xlabel('Time (s)')
legend('Theta (deg)','Thetadot (deg/s)')
grid on

figure();
p1 = plot(tout,dwdtout*180/pi);
grid on
xlabel('Time (s)')
ylabel('Angular Acceleration (deg/s^2)')
%legend('Roll','Pitch','Yaw')

axis = {'Wx','Wy','Wz'};

for idx = 1:3
    figure()
    plot(tout,dwdtout(idx,:)*180/pi,'b-')
    hold on
    plot(tout,Idotterm(idx,:)*180/pi,'r--');
    plot(tout,crossterm(idx,:)*180/pi,'g--');
    %plot(tout,(Idotterm(idx,:)+crossterm(idx,:))*180/pi,'k--')
    legend('Total','Idot Term','Cross Product Term')
    title(axis{idx})
    xlabel('Time (s)')
    ylabel('Angular Acceleration (deg/s^2)')
end

figure()
plot(tout,H);
grid on
xlabel('Time (s)')
ylabel('Angular Momentum (N-m-s)')