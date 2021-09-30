function M = BLIMPAGE_rev1()
clc
close all
disp('THE BLIMPAGE')

global WP xcommand ycommand Num
WP = [1;1;1;1;1];
% b1xwp = [-2];
% b2xwp = [4];
% b3xwp = [10];
% b4xwp = [16];
% b5xwp = [22];
% b1ywp = [5-.4827];
bxwp = 11;
b1ywp = [2.5];
% y2wp = [2.5;2.5;2.5;2.5];
% y3wp = [4;4;4;4];

xcommand = bxwp*[1 1 1 1 1];
% xcommand = [b1xwp b2xwp b3xwp b4xwp b5xwp];
ycommand = [b1ywp b1ywp b1ywp b1ywp b1ywp];

tf = 100;                                % Define TIME Variables (final & initial)
ti = 0;
timestep = 0.01;
time = ti:timestep:tf;

Num = 5;                                  % Number of blimps
State = zeros(6,Num);                     % Initial State Vector: [x;y;psi;u;v;r]
StateOUT = zeros(6,length(time),Num);

T1 = zeros(length(State),Num);
T2 = zeros(length(State),Num);
dpsi = zeros(1,Num);
% StateDot = Derivatives(State,T1,T2);   % Initialize StateDot Vector

%%%%%% RK4 %%%%%%
% State(:,1) = [-10;1.4;0*pi/180;0;0;0];         % Hallway
% State(:,2) = [-10;2.5;180*pi/180;0;0;0];
% State(:,3) = [-10;4;10*pi/180;0;0;0];

State(:,1) = [-9.5;4.5;0*pi/180;0;0;0];                  %the 5th gate
if Num > 1
  State(:,2) = [-9;3.5;0*pi/180;0;0;0];
  State(:,3) = [-9.5;2.5;0*pi/180;0;0;0];
  State(:,4) = [-9;1.5;0*pi/180;0;0;0];
  State(:,5) = [-9.5;.5;0;0;0;0];
end

phi = zeros(6,Num);
tnext = 0;
for idx = 1:length(time)            % For each time step
    phi = 0*phi;
    if 100*time(idx)/tf > tnext
      tnext = tnext + 10;
      disp(['Simulation ',num2str(time(idx)/tf*100),' Percent Complete'])    
    end
    for B = 1:Num                     % For each blimp
        StateOUT(:,idx,B) = State(:,B);                     % 3D matrix
        [u1,u2,delpsi] = Control(State(:,B),B,time(idx));

        k1 = Derivatives(State,u1,u2,B,time(idx));                  % Compute k1
        k2 = Derivatives(State + k1*timestep/2,u1,u2,B,time(idx));  % Compute k2
        k3 = Derivatives(State + k2*timestep/2,u1,u2,B,time(idx));  % Compute k3
        k4 = Derivatives(State + k3*timestep,u1,u2,B,time(idx));    % Compute k4
        phi = phi + (1/6)*(k1 + 2*k2 + 2*k3 + k4);        % Compute phi

        T1(idx,B) = u1;
        T2(idx,B) = u2;
        dpsi(idx,B) = delpsi;
    end 
    State = State + phi*timestep;       % Apply to all blimps at once
end

% figure
% set(gca,'FontSize',18);
% xlabel('X(meters)'); ylabel('Y(meters)');
% title('Airships Top Down View'); grid on; hold on;
% axis ([0 10 0 10])
% axis ([-10 10 -6 12])
% vidObj = VideoWriter('Random.avi');      % New animation file name
% open(vidObj);                            % Open file
% colorV = [[1 0 0];[0 1 0];[0 0 1];[.8 0 .8];[.8 .8 0]]; 
% 
% for idx = 1:15:length(time)              % Re-play animation
%     cla;
%     
%     s = 0.4827;                         % Heading vector size
%     plot(-10:.1:10,5,'k.','MarkerSize',12)
%     plot(-10:.1:10,0,'k.','MarkerSize',12)
%     plot(xcommand,ycommand,'k*','MarkerSize',10)
%     
%     for B = 1:Num
%         State(:,B) = StateOUT(:,idx,B);
%         xBd = State(1,B);
%         yBd = State(2,B);
%         psiB = State(3,B);
%         headingB = [cos(psiB);sin(psiB)];
%         circle(xBd,yBd,s,colorV(B,:));
%         plot([xBd xBd+headingB(1)*s],[yBd yBd+headingB(2)*s],'k')
%     end
%     drawnow
%     currFrame = getframe;
%     writeVideo(vidObj,currFrame);
% end
% close(vidObj);      % Close the file.

fig = figure();
set(fig,'color','white');
hold on
colors = {'b-','r-','g-','m-','c-'};
colors2 = {'b--','r--','g--','m--','c--'};
for idx = 1:Num
  x = StateOUT(1,:,idx);
  y = StateOUT(2,:,idx);
  plot(x,y,colors{idx},'LineWidth',2);
end  
set(gca,'FontSize',18);
xlabel('X (m)'); ylabel('Y (m)');
legend('B1','B2','B3','B4','B5');
title('Flight Path'); grid on;
plot(-10:.1:30,5,'k.','MarkerSize',12)
plot(-10:.1:30,0,'k.','MarkerSize',12)
plot(xcommand,ycommand,'k*','MarkerSize',10)

figu = figure; %Velocity plots
set(figu,'color','white')
hold on
for idx = 1:Num  
  u = StateOUT(4,:,idx);
  v = StateOUT(5,:,idx);
  Vtotal = sqrt(u.^2+v.^2);
  plot(time,Vtotal,colors{idx},'LineWidth',2)
  %plot(time,u,colors{idx},'LineWidth',2);
  %plot(time,v,colors2{idx},'LineWidth',2)
end
xlabel('Time (sec)')
ylabel('Velocity (m/s)')
legend('B1','B2','B3','B4','B5');


% figT = figure();
% set(figT,'color','white')
% hold on
% for idx = 1:Num
%     T1B = T1(:,idx);
%     T2B = T2(:,idx);
%     plot(time,T1B,colors{idx},'LineWidth',2)
%     plot(time,T2B,colors2{idx},'LineWidth',2)
% end
% xlabel('Time (sec)')
% ylabel('Thrust (N)')
% 
% figr = figure; %Velocity plots
% set(figr,'color','white')
% hold on
% for idx = 1:Num  
%   r = StateOUT(6,:,idx);
%   plot(time,r,colors{idx},'LineWidth',2);
% end
% xlabel('Time (sec)')
% ylabel('Yaw Rate (rad/s)')
% 
% figD = figure; %Velocity plots
% set(figD,'color','white')
% hold on
% for idx = 1:Num  
%   delpsi = dpsi(:,idx);
%   plot(time,delpsi,colors{idx},'LineWidth',2);
% end
% xlabel('Time (sec)')
% ylabel('Delta Psi (deg)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% legend('Wall 1','Wall 2', 'Waypoint')

% fig = figure();
% set(fig,'color','white');
% plot(time,x1,'c--','LineWidth',2);
% set(gca,'FontSize',18);
% xlabel('Time(sec)'); ylabel('X(m)');
% title('Inertial X1 Position'); grid on;

% fig = figure();
% set(fig,'color','white');
% plot(time,y1,'c--','LineWidth',2);
% set(gca,'FontSize',18);
% xlabel('Time(sec)'); ylabel('Y(m)');
% title('Inertial Y1 Position'); grid on;

% fig = figure();
% set(fig,'color','white');
% plot(time,psi*(180/pi),'m-','LineWidth',2);
% set(gca,'FontSize',18);
% xlabel('Time(sec)'); ylabel('Psi(degrees)');
% title('Yaw(psi) Angle'); grid on;
% 
% fig = figure();
% set(fig,'color','white');
% plot(time,u,'b-','LineWidth',2);
% set(gca,'FontSize',18);
% xlabel('Time(sec)'); ylabel('Velocity in x (m/s)');
% title('X Component of Body Frame Velocity'); grid on;
% 
% fig = figure();
% set(fig,'color','white');
% plot(time,v,'b-','LineWidth',2);
% set(gca,'FontSize',18);
% xlabel('Time(sec)'); ylabel('v');
% title('Y Component of Body Frame Velocity'); grid on;
% 
% fig = figure();
% set(fig,'color','white');
% plot(time,r,'g-','LineWidth',2);
% set(gca,'FontSize',18);
% xlabel('Time(sec)'); ylabel('Angular Velocity(radians/sec)');
% title('Angular Velocity in Body Frame'); grid on;

% fig = figure();
% set(fig,'color','white');
% plot(time,T1,'g--','LineWidth',2);
% set(gca,'FontSize',18);
% xlabel('time(sec)'); grid on;

% hold on;
% plot(time,T2,'b-','LineWidth',2);
% set(gca,'FontSize',18);
% legend('Thrust Left(T1)','Thrust Right(T2)');
% xlabel('time(sec)'); ylabel('Thrust(Newtons)');
% title('Thrust Channels'); grid on;

function StateDot = Derivatives(State,T1,T2,B,time)    
global Num

x = State(1,B);
y = State(2,B);
psi = State(3,B);
u = State(4,B);
v = State(5,B);
r = State(6,B);
V = sqrt(u^2+v^2);
Rbi = [cos(psi) -sin(psi);sin(psi) cos(psi)]; %Rotate body to inertial frames

b = (38/12)/3.28;           % Diameter another direction
radius = b/2;
d = (5/12)/3.28;            % Moment arm
rho = 1.225;
% th=(1+15/16)/3.28;          % Diameter some direction
% R = rho*.0555*b^3*th^2;
% vol = (4/3)*pi*b/2*b/2*th/2;
% m = .1155;

vol = (4/3)*pi*(b/2)^3;
mstructure = 32/1000; %%%%mass of structure
mgondola = 68/1000; %%%mass of gondola and batteries
density_helium = 0.18; %%%%kg/m^3 of helium
m = density_helium*vol + mstructure + mgondola; %%%total mass
I = (2/5)*m*(b/2)^2;     % MOI for solid sphere
k = (8/3)*(b/2)^3;      %apparent mass
k1 = k;
k2 = k;
k3 = 0;

% k = 20;             % Spring Constant
% C = 10;    % Cnr = .029;                  % Coef of drag from yaw moment
% 
% V0 = 2.5;
% rhat = r*b/(2*V0);
% Cnr = Cnr*rhat;
         % Damping Constant
k = 5.95;
C = 5.9;
% Cd = .05188;
Cd = 0.125;
Cnr = 0.023;

xydot = Rbi*[u;v];
psidot = r;

Xa = -(1/2)*rho*V*vol^(2/3)*u*Cd;
Ya = -(1/2)*rho*V*vol^(2/3)*v*Cd;
% Na = -(1/2)*rho*V^2*vol*Cnr;
Na = -(1/2)*rho*vol^(5/3)*abs(psidot)*(Cnr*psidot);

%%% Wall Collision %%%
Fwcont = [0;0];
Fwshr = 0;

if y+radius >= 5 || y-radius <= 0
    ywall1 = [x;5];
    ywall2 = [x;0];
    r1 = [x;y];
    
    r1w1_I = r1 - ywall1;
    r1w1_mag = sqrt(r1w1_I(1)^2+r1w1_I(2)^2);
    
    r1w2_I = r1 - ywall2;
    r1w2_mag = sqrt(r1w2_I(1)^2+r1w2_I(2)^2);
    
% if r1w1_mag <= radius || r1w2_mag <= radius 
    if y+radius >= 5
        ep_hat = [0;-1];
        et_hat = [-1;0];
        rcp1_I = (r1w1_mag)*ep_hat;
        dely = (radius - r1w1_mag);
    elseif y-radius <= 0
        ep_hat = [0;1];
        et_hat = [1;0];
        rcp1_I = (r1w2_mag)*ep_hat;
        dely = (radius - r1w2_mag);
    end
        rcp1_B1 = Rbi'*rcp1_I;
        
        VB1_B1 = [u;v];
        VB1_I = Rbi*VB1_B1;
        Vp1_B1 = [-rcp1_B1(2)*r;rcp1_B1(1)*r];
        Vp1_I = VB1_I - Rbi*Vp1_B1; %% Possible error
        delVp_I = -Vp1_I;
        delydot = dot(delVp_I,ep_hat);
        delxdot = dot(delVp_I,et_hat);
        
        Fwall = -k*dely-C*abs(delydot);      % In the y Inertial direction
        Fwshear = -C*delxdot;
        Fwcontact_I = [ep_hat,et_hat]*[Fwall;Fwshear];
        Fwcont = Rbi'*Fwcontact_I;
        Fwshr = Fwshear*(radius - dely);
end

%%% Blimp Collision %%%
r1 = [x;y];
Fcont = [0;0];
Fspin = 0;

for n = 1:Num
    Fcontact_B = [0;0];
    Fs = 0;
    if n ~= B
        Rbin = [cos(State(3,n)) -sin(State(3,n));sin(State(3,n)) cos(State(3,n))];
        rn = [State(1,n);State(2,n)];
        r1n_I = r1-rn;
        r1n_mag = sqrt(r1n_I(1)^2+r1n_I(2)^2);

        if r1n_mag <= 2*radius
            ep_hat = r1n_I/r1n_mag;         % ep and et are in the inertial frame
            et_vec = [0;0];
            if abs(ep_hat(1)) < 1e-4
                et_vec(1) = 1;
                et_vec(2) = -ep_hat(1)*et_vec(1)/ep_hat(2);
            else
                et_vec(2) = 1;
                et_vec(1) = -ep_hat(2)*et_vec(2)/ep_hat(1);
            end
            et_mag = sqrt(et_vec(1)^2+et_vec(2)^2);
            et_hat = et_vec/et_mag;
            D = ep_hat(1)*et_hat(2)-ep_hat(2)*et_hat(1);
            if D<0
                et_hat = -et_hat;
            end
            
            rcp1_I = (r1n_mag/2)*ep_hat;
            rcpn_I = -rcp1_I;
            delx = (2*radius - r1n_mag)/2;

            VB1_B1 = [u;v];
            VBn_Bn = [State(4,n);State(5,n)];

            VB1_I = Rbi*VB1_B1;
            VBn_I = Rbin*VBn_Bn;

            rcp1_B1 = -Rbi'*rcp1_I;
            rcpn_Bn = Rbin'*rcpn_I;

            Vp1_B1 = [-rcp1_B1(2)*r;rcp1_B1(1)*r];
            Vpn_Bn = [-rcpn_Bn(2)*State(6,n);rcpn_Bn(2)*State(6,n)];

            Vp1_I = VB1_I + Rbi*Vp1_B1;
            Vpn_I = VBn_I + Rbin*Vpn_Bn;  %Vpn_I = [cos(State(3,n)) sin(State(3,n));-sin(State(3,n)) cos(State(3,n))]*Vpn_Bn;
            delVp_I = Vpn_I-Vp1_I;  

            delxdot = dot(delVp_I,ep_hat);
            delydot = dot(delVp_I,et_hat);
            Fn = -k*delx-C*abs(delxdot);
            Ft = -C*delydot;
            Fcontact_I = [ep_hat,et_hat]*[Fn;Ft];
            Fcontact_B = Rbi'*Fcontact_I;
            Fs = Ft*r1n_mag/2;
        end
    end
    Fcont = Fcont + Fcontact_B;
    Fspin = Fspin + Fs;
end

X = Xa + (T1+T2) - Fcont(1) - Fwcont(1);
Y = Ya - Fcont(2) - Fwcont(2);
N = Na + (T2*d-T1*d) + Fspin - Fwshr;
udot = X/(m + k1) + r*v;
vdot = Y/(m + k2) - r*u;
rdot = N/(I + k3);

StateDot = zeros(6,Num);
StateDot(:,B) = [xydot;psidot;udot;vdot;rdot];
% fprintf('t = %f, xdot = %f, udot = %f\n',time,xdot,udot)


function [T1,T2,delpsi] = Control(State,B,time)
global WP xcommand ycommand
x = State(1);
y = State(2);
psi = State(3);
u = State(4);
v = State(5);
psidot = State(6);
d = (5/12)/3.28;

ucontrol = 1;           % Body velocity should be 2.5 m/s
psidotcontrol = 0;

dx = xcommand(WP(B),B)-x;
dy = ycommand(WP(B),B)-y;
psic = atan2(dy,dx);

delWP = (sqrt(dx^2 + dy^2));
if delWP < 1                       % If close to waypoint
    WP(B) = WP(B) + 1;                % Go to next waypoin
    if WP(B) > length(xcommand(:,1))    % After all WPs, start over
        WP(B) = 1;
    end
end

new_controls = 1;
if new_controls == 1
    kp = 0.008; 
    kpN = -0.009;
    kdN = -0.014;
    sig = 0.045;
    if time > 25
        sig = 0.045;
    end
    delpsi = atan2((sin(psi)*cos(psic)-cos(psi)*sin(psic)),(cos(psi)*cos(psic)+sin(psi)*sin(psic)));
    TN = kpN*(delpsi) + kdN*(psidot - psidotcontrol);
    T1 = (pi-abs(delpsi))/pi*kp*(ucontrol-u) - TN/d;
    T2 = (pi-abs(delpsi))/pi*kp*(ucontrol-u) + TN/d;

    if T1 > sig
       T1 = sig;
    end
    if T1 < 0
        T1 = 0;
    end
    if T2 > sig
       T2 = sig;
    end  
    if T2 < 0
        T2 = 0;
    end
else %old controls
    kp = .05;
    kd = 1;
    kpN = .5;
    kdN = .07;
    sig = 0.08;

    delpsi = -atan2((sin(psi)*cos(psic)-cos(psi)*sin(psic)),(cos(psi)*cos(psic)+sin(psi)*sin(psic)));
    TN = kpN*(delpsi) + kdN*(psidotcontrol-psidot);
    T1 = (kp/2)*(delWP) + (kd/2)*(ucontrol-u) - TN/(2*d);
    T2 = (kp/2)*(delWP) + (kd/2)*(ucontrol-u) + TN/(2*d);
    
    if abs(T1)>sig
       T1 = sign(T1)*sig;
    end
    if abs(T2)>sig
       T2 = sign(T2)*sig;
    end  
end