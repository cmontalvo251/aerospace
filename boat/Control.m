function uk = Control(xk,tk)
global WPctr
%%%Extract States
x = xk(1);
y = xk(2);
psi = xk(3);
u = xk(4);
v = xk(5);
r = xk(6);

%%%Compute Control
uc = 5;
kp = 70;
dt = kp*(uc-u);
if dt > 100
    dt = 100;
end
if dt < 0
    dt = 0;
end
%%%%Rudder Angle Command
kpsip = 2.5;
kpsid = 1.25;
yc = [0 0 20 20]*2;
xc = [0 20 20 0]*2;
yci = yc(WPctr);
xci = xc(WPctr);

psic = atan2(yc(WPctr)-y,xc(WPctr)-x);
delpsi = -atan2(sin(psi)*cos(psic)-cos(psi)*sin(psic),cos(psi)*cos(psic)+sin(psi)*sin(psic));

if sqrt((x-xc(WPctr))^2+(y-yc(WPctr))^2) < 10
    WPctr = WPctr + 1;
    if WPctr > length(xc)
        WPctr = 1;
    end
end


rc = 0;
%delpsi = 0; %drive straight
dr = kpsip*delpsi + kpsid*(rc-r);
if dr > 30*pi/180;
    dr = 30*pi/180;
end
if dr < -30*pi/180;
    dr = -30*pi/180;
end

uk = [dt;dr];
%ueq = [50.52;0];




