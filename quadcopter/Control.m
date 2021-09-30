function uk = Control(xk,tk,QMODE,uk)
global WPctr
%%%Extract States
x = xk(1);
z = xk(2);
phi = xk(3);
xdot = xk(4);
zdot = xk(5);
p = xk(6);

switch QMODE
    case 'RATE_HOLD'
        r1 = uk(1);
        r2 = uk(2);
    case 'RATE'
        r1 = uk(1);
        r2 = uk(2);
        ravg = (r1+r2)/2;
        r1 = ravg;
        r2 = ravg;
    case 'STD'
        %%%Compute Control
        phic = uk(3);
        r0 = uk(4);
        m = 1.2; %%%kg
        H = (2/12)/3.28;
        L = 1/3.28;
        Iyy = m/12*(H^2+L^2);
        d = L/2;
        zed = 0.5594;
        wn = 10.2;
        %kp = -Iyy*wn^2/(2*d);
        kp = -1;
        kd = -0.5;
        %kd = -zed*wn*Iyy/d;
        pc = 0;
        rc = kp*(phi-phic) + kd*(p-pc);
        r1 = r0 + rc;
        r2 = r0 - rc;
        %r1 = 0;
        %r2 = 0;
    case 'AUTO'
        %%%Waypoint Controller
        xc = [5];
        zc = [-1];
        if sqrt((x-xc(WPctr))^2+(z-zc(WPctr))^2) < 1
            WPctr = WPctr + 1;
            if WPctr > length(xc)
                WPctr = 1;
            end
        end
        %%%Altitude Hold
        kp = 100; 
        kd = 50;
        zdotc = 0;
        rALT = kp*(z-zc(WPctr)) + kd*(zdot-zdotc);
        r1 = rALT;
        r2 = rALT;
        %%%%%Lateral Control
        xdotc = 0;
        kpx = -0.01;
        kdx = -0.02;
        phic = kpx*(x-xc(WPctr)) + kdx*(xdot-xdotc);
        if abs(phic) > 45*pi/180
            phic = sign(phic)*45*pi/180;
        end
        %%%Roll Controller
        kp = -1;
        kd = -0.5;
        pc = 0;
        rc = kp*(phi-phic) + kd*(p-pc);
        r1 = r1 + rc;
        r2 = r2 - rc;
end

uk = [r1;r2;uk(3);uk(4)];

if uk(1) > 100
   uk(1) = 100;
end
if uk(2) > 100;
   uk(2) = 100;
end
if uk(4) > 100
    uk(4) = 100;
end
if uk(4) < 0
    uk(4) = 0;
end
if uk(1) < 0
   uk(1) = 0;
end
if uk(2) < 0
   uk(2) = 0;
end
%if uk(3) > pi/3
%    uk(3) = pi/3;
%end
%if uk(3) < -pi/3
%    uk(3) = -pi/2;
%end


    




