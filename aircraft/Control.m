function uk = Control(xk,tk,AMODE,uk)
global WPctr
%%%Extract States
x = xk(1);
z = xk(2);
theta = xk(3);
u = xk(4);
w = xk(5);
q = xk(6);

switch AMODE
    case 'MANUAL'
        dt = uk(1);
        de = uk(2);
    case 'AUTO'
        kp = -10;
        uc = 20;
        dt = kp*(u-uc);
        
        zc = -100;
        wc = 0;
        kp = 0.01;
        kd = 0.02;
        thetac = kp*(z-zc) + kd*(w-wc);
        
        if thetac > 30*pi/180
            thetac = 30*pi/180;
        end
        if thetac < -30*pi/180
            thetac = -30*pi/180;
        end
        %thetac = 20*pi/180;
        
        kp = 0.3;
        kd = 0.1;
        qc = 0;
        de = kp*(theta-thetac) + kd*(q-qc);
end

uk = [dt;de];

if uk(1) > 100
    uk(1) = 100;
end
if uk(1) < 0
    uk(1) = 0;
end
if uk(2) > 30*pi/180
    uk(2) = 30*pi/180;
end
if uk(2) < -30*pi/180
    uk(2) = -30*pi/180;
end




