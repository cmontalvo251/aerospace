clear
clc
close all
 
dataloc = '/Users/rockwellgarrido/Documents/MATLAB/MAAWM-master-58f4f52754d6704e3c9db7e564bdc2e3fdfbab7e/Frequency_Analysis/';
addpath 'supplemental_routines'
UVWfrontend
%% U
Nexp = 100; % number of steps
yvals = linspace(-500,500,Nexp); % varibale that moves from -500 to 500 with N number of steps
Umat_recreate = []; % prealocated to be a matrix
am_mat = []; % prealocated to be a matrix
bm_mat = []; % prealocated to be a matrix
for y = yvals % for loop to run for the values of 'yvals' and then stop
    wmax = 20; % maximum fundemental frequency 
    x = xcoord; % coordinate from north-south along windfield
    xexp = linspace(-500,500,Nexp); % vector variable that moves from -500 to 500 with N number of steps
    yexp = xexp;
    dx = x(2)-x(1);  % change in x coordinate for step size
    u = 0*xcoord; %prealocates u to be a row vector of zeros the length of the x coordinate
    for idx = 1:length(xcoord) % nested for loop running from 1 to length of the x coordinate and stopping
        uvw = uvwout(xcoord(idx),y,200,0,dataloc,0); % takes 200 sets of uvwout(xcoord,y) and saves them 
        u(idx) = uvw(1);  % take from the 200 sets above pulls out the first values in the column vector and set to U
    end % ends nested loop after all 200 values of U are saved
    a0u = mean(u); % take the average of the U's to find the a0 amplitude
    us = u - a0u; % QUESTION!!!
    mvec = 1:wmax;
    am = 0*mvec;
    bm = 0*mvec;
    T = x*pi/500;
    urecreate = 0*xexp;
    for m = mvec
        %%%Trapezoidal rule
        for idx = 1:length(x)-1
            a1 = us(idx)*cos(m*T(idx));
            a2 = us(idx+1)*cos(m*T(idx+1));
            am(m) = am(m) + (1/500)*(0.5*(a1+a2)*dx);
 
            b1 = us(idx).*sin(m.*T(idx));
            b2 = us(idx+1).*sin(m.*T(idx+1));
            bm(m) = bm(m) + (1/500)*(0.5*(b1+b2)*dx);
        end
        %%%Recreate the waveform
        urecreate = urecreate + am(m)*cos(m*pi/500*xexp)+bm(m)*sin(m*pi/500*xexp);
    end
    urecreate = urecreate + a0u;
    Umat_recreate = [Umat_recreate;urecreate];
    am_mat = [am_mat;am];
    bm_mat = [bm_mat;bm];
end
Umat = [];
for y = ycoord
    u = 0*xcoord;
    for idx = 1:length(xcoord)
        uvw = uvwout(xcoord(idx),y,200,0,dataloc,0);
        u(idx) = uvw(1);   
    end
    Umat = [Umat;u];
end
[xx,yy] = meshgrid(xcoord,ycoord);
[xxexp,yyexp] = meshgrid(xexp,yexp);
set(figure,'color','white')
mesh(xx,yy,Umat)
hold on
surf(xxexp,yyexp,Umat_recreate)
grid on
title('Windfield in U (Gathered Data vs. Fourier Recreation)','Fontsize',18)
xlabel('Xcoord (m)')
ylabel('Ycoord (m)')
zlabel('U (m/s)')
legend('Gathered Data','Fourier Recreation')

freq = mvec*pi/500;
[ff,~] = meshgrid(freq,yexp);
set(figure,'color','white')
mesh(ff,yexp,am_mat)
title('a_{n} Frequency of Windfield in U','Fontsize',18)
xlabel('Frequency (rad/m)')
ylabel('Ycoord (m)')
zlabel('Amplitue (m/s)')

set(figure,'color','white')
mesh(ff,yexp,bm_mat)
title('b_{n} Frequency of Winfield in U','Fontsize',18)
xlabel('Frequency (rad/m)')
ylabel('Ycoord (m)')
zlabel('Amplitue (m/s)')

set(figure,'color','white')
mesh(ff,yexp,sqrt(bm_mat.^2+am_mat.^2))
title('Magnitude of a_{n} and b_{n} frequencies','Fontsize',18)
xlabel('Frequency (rad/m)')
ylabel('Ycoord (m)')
zlabel('Amplitue (m/s)')
figure()
plot(freq,abs(am),'b.','markersize',16);xlabel('Frequency (rad/meter)','FontSize',13);
ylabel('Amplitudec (m/s)','FontSize',13);grid on
figure()
plot(freq,abs(bm),'b.','markersize',16);xlabel('Frequency (rad/meter)','FontSize',13);
ylabel('Amplitude (m/s)','FontSize',13);grid on
figure()
plot(x,u,'b*',xexp,urecreate,'r--','linewidth',2);xlabel('X (meters)','Fontsize',13);
ylabel('U (m/s)','Fontsize',13);legend('Wind Data','Fourier Recreation');grid on

%% V

Nexp = 100;
yvals = linspace(-500,500,Nexp);
Vmat_recreate = [];
am_mat = [];
bm_mat = [];
for y = yvals
    wmax = 20;
    x = xcoord;
    xexp = linspace(-500,500,Nexp);
    yexp = xexp;
    dx = x(2)-x(1);
    v = 0*xcoord;
    for idx = 1:length(xcoord)
        uvw = uvwout(xcoord(idx),y,200,0,dataloc,0);
        v(idx) = uvw(2);
    end
    a0v = mean(v);
    vs = v - a0v;
    mvec = 1:wmax;
    am = 0*mvec;
    bm = 0*mvec;
    T = x*pi/500;
    vrecreate = 0*xexp;
    for m = mvec
        %%%Trapezoidal rule
        for idx = 1:length(x)-1
 
            a1 = vs(idx)*cos(m*T(idx));
            a2 = vs(idx+1)*cos(m*T(idx+1));
            am(m) = am(m) + (1/500)*(0.5*(a1+a2)*dx);
 
            b1 = vs(idx)*sin(m*T(idx));
            b2 = vs(idx+1)*sin(m*T(idx+1));
            bm(m) = bm(m) + (1/500)*(0.5*(b1+b2)*dx);
        end
        %%%Recreate the waveform
        vrecreate = vrecreate + am(m)*cos(m*pi/500*xexp)+bm(m)*sin(m*pi/500*xexp);
    end
    vrecreate = vrecreate + a0v;
    Vmat_recreate = [Vmat_recreate;vrecreate];
    am_mat = [am_mat;am];
    bm_mat = [bm_mat;bm];
end
Vmat = [];
for y = ycoord
    v = 0*xcoord;
    for idx = 1:length(xcoord)
        uvw = uvwout(xcoord(idx),y,200,0,dataloc,0);
        v(idx) = uvw(2); 
    end
    Vmat = [Vmat;v];
end
[xx,yy] = meshgrid(xcoord,ycoord);
[xxexp,yyexp] = meshgrid(xexp,yexp);
set(figure,'color','white')
mesh(xx,yy,Vmat)
hold on
surf(xxexp,yyexp,Vmat_recreate)
grid on
title('Windfield in V (Gathered Data vs. Fourier Recreation)','Fontsize',18)
xlabel('Xcoord (m)')
ylabel('Ycoord (m)')
zlabel('V (m/s)')
legend('Gathered Data','Fourier Recreation')

freq = mvec*pi/500;
[ff,~] = meshgrid(freq,yexp);
set(figure,'color','white')
mesh(ff,yexp,am_mat)
title('a_{n} Frequency of Windfield in V','Fontsize',18)
xlabel('Frequency (rad/m)')
ylabel('Ycoord (m)')
zlabel('Amplitue (m/s)')

set(figure,'color','white')
mesh(ff,yexp,bm_mat)
title('b_{n} Frequency of Winfield in V','Fontsize',18)
xlabel('Frequency (rad/m)')
ylabel('Ycoord (m)')
zlabel('Amplitue (m/s)')

set(figure,'color','white')
mesh(ff,yexp,sqrt(bm_mat.^2+am_mat.^2))
title('Magnitude of a_{n} and b_{n} Frequencies','Fontsize',18)
xlabel('Frequency (rad/m)')
ylabel('Ycoord (m)')
zlabel('Amplitue (m/s)')

%% W

Nexp = 100;
yvals = linspace(-500,500,Nexp);
Wmat_recreate = [];
am_mat = [];
bm_mat = [];
for y = yvals
    wmax = 20;
    x = xcoord;
    xexp = linspace(-500,500,Nexp);
    yexp = xexp;
    dx = x(2)-x(1);
    w = 0*xcoord;
    for idx = 1:length(xcoord)
        uvw = uvwout(xcoord(idx),y,200,0,dataloc,0);
        w(idx) = uvw(3);
    end
    a0w = mean(w);
    ws = w - a0w;
    mvec = 1:wmax;
    am = 0*mvec;
    bm = 0*mvec;
    T = x*pi/500;
    wrecreate = 0*xexp;
    for m = mvec
        %%%Trapezoidal rule
        for idx = 1:length(x)-1
 
            a1 = ws(idx)*cos(m*T(idx));
            a2 = ws(idx+1)*cos(m*T(idx+1));
            am(m) = am(m) + (1/500)*(0.5*(a1+a2)*dx);
 
            b1 = ws(idx)*sin(m*T(idx));
            b2 = ws(idx+1)*sin(m*T(idx+1));
            bm(m) = bm(m) + (1/500)*(0.5*(b1+b2)*dx);
        end
        %%%Recreate the waveform
        wrecreate = wrecreate + am(m)*cos(m*pi/500*xexp)+bm(m)*sin(m*pi/500*xexp);
    end
    wrecreate = wrecreate + a0w;
    Wmat_recreate = [Wmat_recreate;wrecreate];
    am_mat = [am_mat;am];
    bm_mat = [bm_mat;bm];
end
Wmat = [];
for y = ycoord
    w = 0*xcoord;
    for idx = 1:length(xcoord)
        uvw = uvwout(xcoord(idx),y,200,0,dataloc,0);
        w(idx) = uvw(3); 
    end
    Wmat = [Wmat;w];
end
[xx,yy] = meshgrid(xcoord,ycoord);
[xxexp,yyexp] = meshgrid(xexp,yexp);
set(figure,'color','white')
mesh(xx,yy,Wmat)
hold on
surf(xxexp,yyexp,Wmat_recreate)
grid on
title('Windfield in W (Gathered Data vs. Fourier Recreation)','Fontsize',18)
xlabel('Xcoord (m)')
ylabel('Ycoord (m)')
zlabel('W (m/s)')
legend('Gathered Data','Fourier Recreation')

freq = mvec*pi/500;
[ff,yy] = meshgrid(freq,yexp);
set(figure,'color','white')
mesh(ff,yexp,am_mat)
title('a_{n} Frequency of Windfield in W','Fontsize',18)
xlabel('Frequency (rad/m)')
ylabel('Ycoord (m)')
zlabel('Amplitue (m/s)')

set(figure,'color','white')
mesh(ff,yexp,bm_mat)
title('b_{n} Frequency of Winfield in W','Fontsize',18)
xlabel('Frequency (rad/m)')
ylabel('Ycoord (m)')
zlabel('Amplitue (m/s)')

set(figure,'color','white')
mesh(ff,yexp,sqrt(bm_mat.^2+am_mat.^2))
title('Magnitude of a_{n} and b_{n} Frequencies','Fontsize',18)
xlabel('Frequency (rad/m)')
ylabel('Ycoord (m)')
zlabel('Amplitue (m/s)')

%% ERROR

