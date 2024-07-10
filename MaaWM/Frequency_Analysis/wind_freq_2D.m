clear
clc
close all
%% File Path Setup
dataloc = '/home/carlos/Files/GitLab_Repos/Research/MaaWM/WRF_Model';
addpath 'supplemental_routines'
UVWfrontend
%% U

for wmax = 1:20 %%Number of fundamental modes
    
    Nexp = 100; % number of steps
    yvals = linspace(-500,500,Nexp); % varibale that moves from -500 to 500 with N number of steps
    Umat_recreate = []; % prealocated to be a matrix
    uam_mat = []; % prealocated to be a matrix
    ubm_mat = []; % prealocated to be a matrix
    
    for y = yvals % for loop to run for the values of 'yvals' and then stop
        % maximum fundemental frequency    
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
        us = u - a0u; % using a0 to shift function so an and bn can be solved for
        mvec = 1:wmax; % makes a row vector the length of wmax
        uam = 0*mvec; % in order to make the math work prealocate am to be the demensions of mvec
        ubm = 0*mvec; % do the same to bm as done for am
        T = x*pi/500; % waveform arguements
        urecreate = 0*xexp; % make recreate dimensions equal to xexp
        for m = mvec
            %%%Trapezoidal rule
            for idx = 1:length(us)-1 % here the length() doesnt matter as long as it has a length of 40 to make matrix dimensions equal
                ua1 = us(idx)*cos(m*T(idx));  
                ua2 = us(idx+1)*cos(m*T(idx+1));
                uam(m) = uam(m) + (1/500)*(0.5*(ua1+ua2)*dx);
 
                ub1 = us(idx).*sin(m.*T(idx));
                ub2 = us(idx+1).*sin(m.*T(idx+1));
                ubm(m) = ubm(m) + (1/500)*(0.5*(ub1+ub2)*dx);
            end
            %%%Recreate the waveform
            urecreate = urecreate + uam(m)*cos(m*pi/500*xexp)+ubm(m)*sin(m*pi/500*xexp);
        end
        urecreate = urecreate + a0u;
        Umat_recreate = [Umat_recreate;urecreate];
        uam_mat = [uam_mat;uam];
        ubm_mat = [ubm_mat;ubm];
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
    xlabel('Xcoord (m)'); ylabel('Ycoord (m)'); zlabel('U (m/s)')
    legend('Fourier Recreation','Gathered Data')
    
end

% Least Square Error (Outside Main Loop)

for wmax = 1:20

    x = xcoord;
    xexp = linspace(-500,500,400);                      % Vector, and sample rate
    dx = x(2)-x(1);                                     % Get step size
    u = 0*xcoord;

    for idx = 1:length(xcoord)
        uvw = uvwout(xcoord(idx),-500,200,0,dataloc,0);
        u(idx) = uvw(1);
    end
 
    a0u = mean(u);
    us = u - a0u;
    mvec = 1:wmax;
    amu = 0*mvec;
    bmu = 0*mvec;
    T = x*pi/500;
    urecreate = 0*xexp;

    for m = mvec
        %%%Trapezoidal rule
        for idx = 1:length(us)-1

            a1u = us(idx)*cos(m*T(idx));
            a2u = us(idx+1)*cos(m*T(idx+1));
            amu(m) = amu(m) + (1/500)*(0.5*(a1u+a2u)*dx);

            b1u = us(idx)*sin(m*T(idx));
            b2u = us(idx+1)*sin(m*T(idx+1));
            bmu(m) = bmu(m) + (1/500)*(0.5*(b1u+b2u)*dx);
        end
        %%%Recreate the waveform
        urecreate = urecreate + amu(m)*cos(m*pi/500*xexp)+bmu(m)*sin(m*pi/500*xexp);
    end
    urecreate = urecreate + a0u;

    uint = interp1(xexp,urecreate,xcoord);
    LSEu(wmax) = sum((uint-u).^2);  

end

freq = mvec*pi/500;
[ff,~] = meshgrid(freq,yexp);
set(figure,'color','white')
mesh(ff,yexp,uam_mat)
title('a_{n} Frequency of Windfield in U','Fontsize',18)
xlabel('Frequency (rad/m)'); ylabel('Ycoord (m)'); zlabel('Amplitude (m/s)')

set(figure,'color','white')
mesh(ff,yexp,ubm_mat)
title('b_{n} Frequency of Winfield in U','Fontsize',18)
xlabel('Frequency (rad/m)'); ylabel('Ycoord (m)'); zlabel('Amplitude')

set(figure,'color','white')
mesh(ff,yexp,sqrt(ubm_mat.^2+uam_mat.^2))
title('Magnitude of a_{n} and b_{n} frequencies in U','Fontsize',18)
xlabel('Frequency (rad/m)'); ylabel('Ycoord (m)'); zlabel('Amplitude')
    
set(figure,'color','white')
plot(LSEu,'linewidth',2);
title('Least Squares Error vs Frequency of U','FontSize',18);
xlabel('w_{max}','FontSize',13);ylabel('Error','FontSize',13);
grid on

%% V

for wmax = 1:20
    Nexp = 100;
    yvals = linspace(-500,500,Nexp);
    Vmat_recreate = [];
    vam_mat = [];
    vbm_mat = [];
    for y = yvals   
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
        vam = 0*mvec;
        vbm = 0*mvec;
        T = x*pi/500;
        vrecreate = 0*xexp;
        for m = mvec
            %%%Trapezoidal rule
            for idx = 1:length(vs)-1
                
                va1 = vs(idx)*cos(m*T(idx));
                va2 = vs(idx+1)*cos(m*T(idx+1));
                vam(m) = vam(m) + (1/500)*(0.5*(va1+va2)*dx);
                
                vb1 = vs(idx)*sin(m*T(idx));
                vb2 = vs(idx+1)*sin(m*T(idx+1));
                vbm(m) = vbm(m) + (1/500)*(0.5*(vb1+vb2)*dx);
            end
            %%%Recreate the waveform
            vrecreate = vrecreate + vam(m)*cos(m*pi/500*xexp)+vbm(m)*sin(m*pi/500*xexp);
        end
        vrecreate = vrecreate + a0v;
        Vmat_recreate = [Vmat_recreate;vrecreate];
        vam_mat = [vam_mat;vam];
        vbm_mat = [vbm_mat;vbm];
        
%         % Least Squares Error
% 
%         vint = interp1(xexp,vrecreate,xcoord);
%         LSEv(wmax) = sum((vint-v).^2);
    
    end
    
    Vmat = [];
    for y = ycoord % Repeat the first half of the first for loop in order to handle the UVW wind map.
        v = 0*xcoord;
        for idx = 1:length(xcoord)
            uvw = uvwout(xcoord(idx),y,200,0,dataloc,0);
            v(idx) = uvw(2); 
        end
        Vmat = [Vmat;v];
    end
    [xx,yy] = meshgrid(xcoord,ycoord);% meshgrid xcoord(1x40) and ycoord(1x40) in order to make them 40x40
    [xxexp,yyexp] = meshgrid(xexp,yexp); % meshgrid xexp(1x100) and yexp(1x100) in order to make them 100x100
    set(figure,'color','white') %set the figure to have a white background
    mesh(xx,yy,Vmat) 
    hold on
    surf(xxexp,yyexp,Vmat_recreate) %surface plot of Fourier Recreation
    grid on
    title('Windfield in V (Gathered Data vs. Fourier Recreation)','Fontsize',18)
    xlabel('Xcoord (m)'); ylabel('Ycoord (m)'); zlabel('V (m/s)')
    legend('Gathered Data','Fourier Recreation')
    
end

% Least Square Error (Outside Main Loop)

for wmax = 1:20

    x = xcoord;
    xexp = linspace(-500,500,400);                      % Vector, and sample rate
    dx = x(2)-x(1);                                     % Get step size
    v = 0*xcoord;

    for idx = 1:length(xcoord)
        uvw = uvwout(xcoord(idx),-500,200,0,dataloc,0);
        v(idx) = uvw(2);
    end
 

    a0v = mean(v);
    vs = v - a0v;

    mvec = 1:wmax;
    amv = 0*mvec;
    bmv = 0*mvec;
    T = x*pi/500;
    vrecreate = 0*xexp;

    for m = mvec
        %%%Trapezoidal rule
        for idx = 1:length(vs)-1

            a1v = vs(idx)*cos(m*T(idx));
            a2v = vs(idx+1)*cos(m*T(idx+1));
            amv(m) = amv(m) + (1/500)*(0.5*(a1v+a2v)*dx);

            b1v = vs(idx)*sin(m*T(idx));
            b2v = vs(idx+1)*sin(m*T(idx+1));
            bmv(m) = bmv(m) + (1/500)*(0.5*(b1v+b2v)*dx);
        end
        %%%Recreate the waveform
        vrecreate = vrecreate + amv(m)*cos(m*pi/500*xexp)+bmv(m)*sin(m*pi/500*xexp);
    end
    vrecreate = vrecreate + a0v;

    vint = interp1(xexp,vrecreate,xcoord);
    LSEv(wmax) = sum((vint-v).^2);  

end

freq = mvec*pi/500; % waveform arguement set to equal x*pi/500; where x increases by mvec
[ff,~] = meshgrid(freq,yexp); 
set(figure,'color','white')
mesh(ff,yexp,vam_mat)
title('a_{n} Frequency of Windfield in V','Fontsize',18)
xlabel('Frequency (rad/m)'); ylabel('Ycoord (m)'); zlabel('Amplitude')

set(figure,'color','white')
mesh(ff,yexp,vbm_mat)
title('b_{n} Frequency of Winfield in V','Fontsize',18)
xlabel('Frequency (rad/m)'); ylabel('Ycoord (m)'); zlabel('Amplitude')

set(figure,'color','white')
mesh(ff,yexp,sqrt(vbm_mat.^2+vam_mat.^2))
title('Magnitude of a_{n} and b_{n} Frequencies in V','Fontsize',18)
xlabel('Frequency (rad/m)'); ylabel('Ycoord (m)'); zlabel('Amplitude')

set(figure,'color','white')
plot(LSEv,'linewidth',2);
title('Least Squares Error vs Frequency of V','FontSize',18);
xlabel('w_{max}','FontSize',13);ylabel('Error','FontSize',13);
grid on

%% W

for wmax = 1:20
    Nexp = 100;
    yvals = linspace(-500,500,Nexp);
    Wmat_recreate = [];
    wam_mat = [];
    wbm_mat = [];
    for y = yvals    
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
        wam = 0*mvec;
        wbm = 0*mvec;
        T = x*pi/500;
        wrecreate = 0*xexp;
        for m = mvec
            %%%Trapezoidal rule
            for idx = 1:length(ws)-1
                
                wa1 = ws(idx)*cos(m*T(idx));
                wa2 = ws(idx+1)*cos(m*T(idx+1));
                wam(m) = wam(m) + (1/500)*(0.5*(wa1+wa2)*dx);
 
                wb1 = ws(idx)*sin(m*T(idx));
                wb2 = ws(idx+1)*sin(m*T(idx+1));
                wbm(m) = wbm(m) + (1/500)*(0.5*(wb1+wb2)*dx);
            end
            %%%Recreate the waveform
            wrecreate = wrecreate + wam(m)*cos(m*pi/500*xexp)+wbm(m)*sin(m*pi/500*xexp);
        end
        wrecreate = wrecreate + a0w;
        Wmat_recreate = [Wmat_recreate;wrecreate];
        wam_mat = [wam_mat;wam];
        wbm_mat = [wbm_mat;wbm];
        
%         % Least Squares Error
% 
%         wint = interp1(xexp,wrecreate,ycoord);
%         LSEw(wmax) = sum((wint-w).^2);
    
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
    xlabel('Xcoord (m)'); ylabel('Ycoord (m)'); zlabel('W (m/s)')
    legend('Gathered Data','Fourier Recreation')

end

% Least Squares Error (Outside Main Loop)

for wmax = 1:20

    x = xcoord;
    xexp = linspace(-500,500,400);                      % Vector, and sample rate
    dx = x(2)-x(1);                                     % Get step size
    w = 0*xcoord;

    for idx = 1:length(xcoord)
        uvw = uvwout(xcoord(idx),-500,200,0,dataloc,0);
        w(idx) = uvw(3);
    end
 
    a0w = mean(w);
    ws = w - a0w;
    mvec = 1:wmax;
    amw = 0*mvec;
    bmw = 0*mvec;
    T = x*pi/500;
    wrecreate = 0*xexp;

    for m = mvec
        %%%Trapezoidal rule
        for idx = 1:length(us)-1

            a1w = ws(idx)*cos(m*T(idx));
            a2w = ws(idx+1)*cos(m*T(idx+1));
            amw(m) = amw(m) + (1/500)*(0.5*(a1w+a2w)*dx);

            b1w = ws(idx)*sin(m*T(idx));
            b2w = ws(idx+1)*sin(m*T(idx+1));
            bmw(m) = bmw(m) + (1/500)*(0.5*(b1w+b2w)*dx);
        end
        %%%Recreate the waveform
        wrecreate = wrecreate + amw(m)*cos(m*pi/500*xexp)+bmw(m)*sin(m*pi/500*xexp);
    end
    wrecreate = wrecreate + a0w;
    wint = interp1(xexp,wrecreate,xcoord);
    LSEw(wmax) = sum((wint-w).^2);  

end
    
freq = mvec*pi/500;
[ff,~] = meshgrid(freq,yexp);
set(figure,'color','white')
mesh(ff,yexp,wam_mat)
title('a_{n} Frequency of Windfield in W','Fontsize',18)
xlabel('Frequency (rad/m)'); ylabel('Ycoord (m)'); zlabel('Amplitude')

set(figure,'color','white')
mesh(ff,yexp,wbm_mat)
title('b_{n} Frequency of Winfield in W','Fontsize',18)
xlabel('Frequency (rad/m)'); ylabel('Ycoord (m)'); zlabel('Amplitude')

set(figure,'color','white')
mesh(ff,yexp,sqrt(wbm_mat.^2+wam_mat.^2))
title('Magnitude of a_{n} and b_{n} Frequencies in W','Fontsize',18)
xlabel('Frequency (rad/m)'); ylabel('Ycoord (m)'); zlabel('Amplitude')

set(figure,'colr','white')
plot(LSEw,'linewidth',2);
title('Least Squares Error vs Frequency of W','FontSize',18);
xlabel('w_{max}','FontSize',13);ylabel('Error','FontSize',13);
grid on

%% Averages

um_mat = sqrt(uam_mat.^2+ubm_mat.^2);
vm_mat = sqrt(vam_mat.^2+vbm_mat.^2);
wm_mat = sqrt(wam_mat.^2+vbm_mat.^2);
Mag_MAT = sqrt(um_mat.^2+vm_mat.^2+wm_mat.^2);
distance = freq.^(-1);
[dd,yy] = meshgrid(distance,yexp);

set(figure,'color','white')
mesh(ff,yexp,MagMAT)
title('Magnitude of a_{n} and b_{n} Frequencies in UVW','Fontsize',18)
xlabel('Frequency (rad/m)'); ylabel('Ycoord (m)'); zlabel('Amplitude')
set(figure,'color','white')
mesh(dd,yexp,MagMAT)
title('Magnitude of a_{n} and b_{n} Frequencies UVW','Fontsize',18)
xlabel('Distance Between Points of Reference (m)'); ylabel('Ycoord (m)'); zlabel('Amplitude')
