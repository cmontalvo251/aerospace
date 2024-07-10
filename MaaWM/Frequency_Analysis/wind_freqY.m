function [wavelength,Euvw_avg] = wind_freqY(dataloc,iplot)

%%%%Pull in Data
UVWfrontend

Nmax = 20;
x = 0;

Euu = zeros(Nmax,length(xcoord));
Evv = Euu;
Eww = Evv;

xctr = 0;

for x = [xcoord 0]
    
    x
    
    xctr = xctr + 1;
    
    Eu = zeros(Nmax,1);
    Ev = zeros(Nmax,1);
    Ew = zeros(Nmax,1);

    for N = 1:Nmax
        
        u = 0*ycoord;
        v = 0*ycoord;
        w = 0*ycoord;
        z = 200; %%%AGL
        
        for idx = 1:length(ycoord)
            uvw = uvwout(x,ycoord(idx),200,0,dataloc,0);
            u(idx) = uvw(1);
            v(idx) = uvw(2);
            w(idx) = uvw(3);
        end
        %%%n*w0 = w rad/m or angular frequency
        %%%w = 2*pi*f thus f = w/(2*pi) which is frequency in 1/m
        %%%lambda = 1/f which is period or in this case wavelength or meter
        %%%lambda = 2*pi/(n*w0)
        %N = 20; %%Maximum number of frequencies
        tau = 500;
        Nf = 1000;
        [a0u,anu,bnu,yfit,urecreate] = fourier(ycoord,u,N,tau,Nf);
        [a0v,anv,bnv,yfit,vrecreate] = fourier(ycoord,v,N,tau,Nf);
        [a0w,anw,bnw,yfit,wrecreate] = fourier(ycoord,w,N,tau,Nf);
        
        w0 = pi/tau;
        lambda = 2*pi./((1:N).*w0);
        
        if N == Nmax && x == 0 && iplot
            plottool(1,'Am Vals',12,'Wavelength (m)','Amplitudec (m/s)','FontSize');
            plot(lambda,abs(anu),'b-','markersize',16);
            plot(lambda,abs(anv),'r-','markersize',16);
            plot(lambda,abs(anw),'g-','markersize',16);
            
            plottool(1,'Bm Vals',12,'Wavelength (m)','Amplitude (m/s)');
            plot(lambda,abs(bnu),'b-','markersize',16);
            plot(lambda,abs(bnv),'r-','markersize',16);
            plot(lambda,abs(bnw),'g-','markersize',16);
            
            cnu = sqrt(anu.^2 + bnu.^2);
            cnv = sqrt(anv.^2 + bnv.^2);
            cnw = sqrt(anw.^2 + bnw.^2);
            
            plottool(1,'Cn,(u,v,w)',12,'Wavelength (m)','C_n (m/s)');
            plot(lambda,cnu,'b-','LineWidth',2)
            plot(lambda,cnv,'r-','LineWidth',2)
            plot(lambda,cnw,'g-','LineWidth',2)
            legend('C_{n,u}','C_{n,v}','C_{n,w}')
            
            plottool(1,'Example Wave Form',12,'X (meters)','Windspeed (m/s)');
            plot(yfit,urecreate,'b-','linewidth',2);
            plot(yfit,vrecreate,'r-','linewidth',2);
            plot(yfit,wrecreate,'g-','linewidth',2);
            legend('U','V','W')
            plot(ycoord,u,'b*','linewidth',2);
            plot(ycoord,v,'r*','linewidth',2);
            plot(ycoord,w,'g*','linewidth',2);
        end
        
        %%%Compute the error between u and urecreate
        for jdx = 1:length(ycoord)
            yj = ycoord(jdx);
            ubar = interp1(yfit,urecreate,yj);
            Eu(N) = Eu(N) + (u(jdx)-ubar)^2;
            vbar = interp1(yfit,vrecreate,yj);
            Ev(N) = Ev(N) + (v(jdx)-vbar)^2;
            wbar = interp1(yfit,wrecreate,yj);
            Ew(N) = Ew(N) + (w(jdx)-wbar)^2;
        end
    end
    Eu = sqrt(Eu./Eu(1));
    Ev = sqrt(Ev./Ev(1));
    Ew = sqrt(Ew./Ew(1));
    
    N_vec = 1:Nmax;
    
    if x == 0 && iplot
        plottool(1,'Error',12,'N','Normalized E_N');
        plot(N_vec,Eu,'b-','LineWidth',2)
        plot(N_vec,Ev,'r-','LineWidth',2)
        plot(N_vec,Ew,'g-','LineWidth',2)
        legend('U','V','W')
    else
        Euu(:,xctr) = Eu;
        Evv(:,xctr) = Ev;
        Eww(:,xctr) = Ew;
    end
    
end

%%%%Create Mesh Plots
Euvw = (Euu+Evv+Eww)/3;
Euvw_avg = sum(Euvw')/40;
wavelength = 2*pi./(w0*N_vec);
if iplot
[xx,yy] = meshgrid(xcoord,N_vec);

plottool(1,'Error Mesh',12,'N','X coordinate (m)','Normalized E_{N,u}');
mesh(yy,xx,Euu)
plottool(1,'Error Mesh',12,'N','X coordinate (m)','Normalized E_{N,v}');
mesh(yy,xx,Evv)
plottool(1,'Error Mesh',12,'N','X coordinate (m)','Normalized E_{N,w}');
mesh(yy,xx,Eww)

%%%Plot Combination of all of them
plottool(1,'Error Mesh ALL',12,'N','X coordinate (m)','Normalized E_{N}');
mesh(yy,xx,Euvw)


%%%%Now Average the whole thing
plottool(1,'Error ALL',12,'N','Normalized E_{N}')
plot(N_vec,Euvw_avg,'LineWidth',2)


plottool(1,'Error ALL',12,'Wavelength (m)','Normalized E_{N}')
plot(wavelength,Euvw_avg,'LineWidth',2)
end
