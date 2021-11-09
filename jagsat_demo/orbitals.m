close all 
L=input('Enter Angular Momentum Quantum Number:');
M=input('Enter Magnetic Quantum Number:');
R=input('Enter Resolution: (Leave blank for 100)');
if isempty(R)
    R=100;
end
RES=[R R];
if L<M, error('ml must be less than or eqaul to l!'); end
THETA=linspace(0,2*pi,RES(1));  % Azimuthal/Longitude/Circumferential
PHI  =linspace(0,  pi,RES(2));  % Altitude /Latitude /Elevation
[THETA,PHI]=meshgrid(THETA,PHI);
Lmn=legendre(L,cos(PHI));
if L~=0
  Lmn=squeeze(Lmn(M+1,:,:));
end
a1=((2*L+1)/(4*pi));
a2=factorial(L-M)/factorial(L+M);
C=sqrt(a1*a2);
Ymn=C*Lmn.*exp(1i*M*THETA);
[Xm,Ym,Zm]=sph2cart(THETA,PHI-pi/2,abs(Ymn).^2);
[Xr,Yr,Zr]=sph2cart(THETA,PHI-pi/2,real(Ymn).^2);
[Xi,Yi,Zi]=sph2cart(THETA,PHI-pi/2,imag(Ymn).^2);
% [Xp,Yp,Zp]=sph2cart(THETA,PHI-pi/2,angle(Ymn).^2);
f=figure; axis off; hold on;
  axes('position',[0.0500 0 0.2666 1]); 
    surf(Xm,Ym,Zm); 
    axis equal off; %rot3d;
    light; lighting phong; camzoom(1.3);
  axes('position',[0.3666 0 0.2666 1]); 
    surf(Xr,Yr,Zr); 
    axis equal off; %rot3d;
    light; lighting phong; camzoom(1.3);
  axes('position',[0.6833 0 0.2666 1]); 
    surf(Xi,Yi,Zi); 
    axis equal off; %rot3d;
    light; lighting phong; camzoom(1.3);
  axes('position',[0 0.9 1 0.1]); axis off;
  text(0.50,0.25,[num2str(L),', ',num2str(M),' Orbital'],'HorizontalAlignment','Center');
  axes('position',[0 0.85 1 0.1]); axis off;
    text(0.50,0.25,['95% Chance of Finding Electron'],'HorizontalAlignment','Center');
  axes('position',[0 0.82 1 0.1]); axis off;
    text(0.50,0.25,['in this Region for l=',num2str(L),' and ml=',num2str(M)],'HorizontalAlignment','Center');
  axes('position',[0 0.3 1 0.1]); axis off;
    text(0.20,0.25,['|Y^',num2str(M),'_',num2str(L),'|^2'],'HorizontalAlignment','Center');
    text(0.50,0.25,['Real(Y^',num2str(M),'_',num2str(L),')^2'],'HorizontalAlignment','Center');
    text(0.80,0.25,['Imag(Y^',num2str(M),'_',num2str(L),')^2'],'HorizontalAlignment','Center');
  %setfig(gcf,10,5,12);
return