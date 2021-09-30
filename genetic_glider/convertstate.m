function [m,I,xbody] = convertstate(x);
%x = [1,1,5,pi/4,1,1,0,pi/4,5];
%l1,w1,d2,g1

numpanels = 2;
rodflag = 0;

lrod = x(end);
m = 0;
xcgR = 0;
ycgR = 0;
zcgR = 0;

for ii = 1:numpanels
  iistr = num2str(ii);

  %%get panels out
  eval(['l',iistr,'=x(',num2str((ii-1)*4+1),');'])
  eval(['w',iistr,'=x(',num2str((ii-1)*4+2),');'])
  eval(['d',iistr,'=x(',num2str((ii-1)*4+3),');'])
  eval(['g',iistr,'=x(',num2str((ii-1)*4+4),');'])

  %Get orientation of panels
  eval(['phi',iistr,'=g',iistr,';'])
  eval(['theta',iistr,'=0;'])
  eval(['psi',iistr,'=0;'])

  %%Get position of panels
  eval(['x',iistr,'= d',iistr,';'])
  eval(['y',iistr,'= l',iistr,'/2*cos(phi',iistr,');'])
  eval(['z',iistr,'= l',iistr,'/2*sin(phi',iistr,');'])

  %%Mirror panel to opposite side and create panelii+numpanels
  eval(['w',num2str(ii+numpanels),'=w',iistr,';'])
  eval(['l',num2str(ii+numpanels),'=l',iistr,';'])
  eval(['x',num2str(ii+numpanels),'=x',iistr,';'])
  eval(['y',num2str(ii+numpanels),'=-y',iistr,';'])
  eval(['z',num2str(ii+numpanels),'=z',iistr,';'])  
  eval(['phi',num2str(ii+numpanels),'=pi/2-phi',iistr,';'])
  eval(['theta',num2str(ii+numpanels),'=0;'])
  eval(['psi',num2str(ii+numpanels),'=0;'])

  %%Get weight and cg location of aircraft
  eval(['A=w',num2str(ii),'*l',num2str(ii),';'])
  m = m + 2*A;
  eval(['xA=2*x',num2str(ii),'*A;'])
  eval(['zA=2*z',num2str(ii),'*A;'])
  xcgR = xcgR + xA;
  zcgR = zcgR + zA;
  
end
m = m + rodflag.*lrod; %%add mass of rod
xcgR = xcgR/m;
zcgR = zcgR/m;

%%Compute Inertia of aircraft
%%Add moment of inertia of rod
I = rodflag.*[0 0 0;0 (lrod^3)/12 0;0 0 (lrod^3)/12];
for ii = 1:numpanels*2
  eval(['A=w',num2str(ii),'*l',num2str(ii),';'])
  %%Inertia of Panel
  eval(['Ixx = (1/12)*A*l',num2str(ii),'^2;'])
  eval(['Iyy = (1/12)*A*w',num2str(ii),'^2;'])
  eval(['Izz = (1/12)*A*(l',num2str(ii),'^2+w',num2str(ii),'^2);'])
  %%Rotate to Body Frame
  eval(['cphi = cos(phi',num2str(ii),');'])
  eval(['sphi = sin(phi',num2str(ii),');'])
  R = [1 0 0;0 cphi sphi;0 -sphi cphi]; %%transformation to go from
                                        %body frame to panel frame
                                        %thus vpanel = R'*vbody
  IpanelP = [Ixx 0 0;0 Iyy 0;0 0 Izz];
  IbodyP = R*IpanelP*R';
  
  %%Parrallel Axis Theorem
  iistr = num2str(ii);
  eval(['rBPvec = [x',iistr,',y',iistr,',z',iistr,'];'])
  rBP = skew(rBPvec);
  IbodyB = IbodyP + A*rBP*rBP';
  
  %%Add moments of inertia
  I = I + IbodyB;
end


%%Load Variables into panels
xbody = zeros(8,numpanels*2);
for ii = 1:numpanels*2
  iistr = num2str(ii);
  eval(['panel = [x',iistr,'-xcgR,y',iistr,',z',iistr,'-zcgR,phi',iistr,',theta',iistr,',psi',iistr,',l',iistr,',w',iistr,'];'])
  xbody(:,ii) = panel';
end

function out = skew(vec)

out = zeros(3,3);

out(1,3) = vec(2);
out(1,2) = -vec(3);
out(2,1) = vec(3);
out(2,3) = -vec(1);
out(3,1) = -vec(2);
out(3,2) = vec(1);



