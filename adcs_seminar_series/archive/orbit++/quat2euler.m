function ptp = quat2euler(q0123,iplot,t,BodyName)
%%ptp = quat2euler(q0123)
%input is a Nx4 vector with quaternions.
%output is a Nx3 vector of 3-2-1 euler angles

if nargin < 4
    BodyName = '';
end

q0 = q0123(:,1);
q1 = q0123(:,2);
q2 = q0123(:,3);
q3 = q0123(:,4);

ptp(:,1) = (atan2(2.*(q0.*q1 + q2.*q3),1-2.*(q1.^2 + q2.^2))); %phi
ptp(:,2) = asin(2.*(q0.*q2-q3.*q1)); %theta
ptp(:,3) = atan2(2.*(q0.*q3 + q1.*q2),1-2.*(q2.^2 + q3.^2)); %psi

% q0 = cos(phi/2).*cos(theta/2).*cos(psi/2) + sin(phi/2).*sin(theta/2).*sin(psi/2)
% q1 = sin(phi/2).*cos(theta/2).*cos(psi/2) - cos(phi/2).*sin(theta/2).*sin(psi/2)
% q2= cos(phi/2).*sin(theta/2).*cos(psi/2) + sin(phi/2).*cos(theta/2).*sin(psi/2)
% q3 = cos(phi/2).*cos(theta/2).*sin(psi/2) - sin(phi/2).*sin(theta/2).*cos(psi/2)

ptp = real(ptp);

if nargin < 2
    iplot = 0;
end
if iplot
    Names = {['\phi(deg) ',BodyName],['\theta(deg) ',BodyName],['\psi(deg) ',BodyName]};
    for ii = 1:3
        plottool(1,'PTP',18,'Time(sec)',Names{ii});
        f = 180./pi;
        plot(t,f.*ptp(:,ii),'k-','LineWidth',2)
    end
end

% Copyright - Carlos Montalvo 2015
% You may freely distribute this file but please keep my name in here
% as the original owner
