clear
clc
close all

%%%%%Ok this works

%%%%Let's make some random roll pitch and yaw values
phi = -pi + rand*2*pi;
theta = -pi + rand*2*pi;
psi = -pi + rand*2*pi;
%phi = 1.7029
%theta = 2.3670
%psi = 2.8373
% phi = pi/4;
% theta = pi/4;
% psi = pi/4;
ptp = [phi,theta,psi]
ptp_deg = ptp*180/pi

%%%Create the Rotation Matrix
TIB = R123(phi,theta,psi);
TBI = TIB';

[phi_,theta_,psi_] = extract_Euler(TBI);

%%%Now let's convert the Euler Angles to Quaternions
q0123_e2q = euler2quat([phi,theta,psi])

%%%Let's test the norm component
n = norm(q0123_e2q);

%%%Now let's make the transformation matrix using quaternions
TIB_quat = RQUAT(q0123_e2q);
TBI_quat = TIB_quat';

%%%Let's extract the quaternion
%%I guess theoretically I could just extract the Euler Angles and
%then convert to quaternions?
[phi_q,theta_q,psi_q] = extract_Euler(TBI_quat);
q0123_extract_Euler = euler2quat([phi_q,theta_q,psi_q])
%%But that just seems roundabout.

%%This method below will give you one of the possible quaternion solutions
%%The reason is because Euler angles are non-unique 
%%For example take psi = pi and psi = -pi
%They are both technically the same direction 
%There is a discontinuity in the Euler angles causing the issue
%here in 4Dimensional space. 
%%What you are doing is mapping 4D space into 3D space and
%unfortunately you lose information
q0123_e = extract_quaternion(TBI_quat)

ptp_e = quat2euler(q0123_e)
ptp_e_deg = ptp_e*180/pi

