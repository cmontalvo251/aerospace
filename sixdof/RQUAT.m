function R = RQUAT(q0,q1,q2,q3)
%compute R such that v(inertial) = R v(body)
%Using Mark Costello's notation this would be TIB

if length(q0) == 4
  arg = q0;
  q0 = arg(1);
  q1 = arg(2);
  q2 = arg(3);
  q3 = arg(4);
end

%%%Taken fro Boom 2010
%xcgdot = (q0**2+q1**2-q2**2-q3**2)*ub + 2*(q1*q2-q0*q3)*vb
%     &                                      + 2*(q0*q2+q1*q3)*wb 
%      ycgdot = 2*(q1*q2+q0*q3)*ub + (q0**2-q1**2+q2**2-q3**2)*vb
%     &                            + 2*(q2*q3-q0*q1)*wb 
%      zcgdot = 2*(q1*q3-q0*q2)*ub + 2*(q0*q1+q2*q3)*vb 
%     &                            + (q0**2-q1**2-q2**2+q3**2)*wb

R = [(q0^2+q1^2-q2^2-q3^2) 2*(q1*q2-q0*q3)  2*(q0*q2+q1*q3);
      2*(q1*q2+q0*q3) (q0^2-q1^2+q2^2-q3^2) 2*(q2*q3-q0*q1);
      2*(q1*q3-q0*q2) 2*(q0*q1+q2*q3) (q0^2-q1^2-q2^2+q3^2)];

% Copyright - Carlos Montalvo 2015
% You may freely distribute this file but please keep my name in here
% as the original owner
