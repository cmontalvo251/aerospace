function q0123 = extract_quaternion(TBI)
%%%Assuming T is a 3x3 matrix, extract the quaternion vector (q0,q1,q2,q3)
%%%assuming a 3-2-1 transformation sequence
%%%Let T be defined such that v(body) = T v(inertial)
%%In Mark Costello's notation, T would TBI
alfa = TBI(1,2) + TBI(2,1);
%alfai = 4*q1i*q2i
bata = TBI(3,1) + TBI(1,3);
%batai = 4*q1i*q3i
gama = TBI(3,2) + TBI(2,3);
%gamai = 4*q2i*q3i
q1squared = alfa*bata/(4*gama);

%%%Solutions split into two different solutions here
q1a = sqrt(q1squared);
q1b = -sqrt(q1squared);
q2a = alfa/(4*q1a);
q2b = alfa/(4*q1b);
q3a = bata/(4*q1a);
q3b = bata/(4*q1b);

%%%These are the same so just pick one
%q0asquared = TBI(1,1) + q2a^2 + q3a^2 - q1a^2
%q0bsquared = TBI(1,1) + q2b^2 + q3b^2 - q1b^2
q0squared = TBI(1,1) + q2a^2 + q3a^2 - q1a^2;

%%%Solution However still splits into 4 possible solutions
q0a = sqrt(q0squared);
q0b = -sqrt(q0squared);

%%%Here are my 4 possible solutions
q0123aa = [q0a,q1a,q2a,q3a];
q0123ba = [q0b,q1a,q2a,q3a];
q0123ab = [q0a,q1b,q2b,q3b];
q0123bb = [q0b,q1b,q2b,q3b];

%According to this website there are multiple solutions 
%that yield the same euler angles. 
%http://planning.cs.uiuc.edu/node151.html
%So what we want to do is compute the euler angles from the matrix
[phi,theta,psi] = extract_Euler(TBI);
ptp = [phi,theta,psi];
%Then check and see which have the same euler angles
quats = [q0123aa;q0123ba;q0123ab;q0123bb];
valid_quats = [];
for j = 1:4
  q0123j = quats(j,:);
  ptpj = quat2euler(q0123j);
  val = abs(sum(ptpj-ptp));
  if val < 1e-10
    valid_quats = [valid_quats;q0123j];
  end
end
%%http://planning.cs.uiuc.edu/node151.html
%After this loop the number of valid quats will only be 2
%Again this article above states that there are 2 quaternions
%for every set of Euler Angles since a quaternion E(q) = E(-q)
%The article then states that to get around this one simply needs
%to adopt a sign convention and force q0 > 0 
for j = 1:2
  if valid_quats(j,1) > 0 
    q0123 = valid_quats(j,:);
    return
  end
end


% Copyright - Carlos Montalvo 2019
% You may freely distribute this file but please keep my name in here
% as the original owner