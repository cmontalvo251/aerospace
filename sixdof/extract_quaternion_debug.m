function valid_quats = extract_quaternion_debug(TBI,q0123i)
%%%Assuming T is a 3x3 matrix, extract the quaternion vector (q0,q1,q2,q3)
%%%assuming a 3-2-1 transformation sequence
%%%Let T be defined such that v(body) = T v(inertial)
%%In Mark Costello's notation, T would TBI
TIB = TBI';

q0i = q0123i(1);
q1i = q0123i(2);
q2i = q0123i(3);
q3i = q0123i(4);

alfa = TBI(1,2) + TBI(2,1);
%alfai = 4*q1i*q2i
bata = TBI(3,1) + TBI(1,3);
%batai = 4*q1i*q3i
gama = TBI(3,2) + TBI(2,3);
%gamai = 4*q2i*q3i
q1squared = alfa*bata/(4*gama);
%q1isquared = q1i^2

%%%Solutions split into two different solutions here
q1a = sqrt(q1squared);
q1b = -sqrt(q1squared);
%q1i

q2a = alfa/(4*q1a);
q2b = alfa/(4*q1b);
%q2i

q3a = bata/(4*q1a);
q3b = bata/(4*q1b);
%q3i

%%%These are the same so just pick one
%q0asquared = TBI(1,1) + q2a^2 + q3a^2 - q1a^2
%q0bsquared = TBI(1,1) + q2b^2 + q3b^2 - q1b^2
q0squared = TBI(1,1) + q2a^2 + q3a^2 - q1a^2;
%q0isquared = q0i^2

%%%Solution However still splits into 4 possible solutions
q0a = sqrt(q0squared);
q0b = -sqrt(q0squared);
%q0i

%%%Here are my 4 possible solutions
q0123aa = [q0a,q1a,q2a,q3a];
q0123ba = [q0b,q1a,q2a,q3a];
q0123ab = [q0a,q1b,q2b,q3b];
q0123bb = [q0b,q1b,q2b,q3b];
%q0123i 

%%%So which one of the 4 is it? 
%%%All of these solutions satisfy the norm component that is
%norm(q0123ij) = 1 for all ij = ab
%What if we compute the TIB matrix?
%%Ok we're getting somewhere
%TIBaa = RQUAT(q0123aa);
%TIBab = RQUAT(q0123ab);
%TIBba = RQUAT(q0123ba);
%TIBbb = RQUAT(q0123bb);
%TIB = TBI'
%Let's check the sign of alfa
%Remember there are only 2 solutions for q1,q2,q3
% alfa_a = 4*q1a*q2a;
% alfa_b = 4*q1b*q2b;
% VALID_ALFAS = (sign(alfa) == sign([alfa_a,alfa_b]))
% %bata
% bata_a = 4*q1a*q3a;
% bata_b = 4*q1b*q3b;
% VALID_BATAS = (sign(bata) == sign([bata_a,bata_b]));
% %gama
% gama_a = 4*q2a*q3a;
% gama_b = 4*q2b*q3b;
% VALID_GAMAAS = (sign(gama) == sign([gama_a,gama_b]));
% VALID_SOLUTIONS = (VALID_ALFAS + VALID_BATAS + VALID_GAMAAS) == 3
% Ok this isn't it. So let's go back to the TIB matrix
TIBaa_s = sign(RQUAT(q0123aa));
TIBab_s = sign(RQUAT(q0123ab));
TIBba_s = sign(RQUAT(q0123ba));
TIBbb_s = sign(RQUAT(q0123bb));
TIB_s = sign(TIB);

aa = sum(sum(TIB_s == TIBaa_s));
ab = sum(sum(TIB_s == TIBab_s));
ba = sum(sum(TIB_s == TIBba_s));
bb = sum(sum(TIB_s == TIBbb_s));

valid_quats = [];

if aa == 9
  q0123o = q0123aa;
  valid_quats = [valid_quats;q0123o];
end
if ab == 9
  q0123o = q0123ab;
  valid_quats = [valid_quats;q0123o];
end
if ba == 9
  q0123o = q0123ba;
  valid_quats = [valid_quats;q0123o];
end
if bb == 9
  q0123o = q0123bb;
  valid_quats = [valid_quats;q0123o];
end

%%%This will yield anywhere from 2 to 4 solutions
%valid_quats
%%What I don't understand is how multiple quaternions will
%%yield the same transformation matrix
%%I guess I'm curious what the corresponding Euler angles are
%[r,c] = size(valid_quats);
%TIB
%for j = 1:r
  %q0123j = valid_quats(j,:)
  %TIBj = RQUAT(q0123j)
  %ptpj = quat2euler(q0123j)
%end

%%An idea I had is to see if the quaternions produce the same euler
%angles and if they do use those quaternions since there are
%multiple solutions.



% Copyright - Carlos Montalvo 2019
% You may freely distribute this file but please keep my name in here
% as the original owner