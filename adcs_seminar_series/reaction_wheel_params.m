%%%Reaction wheel parameters
mr = .13;
rr = 42/1000;
hr = 19/1000;

%%%Maximum Speed
rpm = 8000; %%%rpm
maxSpeed = rpm*2*pi/60;

%%%Maximum Torque
maxTorque = 0.004; %%N-m

%%%Orientation of Reaction Wheel
n3 = [1;0;0];
n1 = [0;1;0];
%n3 = [0;0;1];
n2 = [1;1;1];
n2 = n2/norm(n2);

%%%Move the reaction wheels a bit from the cg
r1 = [4;0;0]/1000;
r2 = [0;4;0]/1000;
r3 = [0;0;4]/1000;

%%%Compute the inertia of the reaction wheels
Idisk = (1/12)*(3*rr^2+hr^2);
IrR = mr*[(1/2)*rr^2 0 0 ;0 Idisk 0; 0 0 Idisk];

%%%Compute maximum angular acceleration
maxAlpha = maxTorque/IrR(1,1);

%%%Transofrmation from RR frame to body frame of satellite
T1 = Rscrew(n1);
T2 = Rscrew(n2);
T3 = Rscrew(n3);
%%%Compute inertia of reaction wheel in body frame of satellite
Ir1B = T1'*IrR*T1;
Ir2B = T2'*IrR*T2;
Ir3B = T3'*IrR*T3;

%%%Compute J for control
J = [Ir1B*n1,Ir2B*n2,Ir3B*n3];
%Jinv = inv(J'*J)*J'; %%%if using 3 reaction wheels this one will work
Jinv = J'*(inv(J*J')); %%Use this if using more than 3 reaction wheels

%%%Compute the inertia of the reaciotn wheel in the body frame in ref to the cg of the satellite
%%%which means we need to use the parallel axis theorem
sr1 = skew(r1);
Ir1Bcg = Ir1B + mr*(sr1')*sr1;
sr2 = skew(r2);
Ir2Bcg = Ir2B + mr*(sr2')*sr2;
sr3 = skew(r3);
Ir3Bcg = Ir3B + mr*(sr3')*sr3;