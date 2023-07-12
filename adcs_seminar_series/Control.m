function [current,rwalphas] = Control(BfieldNav,pqrNav,ptpNav)
global IrR Jinv

%%%%%%%%%%%%%%%%BDOT CONTROLLER%%%%%%%%%%%%%
k = 67200;
%%%Bfield is in teslas - 40000 nT = 4e4e-9 = 4e-5 ~= 1e-5
%%%pqr is in rad/s --    0.1 rad/s  = 1e-1
%%% pqr*Bfield = 1e-1*1e-5 = 1e-6
%%% pqr*Bfield / (n*A) = 6e=7
%%% muB = n*i*A
magtorquer_params
current = k*cross(pqrNav,BfieldNav)/(n*A);
%%%curent to be in amps ~= 40mA = 40e-3 = 4e-2 

%%%%%%%%%%%%%%RW CONTROLLER%%%%%%%%%%%%%%%
%reaction_wheel_params
if sum(abs(pqrNav)) < 0.1
    KP = eye(3)*1.0*IrR(1,1);
    ptpcommand = [1;0;0];
    pqrcommand = [0;0;0];
    KD = eye(3)*45*IrR(1,1);
    Mdesired = -KD*(pqrcommand - pqrNav) - KP*(ptpcommand-ptpNav);
else
    Mdesired = [0;0;0];
end

%%%Invert Mdesired
%rwalphas = Mdesired/IrR(1,1);
rwalphas = Jinv*Mdesired;


