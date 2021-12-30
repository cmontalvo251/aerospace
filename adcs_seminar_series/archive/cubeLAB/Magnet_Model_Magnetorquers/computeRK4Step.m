function [stateOUT,xyzcurrent] = computeRK4Step(fhandle,stateIN,timestep)
%%%FIXED STEP INTEGRATION RK4 LOOP

%%%Only comput control Once per timestep
[LMN,xyzcurrent] = ComputeControl(stateIN);

k1 = feval(fhandle,stateIN,LMN);
k2 = feval(fhandle,stateIN+k1*timestep/2,LMN);
k3 = feval(fhandle,stateIN+k2*timestep/2,LMN);
k4 = feval(fhandle,stateIN+k3*timestep,LMN);

phidot = (1/6)*(k1 + 2*k2 + 2*k3 + k4);

stateOUT = stateIN + phidot*timestep;


