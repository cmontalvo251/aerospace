function state = computeRK4Step(fhandle,state,LMNcontrol)
%%%FIXED STEP INTEGRATION RK4 LOOP
parameters

if nargin == 2
  LMNcontrol = [];
end

k1 = feval(fhandle,state,LMNcontrol);
k2 = feval(fhandle,state+k1*timestep/2,LMNcontrol);
k3 = feval(fhandle,state+k2*timestep/2,LMNcontrol);
k4 = feval(fhandle,state+k3*timestep,LMNcontrol);

phidot = (1/6)*(k1 + 2*k2 + 2*k3 + k4);

state = state + phidot*timestep;


