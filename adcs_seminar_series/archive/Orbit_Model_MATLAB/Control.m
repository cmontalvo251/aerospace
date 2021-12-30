% PID Controller
function u = Control(State,time)
kp = 180; %%% proportional gain or propulsive force
kd = 15;  %%% derivative gain or resistive force
xcommand = 3;
xdotcommand = 0;
u = 0;
if time > 5
    u = -kp*(State(1)-xcommand) - kd*(State(2)-xdotcommand);
end
