uk = [0;0];
fig = figure('Name','Boat_Video_Game','keypressfcn',@keypress,'KeyReleaseFcn',@keyrelease);
xlabel('X (m)')
ylabel('Y (m)')
grid on
%%%Kick off the timer
tic
%%%Reverse the ydirection for aircraft convention
set(gca,'YDir','reverse')