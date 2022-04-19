uk = [0;0;0;0];
fig = figure('Name','Aircraft_Video_Game','keypressfcn',@keypress,'KeyReleaseFcn',@keyrelease);
xlabel('X (m)')
ylabel('Z (m)')
grid on
%%%Kick off the timer
tic
%%%Reverse the ydirection for aircraft convention
set(gca,'YDir','reverse')