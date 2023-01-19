purge

% Plot throttle curve
tmid = 1500;
tin = linspace(1100,1900,100);
plottool(1,'Throttle Curve',18,'Throttle In (us)','Throttle Out (us)')
for texp = [1:0.1:2]
    tout = throttle_curve_FINAL(tin,tmid,texp);
    plot(tin,tout)
end

throttle = 1800;
%pitch = 1500;
yaw = 1500;
gear = 0; % not needed right now
tmid = 1500;
texp = 1;
rate = 1.0;
plottool(1,'ESCs',18,'Roll Command (us)','Pitch Command (us)','ESC Value')
color = {'b*','r*','g*','m*'};
roll_vec = 1100:10:1900;
pitch_vec = 1100:10:1900;
[rr,pp] = meshgrid(roll_vec,pitch_vec);
esc_mesh = zeros(length(roll_vec),length(pitch_vec),4);
rctr = 0;
pctr = 0;
for roll = roll_vec
    rctr = rctr + 1;
    pctr = 0;
    for pitch = pitch_vec
        pctr = pctr + 1;
        esc_commands = acro_mode(throttle,roll,pitch,yaw,gear,tmid,texp,rate);
        for jdx = 1:4
            esc_mesh(rctr,pctr,jdx) = esc_commands(jdx);
        %    plot3(roll,pitch,esc_commands(jdx),color{jdx})
        end
    end
end
for idx = 1:4
    % plottool(1,'ESCs',18,'Roll Command (us)','Pitch Command (us)','ESC Value')
    mesh(rr,pp,esc_mesh(:,:,idx))
end