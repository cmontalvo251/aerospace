function esc_commands = acro_mode(throttle,roll,pitch,yaw,gear,tmid,texp,rate)

% Rates? - Need to figure out max rate
%%%Rate_Max if throttle = 1500 and roll = 1900 we want escL = 1900 and escR = 1100
% escL = throttle + rate_max*roll
% 1900 = 1500 + rate_max*(1900-1500)
% 400/400 = rate_max
% rate_max = 1
if rate > 1
    rate = 1;
end

% Assume roll and pitch rates are the same but yaw_rate is half
rate_yaw = rate/2;

%%%This takes 5 channel signals and returns the ESC vals
throttle_out = throttle_curve_FINAL(throttle,tmid,texp);

% Now that we have throttle_out, in acro mode we just feed through our
% signals
esc_commands = [0;0;0;0];
% Top Left
esc_commands(1) = throttle_out + rate*(roll-1500) + rate*(pitch-1500) + rate_yaw*(yaw-1500);
% Bottom Left
esc_commands(2) = throttle_out + rate*(roll-1500) - rate*(pitch-1500) - rate_yaw*(yaw-1500);
% Top Right
esc_commands(3) = throttle_out - rate*(roll-1500) + rate*(pitch-1500) - rate_yaw*(yaw-1500);
% Bottom Right
esc_commands(4) = throttle_out - rate*(roll-1500) - rate*(pitch-1500) + rate_yaw*(yaw-1500);

% Saturation Block
for idx = 1:4
    if esc_commands(idx) < 1100
        esc_commands(idx) = 1100;
    end
    if esc_commands(idx) > 1900
        esc_commands(idx) = 1900;
    end
end
