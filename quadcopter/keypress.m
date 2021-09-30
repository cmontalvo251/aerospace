function keypress(varargin)
global uk

%%%uk = [dt,dr];

key = get(gcbf,'CurrentKey');

switch key
    %%%Increase Thrust
 case 'uparrow'
  uk(1) = uk(1) + 1;
  uk(2) = uk(2) + 1;
  uk(4) = uk(4) + 1;
  %%%%Decrease Thrust
 case 'downarrow'
  uk(1) = uk(1) - 1;
  uk(2) = uk(2) - 1;
  uk(4) = uk(4) - 1;
  %%%%Turn Left? 
 case 'leftarrow'
  uk(1) = uk(1) - 0.1;
  uk(2) = uk(2) + 0.1;
  uk(3) = uk(3) - 0.01;
  %%%%Turn Right?
 case 'rightarrow'
  uk(1) = uk(1) + 0.1;
  uk(2) = uk(2) - 0.1;
  uk(3) = uk(3) + 0.01;
end

if uk(1) > 100
   uk(1) = 100;
end
if uk(2) > 100;
   uk(2) = 100;
end
if uk(4) > 100
    uk(4) = 100;
end
if uk(4) < 0
    uk(4) = 0;
end
if uk(1) < 0
   uk(1) = 0;
end
if uk(2) < 0
   uk(2) = 0;
end
if uk(3) > pi/3
    uk(3) = pi/3;
end
if uk(3) < -pi/3
    uk(3) = -pi/2;
end
