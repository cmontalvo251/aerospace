function keypress(varargin)
global uk

%%%uk = [dt,dr];

key = get(gcbf,'CurrentKey');

switch key
 %%%Increase Thrust
 case 'rightarrow'
  uk(1) = uk(1) + 1;
  %%%%Decrease Thrust
 case 'leftarrow'
  uk(1) = uk(1) - 1;
  %%%%Pull up 
 case 'downarrow'
  uk(2) = uk(2) - 0.01;
  %%%%Push over?
 case 'uparrow'
  uk(2) = uk(2) + 0.01;
end

if uk(1) > 100
   uk(1) = 100;
end
if uk(1) < 0
   uk(1) = 0;
end
if uk(2) > pi/3
    uk(2) = pi/3;
end
if uk(2) < -pi/3
    uk(2) = -pi/2;
end
