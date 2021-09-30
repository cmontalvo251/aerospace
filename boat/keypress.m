function keypress(varargin)
global uk

%%%uk = [dt,dr];

key = get(gcbf,'CurrentKey');

switch key
    %%%Increase Thrust
 case 'uparrow'
  uk(1) = uk(1) + 1;
  if uk(1) > 100
      uk(1) = 100;
  end
  %%%%Decrease Thrust
 case 'downarrow'
  uk(1) = uk(1) - 1;
  if uk(1) < 0
      uk(1) = 0;
  end
  %%%%Turn Left? 
 case 'leftarrow'
  uk(2) = uk(2) + 0.1;
  if uk(2) > 30*pi/180
      uk(2) = 30*pi/180;
  end
  %%%%Turn Right?
 case 'rightarrow'
  uk(2) = uk(2) - 0.1;
  if uk(2) < -30*pi/180;
      uk(2) = -30*pi/180;
  end
end



