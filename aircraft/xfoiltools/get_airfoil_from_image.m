clear
clc
close all

LOADIMAGE = 0;

if LOADIMAGE
img = imread('number3.jpg');
imshow(img)
title('Crop image and then hit enter')
pause
title('Now select upper surface points (LE to TE) then press enter')
[xu,yu] = ginput
title('Now select lower surface points (LE to TE) then press enter')
[xl,yl] = ginput
%%%For backup
save airfoil.mat
else
    load airfoil.mat
    %%%%Translate to cartesian coordinates
    x0 = (yu(1) + yl(1))/2;
    y0 = (xu(1) + xl(1))/2;
    %%%Shift and rotate
    xupper = -(yu-x0);
    yupper = -(xu-y0);
    xlower = -(yl-x0);
    ylower = -(xl-y0);
    %%%Now scale
    xlen_upper = xupper(end)-xupper(1);
    xlen_lower = xlower(end)-xlower(1);
    chord = (xlen_upper + xlen_lower)/2;
    xupper = xupper/chord;
    yupper = yupper/chord;
    xlower = xlower/chord;
    ylower = ylower/chord;
    %%%Average LE and TE
    xLE = (xupper(1) + xlower(1))/2;
    yLE = (yupper(1) + ylower(1))/2;
    xupper(1) = xLE;
    xlower(1) = xLE;
    yupper(1) = yLE;
    ylower(1) = yLE;
    xTE = (xupper(end) + xlower(end))/2;
    yTE = (yupper(end) + ylower(end))/2;
    xupper(end) = xTE;
    xlower(end) = xTE;
    yupper(end) = yTE;
    ylower(end) = yTE;
    %%%Now shift everything one more time to xLE and yLE are zero
    xupper = xupper - xupper(1);
    xlower = xlower - xlower(1);
    yupper = yupper - yupper(1);
    ylower = ylower - ylower(1);
    %%%Plot
    p0 = plot(xupper,yupper,'b*');
    hold on
    cubicsplines(xupper,yupper,p0)
    plot(xlower,ylower,'r*')
    cubicsplines(xlower,ylower,p0)
    axis equal
    %%%Now we need to smooth upper and lower surfaces
    %%%Use splines
    
    %%%%Fit to Nth order polynomials
    %%%%Assume that y(x=xLE) = yLE and y(x=xTE) = xTE
    %%%Derivation is in Rocketbook pdf same folder
    bigY = yupper - yTE/xTE*xupper;
    N = 2; %%Order of polynomial fit
    H = zeros(length(bigY),N-1);
    for idx = 1:N-1
        H(:,idx) = xupper.^(idx+1) - (xTE^idx)*(xupper.^idx);
    end
    %%Solve for coefficients
    bigA = inv(H'*H)*H'*bigY;
    %%Solve for a1 (we know that a0 = 0 by shifting everything to the
    %%%origin)
    a1 = yTE/xTE;
    for idx = 1:length(bigA)
       a1 = a1 - bigA(idx)*xTE^idx; 
    end
    bigA = [a1;bigA];
    %%%Create smooth upper surface
    xupper_smooth = linspace(xupper(1),xupper(end),1000);
    yupper_smooth = 0*xupper_smooth;
    for idx = 1:length(bigA)
        yupper_smooth = yupper_smooth + bigA(idx)*xupper_smooth.^idx;
    end
    %plot(xupper_smooth,yupper_smooth,'k-')
end