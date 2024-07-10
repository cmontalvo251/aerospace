clear
clc
close all

dataloc = '/home/carlos/Documents/USA/Research/Atmospheric_Windmapper/MaaWM/WRF_Wind_Data/';
UVWfrontend

%%%Coordinates of X is xcoord
%%%Coordinates of Y is ycoord
%%%Coordinates of Z is zcoord

%let 
figure()
for y = ycoord
    z = 200;

    uvec = 0*xcoord;
    vvec = uvec;
    wvec = uvec;
    
    for idx = 1:length(xcoord)
        uvw = uvwout(xcoord(idx),y,z,0,dataloc,0);
        uvec(idx) = uvw(1);
        vvec(idx) = uvw(2);
        wvec(idx) = uvw(3);
    end
    
    plot3(xcoord,y*ones(length(xcoord),1),uvec,'b*')
    hold on
    plot3(xcoord,y*ones(length(xcoord),1),wvec,'r*')
    plot3(xcoord,y*ones(length(xcoord),1),vvec,'g*')
    drawnow
    pause(0.5)
end
