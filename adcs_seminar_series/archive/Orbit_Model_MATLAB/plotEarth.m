purge
% Simulation
fig = figure();
A = imread('earth.png');
set(fig,'color','white');
sky = imread('world-map.jpg');
anakin = 1;
a = 1.2;
[X,Y,Z] = ellipsoid(0,0,0,anakin,anakin,anakin);
for n = 1
    cla 
    hold on 
    surf(X,Y,Z);
    h = findobj('Type','surface');
    %set(h,'CData',flipdim(sky,1),'FaceColor','texturemap','edgecolor','non
    %e','FaceLighting','Gouraud','Clipping','off')
    set(h,'CData',flipdim(sky,1),'FaceColor','texturemap','FaceLighting','Gouraud','Clipping','off')
    view(-27,70)
    axis equal 
    axis ([-a a -a a])
    set(gca,'Fontsize',14)
    set(gca,'ZTickLabel',{'-0.1' '' '0.1'})
    xlabel( '$\hat{X}$(nd)','interpreter','LaTeX')
    ylabel( '$\hat{Y}$(nd)','interpreter','LaTeX')
    title('Unit Position')
    grid on
end 
