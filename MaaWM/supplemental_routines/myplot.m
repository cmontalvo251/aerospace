function myplot(myname,myfont,myx,myy,myz,myview)

fig = figure('Name',myname);
set(fig,'color','white')
set(axes,'FontSize',myfont)
xlabel(myx)
ylabel(myy)
zlabel(myz)
view(myview);
hold on
grid on
