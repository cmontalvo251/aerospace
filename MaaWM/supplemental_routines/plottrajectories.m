function plottrajectories(filename,IMOVIE,PLOTENTIRETRAJECTORY,SAVEFIG)

%%%PLOTTING ROUTINE
data = dlmread(filename);
state = data(:,2:end);
[r,c] = size(state);
NSTATES = 12; %%%The problem here is that hmmmm.. Idea. I think I fixed it....
num_ac = round(c/NSTATES);

colors = {'b','r','g','m','c','k','y'};
linetype = {'-','--','-.','-','--','-.'};
LineWidth = 2;

if num_ac == 2
    plottool(1,'CrissCross',18,'Time (sec)','Distance (ft)')
    dx = zeros(r,1);
    dy = zeros(r,1);
    f = 3;
    for jj = 1:num_ac
        f = f - 2;
        ii = (jj-1)*NSTATES+1;
        %plot3(state(:,ii),state(:,ii+1),state(:,ii+2),[colors{ctr1},linetype{ctr2}],'LineWidth',LineWidth)
        dx = dx + f*state(:,ii);
        dy = dy + f*state(:,ii+1);
    end
    plot(data(:,1),sqrt(dx.^2+dy.^2),'b-','LineWidth',2)
end

%%%%Plot Trajectories
h1 = plottool(1,'Trajectories',12,'X(ft)','Y(ft)','Z(ft)','',[-27 30]);

if IMOVIE
    if IMOVIE == 1
        tend = r;
    else
        tend = IMOVIE;
    end
    skip = 5;
    [DYNSKIP,BOUNDARY,IKIND,K35] = read_boundary('Input_Files/MaaWM.SIM');
else
    tend = 1;
    skip = 1;
end
ctr = 1;
maxdigits = length(num2str(tend));
for tt = 1:skip:tend
    cla;
    ctr1 = 0;
    ctr2 = 1;
    for jj = 1:num_ac
        ctr1 = ctr1 + 1;
        if ctr1 > 7
            ctr1 = 1;
            ctr2 = ctr2 + 1;
        end
        ii = (jj-1)*NSTATES+1;
        if PLOTENTIRETRAJECTORY
            plot3(state(:,ii),state(:,ii+1),state(:,ii+2),[colors{ctr1},linetype{ctr2}],'LineWidth',LineWidth)
        end
        if IMOVIE
            try
                draw_plane(state(tt,ii:(ii+5)),colors{ctr1},[100 2])
            catch me
                x = state(tt,ii);
                y = state(tt,ii+1);
                z = state(tt,ii+2);
                mycube(0.1*max(x),0.1*max(y),0.1*max(z),x,y,z,0,0,0,colors{ctr1})
            end
        end
    end
    if IMOVIE
        view(0,90)
        rectangle('Position',[-BOUNDARY -BOUNDARY 2*BOUNDARY 2*BOUNDARY])
        xlim([-BOUNDARY*1.2 BOUNDARY*1.2])
        ylim([-BOUNDARY*1.2 BOUNDARY*1.2])
    else
        axis equal
    end
    title(num2str(data(tt,1)))
    reverse(['y','z'])
    drawnow
    if SAVEFIG
        filename = getfilename(ctr,maxdigits);
        fullfilename = ['Frames/File',filename,'.jpg']
        saveas(gcf,fullfilename, 'jpg')
        ctr = ctr + 1;
    end
end
