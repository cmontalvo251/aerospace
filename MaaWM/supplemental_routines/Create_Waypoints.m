function Create_Waypoints(num_Aircraft,Boundary,outfile,iplot,freq,pathtype,mult,ALT)
close all

if nargin < 7
    mult = 0;
end

fh = fopen( outfile, 'wt');

names = {' !X(ft)',' !Y(ft)',' !Z(ft)',' Flag(0 or 1)'};
for ii = 1:length(names)
    fprintf(fh,names{ii});
end
fprintf(fh,'\n');

if iplot
    plottool(1,'Trajectories',12,'X(ft)','Y(ft)','Z(ft)','',[-27 30]);
end

colors = {'b','r','g','k','y','m','b','r','g','k','y','m'};

for ii = 1:num_Aircraft

    switch pathtype
        case 'sine' %%%Sinusoid
            ph = (ii-1)*pi/2;
            inc = (pi/freq)/2;
            T = 0:inc:(4*pi);
            X = 1.0*Boundary*cos(freq*T);
            Y = linspace(-Boundary,Boundary,length(T))';
            Z = (-200 - ALT*(ii-1))*ones(length(Y),1)*3.28;
        case 'spiral'
            GridDistance = Boundary*2;
            Spiral_Reduction = GridDistance/freq;
            xy = [Boundary;-Boundary];
            XY = [[-Boundary;-Boundary],xy];
            Z = [-200;-200]*3.28;
            Arot = [0 -1;1 0];
            nhat = Arot*[1;0];
            ctr = 1;
            while GridDistance >= 0
                xy = xy + GridDistance*nhat;
                nhat = Arot*nhat;
                ctr = ctr + 1;
                if ctr > 2
                    ctr = 1;
                    GridDistance = GridDistance - Spiral_Reduction;
                end
                XY = [XY,xy];
                Z = [Z;-200*3.28];
            end
            X = XY(1,:);
            Y = XY(2,:);
%             S = spiral(round(freq));
%             [r,c] = size(S)
%             m = max(max(S));
%             X = zeros(m,1);
%             Y = X;
%             Z = -200*ones(m,1);
%             for idx = m:-1:1
%                 [row,col] = ind2sub(size(S),find(S==idx));
%                 X(idx) = -Boundary + (row-1)*2*Boundary/(r-1);
%                 Y(idx) = -Boundary + (col-1)*2*Boundary/(c-1);
%             end
%             Lremove = [];
%             Xi = X(1);
%             for idx = 2:length(X)
%                 if X(idx) == Xi
%                     Lremove = [Lremove;idx];
%                 end
%                 Xi = X(idx);
%             end
%             Lremove
        otherwise
            X = 0;
            Y = 0;
            Z = -200*3.28;
    end
    
    if mult
       se_pts = round(linspace(1,length(T),num_Aircraft+1));
       s = se_pts(ii);
       e = se_pts(ii+1);
       X = X(s:e);
       Y = Y(s:e);
       Z = Z(s:e);
    end
    
    if iplot
        ctr = ii;
        while ctr > length(colors)
            ctr = ctr-length(colors);
        end
        plot3(X,Y,Z,[colors{ctr},'*-'],'LineWidth',2)
        reverse('y')
    end
    
    flag = ones(length(Z),1);
    
    flag(end) = 0; %%%Tells the code to stop and move to the next aircraft
    
    fprintf(fh,[num2str(ii),' !Aircraft Number \n']);

    for idx = 1:length(X)
        fprintf(fh,[num2str(X(idx)),' ',num2str(Y(idx)),' ',num2str(Z(idx)),' ',num2str(flag(idx)) '\n']);
    end
end

fclose(fh);

