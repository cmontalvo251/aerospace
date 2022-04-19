clear all
close all
clc

% Re = 100,000
% clcdfiles = {'NACA64AclcdRe100k.txt'};
% clcmfiles = {'NACA64AclcmRe100k.txt'};
% clcdfiles = {'NACA6409clcdRe100k.txt'};
% clcmfiles = {'NACA6409clcmRe100k.txt'};
clcdfiles = {'A18cleanclcdRe100k.txt', 'BW3clcdRe100k.txt', 'ClarkYclcdRe100k.txt',...
    'E387clcdRe100k.txt', 'NACA6409clcdRe100k.txt', 'S1223clcdRe100k.txt', 'S7075clcdRe100k.txt'};
clcmfiles = {'A18cleanclcmRe100k.txt', 'BW3clcmRe100k.txt', 'ClarkYclcmRe100k.txt',...
    'E387clcmRe100k.txt', 'NACA6409clcmRe100k.txt', 'S1223clcmRe100k.txt', 'S7075clcmRe100k.txt'};


for i = 1: length(clcdfiles)

    aclcd = dlmread(clcdfiles{i});
    aclcm = dlmread(clcmfiles{i});

    %Airfoils w/ lift, drag, and moment data

    alpha1 = aclcm(:,1)*pi/180;
    cl1 = aclcm(:,2);
    cm = aclcm(:,3);

    alpha2 = aclcd(:,1)*pi/180;
    cl2 = aclcd(:,2);
    cd = aclcd(:,3);
    ltod = cl2./cd;
    powerrat = cl2.^(3/2)./cd;
    
    figure(1)
    hold all
    plot(alpha1*180/pi, cl1, '.')
    figure(2)
    hold all
    plot(alpha1*180/pi, cm, '.')
    figure(3)
    hold all
    plot(alpha2*180/pi, cd, '-o')
    figure(4)
    hold all
    plot(alpha2*180/pi, ltod, '-o')
    figure(5)
    hold all
    plot(alpha2*180/pi, powerrat, '-o')
end

figure(1)
xlabel('Angle of Attack [deg]')
ylabel('Lift Coeff')
legend('A18', 'BW3', 'ClarkY', 'E387', 'NACA6409', 'S1223', 'S7075')
figure(2)
xlabel('Angle of Attack [deg]')
ylabel('Pitching Mom. Coeff')
legend('A18', 'BW3', 'ClarkY', 'E387', 'NACA6409', 'S1223', 'S7075')
figure(3)
xlabel('Angle of Attack [deg]')
ylabel('Drag Coeff')
legend('A18', 'BW3', 'ClarkY', 'E387', 'NACA6409', 'S1223', 'S7075')
figure(4)
xlabel('Angle of Attack [deg]')
ylabel('Lift to Drag Ratio')
legend('A18', 'BW3', 'ClarkY', 'E387', 'NACA6409', 'S1223', 'S7075')
figure(5)
xlabel('Angle of Attack [deg]')
ylabel('Power Ratio')
legend('A18', 'BW3', 'ClarkY', 'E387', 'NACA6409', 'S1223', 'S7075')

%Re = 200,000
% clcdfiles = {'NACA64AclcdRe200k.txt'};
% clcmfiles = {'NACA64AclcmRe200k.txt'};
% clcdfiles = {'NACA6409clcdRe200k.txt'};
% clcmfiles = {'NACA6409clcmRe200k.txt'};
clcdfiles = {'A18cleanclcdRe200k.txt', 'BW3clcdRe200k.txt', 'ClarkYclcdRe200k.txt',...
    'E387clcdRe200k.txt', 'NACA6409clcdRe200k.txt','S1223clcdRe200k.txt', 'S7075clcdRe200k.txt'};
clcmfiles = {'A18cleanclcmRe200k.txt', 'BW3clcmRe200k.txt', 'ClarkYclcmRe200k.txt',...
    'E387clcmRe200k.txt', 'NACA6409clcmRe200k.txt', 'S1223clcmRe200k.txt', 'S7075clcmRe200k.txt'};

for i = 1: length(clcdfiles)

    aclcd = dlmread(clcdfiles{i});
    aclcm = dlmread(clcmfiles{i});

    %Airfoils w/ lift, drag, and moment data

    alpha1 = aclcm(:,1)*pi/180;
    cl1 = aclcm(:,2);
    cm = aclcm(:,3);

    alpha2 = aclcd(:,1)*pi/180;
    cl2 = aclcd(:,2);
    cd = aclcd(:,3);
    ltod = cl2./cd; 
    powerrat = cl2.^(3/2)./cd;
    
    figure(6)
    hold all
    plot(alpha1*180/pi, cl1, '.')
    figure(7)
    hold all
    plot(alpha1*180/pi, cm, '.')
    figure(8)
    hold all
    plot(alpha2*180/pi, cd, '-o')
    figure(9)
    hold all
    plot(alpha2*180/pi, ltod, '-o')
    figure(10)
    hold all
    plot(alpha2*180/pi, powerrat, '-o')

end

figure(6)
xlabel('Angle of Attack [deg]')
ylabel('Lift Coeff')
legend('A18', 'BW3', 'ClarkY', 'E387', 'NACA6409', 'S1223', 'S7075')
figure(7)
xlabel('Angle of Attack [deg]')
ylabel('Pitching Mom. Coeff')
legend('A18', 'BW3', 'ClarkY', 'E387', 'NACA6409', 'S1223', 'S7075')
figure(8)
xlabel('Angle of Attack [deg]')
ylabel('Drag Coeff')
legend('A18', 'BW3', 'ClarkY', 'E387', 'NACA6409', 'S1223', 'S7075')
figure(9)
xlabel('Angle of Attack [deg]')
ylabel('Lift to Drag Ratio')
legend('A18', 'BW3', 'ClarkY', 'E387', 'NACA6409', 'S1223', 'S7075')
figure(10)
xlabel('Angle of Attack [deg]')
ylabel('Power Ratio')
legend('A18', 'BW3', 'ClarkY', 'E387', 'NACA6409', 'S1223', 'S7075')


% 
% Cl_alpha = polyfit(alpha1(1:17), cl1(1:17), 1);
% a0 = Cl_alpha(1)
% alpha0 = interp1(cl1(1:17),alpha1(1:17),  0)
% 
% Cd_alpha = polyfit(alpha2, cd, 2)
% Cd0 = Cd_alpha(3)
% % cd = Cd_alpha(1)*(3*pi/180)^2 + Cd_alpha(2)*(3*pi/180) + Cd_alpha(3)
% 
% cm_ac = mean(cm)


