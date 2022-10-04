clear
clc
close all

A = imread('file.png');

A = A(37:615,1:1300,:);

[r,c,d] = size(A);

for ii = 1:r
    for jj = 1:c
        red = A(ii,jj,1);
        green = A(ii,jj,2);
        blue = A(ii,jj,3);
        if green == red && red == blue
            A(ii,jj,1) = 0;
            A(ii,jj,3) = 0;
            A(ii,jj,2) = 0;
        end
        if green > red && green > blue
            A(ii,jj,1) = 255;
            A(ii,jj,3) = 255;
            A(ii,jj,2) = 255;
        end
    end
end

imshow(A);

Agrey = rgb2gray(A);

figure()
imshow(Agrey)

[r,c] = size(Agrey);
xcoordinates = 1:c;
ycoordinates = zeros(1,length(xcoordinates));
for x = 1:c
    %%%At this current x coordinate search for my y coordinate
    slice = Agrey(:,x);
    y = find(slice == 255,1);
    if ~isempty(y)
        ycoordinates(x) = y;
    end
end
L = ycoordinates == 0;
ycoordinates(L) = [];
xcoordinates(L) = [];

figure()
plot(xcoordinates,ycoordinates,'b*')

xval = 97;
scalex = 11.8/(1244-xval);

yval = 423;
scaley = 90/(25-yval);

xshift = (xcoordinates - xval)*scalex;
yshift = (ycoordinates - yval)*scaley;

fig = figure();
plot(xshift,yshift,'g-')
grid on
xlabel('Time (s)')
ylabel('Speed (m/s)')
title('NO CALL - 1254/8')
legend('Speed')
set(fig,'color','white')



    



