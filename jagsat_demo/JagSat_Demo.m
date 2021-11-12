%Matlab has to re-define the serial port every time it runs for some reason
if exist('s','var') == 1
    delete(s)
end

s = 1

s = serialport("COM5",38400,"FlowControl","none"); %Define the serial port to pull gyro/magnetometer measurements
pause(3)

t = 0; %Initialize all parameters
thetax = 0;
thetay = 0;
thetaz = 0;
i = 0;

t = 0;

tic; %Start clock
while 1
    writeline(s,'a') %The gyro is set up to print the gyro measurements when sent a character
    
    chdata = char(readline(s)); %Take the gyro measurements in as a character string
    gyro = str2num(chdata); % Convert that string to numbers (double for some reason won't work here)
    
    if size(gyro) > 0 %The first iteration of this code yields a size zero array and crashes the code, so this "if" takes care of the nonzero first array
        toc; %Stop clock
        t1 = toc; %Define a time unit based on how long the clock was running
        t = t + t1;
        tic; %Start the clock again for the next loop
        
        %Rounds all measurements down
        gyro(1) = round(gyro(1),1); 
        gyro(2) = round(gyro(2),1);
        gyro(3) = round(gyro(3),1);

        gyro = gyro*180./pi;

        %%%ADD JITTER
        %gyro(1) = gyro(1) + 30.0*sin(40*t);
        %gyro(2) = gyro(2) + 30.0*sin(40*t);
        %gyro(3) = gyro(3) + 30.0*sin(40*t);

        %gyro(1) = gyro(1) + 10;
        %gyro(2) = gyro(2) - 20;
        %gyro(3) = gyro(3) + 30.*sin(t);

        thetax = thetax + (gyro(1)*t1); %Multiply the processed gyro data by the time unit, converts it to degrees, and adds it to the previous angle      
        thetay = thetay + (gyro(2)*t1);     
        thetaz = thetaz + (gyro(3)*t1);
        
        if i >= 20 %Plot the shape every twenty iterations--without this step, it tries to draw the shape too quickly and lags severely

            h1 = trimesh(model);
            h1.FaceColor = 'b';
            h1.EdgeColor = [0.0 0.0 0.1];
            title(['DEGREES - X: ', num2str(thetax), ' Y: ', num2str(thetay), ' Z: ', num2str(thetaz),' DEGREES/SEC X: ',num2str(gyro(1)),' Y: ',num2str(gyro(2)),' Z: ',num2str(gyro(3))]);
            rotate(h1,[1 0 0],thetaz);
            rotate(h1,[0 1 0],-thetay);
            rotate(h1,[0 0 1],thetax);
            axis equal

            i = 0; %Reset the number of iterations
        else
            i = i + 1; %Count the number of iterations up
        end
    end
end