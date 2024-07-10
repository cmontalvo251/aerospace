function Create_Initial_Conditions(aircraft_types, Boundary_Size,outfile,ALT)

%%%The variable aircraft_types will be a vector with two numbers. The first number
%%%will be the number of airplanes and the second number will be the number of quadcopters
num_Aircraft = sum(aircraft_types);
num_Airplanes = aircraft_types(1);
num_Quads = aircraft_types(2);

%%%Position
AircraftXY = Boundary_Size.*(-0.5+rand(2,num_Aircraft));
AircraftZ = -200*3.28*ones(1,num_Aircraft) - ALT*(0:(num_Aircraft-1));

%%%Attitude
phi = zeros(1,num_Aircraft);
theta = phi;
Heading = (2*pi).*rand(1,num_Aircraft);

%%%Velocity
u = [20*3.28*ones(1,num_Airplanes),zeros(1,num_Quads)];
v = 0*u;
w = 0*u;

%%%%Angular rates
p = 0*u;
q = p;
r = p;

%%%States
states = [AircraftXY;AircraftZ;phi;theta;Heading;u;v;w;p;q;r];

names = {' !X(ft)',' !Y(ft)',' !Z(ft)',' !Phi(rad)',' !Theta(rad)',' !Psi(rad)', ...
	 ' !U(ft/s)',' !V(ft/s)',' !W(ft/s)',' !P(rad/s)',' !Q(rad/s)',' !R(rad/s)'};

fh = fopen( outfile, 'wt');

fprintf(fh,[num2str(num_Aircraft),' !Number of Aircraft \n']);
fprintf(fh,[num2str(num_Airplanes),' !Number of Airplanes \n']);

for idx = 1:num_Aircraft
	if idx <= num_Airplanes
    	fprintf(fh, 'Aircraft ');
    else
    	fprintf(fh, 'Quad ');
    end
    %fprintf(fh, num2str(idx));    
    fprintf(fh, '\n');
    for jdx = 1:12
      fprintf(fh, [num2str(states(jdx,idx)),names{jdx}]);
      fprintf(fh, '\n');
    end
end

fclose(fh);