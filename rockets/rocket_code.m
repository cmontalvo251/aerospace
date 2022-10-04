function rocket_code() %%%This is the name of the file. Don't change this
global Force %%We're going to make Force a global variable so we can use it between functions

%%%What is the initial velocity of the rocket?
velocity_initial = 0; %%Let's assume it's in space and not moving
%%%What is velocity as a function of time?
%%%First we need a time vector
dt = 0.1; %%%This is our numerical timestep. The smaller we make it. 
%%The more accurate our simulation is but also the longer it takes to
%%simulate
time_vector = 0:dt:1000; %%%Here we are starting at 0 and stopping at 1000 seconds

%%%We're going to use the built in toolbox ode45 to compute velocity as a
%%%function of time. You have to give it the Equations of Motion, the time
%%%vector and the initial velocity.
[t_out,velocity_out] = ode45(@Equations_of_Motion,time_vector,velocity_initial);

%%%Then let's plot the velocity as a function of time
close all %%%Close all figures just in case
figure() %%Make a new figure
plot(t_out,velocity_out,'b-','LineWidth',2)
xlabel('Time (sec)')
ylabel('Velocity (m/s)')
grid on

%%%Plotting the Force as a function of time is a bit more complex
%%First we have to make a vector of zeros
force_out = 0*t_out; %%This is the easiest way to do this
%%%Then we need to compute force at every time
for idx = 1:length(t_out) %%THis is a loop that loops through every time
    %%%This computes the equations of motion at every point in time
    Equations_of_Motion(t_out(idx),velocity_out(idx));
    %%%We then save the force into our force vector
    force_out(idx) = Force;
end

%%%Now we can plot Force
figure() %%%This makes a new figure
plot(t_out,force_out,'b-','LineWidth',2)
xlabel('Time (sec)')
ylabel('Force (N)')
grid on


%%%These are the equations of motion of the rocket
%%%The output of the function if acceleration and the inputs are time and
%%%velocity
function acceleration = Equations_of_Motion(time,velocity) 
global Force

%%%This says force is 100 N when the time is less than 510 secondss
if time < 510
    Force = 100; %%For now let's just assume we're in space and we fire a rocket for 
else
    %%%And zero when time is greater than 510 seconds
    Force = 0;
end

%%%Mass of the rocket
mass = 1000; %%%1000 kilograms for now

%%%Acceleration of the rocket using Newton's second law
acceleration = Force/mass;