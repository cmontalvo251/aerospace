addpath 'supplemental_routines/'

purge

global NSTATES

montefile = fopen('MonteCarlo.txt','wt');
%fprintf(montefile,'%s \n','Freq Su Sv Sw Avg_ErrorU Avg_ErrorV Avg_ErrorW MaxUE MaxVE MaxWE DataPts');
fprintf(montefile,'%s \n','NUM_AIRCRAFT Su Sv Sw Avg_ErrorU Avg_ErrorV Avg_ErrorW MaxUE MaxVE MaxWE DataPts Tend(sec)');
system('rm Output_Files/RCOND.txt');
RCONDFILE = fopen('Output_Files/RCOND.txt','wt');
fprintf(RCONDFILE,' ');

%%%%Number of turns for spiral = freq*2 so if you want 10 turns freq = 5
%freq = [4 28 52 76]/2; %%%Make frequency a vector and it will do multiple runs
freq = 52/2;
%freq = 35; %%%for grid space
%freq = inf;
%data_pts = [10:1:20]; %%%This is the number of data points along each axis
%data_pts = ceil(sqrt(([16:2:64].^2)*1.4));
%data_pts = [4:8:76]/2;
%data_pts = 35

AIRCRAFT_VEC = [1 0]; %%%The first number is the number of aircraft and the second number is the number of quadcopters
NSTATES = 12;

[num_Runs,c] = size(AIRCRAFT_VEC);

for idx = 1:num_Runs

  AIRCRAFT_TYPES = AIRCRAFT_VEC(idx,:);
  NUM_AIRCRAFT = sum(AIRCRAFT_TYPES);


rs = [0 0 0];

%for input_arg = 1:10
%while sum(abs(rs)) < 2.7
%for input_arg = data_pts %%%make sure to change this to either freq or data_pts depending on if you are doing a sine/spiral or grid points
for input_arg = freq

close all

MULT = 0; %%%I think this makes the waypoints for every aircraft non-intersecting
ALT = 0; %%DIFFERENCE IN ALTITUDE OF AIRCRAFT
RUNCODE = 1;

CREATE_ICS = 0; %%This will create new initial conditions for the aircraft (Does this work for quad too?)

CREATE_WAYPOINTS = 0; %%%This is to create new waypoints
PATHTYPE = 'sine'; %%%Other options are spiral

CHANGE_UPDATE_RATE = 0; %%%Change the sampling rate of the aircraft

CREATE_GRID_SPACE = 0;
GRIDTYPE = 'WRF'; %%%THIS IS GRIDTYPE OF WIND MODEL. SET TO WRF to simply have the routine grab the .ATM file dataloc

%%%what do you want to plot?
PLOTSTATES = 0;
PLOTCONTROLS = 0;
PLOTWINDS = 0;
PLOTTRAJS = 0; %%%This will plot the entire trajectory (However for some reason you need to set PLOTENTIRETRAJECTORY = 1)

IMOVIE = 0; %%%This should play a movie of the aircraft flying. Only set to 1 if you are using MATLAB as opposed to Octave
PLOTENTIRETRAJECTORY = 0; %%%This will plot the trajectory of the aircraft (IMOVIE must be set to 1 for this to work)

SAVEFIG = 0;

COMPUTEFIT = 1; %%plot the RBF results
IMESH = 1; %%%plot meshes as opposed to single data points
IMASK = 0; %%%This mask is to get the number of data points for a certain square area
OFFSET = 100*3.28; %%%ADD A SAFETY FACTOR FOR FITS meters to feet
%%%(Typically set to 100. But if you want to characterize some error as a 
%%function of spatial distribution set this to zero or even -50) 
%%%Actually set it to -50

COMPUTE_EIGS = 0; %%%compute eigenvalues of aircraft 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%YOU SHOULDN'T HAVE TO CHANGE ANYTHING BELOW HERE (UNLESS THERE'S A BUG)%%%%%%%%%%%%%%%%%%%%%

%%%SOME CHECKS
if IMOVIE
    PLOTTRAJS = 1;
end

%%%Read the boundary 
[DYNSKIP,BOUNDARY,IKIND,K35] = read_boundary('Input_Files/MaaWM.SIM');

dataloc = getdataloc('Input_Files/MaaWM.ifiles');
UVWfrontend

%%%Create a Grid Space of Points
if CREATE_GRID_SPACE
    %%%Note this routine will auto generate the  Winds.OUT file
    %%%This is the same file that the RBF routine will read in
    %%%If you look in the .SIM file there is a variable that says skip
    %%%dynamic model. Set that to 1 if you want to only run the RBF model
    %%%%Note that if you want to define your own dataloc you need to pass
    %%%%it to this routine as an argument and set GRIDTYPE to USER
    %%%%NVM. YOU DO NOT NEED TO DO THIS ANYMORE. I have created a routine
    %%%%called getdataloc.m which will get the dataloc automatically.
  freqx = grid_space(input_arg,BOUNDARY,GRIDTYPE,K35,IKIND);
end

%%%%Change UPDATE RATE
if CHANGE_UPDATE_RATE 
    write_to_csfile(input_arg)
end

if CREATE_ICS && RUNCODE
    Create_Initial_Conditions(AIRCRAFT_TYPES,BOUNDARY,'Input_Files/Aircraft.ICS',ALT);
end

if CREATE_WAYPOINTS && RUNCODE
    Create_Waypoints(NUM_AIRCRAFT,BOUNDARY,'Input_Files/Waypoints.WAY',1,input_arg,PATHTYPE,MULT,ALT);
end

if RUNCODE
  %system('make'); %%This won't work on Windows. You need to simulate ahead of time
  if ~DYNSKIP
      system('rm Output_Files/*.OUT');
  end
  system('./Run.exe Input_Files/MaaWM.ifiles'); %%%This won't work on windows either
  %system('Run.exe Input_Files/MaaWM.ifiles'); %%%Try this for WINDOWS
end

%%%%Plot data
filename = 'Output_Files/State.OUT';
if PLOTTRAJS
  plottrajectories(filename,IMOVIE,PLOTENTIRETRAJECTORY,SAVEFIG);
end
if PLOTSTATES
  plotstates(1:12,filename);
end

%%%Verify Eigenvalues of A matrix
if COMPUTE_EIGS
  A = dlmread('Output_Files/A.txt');
  A = A(1:12,1:12);
  eig(A)
end

if PLOTWINDS
    windname = 'Output_Files/Winds.OUT';
    plotwinds(windname);
end

if PLOTCONTROLS
    controlname = 'Output_Files/Controls.OUT';
    plotcontrols(controlname);
end

if COMPUTEFIT
    Analyze_Fit;
end

end
end

fclose(montefile);

system('cat MonteCarlo.txt');
