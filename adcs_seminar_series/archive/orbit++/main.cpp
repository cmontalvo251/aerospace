// Multiple CubeSat Simulation Software 
// Copyright - Carlos Montalvo 2015
// You may freely distribute this file but please keep my name in here
// as the original owner

//To compile this on linux just type 
//$ make
//On windows you need to get MinGW then type
//C:\ mingw32-make.exe

#include "main.h"
#include "satellite_physics.h"
#include "visualizer.h"
#ifndef __linux__ //This header file is not needed for Linux
#include "define.h" //this is used mainly for Visual Studio, as the makefile doesn't affect VS projects
#endif 
#include <fstream>

using namespace std;

#ifdef BOOST_FLAG
//In order to get this to work on linux you need to run sudo apt-get install libboost-all-dev
#include <boost/thread.hpp> 
using namespace boost;
#endif

Satellite* p_USASat[10]; //Pointer to satellite physics engine - Allocate for a maximum of 10 satellites
Satellite* p_USASatModel[10]; //Pointer to model physics engine
Controller* p_Control[10]; //pointer to control routine

int main(int argc,char** argv) {

  //Seed random number generator to get different results everytime
  srand(time(NULL));

  //This will also read in Satellite.SIM and pull in the number of satellites
  //To simulate
  Engine* p_Engine = new Engine(); //Pointer to physics engine
  printf("Physics Engine Online \n");

  #ifdef OPENGL_FLAG
  if (p_Engine->NSATS != 1) {
    p_Engine->NSATS = 1;
    printf("Sorry OpenGL Only works with 1 satellite at the moment. Make sure to fix this in Satellite.SIM \n");
    printf("For now your first satellite will be created and plotted. \n");
    PAUSE();
  }
  #endif

  //This loop creates all of the satellites and controllers
  for (int idx = 0;idx<p_Engine->NSATS;idx++) {
    //Create a new control variable
    p_Control[idx] = new Controller(p_Engine->Magnetic_Flag); 
    printf("Controller %d Created \n",idx);
    //Here we will create each satellite object and pass it to the pointer
    MATLAB state_XYZ0,state_PTP0;
    Create_Initial_Conditions(&state_XYZ0,&state_PTP0,idx);
    printf("Creating Satellite %d \n",idx);
    state_XYZ0.disp();
    state_PTP0.disp();
    p_USASat[idx] = new Satellite(0,state_XYZ0,state_PTP0,idx); //the zero means leave errors off since this is the real satellite
    //Pass the magnetic and gravity parameters
    p_USASat[idx]->Gravity_Flag = p_Engine->Gravity_Flag;
    p_USASat[idx]->Magnetic_Flag = p_Engine->Magnetic_Flag;
    p_USASat[idx]->yr = p_Engine->yr;
    printf("Satellite %d Created \n",idx);
  }

  //Now we kick off all of the threads
  #ifdef BOOST_FLAG
  boost::thread runSatellite(SatellitePhysicsThread,p_USASat,p_Control,p_Engine);  
  #else
  SatellitePhysicsThread(p_USASat,p_Control,p_Engine);
  #endif

  //We are going to kick off an openGL plotting routine
  #ifdef OPENGL_FLAG
  boost::thread runVisualizer(Visualizer,p_USASat,argc,argv);
  #endif

  //We will then kick off the Model with the ierror flag set to 1.
  //Since I moved the physics engine to the class Engine I do not think
  //That this will work anymore so hold on on KALMAN or SENSOR_FLAG
  //#ifdef KALMAN_FLAG
  //Before we kick off the model we need to sample model using sensors. Probably need to make
  //a sensor class and maybe a sensor.cpp file - REVISIT
  //#ifdef SENSOR_FLAG
  //GPS(); //- GPS
  //SS(); //- Solar Sensor
  //RG(); //- Rate Gyro - p,q,r
  //MAG(); //- magnetometer
  //ACCEL(); //- accelerometer
  //#endif 
  //p_USASatModel = new Satellite(1,state_XYZ0,state_PTP0);
  //#ifdef BOOST_FLAG
  //boost::thread KalmanFilter(SatellitePhysicsThread,p_USASatModel,p_Control,p_Engine);  
  //#else
  //SatellitePhysicsThread(p_USASatModel,p_Control,p_Engine);
  //#endif
  //#endif

  #ifdef BOOST_FLAG
  //Wait two seconds to make sure all threads have kicked off properly
  cross_sleep(2);
  #endif
    
  //Main loop will run continuously unless someone hits the 'q' button
  #ifdef SIMFOREVER_FLAG
  char key = ' ';
  while (key != 'q') 
    {
       #ifdef OPENGL_FLAG
       cout << "Press ESC to quit" << endl;
       #else
       #ifdef BOOST_FLAG
       cout << "Press q then enter to quit\n ";
       #endif
       #endif
       cin >> key;
    }
  #else
  while (p_Engine->simOK()){
    cross_sleep(1);
  }
  #endif

  #ifdef BOOST_FLAG
  //Destroy all threads
  printf("Destroying All Threads...\n");
  pthread_cancel(runSatellite.native_handle());
  printf("Function Quit \n");
  #endif

  #ifdef OPENGL_FLAG
  pthread_cancel(runVisualizer.native_handle());
  #endif

  return 0;

}

void Create_Initial_Conditions(MATLAB* state_XYZ0,MATLAB* state_PTP0,int satIDX) {

  //This creates two empty vectors for us to put stuff in here.
  //But only if we haven't initialized these yet
  if (state_XYZ0->getName() == 0) {
    state_XYZ0->zeros(NUMXYZ,1,"state_XYZ0");
    state_PTP0->zeros(NUMPTP,1,"state_PTP0");
  }

  //Alright instead of hardcoding all of the stuffs here we will import a file.
  //What do we want to call the file?
  //Satellite#.INIT --- sure. How do we deal with multiple satellites? Throw a number in there
  //For now just hard code the zero.
  ifstream initfile;
  char filename[256];
  sprintf(filename,"%s%d%s","Input_Files/Satellite",satIDX,".INIT");
  initfile.open(filename);
  if (initfile.is_open()) {
    printf("Initial Conditions File Found: %s \n",filename);
    //Alright let's grab everything we need
    //If you look in mathp.h the values MEARTH, REARTH and some other Earth parameters are defined as globals which means
    //we can use them anywhere
    //So to set the initial condition we will write
    //Make some dummy variables
    string input;
    double Distance_of_Orbit_from_Surface; //So 600 * 1000 is 600 km above the surface
    double floatvalue;   
    //Grab first line
    getline(initfile,input);
    Distance_of_Orbit_from_Surface = atof(input.c_str());
    if (Distance_of_Orbit_from_Surface == -99) {
      //This means that this is simply a satellite model rather than an orbit model
      getline(initfile,input);
      state_XYZ0->set(1,1,atof(input.c_str()));
    }
    else {
      //Then we will say the initial condition of the orbit is just REARTH + Distance_of_Orbit_from_Surface
      state_XYZ0->set(1,1,REARTH+Distance_of_Orbit_from_Surface*1000); //x0 meters
    }
    //Grab second line
    getline(initfile,input);
    floatvalue = atof(input.c_str());
    state_XYZ0->set(2,1,floatvalue);           //y0  
    //Grab third line
    getline(initfile,input);
    floatvalue = atof(input.c_str());
    state_XYZ0->set(3,1,floatvalue);           //z0
    //Grab fourth line
    getline(initfile,input);
    floatvalue = atof(input.c_str());
    state_XYZ0->set(4,1,floatvalue);           //%%%meters/second  
    if (Distance_of_Orbit_from_Surface == -99) {
      getline(initfile,input);
      state_XYZ0->set(5,1,atof(input.c_str())); //ydot
      getline(initfile,input);
      state_XYZ0->set(6,1,atof(input.c_str())); //zdot
    } else {
      //Compute Velocity based on distance from Surface
      double Vsat = sqrt(MEARTH*GSPACE/state_XYZ0->get(1,1));
      //need to figure out how to fix this or rather how to get this to work
      //Ok well this works now. 
      double inclination_angle; 
      getline(initfile,input);
      inclination_angle = atof(input.c_str());
      state_XYZ0->set(5,1,Vsat*cos(inclination_angle*PI/180.0)); //ydot - this assumes that the z axis of the earth is going through
      state_XYZ0->set(6,1,Vsat*sin(inclination_angle*PI/180.0)); //zdot - the north pole
    }
    //Grab Euler angles 
    MATLAB q0123,ptp;
    q0123.zeros(4,1,"q0123");
    ptp.zeros(3,1,"ptp");
    //Roll 
    getline(initfile,input);
    floatvalue = atof(input.c_str());    
    ptp.set(1,1,floatvalue*PI/180.0); //phi0 rad
    //Pitch
    getline(initfile,input);
    floatvalue = atof(input.c_str());    
    ptp.set(2,1,floatvalue*PI/180.0); //theta0 rad
    //Yaw
    getline(initfile,input);
    floatvalue = atof(input.c_str());    
    ptp.set(3,1,floatvalue*PI/180.0); //psi0 rad
    //Convert to quaternions
    q0123.euler2quat(ptp);  
    //If debug mode is on print the quaternions
    #ifdef DEBUG_MODE
    q0123.disp();
    PAUSE();
    #endif
    //Send the quaternion state to the state_PTP0 vector
    for (int idx = 1;idx<=4;idx++)
    {
      state_PTP0->set(idx,1,q0123.get(idx,1));
    }
    //
    //#ifdef DEBUG_MODE -- Random or zero rates
    //srand(time(NULL)); //seeds the RNG
    // state_PTP0->set(5, 1, rand() % 7); //sets PQR to random integers between and including 0 and 6 The percent sign is the modulo function
    // state_PTP0->set(6, 1, rand() % 7); //If you want to set this to zero comment out the lines below  
    // state_PTP0->set(7, 1, rand() % 7);
    //state_PTP0->set(5,1,0);
    //state_PTP0->set(6,1,0);
    //state_PTP0->set(7,1,0);
    //state_PTP0->disp();
    //PAUSE();
    //#endif

    //Get rates
    //Roll Rate
    getline(initfile,input);
    floatvalue = atof(input.c_str());
    state_PTP0->set(5,1,floatvalue);
    //Pitch Rate
    getline(initfile,input);
    floatvalue = atof(input.c_str());
    state_PTP0->set(6,1,floatvalue);
    //Yaw Rate
    getline(initfile,input);
    floatvalue = atof(input.c_str());    
    state_PTP0->set(7,1,floatvalue);
    initfile.close();
    //worst-case scenario is 5.6 rad/s so I've got them set to 2 rad/s here
    //state_PTP0->set(5,1,2); //state_PTP0 is q0,q1,q2,q3,P,Q,R
    //state_PTP0->set(6,1,2); //order is 1-7
    //state_PTP0->set(7,1,2); //P,Q,R is 5-7
  }
  else {
    printf("%s%s%s","Sorry ",filename," could not be found \n");
    exit(1);
  }


}
