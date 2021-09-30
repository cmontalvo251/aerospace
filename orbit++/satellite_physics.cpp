// Copyright - Carlos Montalvo 2015
// You may freely distribute this file but please keep my name in here
// as the original owner

#include "main.h"
#include "satellite_physics.h"
#include "../timer.h" 

#ifndef __linux__ //This define.h header is not needed for linux
#include "define.h" //this is used mainly for Visual Studio, as the makefile doesn't affect VS projects
#endif

using namespace std;
using namespace GeographicLib; //Necessary for Grav and Magnet models

//Controller Class
Controller::Controller(int Magnetic_Flag) 
{
  //Initialize Variables
  LMN.zeros(3,1,"LMN");
  state_XYZ_.zeros(NUMXYZ,1,"Controller state_XYZ");
  state_PTP_.zeros(NUMPTP,1,"Controller state_PTP");
  //All of this is archaic code from a much different control scheme
  //q0123.zeros(4,1,"Controller q0123");
  pqr.zeros(3,1,"Controller pqr");
  //ptp.zeros(3,1,"Controller ptp");
  //Kpptp.zeros(3,1,"Kpptp");
  //Kdpqr.zeros(3,1,"Kdpqr");
  //pqrc.zeros(3,1,"pqrc");
  //pqrc_pqr.zeros(3,1,"pqrc_pqr");
  //Kp = diag([1.0,1.0,1.0]);
  //Kp.zeros(3,3,"Kp");
  //Kp.set(1,1,1.0);
  //Kp.set(2,2,1.0);
  //Kp.set(3,3,1.0);
  //Kd = diag([0.5,0.5,0.5]);
  //Kd.zeros(3,3,"Kd");
  //Kd.set(1,1,0.5);
  //Kd.set(2,2,0.5);
  //Kd.set(3,3,0.5);
  //Km = 0.01*diag([1.0,1.0,0.4]);
  //Km.zeros(3,3,"Kd");
  //Km.set(1,1,1.0);
  //Km.set(2,2,1.0);
  //Km.set(3,3,0.4);
  //Km.mult_eq(0.02);
  BVEC_.zeros(3, 1, "Magnetic Field in the Body frame");
  //BVEC.set(1, 1, 300);
  //BVEC.set(2, 1, 400);
  //BVEC.set(3, 1, 500);
  BVEC_Tesla.zeros(3, 1, "Magnetic Field in Body Frame in Teslas");
  BVEC_Tesla_OLD.zeros(3, 1, "Old Magnetic Field in Body Frame in Teslas");
  MMTVEC.zeros(3, 1, "Magnetic moment");
  XMMTVEC.zeros(3, 1, "X Magnetic moment");
  YMMTVEC.zeros(3, 1, "Y Magnetic moment");
  ZMMTVEC.zeros(3, 1, "Z Magnetic moment");
  MTVEC.zeros(3, 1, "Magnetorque vector");

  //These are unit vectors to tell us what direction our magnetorque vector is in
  XMTVEC.zeros(3, 1, "X Magnetorque unit vector");
  YMTVEC.zeros(3, 1, "Y Magnetorque unit vector");
  ZMTVEC.zeros(3, 1, "Z Magnetorque unit vector");

  XMTVEC.set(1,1,1); //Set the first one to 1 and then that's it
  YMTVEC.set(2,1,1); //set this to 0,1,0
  ZMTVEC.set(3,1,1); //set this to 0,0,1

  CURRENT_VEC.zeros(3,1,"Magnetorquer Current");
  
  //PQR.zeros(3, 1, "Angular velocity");
  //P.zeros(3, 1, "P vector");
  //Q.zeros(3, 1, "Q vector");
  //R.zeros(3, 1, "R vector");

  ///Alright. All of this needs to get moved to an input file
  //What do we call the file? -- I tihnk we need two actually.
  //We need a RW file and a MAG file as well as a CS file
  //Satellite.RW, Satellite.MAG and Satellite.CS

  //Dummy variables
  string input;
  double floatvalue;

  ifstream CS_file;
  CS_file.open("Input_Files/Satellite.CS");
  if (CS_file.is_open()) {    
    //Grab everything we need
    //Update Rate of Controller
    getline(CS_file,input);
    floatvalue = atof(input.c_str());
    UpdateRate_ = floatvalue; // in seconds
    timeNext_ = 0; //initialize at zero
    
    //Get magnetorquer gain
    getline(CS_file,input);
    floatvalue = atof(input.c_str());    
    //K = -1000; //this is the gain for the magnetorquers - THIS WORKS NOW - Ok this only works
    K = floatvalue;
    //if you are using MAGNETORQUER_CONTROLLER_TYPE = 1 or 2 this will work
    //K is irrelevant for the Bdot controller

    //Get magnetorquer controller type
    //MAGNETORQUER_CONTROLLER_TYPE = 2; //Controller Type
    getline(CS_file,input);
    MAGNETORQUER_CONTROLLER_TYPE = atoi(input.c_str());
    //0 = off
    //1 = Carlos Montalvo wierd RREF controller
    //2 = k*(BVEC cross omega)
    //3 = Bdot controller
    //Need to check if magnetic field flag is even on
    if (Magnetic_Flag == 0) {
      MAGNETORQUER_CONTROLLER_TYPE = 0;
    }
    
  }
  else {
    printf("File: Input_Files/Satellite.CS Not found \n");
    exit(1);
  }

  //Alright let's get the .MAG file opened
  ifstream magfile;
  magfile.open("Input_Files/Satellite.MAG");
  if (magfile.is_open()) {
    //Grab everything
    //Moment Limits
    //Moment_Limit = 0.02; //N-meters general -- works in 1 orbit
    //Moment_Limit = 1e-6; //saturation for magnetorquers found by Dr.Spencer - Not used at the moment 3/26/2016

    ///Given to me by Matthew Monkel on 3/27/2017
    getline(magfile,input);
    floatvalue = atof(input.c_str());
    magnnumturn = floatvalue;

    //Assuming a 10 cm x 20 cm magnetoruer wrapped around each solar panel
    //areamagnetorque = (10.0/100.0)*(20.0/100.0); //area of each magnetorque, m^2
    getline(magfile,input);
    floatvalue = atof(input.c_str());
    areamagnetorque = floatvalue;

    //maxcurrent = 0.01; //maximum current allowed to flow through each magnetorquer, A
    getline(magfile,input);
    floatvalue = atof(input.c_str());    
    maxcurrent = floatvalue;
    //this is with the safety factor included.
  }
  else {
    printf("File: Input_Files/Satellite.MAG not found \n");
    exit(1);
  }

  
}

void Controller::setState(MATLAB XYZin,MATLAB PTPin) 
{
  state_XYZ_.overwrite(XYZin);
  state_PTP_.overwrite(PTPin);
}

void Controller::setMagVec(MATLAB BVECin) 
{
  //Put Tesla_OLD here?
  BVEC_Tesla_OLD.overwrite(BVEC_Tesla);
  
  //Takes the magnetic field in the body frame and gives it to the
  //Controller class
  BVEC_.overwrite(BVECin);

  //We need to update BVEC_Tesla as well here I think
  //Convert Magnetic field to Teslas
  BVEC_Tesla.mult(BVEC_,1e-9);

  //How do we update BVEC_Tesla_OLD?
}


void Controller::Compute(double simtime) {
  //If it's time to compute control go ahead and do it
  if (simtime >= timeNext_) {
    //Reset timer
    timeNext_ = simtime + UpdateRate_;
    //q0123 = state_PTP_tilde(1:4);
    // q0123.vecset(1, 4, state_PTP_, 1);
    //pqr = state_PTP_tilde(5:7); %%rad/s
    // pqr.vecset(1, 3, state_PTP_, 5);
    //ptp = quat2euler(q0123')';
    // ptp.quat2euler(q0123);
    //pqrc = Kp*-ptp + Kd*-pqr;
    //ptp = -ptp
    // ptp.mult_eq(-1.0);
    // Kpptp.mult(Kp, ptp);
    //pqr = -pqr;
    // pqr.mult_eq(-1.0);
    // Kdpqr.mult(Kd, pqr);
    // pqrc.plus(Kpptp, Kdpqr);
    //LMNcontrol = Km*(pqrc+pqr);
    // pqrc_pqr.plus(pqrc, pqr);
    
    //LMN.mult(Km, pqrc_pqr);
    //Check for saturation
    LMN.mult_eq(0);

    //Magnetorquers! How do they work?
    //LMN_Magnet_Torquers = mu x magnetic field
    //Both mu and Magnetic field need to be in the body frame
    //Thankfully the magnetic field is already in the controller class
    //it is called BVEC_.

    //No need to call the magnetorquers if the magnetic field is on
    if (MAGNETORQUER_CONTROLLER_TYPE > 0) {
      //mu is the magnetic moment of the satellite the equation for magnetic moment is
      // number of turns * current * Area * normal vector
      // The magnetic moments are computed in this routine below
      computeMagneticMoment();

      //Once that routine runs you should be able to call
      //MMTVEC.disp() and see the magnetic moment
      //MMTVEC.disp();

      //Once you have the magnetic moment you can now cross it with the magnetic field vector
      //to do that we call
      computeMagneticTorque();
    }
    
    //BVEC.disp();
    //checkCurrent();
    //getMagnetorque();
    //for (int idx = 1; idx <= 3; idx++) {
    //if (fabs(LMN.get(idx, 1)) > Moment_Limit) {
    //LMN.set(idx, 1, copysign(Moment_Limit, LMN.get(idx, 1)));
    //}
    //}
  }

  //LMN.disp();
  
}

//Again this routine needs to be renamed properly
void Controller::computeMagneticTorque() 
{
  //These lines aren't needed as they have been moved to computeMagneticMoment()
  //MMTVEC.set(1, 1, (magnnumturn*xcurrent*areamagnetorque));
  //MMTVEC.set(2, 1, (magnnumturn*ycurrent*areamagnetorque));
  //MMTVEC.set(3, 1, (magnnumturn*zcurrent*areamagnetorque));

  //This is still needed but BVECCON has been changed to BVEC_Tesla to denote the fact
  //that this is simply the magnetic field vector in the body frame in Teslas rather
  //than nanoTesla. Also the routine has been replaced with the mult function in MATLAB.cpp
  //command
  //BVECCON.set(1, 1, (BVEC.get(1, 1))*1e-9); //this converts the values in BVEC from nT to T
  //BVECCON.set(2, 1, (BVEC.get(2, 1))*1e-9);
  //BVECCON.set(3, 1, (BVEC.get(3, 1))*1e-9);
  //BVEC_Tesla.mult(BVEC_,1e-9); - Moved to computeMagneticCurrent Routine


  //Now we run the cross product
  //BVEC_Tesla.disp();
  MTVEC.cross(MMTVEC, BVEC_Tesla);

  //MTVEC.disp();
  
  LMN.plus_eq(MTVEC);

  //LMN.set(2,1,-0.01*state_PTP_.get(6,1));

  //LMN.disp();
  
}

void Controller::computeMagneticMoment() //Check current is misleading. Better function names are required here
{
  //Spencer. What is this?
  //P.set(1, 1, state_PTP_.get(5, 1));
  //Q.set(2, 1, state_PTP_.get(6, 1));
  //R.set(3, 1, state_PTP_.get(7, 1));
  //XMMTVEC.cross(BVECCON, P);
  //XMMTVEC.mult_eq(-1 * gain);
  //YMMTVEC.cross(BVECCON, Q);
  //YMMTVEC.mult_eq(-1 * gain);
  //ZMMTVEC.cross(BVECCON, R);
  //ZMMTVEC.mult_eq(-1 * gain);

  //Ok not sure what's going on above but to compute magnetic moment we need to do
  //the following

  // number of turns * current * Area * normal vector

  // First things first - number of turns?
  // This is set in the Controller::Controller routine. The name of the variable
  // is magnnumturn

  // Ok now the area? This again is set in Controller::Controller and is called
  // areamagnetorque

  // Then we have the three normal vectors (XMTVEC, YMTVEC, ZMTVEC)

  // Finally we need to compute the current.
  //To do that we run the function computeMagneticCurrent()
  //This routine will determine what currents to apply in order to
  //slow down the satellite or spin it up. Right now only de-tumbling
  //is coded in.
  computeMagneticCurrent();

  //Now we can compute the magnetic moments (XMMTVEC,YMMTVEC,ZMMTVEC)
  //printdouble(magnnumturn,"NTurns");
  //printdouble(areamagnetorque,"Area");
  //printdouble(xcurrent,"xccurent");
  XMMTVEC.mult(XMTVEC,magnnumturn*areamagnetorque*xcurrent);
  //XMMTVEC.disp();
  YMMTVEC.mult(YMTVEC,magnnumturn*areamagnetorque*ycurrent);
  //Even though this one is zero we leave it here for consistency
  ZMMTVEC.mult(ZMTVEC,magnnumturn*areamagnetorque*zcurrent);

  //Ok and that's it. So I guess we can just add them all together now
  MMTVEC.mult_eq(0);
  MMTVEC.plus_eq(XMMTVEC);
  MMTVEC.plus_eq(YMMTVEC);
  MMTVEC.plus_eq(ZMMTVEC);

  //Not sure what this is down here
  // xcurrent = ((MMTVEC.get(1, 1)) / (magnnumturn * areamagnetorque));
  // ycurrent = ((MMTVEC.get(2, 1)) / (magnnumturn * areamagnetorque));
  // zcurrent = 0; //this compensates for the fact that there is no magnetorquer on the z axis
  // //zcurrent = ((MMTVEC.get(3, 1)) / (magnnumturn * areamagnetorque));
  // MMTVEC.mult_eq(0);
  // if (abs(xcurrent) > maxcurrent)
  //   {
  //     xcurrent = copysign(maxcurrent, xcurrent);
  //   }
  // if (abs(ycurrent) > maxcurrent)
  //   {
  //     ycurrent = copysign(maxcurrent, ycurrent);
  //   }
  // if (abs(zcurrent) > maxcurrent)
  //   {
  //     zcurrent = copysign(maxcurrent, zcurrent);
  //   }
}

void Controller::computeMagneticCurrent() {
  //In order to check and make sure the magnetorquers will actually
  //slow us down we need to compute the sign of the moment from current. Unless
  //Hmmmm. Can current be negative at which point we could run the current backwards
  //and get the opposite moment? Yea. Let's get the info we need and then
  //email Dr. Spencer.

  // //Here we can just cross XMTVEC and BVEC since units don't matter
  // XMMTVEC.cross(XMTVEC,BVEC_);
  // //Now we need to determine the largest component and it's sign
  // double max_val = XMMTVEC.max_abs();

  // //Finally we apply maximum current such that maxcurrent*max_val is negative
  // xcurrent = -copysign(maxcurrent,max_val); //This set xcurrent to maxcurrent with the
  // //opposite sign of max_val. If max_val is negative we make xcurrent a positive
  // //ok wait I think this is way more complex that I thought

  // // XMMTVEC.disp();
  // // printdouble(max_absolute,"Max Absolute");
  
  // YMMTVEC.cross(XMTVEC,BVEC_);
  // max_val = YMMTVEC.max_abs();
  // ycurrent = -copysign(maxcurrent,max_val);

  // //AWW SHIT. I somehow assumed that the moments should always me negative but that's
  // //not always true. I want the moment to be opposite to the current angular rate.
  // //Dammit. Gonna need to look at this on paper. Gonna wait till CubeSat meeting on Thurs.
  // //Hopefully this will be fixed before 3/9/2017

  // //Not needed but gonna leave it for consistency
  // ZMMTVEC.cross(ZMTVEC,BVEC_);
  // zcurrent = 0; //This is set to zero because we do not have a zaxis magnetorquer

  //Convert Magnetic field to Teslas
  //BVEC_Tesla.mult(BVEC_,1e-9);-- THis was moved to setMagVec

  //Ok so here is where we compute the three types of controllers
  if (MAGNETORQUER_CONTROLLER_TYPE == 0) {
    xcurrent = 0;
    ycurrent = 0;
    zcurrent = 0;
  }
  else if (MAGNETORQUER_CONTROLLER_TYPE == 1) {
      //Carlos Montalvo Derivation
    float Mx = K*state_PTP_.get(5,1);
    float Mz = K*state_PTP_.get(7,1);
    float Mxbar = Mx/(magnnumturn*areamagnetorque);
    float Mzbar = Mz/(magnnumturn*areamagnetorque);
    //If zcurrent is set to zero Mybar isn't even used
    float Mybar = -Mxbar*BVEC_Tesla.get(1,1)/BVEC_Tesla.get(2,1) - Mzbar*BVEC_Tesla.get(3,1)/BVEC_Tesla.get(2,1);
    zcurrent = 0; //Make zcurrent zero for simplicity
    //Solve for iy and ix
    
    //It's possible that BVEC_Tesla.get(3,1) might be really small. Need a fix for that
    if (abs(BVEC_Tesla.get(3,1)) < 1e-25) {
      ycurrent = 0;
    }
    else {
      ycurrent = (Mxbar + BVEC_Tesla.get(2,1)*zcurrent)/BVEC_Tesla.get(3,1);
    }
    if (abs(BVEC_Tesla.get(2,1)) < 1e-25) {
      xcurrent = 0;
    }
    else {
      xcurrent = (Mzbar + BVEC_Tesla.get(1,1)*ycurrent)/BVEC_Tesla.get(2,1);
    }
  }
  else if (MAGNETORQUER_CONTROLLER_TYPE == 2) {
    //mu_ideal = k*(BVEC cross omega)
    pqr.vecset(1, 3, state_PTP_, 5);
    //pqr.disp();
    MMTVEC.cross(pqr,BVEC_Tesla);
    //MMTVEC.disp();
    MMTVEC.mult_eq(K);
    //MMTVEC.disp();
    //No from here we need to compute x,y,z current based on this
    xcurrent = MMTVEC.get(1,1)/(magnnumturn*areamagnetorque);
    ycurrent = MMTVEC.get(2,1)/(magnnumturn*areamagnetorque);
    zcurrent = MMTVEC.get(3,1)/(magnnumturn*areamagnetorque);
  }
  else if (MAGNETORQUER_CONTROLLER_TYPE == 3) {
    // printf("Running \n");
    
    //This is the Bdot Controller
    MMTVEC.minus(BVEC_Tesla,BVEC_Tesla_OLD);
    MMTVEC.mult_eq(UpdateRate_);

    //Reset BVEC_Tesla_OLD
    //BVEC_Tesla_OLD.overwrite(BVEC_Tesla); //Do I really need this??
    //Ok. We need to do this only when we update the magnetic model

    // BVEC_.disp();
    // BVEC_Tesla.disp();
    // BVEC_Tesla_OLD.disp();
    // printdouble(UpdateRate_,"UpdateRate");
    // MMTVEC.disp();
    
    //Now from here we need to compute x,y,z current based on this
    xcurrent = K*MMTVEC.get(1,1)/(magnnumturn*areamagnetorque);
    ycurrent = K*MMTVEC.get(2,1)/(magnnumturn*areamagnetorque);
    zcurrent = K*MMTVEC.get(3,1)/(magnnumturn*areamagnetorque);
  }

  //Saturation Blocks - max current is set in the Controller::Controller routine
  double norm_current = fabs(xcurrent) + fabs(ycurrent) + fabs(zcurrent);
  if (norm_current > maxcurrent) {
    xcurrent = xcurrent * maxcurrent / norm_current;
    ycurrent = ycurrent * maxcurrent / norm_current;
    zcurrent = zcurrent * maxcurrent / norm_current;
  }

  //printf("Current = %lf %lf %lf \n",xcurrent,ycurrent,zcurrent);
  
  CURRENT_VEC.set(1,1,xcurrent);
  CURRENT_VEC.set(2,1,ycurrent);
  CURRENT_VEC.set(3,1,zcurrent);

  //Debug Prints
  // printdouble(xcurrent,"xcurrent");
  // CURRENT_VEC.disp();
  
}

double Controller::getCurrent(int row) {
  return CURRENT_VEC.get(row,1);
}

void Controller::addControl(MATLAB LMNin)
{
  LMN.plus_eq(LMNin);
}

void Controller::getControl(MATLAB LMNout) 
{
  LMNout.overwrite(LMN);
}

Engine::Engine() {
  //Alright all of this needs to go into an input file
  //Let's call this one Satellite.SIM
  //TIME VECTOR (units are seconds)%%
  simtime = 0;
  itime = 0; 
  tprint = 0; //print to file flag

  simulate = 1; //so far everything is ok

  fstream simfile;
  string input;
  double floatvalue;
  int intvalue;
  simfile.open("Input_Files/Satellite.SIM");

  if (simfile.is_open()) {
    //Grab everything
    //SIMXYZ
    getline(simfile,input);
    SIMXYZ = atoi(input.c_str()); 
    //SIMPTP
    getline(simfile,input);
    SIMPTP = atoi(input.c_str());
    //Magnet and Gravity Model Stuff as well as Magnetorquer stuff
    getline(simfile,input);
    Gravity_Flag = atoi(input.c_str());

    getline(simfile,input);
    yr = atoi(input.c_str());
    
    getline(simfile,input);
    Magnetic_Flag = atoi(input.c_str());
    
    //time_gravity_next = 1; //These have been removed because we have
    //to call the gravity model every timestep now
    getline(simfile,input);
    time_magnet_next = atof(input.c_str()); //This says only call the magnetic model every 1 second
    //time_gravity = 0; //same with this one

    //ALL OF THIS BELOW IS ARCHAIC CODE. I DON"T WANT TO DELETE IT JUST BECAUSE
    //I MIGHT WANT TO USE THESE PARAMTERS
    getline(simfile,input);
    timestep = atof(input.c_str());

    getline(simfile,input);
    tfinal = atof(input.c_str());

    getline(simfile,input);
    tprintnext = atof(input.c_str());

    //////////////////////////////////

    //if (SIMPTP)
    //  {
    ////FOR ROTATIONAL DYNAMICS 
    //Limit for my laptop is 4e-6 = 4000 kHz. I have an i7 quad core processor @ 2.9 GHz 
    //Only running attitude dynamics
      
    // 4e-5 = 400 kHz
    // 10e-5 = 1000 kHz = 1 Mhz

    //timestep = 0.000004; 
    //timestep = 1e-2;
    //tfinal = 200; //just for debug purposes
    //tprintnext = 0.1;
    //tfinal = 40 * 60; //minutes of orbit
    //tfinal = 1.5*60*60; //approx 1 orbit
    //tfinal = 60 * 60 * 24 * 0.5; //4 days
    //tfinal = 60;
    //tprintnext = 10; //rate to print to file
    // }
    //else
    // {
    ////FOR ORBIT DYNAMICS
    //   timestep = 0.01*3600;
    //timestep = 0.001*3600;
    //   tfinal = 3.0*3600;
    //   tprintnext = 10; //rate to print to file
    // }

    //Number of satellits
    getline(simfile,input);
    NSATS = atoi(input.c_str());
    
  }
  else {
    printf("File: Input_Files/Satellite.SIM not found \n");
    exit(1);
  }

  //If opengl is on you don't have a choice you need to simulate both
  ///DO NOT CHANGE THESE TWO LINES OF CODE!!!!! IF YOU WANT TO CHANGE SIMXYZ OR SIMPTP
  //SCROLL UP A BIT AND CHANGE IT THERE!!!!
  #ifdef OPENGL_FLAG
  SIMXYZ = 1;
  SIMPTP = 1;
  #endif

}

//Constructor
Satellite::Satellite(int ierror,MATLAB state_XYZ0,MATLAB state_PTP0,int satnum) 
{
  Initialize(ierror,state_XYZ0,state_PTP0,satnum);
}

//Satellite Threads
void SatellitePhysicsThread(Satellite** p_USASat,Controller** p_Control,Engine* p_Engine) 
{
  //Connect all control classes to the satellite classes
  for (int idx=0;idx<p_Engine->NSATS;idx++) {
    p_USASat[idx]->setControl(p_Control[idx]);
  }
  //Now comes the hard part because we can't just run physics for each one since we need to run all of them
  //Simultaneously. So this needs to be an Engine routine
  p_Engine->run_satellite_physics(p_USASat); 
}

//Get State (REVISIT THIS NEEDS A MUTEX)
void Satellite::getState(MATLAB state_PTP,MATLAB state_XYZ) 
{
  state_PTP.overwrite(state_PTP_);
  state_XYZ.overwrite(state_XYZ_);
}

void Satellite::setControl(Controller* p_Control) 
{
  control = p_Control;
}

double Engine::getTime() 
{
  return simtime;
}

int Engine::simOK() 
{
  return simulate;
}

void Satellite::Initialize(int ierror,MATLAB state_XYZ0,MATLAB state_PTP0,int satnum) 
{
  if (ierror) {
    //Model Satellite
    name_ = "Model Satellite";
  }
  else {
    name_ = "Satellite";
  }
  imodel = ierror;
  time_magnet = 0;
  satnum_ = satnum;

  //When you invoke this command the computer is going to search for
  //these models. If it can't find them it will throw an error like this
  //Cannot open /usr/local/share/GeographicLib/gravity/egm2008.egm
  //Note it's possible to send the code a different path so we will fix this by putting another
  //#ifndef. If you're on windows you need to download the model otherwise
  //we will simply tell the computer to look in our user specified directory
  //Notice that in satellite_physics.h we created a pointer to grav and mag
  //In order to initialize them we need to use the "new" function much like
  //we use the new function in main.cpp when we create a new Satellite
  #ifndef __linux__
  grav = new GravityModel("egm2008"); //Initializing gravity model
  mag = new MagneticModel("emm2015"); //Initializing magnetic model
  #else
  grav = new GravityModel("egm2008","Model_Resources/EGM_EMM"); 
  mag = new MagneticModel("emm2015","Model_Resources/EGM_EMM"); 
  #endif
  printf("Gravity and Magnetic Models Imported \n");
  //Note. Since these are variables, it would be much nicer to put these in the
  //Satellite class rather than these massive global variables. We can let it slide
  //for now since it works but just a thought.

  //Mass and Inertia of Satellite
  //Sam thing. Let's move this to Satellite.MASS
  //mass = 2.6*(1.0+ierror*randnum(0,0.1)); //kg
  //a_x = 10.0*(1.0+ierror*randnum(0,0.1))/100.0; //meters
  //b_y = 10.0*(1.0+ierror*randnum(0,0.1))/100.0; //meters
  //c_z = 20.0*(1.0+ierror*randnum(0,0.1))/100.0; //meters
  //For now let's just compute inertia based on mass and geometry
  //Ixx = (1.0/12.0)*mass*(SQUARE(b_y)+SQUARE(c_z));
  //Iyy = (1.0/12.0)*mass*(SQUARE(a_x)+SQUARE(c_z));
  //Izz = (1.0/12.0)*mass*(SQUARE(a_x)+SQUARE(b_y));

  //Here we are assuming that every satellite has the same mass
  fstream massfile;
  string input;
  massfile.open("Input_Files/Satellite.MASS");
  if (massfile.is_open()) {
    printf("Mass File Found: Input_Files/Satellite.MASS \n");
    //Mass
    getline(massfile,input);
    mass = atof(input.c_str());
    //Inertias
    getline(massfile,input);
    Ixx = atof(input.c_str());
    getline(massfile,input);
    Iyy = atof(input.c_str());
    getline(massfile,input);
    Izz = atof(input.c_str());
    massfile.close();
  }
  else {
    printf("File: Input_Files/Satellite.MASS not found \n");
    exit(1);
  }

  printf("Initializing Matrices and Vectors for Satellite %d \n",satnum);
  Ivec.zeros(3,1,"Ivec");
  Ivec.set(1,1,Ixx);
  Ivec.set(2,1,Iyy);
  Ivec.set(3,1,Izz);
  I.diag(Ivec,"I");  
  Iinv.overwrite(I,"Iinv");
  Iinv.inverse();
  //Copy Initial Condition Vectors
  state_XYZ_.copy_init(state_XYZ0,"state_XYZ");
  state_PTP_.copy_init(state_PTP0,"state_PTP");
  //Initialize Vectors
  q0123.zeros(4,1,"q0123");
  ptp.zeros(3,1,"q0123");
  pqr.zeros(3,1,"pqr");
  LMNdrag.zeros(3,1,"LMNdrag");
  LMNcontrol.zeros(3,1,"LMNcontrol");
  LMN.zeros(3,1,"LMN");
  I_pqr.zeros(3,1,"I_pqr");
  pqrskew_I_pqr.zeros(3,1,"pqrskew_I_pqr");
  pqrdot.zeros(3,1,"pqrdot");
  rk4_PTP.zeros(NUMPTP,1,"rk4_PTP");
  ktimestep.zeros(NUMPTP,1,"ktimestep");
  k1.zeros(NUMPTP,1,"k1");
  k2.zeros(NUMPTP,1,"k2");
  k3.zeros(NUMPTP,1,"k3");
  k4.zeros(NUMPTP,1,"k4");
  rk4_XYZ.zeros(NUMXYZ,1,"rk4_PTP");
  kxyztimestep.zeros(NUMXYZ,1,"kxyztimestep");
  kxyz1.zeros(NUMXYZ,1,"kxyz1");
  kxyz2.zeros(NUMXYZ,1,"kxyz2");
  kxyz3.zeros(NUMXYZ,1,"kxyz3");
  kxyz4.zeros(NUMXYZ,1,"kxyz4");
  xyz.zeros(3,1,"xyz");
  xyzdot.zeros(3,1,"xyzdot");
  GVEC.zeros(NUMGRAV, 1, "Gravity vector");
  BVEC.zeros(NUMMAGN, 1, "Magnetic vector");
  BVECSPH.zeros(NUMMAGN, 1, "Magnetic vector (SPH)");
  BVECINE.zeros(NUMMAGN, 1, "Magnetic vector (INE)");
  //BVEC_Tesla.zeros(NUMMAGN, 1, "Magnetic vector in Tesla");
  MTVEC.zeros(3, 1, "Magnetorque Vector");
  MMTVEC.zeros(3, 1, "Magnetic Moment Vector");
  sph_coord.zeros(3,1,"Spherical Coordinate (Phi and Theta)");
  ine_coord.zeros(4, 1, "Inertial Coordinate (Quaternions)");
  printf("Matrices and Vectors for Satellite %d Initialized \n",satnum);
  
  //Open File for Writing
  #ifdef PRINT_TO_FILE_FLAG
  char filename[256];
  if (ierror) {
    //Model Satellite
    sprintf(filename,"%s%d%s","Output_Files/State",satnum,"Model.OUT");
  }
  
  else {
    sprintf(filename,"%s%d%s","Output_Files/State",satnum,".OUT");
  }
  //Before we open the file for writing I think it makes sense to try and
  //remove all files that way we gaurantee that we don't plot files that
  //don't exist. We can do this in our MATLAB wrapper but I'd rather do it here
  //int remove ( const char * filename )
  if (!remove(filename)){
    printf("Old .OUT file found. Removing before we open for writing \n");
  }
  stateoutfile = fopen(filename,"wb");
  if (!stateoutfile)
    {
      printf("Output File Not Created in Output_Files/ for %s \n",name_);
      printf("File %s not created \n",filename);
      //simulate = 0;
      exit(1);
    }
  else {
    printf("Output File Created in: %s \n",filename);
  }
  #endif
  #ifdef OPENGL_FLAG
  //Initialize VISOBJECT
  //The first variable is time which is always zero or should be
  statetime.Initialize(0,2);//1 for the satellite and the next for earth
  //Initialize EARTH
  XYZ_EARTH.zeros(3,1,"XYZ_EARTH");
  PTP_EARTH.zeros(3,1,"PTP_EARTH");
  #endif
}

void Satellite::print_to_file(double simtime)
{
  fprintf(stateoutfile,"%lf ",simtime);
  for (int idx = 1;idx<=NUMXYZ;idx++)
    {
      fprintf(stateoutfile,"%.8e ",state_XYZ_.get(idx,1));
    }
  for (int idx = 1;idx<=NUMPTP;idx++)
    {
      fprintf(stateoutfile,"%.8e ",state_PTP_.get(idx,1));
    }
  //if (abs(remainder(simtime, timei)) < .1)
  //{
  //I don't think we want this if statement. 
  for (int idx = 1; idx <= NUMGRAV; idx++)
    {
      fprintf(stateoutfile, "%.8e ", GVEC.get(idx, 1));
    }
  for (int idx = 1; idx <= NUMMAGN; idx++)
    {
      fprintf(stateoutfile, "%.8e ", BVECSPH.get(idx, 1)); //these are spherical coordinates
    }
  for (int idx = 1; idx <= NUMMAGN; idx++)
    {
      fprintf(stateoutfile, "%.8e ", BVECINE.get(idx, 1)); //these are inertial frame components
    }
  for (int idx = 1; idx <= NUMMAGN; idx++)
    {
      fprintf(stateoutfile, "%.8e ", BVEC.get(idx, 1)); //These are body frame components
    }
  for (int idx = 1;idx <= 3; idx++)
    {
      fprintf(stateoutfile,"%.8e ",control->getCurrent(idx)); //This gets the magnetorquer current values
    }
  //}
  fprintf(stateoutfile,"\n");
}

void Satellite::computeXYZDerivs(MATLAB deriv,MATLAB rk4_XYZ)
{
  //Extract State
  xyz.vecset(1,3,rk4_XYZ,1);
  xyzdot.vecset(1,3,rk4_XYZ,4);
  double rSat = xyz.norm();
  //Kinematics
  deriv.vecset(1,3,xyzdot,1);
  //Dynamics
  //Wow. Looks like I wasn't even using the forces and moments in the first place
  //double muSat = -GSPACE*(mass+MEARTH);
  //xyz.mult_eq(muSat/pow(rSat,3));
  //This should set xyzdot to gx,gy,gz. Is that right? Yea. I think so.
  //Problem is since this derivatives routine is running inside the RK4 loop we need to call
  //the gravity model every single timestep. Holy shit....
  getCurrentGravity(rk4_XYZ);
  deriv.vecset(4,6,GVEC,1);
}

void Satellite::computePTPDerivs(MATLAB deriv,MATLAB rk4_PTP)
{
  //Extract State
  q0123.vecset(1,4,rk4_PTP,1);
  pqr.vecset(1,3,rk4_PTP,5);
  double q0 = q0123.get(1,1);
  double q1 = q0123.get(2,1);
  double q2 = q0123.get(3,1);
  double q3 = q0123.get(4,1);
  double p = pqr.get(1,1);
  double q = pqr.get(2,1);
  double r = pqr.get(3,1);

  ///KINEMATICS
  deriv.set(1,1,(-p*q1-q*q2-r*q3)/2.0);
  deriv.set(2,1,(p*q0+r*q2-q*q3)/2.0);
  deriv.set(3,1,(q*q0-r*q1+p*q3)/2.0);
  deriv.set(4,1,(r*q0+q*q1-p*q2)/2.0);

  ///////////DYNAMICS/////////////

  //LMNdrag
  LMNdrag.mult_eq(0); //Set to zero for now

  //LMN = LMNdrag + LMNcontrol;
  control->getControl(LMNcontrol);
  //LMNcontrol.disp();
  LMN.plus(LMNdrag,LMNcontrol);

  //LMN.disp();
  //PAUSE();

  //pqrskew = [0 -r q;r 0 -p;-q p 0];  
  //pqrdot = Iinv*(LMN-pqrskew*I*pqr);
  //pqr.disp();
  I_pqr.mult(I,pqr);
  //I.disp();
  //I_pqr.disp();
  pqrskew_I_pqr.cross(pqr,I_pqr);
  //pqrskew_I_pqr.disp();
  LMN.minus_eq(pqrskew_I_pqr);
  //LMN.disp();
  //Iinv.disp();
  pqrdot.mult(Iinv,LMN);
  
  //Save state derivatives
  deriv.vecset(5,NUMPTP,pqrdot,1);

  //pqrdot.disp();
  //PAUSE();

  ///////////////////////////////
}

void Satellite::normalize_quats(int SIMPTP) {
  //Even if we aren't integrating the rotational dynamics we still need to set q0123
  if (!SIMPTP) {
    q0123.vecset(1,4,state_PTP_,1);
  } else {
    q0123.vecset(1,4,state_PTP_,1);
    q0123.mult_eq(1.0/q0123.norm());
    state_PTP_.vecset(1,4,q0123,1);
  }
}

void Satellite::runge_kutta_PTP(double timestep) {
  rk4_PTP.overwrite(state_PTP_);
  computePTPDerivs(k1,rk4_PTP);

  ktimestep.mult(k1,timestep/2.0);
  rk4_PTP.plus(state_PTP_,ktimestep);
  computePTPDerivs(k2,rk4_PTP);

  ktimestep.mult(k2,timestep/2.0);
  rk4_PTP.plus(state_PTP_,ktimestep);
  computePTPDerivs(k3,rk4_PTP);

  ktimestep.mult(k3,timestep);
  rk4_PTP.plus(state_PTP_,ktimestep);
  computePTPDerivs(k4,rk4_PTP);

  //Add ktimestep = timestep*(1/6)*(k1 + 2*k2 + 2*k3 + k4)
  k2.mult_eq(2.0);
  k3.mult_eq(2.0);
  ktimestep.plus(k1,k2);
  ktimestep.plus_eq(k3);
  ktimestep.plus_eq(k4);
  ktimestep.mult_eq(timestep/6.0);
  //state_PTP_ = state_PTP_ + ktimestep);
  state_PTP_.plus_eq(ktimestep);
}

void Satellite::runge_kutta_XYZ(double timestep) {
  rk4_XYZ.overwrite(state_XYZ_);
  computeXYZDerivs(kxyz1,rk4_XYZ);

  kxyztimestep.mult(kxyz1,timestep/2.0);
  rk4_XYZ.plus(state_XYZ_,kxyztimestep);
  computeXYZDerivs(kxyz2,rk4_XYZ);

  kxyztimestep.mult(kxyz2,timestep/2.0);
  rk4_XYZ.plus(state_XYZ_,kxyztimestep);
  computeXYZDerivs(kxyz3,rk4_XYZ);

  kxyztimestep.mult(kxyz3,timestep);
  rk4_XYZ.plus(state_XYZ_,kxyztimestep);
  computeXYZDerivs(kxyz4,rk4_XYZ);

  //Add ktimestep = timestep*(1/6)*(k1 + 2*k2 + 2*k3 + k4)
  kxyz2.mult_eq(2.0);
  kxyz3.mult_eq(2.0);
  kxyztimestep.plus(kxyz1,kxyz2);
  kxyztimestep.plus_eq(kxyz3);
  kxyztimestep.plus_eq(kxyz4);
  kxyztimestep.mult_eq(timestep/6.0);
  //state_XYZ_ = state_XYZ_ + kxyztimestep);
  state_XYZ_.plus_eq(kxyztimestep);
}

#ifdef OPENGL_FLAG
void Satellite::UpdateOPENGL(double simtime) {
  ///When we aren't simulating the model we want to pass information to
  //the following variable. The visualizer routine has a copy of this
  //so if we change this variable the visualizer will update as well.
  if (!imodel) {
    q0123.vecset(1,4,state_PTP_,1);
    ptp.quat2euler(q0123);
    statetime.setState(state_XYZ_,ptp,1); //1 is for the first object
    //the second object is the earth
    PTP_EARTH.set(3,1,simtime*ANGULAREARTH);  
    statetime.setState(XYZ_EARTH,PTP_EARTH,2);
    //printf("Setting Time of statetime to %lf \n",simtime);
    statetime.setTime(simtime);
  }
}
#endif

void Satellite::DebugPrint() {
  printdouble(time_magnet,"Time Magnet = ");
  state_XYZ_.disp();
  state_PTP_.disp();
  GVEC.disp();
  BVECSPH.disp();
  BVECINE.disp();
  BVEC.disp();
}

void Engine::run_satellite_physics(Satellite** p_USASat) {
  
  TIMER track_time; //Timer to track time this TIMER class has a constructor
  //function that initializes time

  //printf("%s Simulation Initialized and Running...\n",name_);
  printf("Simulation Initialized and Running...\n");
  #ifdef ADCS_FLAG || OPENGL_FLAG || SENSOR_FLAG
  printf("Simulation will Simulate Forever \n");
  #endif

  while (simtime <= tfinal+1e-10) {

    //Ok I think this is where we need to loop on our Satellites
    for (int idx = 0;idx<NSATS;idx++) {

      //Kalman Filter Loop
      if (p_USASat[idx]->imodel) { //this means we are simulating the model
	
	//So in here we need to do a few things
	//First we need to make sensor measurements if needed
        #ifdef SENSOR_FLAG
	//GPS(); //- GPS - 10 Hz accuracy 2 meters
	//SS(); //- Solar Sensor

	//9DOF Sensor 100kHz (4e-6 sec, by the way this is the limit that my computer can simulate) 
	// 100 kHz = 100,000 Hz = 10e4 Hz = 1e5 Hz = 1e-5 sec
	//RG(); //- Rate Gyro - p,q,r - //Noise = 0.2 deg/s for rate gyro
	//MAG(); //- magnetometer - (Nx,Ny,Nz) = 3.2,3.2,4.1 mGauss
	//ACCEL(); //- accelerometer
    
	//In here we need to set iGPS,iSolar,iRG,iMAG and iACCEL 
	//to 1 if a measurement is obtained. I think we should make
	//each of these sensors have their own class
	//I'm going class crazy here. Each method will have the following
	//Method
	//int Update(simtime) - this will update if it's time and return true if it updated, everytime you update to write to a file
	//MATLAB getRaw() - get last measurement of sensor
	//attributes will be
	//double timeNext
	//double UpdateRate
	//MATLAB raw
	//int iupdate - update succeeded or not

	//Ok after we sample all the measure we then run the KalmanUpdate depending on whether or 
	//not we received a new measurement
	
	//KalmanUpdate();
	
        #endif //SENSOR_FLAG

	//If we aren't using sensors it means we are just using the model
	//as the sensor measurement so we pass the current
	//state to the controller - there is no need to run the kalman update when you
	//are simply believing the model because there is no other measurement to average between
	//there also isn't an update rate or anything so we just grab the measurment everytime
	//we can always come back and add in an updaterate of the sensors but you may as well
	//simulate the sensors if that's what you want. 
	//ok so this means we need to send the actual value of the model to the controller
	//but we're going to do that below here so we actually don't do anything 
	
	//At this point we have sampled and potentially run a KalmanUpdate.
	//the result is that the model is currently up to date and is the best 
	//guess we have for control computation. Ok so at this point we send
	//the current state of the model to the controller for computation
	
	//Again we will have another class called CONTROL which will have attributes and stuff to compute control

	//control->setState();
	//control->Compute(simtime); //this will compute the controller values based on an update rate and the current time

      } //End imode loop
      #ifndef KALMAN_FLAG
      //This means I am not the model and the Kalman filter is off
      else  {
	//So i figured it'd be good functionality to run the Controller even if the kalman 
	//filter is not running. So if the Kalman Filter is off or rather the KALMAN FLAG is not defined
	//then when imodel = 0 we send the actual state of the system to controller
	p_USASat[idx]->UpdateMagneticFieldANDController(time_magnet_next,simtime);
      }
      #endif

      
      #ifdef OPENGL_FLAG
      p_USASat[idx]->UpdateOPENGL(simtime);
      #endif 

      #ifdef DEBUG_MODE
      //Print state and pause if DEBUG_MODE IS ON
      printdouble(simtime,"Time");
      pUSASat[idx]->DebugPrint();
      PAUSE();
      #endif

      //Print state to file ever tskip
      //This is not how we want to get gravity. I moved some stuff around
      //Make sure this still works.

      //Get gravity and Magnet values
      //timei right now is set in the Satellite Physics Initialize routine
      //It dictates how often we compute the magnet and gravity models. I
      //changed the code however to have a time_magnet and time_gravity so 
      //that the models are called independently    
      
      #ifdef PRINT_TO_FILE_FLAG
      if (simtime > tprint)
	{
	  p_USASat[idx]->print_to_file(simtime);
	}
      #endif

      //Integrate State Using RK4
      if (SIMPTP) {
	//PTP_RK4
	p_USASat[idx]->runge_kutta_PTP(timestep);

      }
      //Normalize Quaternions
      //Even if we aren't integrating the rotational dynamics we still need to
      //set q0123 and normalize quats
      p_USASat[idx]->normalize_quats(SIMPTP);
      
      if (SIMXYZ) {
	//XYZ_RK4
	p_USASat[idx]->runge_kutta_XYZ(timestep);
      }
    } //End of looping through satellites

    //Step print next
    #ifdef PRINT_TO_FILE_FLAG
    if (simtime > tprint) {
      tprint += tprintnext;
    }
    #endif
    
    //Step Time
    simtime+= (timestep);

    //if Hardware is in the loop just simulate forever
    //or if the visualizer is on
    #ifdef SIMFOREVER_FLAG
    tfinal = simtime+1;
    #endif
    //Check and see how much time has elapsed. If time > timeelapsed you need to wait and 
    //then continue. If time elapsed is greater than time it means your computer is too slow
    //and you need to increase the timestep.
    #ifdef REALTIME_FLAG
      time_since_start = track_time.getTimeSinceStart();
      if (simtime < time_since_start-1e-2) {
	//printf("Computation too slow, increase your timestep \n");
      }
      while (simtime > time_since_start) {
        // printf("Computation too fast - Sleeping for %lf seconds \n",timestep);
      	// printf("Time = %lf \n",time);
      	// printf("Time since start = %lf \n",time_since_start);
      	cross_sleep(timestep);
      	time_since_start = track_time.getTimeSinceStart();
      }
      
    #endif

    #ifdef PRINT_TIME
      if (simtime > itime) {
	//name_ doesn't work now that I am simulating multiple satellites. Not sure what to do about that.
	 //printf("%s Time = %lf out of %lf \n",name_,time,tfinal);
	printf("Simulation Time = %lf out of %lf \n",simtime,tfinal);
         itime++;
      }
    #endif
  }
  //printf("%s Simulation Complete \n",name_);
  //Again name_ doesn't work now. 
  printf("Simulation Complete \n");
  simulate = 0; //make simOK return 0.
  return;

}

void Satellite::UpdateMagneticFieldANDController(double time_magnet_next,double simtime) {
  control->setState(state_XYZ_,state_PTP_);
  if (time_magnet_next == -99 && BVEC.get(1,1)!=0) { 
    time_magnet = simtime + 100;
  }
  if (simtime > time_magnet)
    {
      getCurrentMagnetic();
      time_magnet += time_magnet_next;
    }
  //No matter what we need to do this every timestep
  ine2bod(BVEC,BVECINE); //Rotate inertial frame magnetic field to body frame
  control->setMagVec(BVEC); //Send magnetic field in the body frame to the controller
  // Unfortunately I think we need to call the gravity model every timestep or this won't work
  // if (time > time_gravity)
  // {
  // 	getCurrentGravity(state_XYZ_);
  // 	time_gravity += time_gravity_next;
  // }
  //Compute Control using Magnetorquers and Reaction Wheels
  control->Compute(simtime);
}

void Satellite::getCurrentGravity(MATLAB XYZ_current)
{
	double gx, gy, gz;
	double x = XYZ_current.get(1, 1);
	double y = XYZ_current.get(2, 1);
	double z = XYZ_current.get(3, 1);
	double muSat = -GSPACE*(mass + MEARTH);
	double rSat = state_XYZ_.norm();
	//printf("Gravity_Flag == %d \n",Gravity_Flag);
	if (Gravity_Flag == 1)
	{
		grav->W(x, y, z, gx, gy, gz);
		GVEC.set(1, 1, gx);
		GVEC.set(2, 1, gy);
		GVEC.set(3, 1, gz);
	}
	else if (Gravity_Flag == 0) //This is the point mass model here
	{
		gx = (muSat / pow(rSat, 3))*x;
		gy = (muSat / pow(rSat, 3))*y;
		gz = (muSat / pow(rSat, 3))*z;
		GVEC.set(1, 1, gx);
		GVEC.set(2, 1, gy);
		GVEC.set(3, 1, gz);
	}
	else if (Gravity_Flag == -1) { //Gravity is off
	  GVEC.mult_eq(0);
	}
	//GVEC.disp();
}

void Satellite::getCurrentMagnetic()
{
  if (Magnetic_Flag == 1)
    {
      double b_east, b_north, b_vertical;
      double bx,by,bz;
      double x = state_XYZ_.get(1, 1);
      double y = state_XYZ_.get(2, 1);
      double z = state_XYZ_.get(3, 1);
      //x = rho*sin(phi)*cos(theta);
      //y = rho*sin(phi)*sin(theta);
      // y/x = tan(theta)
      //z = rho*cos(phi);
      double rho = sqrt((pow(x, 2) + pow(y, 2) + pow(z, 2)));
      double phi = (acos(z / rho));
      double the = atan2(y , x);
      double lat = 90 - phi*(180 / PI);
      double lon = the*(180/PI);
      double h = rho-REARTH; //need to send the model the height above the earth's surface
      #ifdef DEBUG_MODE
      state_XYZ_.disp();
      printdouble(rho,"rho");
      printdouble(phi,"phi");
      printdouble(the,"the");
      printdouble(lat,"lat");
      printdouble(lon,"lon");
      printdouble(h,"h");
      PAUSE();
      #endif
      /**
       * Evaluate the components of the geomagnetic field.
       *
       * @param[in] t the time (years).
       * @param[in] lat latitude of the point (degrees).
       * @param[in] lon longitude of the point (degrees).
       * @param[in] h the height of the point above the ellipsoid (meters).
       * @param[out] Bx the easterly component of the magnetic field (nanotesla).
       * @param[out] By the northerly component of the magnetic field (nanotesla).
       * @param[out] Bz the vertical (up) component of the magnetic field
       *   (nanotesla).
       **********************************************************************/
      mag->operator()(yr, lat, lon, h, b_east, b_north, b_vertical);

      //I'm assuming that vertical is Down to preserve the right hand rule.

      //Our spherical reference frame is as follows - b_east is our y
      //component, b_north is our x component and b_vertical is our z
      //component
      bx = b_north;
      by = b_east;
      bz = b_vertical;
      BVECSPH.set(1, 1, bx); //Btdubs these are all in nT
      BVECSPH.set(2, 1, by);
      BVECSPH.set(3, 1, bz);
      
      // if (BVECSPH.found_nans()) {
      //   printf("Sorry Nans detected. Suggest Stop \n");
      //   printf("%lf %lf %lf \n",yr,lat,lon);
      //   printf("%lf %lf %lf \n",x,y,z);
      //   printf("%lf %lf \n",rho,phi);
      //   exit(1);
      // }

      //Our inertial frame is set such that x goes out the equator at
      //the prime meridian, y is orthogonal to x and z goes through
      //the north pole

      //In order to go from spherical to inertial we need to
      //understand that aircraft convention uses the 3-2-1 Euler angle
      //convention

      //Rotation about the z-axis (psiE) - 3
      //rotation about the y-axis (thetaE) - 2
      //rotation about the x-axis (phiE) - 1
      //However, these are not the same as phi and the from the
      //spherical reference frame. Longitude is measure to the right
      //but this is actually a positive rotation about z. Thus
      double psiE = the;
      //Furthermore, phi is a positive rotation about the y-axis but you need to add pi to make sure
      //the z-component points downwards
      double thetaE = phi+PI;
      //Finally, there is no rotation about the x-axis
      double phiE = 0;

      //With these "Euler" Angles defined we can convert the spherical
      //coordinates to inertial coordinates.

      sph2ine(BVECINE, BVECSPH, phiE, thetaE, psiE);
      BVEC.overwrite(BVECINE);
      //We need to rotate the BVECINE TO BODY EVERY TIMESTEP so this has been moved
      //ine2bod(BVEC, BVECINE); //Why are there 3 zeros here? I fixed this but still need to test      
      //Detect for Nans
    }
  else
    {
      //We need an else here. Maybe just set BVEC to zero? Yea
      //that should be fine
      BVEC.mult_eq(0.0);
    }
  //control->BVECSPH.overwrite(BVECSPH);
  //control->BVECINE.overwrite(BVECINE);
  //BVEC.disp();

  //This line of code sends to the satellite class body frame magnetic
  //field vector to the controller class magnetic field vector
  //This also needs to get sent to the controller everytimestep so this has been moved.
  //control->setMagVec(BVEC);
}

void Satellite::sph2ine(MATLAB vecI,MATLAB vecSPH, double phi, double the, double psi)
{
  sph_coord.set(1,1,phi);
  sph_coord.set(2,1,the);
  sph_coord.set(3,1,psi);
  //sph_coord.disp();
  sph2ine32.L321(sph_coord,0); //0 for Euler Angles
  //sph2ine32.disp();
  //vecSPH.disp();
  sph2ine32.rotateBody2Inertial(vecI,vecSPH);
  //vecI.disp();
}

void Satellite::ine2bod(MATLAB vecB, MATLAB vecI)
{
  #ifdef DEBUG_MODE
  q0123.disp();
  PAUSE();
  #endif
  ine2bod123.L321(q0123, 1);
  ine2bod123.rotateInertial2Body(vecB, vecI);
}
