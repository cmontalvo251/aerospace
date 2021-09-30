// Multiple CubeSat Simulation Software 
// Copyright - Carlos Montalvo 2015
// You may freely distribute this file but please keep my name in here
// as the original owner

#ifndef SAT_H
#define SAT_H

//Number of States
#define NUMXYZ 6
#define NUMPTP 7
#define NUMGRAV 3
#define NUMMAGN 3

#ifdef OPENGL_FLAG //This is set in main.h and is only on if you want the visualizer running
#include "opengl.h" //this is a symbolic link to BlackBox.git/C++ 
#endif

#include "../MATLAB.h" //this is a symbolic link to BlackBox.git/C++ (MATLAB.h requires mathp.h)
#include "../Rotation3.h" //this is a symbolic link to BlackBox.git/C++ (This requires MATLAB.h)

//#includes for header files need to be in here. This way we only need
//to add them once

//The GeographicLib library can be downloaded from the internet by following the README.md locating at the ~/ of this folder
//However for the convenience of the linux user the libraries have been placed in the source folder here so ther
//is really no need to download it. In order to make this work a -I and -L have been added to the Makefile
#ifdef __linux__ //this is really annoying but windows uses \ and linux uses /.
#include "GeographicLib/GravityModel.hpp" //Necessary for gravity model
#include "GeographicLib/MagneticModel.hpp" //Necessary for magnetic model
#else
#include "GeographicLib\GravityModel.hpp" //Necessary for gravity model
#include "GeographicLib\MagneticModel.hpp" //Necessary for magnetic model
#endif

#include <fstream> //For file IO

using namespace GeographicLib; //Necessary for Grav and Magnet models

class Controller {
 private:
  int MAGNETORQUER_CONTROLLER_TYPE;
  double UpdateRate_,timeNext_,Moment_Limit;//, gain;
  double K;
  MATLAB LMN,q0123,pqr,ptp,Kpptp,Kdpqr,pqrc,pqrc_pqr;
  MATLAB Kp, Kd, Km, state_XYZ_, state_PTP_;
 
  double xcurrent, ycurrent, zcurrent, areamagnetorque, magnnumturn, maxcurrent;
  MATLAB BVEC_,XMTVEC, YMTVEC, ZMTVEC, MTVEC, XMMTVEC, YMMTVEC, ZMMTVEC, MMTVEC;
  MATLAB BVEC_Tesla,BVEC_Tesla_OLD;
  MATLAB CURRENT_VEC;

  void computeMagneticMoment();
  void computeMagneticTorque();
  void computeMagneticCurrent();
  
 public:
  void getControl(MATLAB);
  void setState(MATLAB,MATLAB);
  void Compute(double);
  void addControl(MATLAB);
  void setMagVec(MATLAB);
  double getCurrent(int);
  //void Initialize();
  //Constructor
  Controller(int);
};

class Satellite {
 private:
  int satnum_;
  double mass,a_x,b_y,c_z,Ixx,Iyy,Izz,timestep;
  double x0,xd0,y0,yd0,Vsat,inclination_angle;
  double z0,zd0,phi0,theta0,psi0,p0,q0,r0;
  double time_since_start,time_magnet;
  MATLAB Ivec,I,Iinv,ptp,q0123,state_PTP_,state_XYZ_;
  MATLAB k1,k2,k3,k4,quatdot,pqrdot,pqr,I_pqr;
  MATLAB pqrskew,LMN,pqrskew_I_pqr,LMNcontrol;
  MATLAB LMNdrag,rk4_PTP,ktimestep,kxyz1,kxyz2,kxyz3,kxyz4;
  MATLAB kxyztimestep,rk4_XYZ,xyz,xyzdot;
  MATLAB XYZ_EARTH,PTP_EARTH;
  MATLAB GVEC, MTVEC,BVEC;
  MATLAB sph_coord, BVECSPH, ine_coord, BVECINE;
  Rotation3 sph2ine32, ine2bod123;
  Controller* control;
  FILE* stateoutfile=NULL;
  char* name_;
  void Initialize(int,MATLAB,MATLAB,int);
  void computePTPDerivs(MATLAB,MATLAB);
  void computeXYZDerivs(MATLAB,MATLAB);
  void getCurrentGravity(MATLAB);
  void getCurrentMagnetic();
  void ine2bod(MATLAB, MATLAB);
  void sph2ine(MATLAB, MATLAB, double, double, double);
  //For Variables like GravityModel and MagneticModel we need to make a pointer
  //otherwise C++ won't like the declaration here
  //Then in the initialize routine we will actually make the variables and initialize them
  GravityModel* grav;
  MagneticModel* mag;
 public:
  #ifdef OPENGL_FLAG
  VISOBJECT statetime;
  void UpdateOPENGL(double);
  #endif
  MATLAB MMTVEC;
  int imodel;
  int Gravity_Flag,Magnetic_Flag;
  double yr;
  void print_to_file(double);
  void setControl(Controller*);
  void getState(MATLAB,MATLAB);
  void UpdateMagneticFieldANDController(double,double);
  void normalize_quats(int);
  void runge_kutta_PTP(double);
  void runge_kutta_XYZ(double);
  void DebugPrint();
  Satellite(int,MATLAB,MATLAB,int);
  //~Satellite();
};

//Physics Engine Class
class Engine {
 private:
  int SIMXYZ,SIMPTP,itime,simulate;
  double simtime,timestep,tfinal;
  double tprint,tprintnext,time_magnet_next;
  double time_magnet,time_gravity;
  double time_gravity_next;
  
 public:
  int NSATS;
  void run_satellite_physics(Satellite**);
  double getTime();
  int Gravity_Flag,Magnetic_Flag;
  double yr;
  int simOK();
  Engine();
};


//Satellite Threads
void SatellitePhysicsThread(Satellite**,Controller**,Engine*); //This is a double pointer because it is a pointer to a vector of satellites. I hope this works

#endif
