#ifndef VEHICLE_H
#define VEHICLE_H

//Because C++ does not have any native matrix multiplication
//I create my own. They are included here.
#include "MATLAB.h"
#include "Rotation3.h"
#include <fstream>
#include "Motor.h"
#include "Propulsion.h"
#include "params.h"
#include "timer.h"

class Vehicle {
 private:
  float mass;
  double zint=0,theta_int=0,roll_int=0,yaw_int=0,rx,ry,rz,rxpusher,d_pitch,d_yaw,d_roll,thrust;
  double perror,qerror,rerror,d_thrust;
  char input_file_[256];
  MATLAB simdata,massdata,icdata,stateinitial,k1,k2,k3,k4,state,rkstate,phi,F,M,Fgrav,Fbody;
  MATLAB Iinv,I,pqr,uvw,uvwdot,xyzdot,ptp,pqrdot,I_pqr,LMN,ptpdot,Kuvw_pqr,pqrskew_I_pqr,rxyz;
  Propulsion propel;
  void pid_loop(MATLAB,MATLAB,double);
  void feedback_linearization(MATLAB,double t);
 public:
  //Constructor
  Vehicle(char*);
  int ok;
  int ImportFiles();
  int ImportFile(char*,MATLAB*,char*);
  void Integrate();
  void Derivatives(MATLAB,double,MATLAB);
  void ForceandMoment(MATLAB,MATLAB,double);
  Rotation3 ine2bod123;
};

#endif
