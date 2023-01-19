#ifndef PROP_H
#define PROP_H

#include "MATLAB.h"
#include "Rotation3.h"
#include <fstream>
#include "Motor.h"
#include "Propulsion.h"
#include "params.h"
#include "timer.h"
#include <Common/FASTPWM.h>

class Propulsion {
 private:
 public:
  MATLAB Fmotors,TauMotors,M,H,U,HT,HHT,HHT_inv,Q,HT_inv_HHT,CHI;
  MATLAB Hprime,Qprime,HTprime,CHIprime,HT_inv_HHTprime,I_,rxyz_,datapts;
  Motor MOTORS[20]; //Make room for 20 motors just cuz I can
  int NUMMOTORS=0,MOTORSOFF=0; //Motors off initialized to none off
  int MOTORSRUNNING,MINSIGNAL,MAXSIGNAL;
  double Tdatapt_,pwm_datapt_,omegaRPMdatapt_,R_;
  double rx_,ry_,rz_,rxpusher_,mass_;
  void MotorsSetup(double m,MATLAB I,double mins,double max,MATLAB rxyz,double rxp,MATLAB datapts,double R);
  void addMotor(double,double,double,double,double,double,int);
  void computeControl(MATLAB,MATLAB,double);
  void MotorBeep(int,int,MATLAB,double);
  void RemoveMotors(int);
  void computeReconfigurable(double,double,double,double);
  void ForceMomentMotors(MATLAB);
  void print_to_file(FILE*);

  //Constructor
  Propulsion();
};

#endif
