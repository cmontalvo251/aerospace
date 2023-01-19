#ifndef MOTOR_H
#define MOTOR_H

#include "MATLAB.h"

class Motor {
private:
public:
  MATLAB r,n,Tn,Iromegan,datapts;
  double mass_prop,mass_motor,kt,kq,Max_Thrust,Nominal_Thrust;
  //Constructor
  Motor();
  int sign;
  double pwm_signal,thrust,omega,torque,ct,cq,rho,A,R,a0,a1,Ir,MINSIGNAL,HOVERSIGNAL,MAXSIGNAL;
  double pwm_nominal,thrust_filtered=0,torque_filtered=0,a;
  void compute_omega();
  void compute_thrust();
  void compute_torque();
  void compute_aero();
  void MotorsSetup();
  double compute_signal(double);
  double compute_signal_NO_OMEGA(double);
  void compute_thrust_NO_OMEGA();
  void compute_torque_NO_OMEGA();
  void MotorCalcs(double ms,double ma,double T,double p,double w,double R);
};

#endif
