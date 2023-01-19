//#include "Vehicle.h"
#include "Motor.h"
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include "mathp.h"

using namespace std;

//Constructor 
Motor::Motor() {
  r.zeros(3,1,"r");
  n.zeros(3,1,"n");
  Tn.zeros(3,1,"Tn");
  datapts.zeros(3,1,"Motor data points");
  Iromegan.zeros(3,1,"Iromegan");
}

void Motor::MotorCalcs(double MINSIGNAL_IN,double MAXSIGNAL_IN,double Tdatapt,double pwm_datapt,double omegaRPMdatapt,double R_IN) {
  
  //Compute Kt
  float dpwm = pwm_datapt - MINSIGNAL_IN;
  kt = Tdatapt/(dpwm*dpwm);
  
  //Angular Velocity computation
  a = (omegaRPMdatapt*2*PI/60.0)/dpwm;

  //Constants
  R = R_IN;
  rho = 0.00238; //slug/ft^3
  A = PI*R*R;  //ft^2

  //Compute ct
  ct = kt/(rho*A*R*R*a*a); //.0335
  cq = sqrt(ct*ct*ct)/sqrt(2); //lbf-s/ft^2 /.004

  MINSIGNAL = MINSIGNAL_IN;
  MAXSIGNAL = MAXSIGNAL_IN;

  //printf("Ct and Cq = %lf %lf \n",ct,cq);
  
  ////ALL OF THIS IS NO LONGER NECESSARY
  //Nemo Propulsion Parameters and some random ones
  //ct = 0.0335; //lbf/ft^2
  //R = 0.375; //feet
  //mass_prop = ((32.64/1000.0)*KG2LBF)/GEARTHENG;
  //mass_motor = ((47/1000.0)*KG2LBF)/GEARTHENG;

  //Propulsion Parameters from Bumblebee
  //This kt computation came from an experiment
  //R = (5.0/12.0); //5 inch props to feet

  //No need to change these
  //rho = 0.00238; //slug/ft^3
  //pwm_signal = MINSIGNAL;
  
  //printf("MINSIGNAL = %lf \n",MINSIGNAL); //This is the signal when the rotors aren't spinning
  //We don't necessarily want the rotors to get to that value
  //Ir = (2.0/3.0)*(mass_prop+mass_motor)*(2*R)*(2*R); //slugs-ft^2 -- Where did this come from????? - This is the inertia of the motor + prop and some other stuff I called Nemo to clarify
  //compute_aero();
}

double Motor::compute_signal_NO_OMEGA(double thrust_req) {
  //Assume Thrust = kt*(pwm_signal - MINSIGNAL)**2
  double pwm_req = sqrt(thrust_req/kt) + MINSIGNAL;
  return pwm_req;
}

void Motor::compute_thrust_NO_OMEGA() {
  thrust = kt*(pwm_signal - MINSIGNAL)*(pwm_signal - MINSIGNAL);
  Tn.mult(n,thrust);
}

void Motor::compute_torque_NO_OMEGA() {
  // printdouble(thrust,"thrust");
  // printdouble(R,"R");
  // printdouble(cq,"cq");
  // printdouble(ct,"ct");
  torque = thrust*R*cq/ct;
}

// Commented these all out so I don't confuse the crap out of everyone using
// these codes
// void Motor::compute_aero() {
//   compute_omega();
//   compute_thrust();
//   compute_torque();
// }

// void Motor::compute_omega() {
//   double omega_RPM = 0;
//   if (pwm_signal > MINSIGNAL) {
//     omega_RPM = a0*pwm_signal + a1; // this is in RPM
//   }
//   omega = omega_RPM*2*PI/60.0; //This converts to radians/second
//   //3.57 comes from the following
//   //OMEGA_MAX = KV*Voltage_battery
//   //KV = 2300
//   //Voltage_Battery = 14.8
//   //OMEGA_MAX is in RPM and must be converted to rad/s
//   //OMEGA_MAX = w*(MAXSIGNAL - MINSIGNAL) -- solve for w
//   //omega = 3.57*(pwm_signal - MINSIGNAL); -- This doesn't work though
//   Iromegan.mult(n,omega*sign*Ir);
// }

// double Motor::compute_signal(double thrust_req) {
//   double omega_req = sqrt(thrust_req/(ct*rho*A*(R*R))); //this is in rad/s
//   double omega_RPM_req = omega_req*60.0/(2*PI);
//   double pwm_req = (omega_RPM_req - a1)/a0;
//   return pwm_req;
// }

// void Motor::compute_thrust() {
//   thrust = ct*rho*A*(R*R)*omega*omega;
//   //cout << "thrust = " << thrust << endl;
//   Tn.mult(n,thrust);
// }

// void Motor::compute_torque() {
//   torque = cq*rho*A*(R*R*R)*omega*omega;
// }
