#include "Propulsion.h"
#include "MATLAB.h"
#include "Motor.h"
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include "mathp.h"

using namespace std;

//Constructor
Propulsion::Propulsion() {
  Fmotors.zeros(3,1,"Force Motors");
  TauMotors.zeros(3,1,"Moment Motors");

  //Set defaults to Bumblebee
  //Estimate the mass of the vehicle
  mass_ = ((735./1000.)*KG2LBF)/GEARTHENG; //grams to kg to lbf to slugs
  //Aren't these supposed to be slug-ft^2???
  //double Ixx = (0.045*KG2LBF)*METERS2FT*METERS2FT/GEARTHENG; //kg-m^2 to slugs-ft^2
  double Iyy = (0.07413*KG2LBF)*METERS2FT*METERS2FT/GEARTHENG;
  double Ixx = Iyy;
  float a = (6.5/100.)*METERS2FT; //cm to meters to feet
  double Izz = Ixx + Iyy - mass_*a*a/6.0;
  I_.zeros(3,3,"Inertia Matrix");
  I_.set(1,1,Ixx);
  I_.set(2,2,Iyy);
  I_.set(3,3,Izz);
  //Position of motors
  //rx_ = (0.0762*METERS2FT); //Meters converted to Feet
  ry_ = (0.10115*METERS2FT);
  rx_ = ry_;
  rz_ = (0.0508*METERS2FT);
  rxpusher_ = 0.1524*METERS2FT;
  rxyz_.zeros(3,1,"Position of Motors");
  rxyz_.set(1,1,rx_);
  rxyz_.set(2,1,ry_);
  rxyz_.set(3,1,rz_);
  //MIN AND MAX SIGNAL DEFAULTS
  //These are from the receiver
  MINSIGNAL = SERVO_MIN; 
  MAXSIGNAL = SERVO_MAX;
  //These come from data sheets
  //at a signal of
  pwm_datapt_ = 1766.0;
  //thrust is
  Tdatapt_ = (5.0/GEARTH)*KG2LBF; //Newtons to kg to lbf
  //angular velocity is
  omegaRPMdatapt_ = 3500.0;
  //Rotor Size
  R_ = (5.0/12.0); //5 inch props to feet
  ///////////////////////////////////////////////
  //Wrap defaults into matlab vector
  datapts.zeros(3,1,"Data points for motors");
  datapts.set(1,1,Tdatapt_);
  datapts.set(2,1,pwm_datapt_);
  datapts.set(3,1,omegaRPMdatapt_);
}

//Motor Setup
void Propulsion::MotorsSetup(double mass,MATLAB I,double MINSIGNAL_IN,double MAXSIGNAL_IN,MATLAB rxyz,double rxpusher,MATLAB datapts_IN,double R_IN) {

  //Reset vars in the event vehicle is different than default
  rx_ = rxyz.get(1,1);
  ry_ = rxyz.get(2,1);
  rz_ = rxyz.get(3,1);
  rxpusher_ = rxpusher;
  I_.set(1,1,I.get(1,1));
  I_.set(2,2,I.get(2,2));
  I_.set(3,3,I.get(3,3));
  mass_ = mass;
  ////////////////////////////////////////////

  //Bottom Motors
  addMotor(rx_,-ry_,rz_,0.0,0.0,-1.0,-1);  //Add an extra variable here 
  addMotor(rx_,ry_,rz_,0.0,0.0,-1.0,1);  //that's a 0 or 1.
  addMotor(-rx_,ry_,rz_,0.0,0.0,-1.0,-1);/// var = 0, motorOff, var = 1, motorOn
  addMotor(-rx_,-ry_,rz_,0.0,0.0,-1.0,1);
  //Top Motors
  addMotor(rx_,-ry_,-rz_,0.0,0.0,-1.0,1); //top_front_left A ROTATION OF 1 = ccw and -1 = CW
  addMotor(rx_,ry_,-rz_,0.0,0.0,-1.0,-1); //top_front_right
  addMotor(-rx_,ry_,-rz_,0.0,0.0,-1.0,1); //top_back_right
  addMotor(-rx_,-ry_,-rz_,0.0,0.0,-1.0,-1); //top_back_left
  //Pusher
  addMotor(-rxpusher_,0.0,0.0,1.0,0.0,0.0,1);
  NUMMOTORS--;

  //Finalize motor calcs
  MotorBeep(MINSIGNAL_IN,MAXSIGNAL_IN,datapts_IN,R_IN);
}

void Propulsion::addMotor(double rx,double ry,double rz,double nx,double ny, double nz,int direction) {
  //Increment the number of motors by 1
  NUMMOTORS++;
  MOTORS[NUMMOTORS-1].r.set(1,1,rx);
  MOTORS[NUMMOTORS-1].r.set(2,1,ry);
  MOTORS[NUMMOTORS-1].r.set(3,1,rz);
  MOTORS[NUMMOTORS-1].n.set(1,1,nx);
  MOTORS[NUMMOTORS-1].n.set(2,1,ny);
  MOTORS[NUMMOTORS-1].n.set(3,1,nz);
  MOTORS[NUMMOTORS-1].sign = direction;
  cout << "Motor " << NUMMOTORS << " Set" << endl;
  // MOTORS[NUMMOTORS-1].r.disp();
  // MOTORS[NUMMOTORS-1].n.disp();
}

void Propulsion::RemoveMotors(int removemotors) {
  MOTORSOFF = removemotors;
  //With motors removed we need to reinitialize all the control matrices
  MOTORSRUNNING = NUMMOTORS - MOTORSOFF;

  printf("Total Motors = %d, Removing %d motors. Running Motors = %d \n",NUMMOTORS,removemotors,MOTORSRUNNING);

  //This is where we make our configuration matrices
  Hprime.zeros(4,MOTORSRUNNING,"Configuration Matrix Prime");
  HTprime.zeros(MOTORSRUNNING,4,"Configuration Matrix Transposed Prime");
  Qprime.zeros(MOTORSRUNNING,4,"Control Matrix Prime");
  CHIprime.zeros(MOTORSRUNNING,1,"Thrust Required Prime");

  for (int i = 1;i<=MOTORSRUNNING;i++) {
    Hprime.set(1,i,1.0);
    Hprime.set(2,i,-MOTORS[i-1].r.get(2,1));
    Hprime.set(3,i,-MOTORS[i-1].r.get(1,1));
    Hprime.set(4,i,MOTORS[i-1].sign);
  }
  //HT*inv(HHT)*M
  //H = 4xN
  //HT = Nx4
  //HHT = 4x4
  HTprime.transpose_not_square(Hprime);
  HHT.mult(Hprime,HTprime);

  //M.disp();
  //H.disp();
  //HT.disp();
  //HHT.disp();
  
  //Shit. Need a 4x4 inverse routine
  //Alright added crazy shit inverse functions
  HHT_inv.matrix_inverse(HHT,4);

  //HHT_inv.disp();

  //Then we need HT*inv(HHT) == N-1 x 4 - 4 x 4
  HT_inv_HHTprime.zeros(MOTORSRUNNING,4,"HT_inv_HHT");
  HT_inv_HHTprime.mult(HTprime,HHT_inv);
  //HT_inv_HHT.disp();

  //Finally we can get Q
  Qprime.mult(HT_inv_HHTprime,M);
  //Q.disp();
  
}

void Propulsion::MotorBeep(int Servo_Min,int Servo_Max,MATLAB datapts_IN,double R_IN) {
  //Alright we know the weight of the vehicle
  double weight = GRAVITYENG*mass_;
  //Rather than using a fancy calculation here I'm just going to compute the
  //thrust required for each motor
  double nominal_thrust = weight/(double)NUMMOTORS;

  //Pull out values for motors
  Tdatapt_ = datapts_IN.get(1,1);
  pwm_datapt_ = datapts_IN.get(2,1);
  omegaRPMdatapt_ = datapts_IN.get(3,1);
  R_ = R_IN;
  //////////////////////////////

  cout << "weight = " << weight << " mass = " << mass_ << endl;
  //cout << "thrust req = " << thrust_req << endl;

  //Compute Max thrust too
  MOTORSRUNNING = NUMMOTORS;

  //Set Nominal Thrust of all motors
  for (int i = 0;i<NUMMOTORS;i++) {
    //Re run motor calcs based on new values
    MOTORS[i].MotorCalcs(Servo_Min,Servo_Max,Tdatapt_,pwm_datapt_,omegaRPMdatapt_,R_);
    MOTORS[i].Nominal_Thrust = nominal_thrust;
    MOTORS[i].pwm_nominal = MOTORS[i].compute_signal_NO_OMEGA(nominal_thrust);
    //This computes Max Thrust so you have two data points
    MOTORS[i].pwm_signal = Servo_Max;
    MOTORS[i].compute_thrust_NO_OMEGA();
    MOTORS[i].Max_Thrust = MOTORS[i].thrust;
    MOTORS[i].pwm_signal = Servo_Min;
  }
  //Don't forget to run MotorCalc for the pusher
  cout << "NUMMOTORS = " << NUMMOTORS << endl;
  MOTORS[NUMMOTORS].MotorCalcs(Servo_Min,Servo_Max,Tdatapt_,pwm_datapt_,omegaRPMdatapt_,R_);
  cout << "MINSIGNAL = " << MOTORS[0].MINSIGNAL << endl;

  printf("Maximum Thrust at %d = %lf \n",Servo_Max,MOTORS[0].Max_Thrust);
  printf("Nominal Thrust at %lf = %lf \n",MOTORS[0].pwm_nominal,MOTORS[0].Nominal_Thrust);

  //This is where we make our configuration matrices
  H.zeros(4,NUMMOTORS,"Configuration Matrix");
  HT.zeros(NUMMOTORS,4,"Configuration Matrix Transposed");
  HHT.zeros(4,4,"HHT");
  HHT_inv.zeros(4,4,"HHT_inv");
  Q.zeros(NUMMOTORS,4,"Control Matrix");
  M.zeros(4,4,"Mass Matrix");
  M.set(1,1,mass_);
  M.set(2,2,I_.get(1,1));
  M.set(3,3,I_.get(2,2));
  M.set(4,4,I_.get(3,3));
  U.zeros(4,1,"Gamma Matrix");
  CHI.zeros(NUMMOTORS,1,"Thrust Required");

  //The controller works like this
  //M*U = H*CHI -- H is a 4xN thus there is no unique solution
  //thus you must use Lagranges Method -- Q = H'*inv(H*H')*M
  //then CHI = Q*U where U is our desired values
  

  for (int i = 1;i<=NUMMOTORS;i++) {
    H.set(1,i,1.0); //The 1 is because all motors create thrust
    H.set(2,i,-MOTORS[i-1].r.get(2,1));
    H.set(3,i,-MOTORS[i-1].r.get(1,1));
    H.set(4,i,MOTORS[i-1].sign);
  }
  //HT*inv(HHT)*M
  //H = 4xN
  //HT = Nx4
  //HHT = 4x4
  HT.transpose_not_square(H);
  HHT.mult(H,HT);

  //M.disp();
  //H.disp();
  //HT.disp();
  //HHT.disp();
  
  //Shit. Need a 4x4 inverse routine
  //Alright added crazy shit inverse functions
  HHT_inv.matrix_inverse(HHT,4);

  //HHT_inv.disp();

  //Then we need HT*inv(HHT) == N-1 x 4 - 4 x 4
  HT_inv_HHT.zeros(NUMMOTORS,4,"HT_inv_HHT");
  HT_inv_HHT.mult(HT,HHT_inv);
  //HT_inv_HHT.disp();

  //Finally we can get Q
  Q.mult(HT_inv_HHT,M);
  Q.disp();
}

void Propulsion::computeReconfigurable(double thrust,double d_roll,double d_pitch,double d_yaw) {
  //Set U
  U.set(1,1,thrust);
  U.set(2,1,d_roll);
  U.set(3,1,d_pitch);
  U.set(4,1,d_yaw);

    //Then we compute the Thrusts required for control
  //chi = Q*U
  //Q.disp();
  //U.disp();
  double value;
  if (MOTORSRUNNING == NUMMOTORS) {
    CHI.mult(Q,U);
    //Once we know the thrust required we need to compute the pwm signal required
    //Constrain Motors between 0 and Maximum thrust
    for (int i = 0;i<NUMMOTORS;i++) {
      //Need a bunch of saturation filters
      //Need to make sure if any thrust is less than zero we change it using the CONSTRAIN function
      value = CHI.get(i+1,1);
      value = CONSTRAIN(value,0.0,MOTORS[0].Max_Thrust);
      CHI.set(i+1,1,value);
      //Then compute Motor signals required for that amount of thrust
      MOTORS[i].pwm_signal = MOTORS[i].compute_signal_NO_OMEGA(CHI.get(i+1,1));
      //printf("Motor signals = %lf \n",MOTORS[i].pwm_signal);
    }
  } else {
    //Motors out
    //Qprime.disp();
    //U.disp();
    CHIprime.mult(Qprime,U);
    //We need to multiply CHI by a factor otherwise will blow this shit up
    //I don't think we need this anymore
    //CHIprime.mult_eq((double)MOTORSRUNNING/(double)NUMMOTORS);
    ///CHIprime.disp();
    //CHIprime.disp();
    for (int i = 0;i<NUMMOTORS;i++) {
      //Need a bunch of saturation filters
      //Need to make sure if any thrust is less than zero we change it
      if (i < MOTORSRUNNING) {
	value = CHIprime.get(i+1,1);
	value = CONSTRAIN(value,0.0,MOTORS[0].Max_Thrust);
	CHIprime.set(i+1,1,value);
	//Then compute Motor signals required for that amount of thrust
	MOTORS[i].pwm_signal = MOTORS[i].compute_signal_NO_OMEGA(CHIprime.get(i+1,1));
      } else {
	MOTORS[i].pwm_signal = MOTORS[0].MINSIGNAL;
      }
    }
  }
}

void Propulsion::ForceMomentMotors(MATLAB pqr) {
  Fmotors.mult_eq(0);
  TauMotors.mult_eq(0);
  //Then compute all our forces and torques
  for (int i = 0;i<NUMMOTORS+1;i++) { //Add 1 because of pusher
    //printf("Motor signals in force routine = %lf \n",MOTORS[i].pwm_signal);
    MOTORS[i].compute_thrust_NO_OMEGA();
    MOTORS[i].compute_torque_NO_OMEGA();
    //Fmotors = Fmotors + MOTORS[i].thrust*MOTORS[i].n
    //Filter Motors
    MOTORS[i].thrust_filtered += 0.005*(MOTORS[i].thrust - MOTORS[i].thrust_filtered);
    Fmotors.plus_mult_eq(MOTORS[i].n,MOTORS[i].thrust_filtered);
    //MOTORS[i].n.disp();
    //Fmotors.disp();
    //printdouble(MOTORS[i].torque,"torque");
    //Torque comes from three sources
    //First is the torque from the direction
    MOTORS[i].torque_filtered += 0.005*(MOTORS[i].torque - MOTORS[i].torque_filtered);
    TauMotors.plus_mult_eq(MOTORS[i].n,MOTORS[i].torque_filtered*MOTORS[i].sign);
    //The next source comes from rxF terms
    TauMotors.plus_cross_eq(MOTORS[i].r,MOTORS[i].Tn);
    //MOTORS[i].r.disp();
    //MOTORS[i].Tn.disp();
    //TauMotors.disp();
    //The final source comes from gyroscopic terms from angular velocity
    // TauMotors.plus_cross_eq(pqr,MOTORS[i].Iromegan);
    // printf("(%d) = %lf ",i,MOTORS[i].torque*MOTORS[i].sign);
    //printf("(%d) = %lf ",i,MOTORS[i].thrust);
  }
  //printf("i = %d \n",i);
  //printf("\n");
  //Fmotors.disp();
  //TauMotors.disp();
  //PAUSE();
}

void Propulsion::print_to_file(FILE* outfile) {
  for(int i = 0;i<NUMMOTORS+1;i++){
    fprintf(outfile,"%lf ",MOTORS[i].thrust_filtered);
    //    printdouble(MOTORS[i].thrust_filtered,"thrust");
  }
  //Print Motor Signal
  for(int i = 0;i<NUMMOTORS+1;i++){
    fprintf(outfile,"%lf ",MOTORS[i].pwm_signal);
  }
  //Print Motor torque
  for(int i = 0;i<NUMMOTORS+1;i++){
    fprintf(outfile,"%lf ",MOTORS[i].torque_filtered);
  }
}
