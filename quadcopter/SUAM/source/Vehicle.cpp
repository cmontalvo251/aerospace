#include "Vehicle.h"
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>

using namespace std;

//Constructor 
Vehicle::Vehicle(char input_file[]) {
  sprintf(input_file_,input_file);
  //Once we copy the file over it's time to import the input deck
  ok = ImportFiles();
  if (ok) {
    //Once we read in the data we need to do a bit of cleanup
    mass = massdata.get(1,1);
    I.eye(3,"I");
    I.set(1,1,massdata.get(2,1));
    I.set(2,2,massdata.get(3,1));
    I.set(3,3,massdata.get(4,1));
    rx = massdata.get(5,1);
    ry = massdata.get(6,1);
    rz = massdata.get(7,1);
    rxpusher = massdata.get(8,1);
    Iinv.zeros(3,3,"Iinv");
    Iinv.overwrite(I);
    Iinv.inverse();
    pqr.zeros(3,1,"pqr");
    I_pqr.zeros(3,1,"I_pqr");
    ptp.zeros(3,1,"ptp");
    uvw.zeros(3,1,"uvw");
    uvwdot.zeros(3,1,"uvwdot");
    pqrdot.zeros(3,1,"pqrdot");
    xyzdot.zeros(3,1,"xyzdot");
    ptpdot.zeros(3,1,"ptpdot");
    Kuvw_pqr.zeros(3,1,"Kuvw_pqr");
    pqrskew_I_pqr.zeros(3,1,"pqrskew_I_pqr");
    LMN.zeros(3,1,"LMN");
    //Create Some Derivatives Vectors
    F.zeros(3,1,"F");
    Fbody.zeros(3,1,"Fbody");
    Fgrav.zeros(3,1,"Fgrav");
    M.zeros(3,1,"M");
    stateinitial.copy_init(icdata,"Initial State");
    //stateinitial.disp();

    rxyz.zeros(3,1,"position of motors");
    rxyz.set(1,1,rx);
    rxyz.set(2,1,ry);
    rxyz.set(3,1,rz);

    //Setup motors
    propel.MotorsSetup(mass,I,propel.MINSIGNAL,propel.MAXSIGNAL,rxyz,rxpusher,propel.datapts,propel.R_);
    
    //Turn off some motors right away
    //propel.RemoveMotors(4); //This will remove that last 4 motors added
    propel.RemoveMotors(2); //This will remove top back left and top back right
    //propel.RemoveMotors(1);

    //Reset Motors so you can turn it back on
    propel.MOTORSRUNNING = propel.NUMMOTORS;
    
    cout << "Import Complete" << endl;
  }
  else {
    cout << "Import Failed" << endl;
  }
}

int Vehicle::ImportFiles() {
  FILE* mainfile;
  mainfile = fopen(input_file_,"r");
  char simfile[256],massfile[256],icfile[256],motorfile[256];
  //cout << "Reading Input File" << endl;
  if (mainfile) {
    fscanf(mainfile,"%s \n",&simfile);
    fscanf(mainfile,"%s \n",&massfile);
    fscanf(mainfile,"%s \n",&icfile);
    //Read in all file nows
    ok = ImportFile(simfile,&simdata,"simdata");
    if (ok) {
      ok = ImportFile(massfile,&massdata,"massdata");
      if (ok) {
	ok = ImportFile(icfile,&icdata,"icdata");
      } else {
	return 0;
      }
    } else {
      return 0;
    }
    //At this point the data should be read in. Should be able to run .disp() to see the contents.
    simdata.disp();
    massdata.disp();
    icdata.disp();
    fclose(mainfile);

    return 1;
  } else {
    cout << "File not found = " << input_file_ << endl;
    return 0;
  }
}

int Vehicle::ImportFile(char* filename,MATLAB* data,char* name) {
  cout << "Reading " << filename << endl;
  //Gonna open this file the old school way
  FILE *file;
  file = fopen(filename,"r");
  int length=0;
  char dummy[256];
  if (file) {
    while (!feof(file)) {
      fgets(dummy,256,file);
      //cout << dummy << endl;
      length++;
    }
    fclose(file);
    //cout << "File length is " << length << endl;
    //Now that we know how big the file is we need to make a MATLAB vector the same size
    data->zeros(length,1,name); //Because we passed a pointer we have to use the -> instead of .		//data.disp();
    //Then we open the file with FSTREAM and throw it in there. Done. Boom.
    fstream datafile;
    datafile.open(filename);
    if (datafile.is_open()) {
      string input;
      for (int idx = 0;idx<length;idx++) {
	getline(datafile,input);
	data->set(idx+1,1,atof(input.c_str()));
      }
      //For debug purposes let's make sure we read everything correctly. File I/O in C++
      //is always hit or miss for me.
      //data->disp();
    } else {
      cout << "Something went wrong in FSTREAM. Maybe the file wasn't closed properly?" << endl;
      return 0;
    }
  } else {
    cout << "File not found = " << filename << endl;
    return 0;
  }
  //This code will automatically generate a MATLAB vector based on how many rows are in
  //the file. Probably need to do an FSEEK or something and then do a ZEROS call to a MATLAB
  //Array.
  //If everything checks out we just return 1
  return 1;
}

void Vehicle::Integrate() {

  cout << "Begin Integration" << endl;

  //In order to run the RK4 engine we need the simdata
  //block.
  double tinitial = simdata.get(1,1);
  double tfinal = simdata.get(2,1);
  double timestep = simdata.get(3,1);
  //We don't need to make a time vector in C++. It's a waste of resources.
  //Honestly it's probably a waste of resource in all the other languages
  //too. But whatever. 
  double t = tinitial;

  //For c++ we have to open up an outfile. I'll just use the old school way.
  FILE* outfile = fopen("Output_Files/C++.OUT","w");

  //Create RK4 vectors.
  k1.copy_init(stateinitial,"k1");
  k2.copy_init(stateinitial,"k2");
  k3.copy_init(stateinitial,"k3");
  k4.copy_init(stateinitial,"k4");
  state.copy_init(stateinitial,"state");
  rkstate.copy_init(stateinitial,"rkstate");
  phi.copy_init(stateinitial,"phi");

  int MOTORSREMOVED = 0;

  while (t <= tfinal) {
    //cout << "Time = " << t << endl;
    //Output contents to file.
    fprintf(outfile,"%lf ",t);
    state.vecfprintf(outfile);
    //Print Motor Thrust
    propel.print_to_file(outfile);
    fprintf(outfile,"\n");
    //RK4 Calls
    Derivatives(k1,t,state);
    //void plus_mult_eq(MATLAB,double); //This is this = this + MATLAB*double
    rkstate.overwrite(state);
    rkstate.plus_mult_eq(k1,timestep/2.0);
    Derivatives(k2,t+timestep/2.0,rkstate);

    rkstate.overwrite(state);
    rkstate.plus_mult_eq(k2,timestep/2.0);
    Derivatives(k3,t+timestep/2.0,rkstate);

    rkstate.overwrite(state);
    rkstate.plus_mult_eq(k3,timestep);
    Derivatives(k4,t+timestep,rkstate);

    phi.mult_eq(0);
    phi.plus_mult_eq(k1,1.0/6.0);
    phi.plus_mult_eq(k2,2.0/6.0);
    phi.plus_mult_eq(k3,2.0/6.0);
    phi.plus_mult_eq(k4,1.0/6.0);

    //if ((t > 40) && (MOTORSREMOVED == 0)) {
      //MOTORSREMOVED = 1;
      //This hase been moved to before the loop
      //propel.RemoveMotors(5); //Testing the limits -- This won't work. Matrix becomes singular
      //propel.RemoveMotors(4); //This will remove that last 4 motors added
      //propel.RemoveMotors(3);
      //propel.RemoveMotors(2); //This will remove top back left and top back right
      //propel.RemoveMotors(1); 
    //}
    
    // if ((t > 90) && (t<110)) {
    //   propel.MOTORSRUNNING = propel.NUMMOTORS - propel.MOTORSOFF;
    // }
    // if (t > 110) {
    //   propel.MOTORSRUNNING = propel.NUMMOTORS;
    // }
  
    //Step state
    t+=timestep;
    state.plus_mult_eq(phi,timestep);
  }
  cout << "Integration Complete" << endl;
  fclose(outfile);
}

void Vehicle::Derivatives(MATLAB dxdt,double t,MATLAB state) {
  //double x = state.get(1,1);
  //double y = state.get(2,1);
  //double z = state.get(3,1);
  double phi = state.get(4,1);
  double theta = state.get(5,1);
  double psi = state.get(6,1);
  double u = state.get(7,1);
  double v = state.get(8,1);
  double w = state.get(9,1);
  double p = state.get(10,1);
  double q = state.get(11,1);
  double r = state.get(12,1);

  //Set up vectors
  ptp.set(1,1,phi);
  ptp.set(2,1,theta);
  ptp.set(3,1,psi);
  uvw.set(1,1,u);
  uvw.set(2,1,v);
  uvw.set(3,1,w);
  pqr.set(1,1,p);
  pqr.set(2,1,q);
  pqr.set(3,1,r);

  //Kinematics
  ine2bod123.L321(ptp,0);
  ine2bod123.rotateBody2Inertial(xyzdot,uvw);
  dxdt.vecset(1,3,xyzdot,1);
  ptpdot.mult(ine2bod123.H,pqr);
  dxdt.vecset(4,6,ptpdot,1);

  //Force Model
  ForceandMoment(state,dxdt,t);

  //Translational Dynamics
  Kuvw_pqr.cross(pqr,uvw);
  F.mult_eq(1.0/mass); //Will this work? I think so
  uvwdot.minus(F,Kuvw_pqr); ///Hey yo. We need to divide by mass here. I mean F needs to be divided by mass
  dxdt.vecset(7,9,uvwdot,1);

  //Rotational Dynamics
  //pqr.disp();
  I_pqr.mult(I,pqr);
  //I.disp();
  //I_pqr.disp();
  //PAUSE();
  pqrskew_I_pqr.cross(pqr,I_pqr);
  M.minus_eq(pqrskew_I_pqr); 
  pqrdot.mult(Iinv,M);
  dxdt.vecset(10,12,pqrdot,1);

  //Integrations
  dxdt.set(13,1,perror);
  dxdt.set(14,1,qerror);
  dxdt.set(15,1,rerror);
}

void Vehicle::feedback_linearization(MATLAB state,double t) {
  
  // roll_command = (50.0/500.0)*(float(rollrc)-1500.0);
  // pitch_command = (50.0/500.0)*(float(pitchrc)-1500.0);
  // yaw_rx = (135.0/500.0)*(float(yawrc)-1500.0); //this is really a yawRate command
  // thrust_desired = ((8.0*propel.MOTORS[0].Max_Thrust)/((double)THROTTLE_MAX-(double)THROTTLE_MIN))*(throttle-THROTTLE_MIN);
  double roll_command = 0*PI/180;
  double pitch_command = 10*PI/180;
  double yaw_rx = 0; //deg/s
  d_thrust = 1.005*(8.0*propel.MOTORS[0].Nominal_Thrust)/mass;
  //printf("%lf \n",propel.MOTORS[0].Nominal_Thrust);

  //Unwrap States
  double roll = state.get(4,1)*RAD2DEG;
  double pitch = state.get(5,1)*RAD2DEG;
  double yaw = state.get(6,1)*RAD2DEG;
  double roll_rate = state.get(10,1)*RAD2DEG;
  double pitch_rate = state.get(11,1)*RAD2DEG;
  double yaw_rate = state.get(12,1)*RAD2DEG;
  double perrorIntegral = state.get(13,1);
  double qerrorIntegral = state.get(14,1);
  double rerrorIntegral = state.get(15,1);
  
  float kp,kd,ki,kyaw;
  kp = 3.0;
  kd = 0.16;
  ki = 0.0225/1.0;
  //ki = 0.0;
  kyaw = 0.2;
  
  //Compute outer loop commands
  float phidot_command = kp*(roll_command - roll); //deg/s
  float thetadot_command = kp*(pitch_command - pitch); //deg/s

  //Constrain by 250
  phidot_command = CONSTRAIN(phidot_command,-250.0,250.0);
  thetadot_command = CONSTRAIN(thetadot_command,-250.0,250.0);

  //Yaw is a little wierd
  float psidot_command = yaw_rx;

  //These ptpdot commands are Euler angle derivatives and need to be converted to body frame angular rates
  float cos_theta = cos(pitch*DEG2RAD);
  float sin_theta = sin(pitch*DEG2RAD);
  float cos_phi = cos(roll*DEG2RAD);
  float sin_phi = sin(roll*DEG2RAD);
  float p_command = phidot_command - sin_theta*psidot_command;
  float q_command = cos_phi*thetadot_command + sin_phi*cos_theta*psidot_command;
  float r_command = -sin_phi*thetadot_command + cos_phi*cos_theta*psidot_command;

  //We then compute the inner loop commands which is a PI controller except for yaw
  perror = (p_command - roll_rate);
  qerror = (q_command - pitch_rate);
  rerror = (r_command - yaw_rate);

  //Then compute your commands note that yaw_out is just proportional
  d_roll = kd*perror + ki*perrorIntegral;
  d_pitch = kd*qerror + ki*qerrorIntegral;
  d_yaw = kyaw*rerror;

  //Constrain individual commands
  d_roll = CONSTRAIN(d_roll,-500,500);
  d_pitch = -CONSTRAIN(d_pitch,-500,500); //For some reason d_pitch is backwards
  d_yaw = -CONSTRAIN(d_yaw,-500,500); //Not sure why but this is negative too

  //printf("%lf %lf %lf %lf \n",d_thrust,d_roll,d_pitch,d_yaw);
  //PAUSE();
  
}

void Vehicle::pid_loop(MATLAB state,MATLAB dxdt,double t) {
    //Commands
  double zc = -10;
  //double phic = 0; //eventually this will be set by the user from the RC_input
  double v = state.get(8,1); //Feed side velocity to roll angle
  double phic = 0;
  // if (simtime > 50) {
  //   phic = -10*PI/180.0;
  // }
  double thetac = 0; //eventually this will be set by the user from the RC_input
  // if (t > 50) {
  //   thetac = -10*PI/180.0;
  // }
  double psic = 0; //eventually this will be set by the user from the RC_input
  if (t > 80) {
    psic = (t-80)*PI/(400.0);
    //psic = 100*PI/180.0;
  }

  //First compute our motor signals
  //propel.computeControl(state,dxdt,t);
  //Alright what we want is to set the U variable
  double kpz = 10;
  double kdz = 40;
  double z = state.get(3,1);
  double zdot = dxdt.get(3,1);
  double zdotc = 0;
  double zerror = (z-zc);
  zint += 0.01*zerror; //Note this is timestep dependent but I don't care I just need something that works right now
  d_thrust = kpz*zerror + kdz*(zdot-zdotc) + zint;

  double kpphi = 1;
  double kdphi = 2;
  double phi = state.get(4,1);
  double p = state.get(10,1);
  double pc = 0;
  double roll_error = (phi-phic);
  //roll_int += 0.001*roll_error;
  d_roll = -kpphi*(phi-phic) - kdphi*(p-pc) - state.get(13,1);

  double kptheta = -1;
  double kdtheta = -2;
  double theta = state.get(5,1);
  double q = state.get(11,1);
  double qc = 0;
  double theta_error = (theta-thetac);
  //theta_int += -0.001*theta_error;
  d_pitch = -kptheta*theta_error - kdtheta*(q-qc) - state.get(14,1);

  //printf("Thetac,d_pitch = %lf,%lf \n",thetac,d_pitch);

  double kppsi = -8;
  double kdpsi = -30;
  double psi = state.get(6,1);
  double r = state.get(12,1);
  double rc = 0;
  double yaw_error = (psi-psic);
  //yaw_int += 0.01*yaw_error;
  d_yaw = -kppsi*yaw_error - kdpsi*(r-rc) + state.get(15,1);
}

void Vehicle::ForceandMoment(MATLAB state,MATLAB dxdt,double t) {
  //Gravity in the inertial frame
  Fgrav.set(1,1,0);
  Fgrav.set(2,1,0);
  Fgrav.set(3,1,GRAVITYENG*mass);

  //Rotate Gravity to Body Frame
  ine2bod123.rotateInertial2Body(Fbody,Fgrav);

  //pid_loop(state,dxdt,t);
  feedback_linearization(state,t);

  //Then call the computeReconfigurable command
  propel.computeReconfigurable(d_thrust,d_roll,d_pitch,d_yaw);

  //Need Control for the pusher
  double uc = 0;
  // if (t > 40) {
  //   uc = 15;
  // }
  double u = state.get(7,1);
  double pwm_pusher = 10*(uc-u) + propel.MOTORS[8].MINSIGNAL;
  pwm_pusher = 0;
  pwm_pusher = CONSTRAIN(pwm_pusher,propel.MOTORS[8].MINSIGNAL,propel.MOTORS[8].MAXSIGNAL);
  //Send control to pusher
  propel.MOTORS[8].pwm_signal = pwm_pusher;
  // printdouble(propel.MOTORS[8].pwm_signal,"pusher signal");
  // printdouble(propel.MOTORS[8].MINSIGNAL,"MINSIGNAL");

  for (int i = 0;i<=8;i++) {
    if (i < propel.MOTORSRUNNING) {
      propel.MOTORS[i].pwm_signal = CONSTRAIN(propel.MOTORS[i].pwm_signal,IDLE,SERVO_MAX);
    } else {
      propel.MOTORS[i].pwm_signal;
    }
  }  
  
  //Then compute Force and Moment from Motors
  propel.ForceMomentMotors(pqr);

  //Fgrav.disp();
  //Fbody.disp();

  //Add everything together
  F.mult_eq(0);
  F.plus_eq(Fbody);
  F.plus_eq(propel.Fmotors);
  M.mult_eq(0);
  M.plus_eq(propel.TauMotors);

  //F.disp();
  //M.disp();
  //PAUSE();
  //printf("%lf \n",F.get(3,1));
  
}
