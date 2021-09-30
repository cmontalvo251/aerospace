// Multiple CubeSat Simulation Software 
// Copyright - Carlos Montalvo 2015
// You may freely distribute this file but please keep my name in here
// as the original owner

#ifndef MAIN_H
#define MAIN_H

#include "../MATLAB.h"
#include "../timer.h" 

///////////////Sanity checks//////////

//If openGL flag is on you need to include boost since they run in parallel
#ifdef OPENGL_FLAG
#define BOOST_FLAG
#endif
//Note, COMMS FLAG is ON YOU MUST SIMULATE IN REALTIME
#ifdef COMMS_FLAG
#define REALTIME_FLAG
#endif
//If OPENGL or SENSOR FLAG IS ONE YOU MUST SIMULATE FOREVER
#ifdef OPENGL_FLAG
#define SIMFOREVER_FLAG
#endif
#ifdef SENSOR_FLAG
#define SIMFOREVER_FLAG
#endif

void Create_Initial_Conditions(MATLAB*,MATLAB*,int);

#endif
