// Multiple CubeSat Simulation Software 
// Copyright - Carlos Montalvo 2015
// You may freely distribute this file but please keep my name in here
// as the original owner

#include "../MATLAB.h"
#include "../timer.h" //this is a symbolic links to BlackBox.git/C++ //this is a symbolic links to BlackBox.git/C++ 
#include "satellite_physics.h"
#include "../opengl.h" //this is a symbolic links to BlackBox.git/C++ //this is a symbolic links to BlackBox.git/C++ 
//This creates a variable called glhandle_g (_g for global)
//You need opengl for this
// sudo apt-get install freeglut3-dev
// For windows
// download glut

void Visualizer(Satellite** p_USASat,int argc,char** argv) {

  printf("Initializing Visualizer\n");

  //OPENGL glhandle_g; -- This variable is created in opengl.h
  //Once you call this the entire routine stops and nothing can run after this
  //So you need to give it a pointer to an object file with a method called
  //getState(MATLAB,MATLAB) and getTime() There is a VISOBJECT class prototype
  //in opengl.h that you can use to do what you need. 
  #ifdef OPENGL_FLAG
  //Hardcode p_USASat[0] since this only works with 1 satellite right now
  glhandle_g.Initialize(&p_USASat[0]->statetime,argc,argv,10000,600,600,0);
  #endif

  return;

}

