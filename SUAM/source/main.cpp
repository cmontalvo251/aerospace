//Alright. Time for some c++ code.
#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include "Vehicle.h"

using namespace std;

char input_file[256];

//First thing to do is grab the input argument if it exists.
int main(int argc,char* argv[]) {
	if (argc > 1) {
		sprintf(input_file,argv[1]);
	}
	else {
		cout << "Sorry no input file given" << endl;
		cout << "Defaulting to Input_Files/NAH.files" << endl;\
		sprintf(input_file,"Input_Files/NAH.files");
	}

	//Then we can create the vehicle class which also imports the input deck
	Vehicle nah = Vehicle(input_file);

	//Then we can integrate the equations of motion
	if (nah.ok) {
		nah.Integrate();
	}
}