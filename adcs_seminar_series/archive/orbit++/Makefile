#DIRECTIVE DEFINITIONS

#0 and 1 means defined or not
#if the variable is defined it's a 1 and if not
#it's a zero

#OPENGL_FLAG 
#0 = Do not open an openGL window
#1 = open an openGL window, for now this means the code will simulate 
#forever until you type in q for quit

#KALMAN_FLAG
#0 = Simulate Kalman Filter
#1 = Do not simulate kalman filter

#SENSOR_FLAG
#0 = simulate sensors
#1 = get sensor information from actual hardware

#BOOST_FLAG
#0 = Run everything serially
#1 = Run Everything in parrallel with multiple threads. You must have boost installed
#In order to get this to work on linux you need to run sudo apt-get install libboost-all-dev

#COMMS_FLAG
#0 = No communication required. Functions called internally using mutexes
#1 = Communication over SERIAL

#REALTIME_FLAG
#0 = Simulate as fast as possible
#1 = Simulate in Real Time

#PRINT_TIME
#0 = do not print time
#1 = print time of simulation

CC=g++
EXECUTABLE=Orbit.exe
#MATLAB and mathp are symbolic links to BlackBox.git/C++
#In order for the gravity/magnetic model to compile you need to include
#all cpp files in Model_Resources/Resources
SOURCES=$(wildcard Model_Resources/Resources/*.cpp) main.cpp satellite_physics.cpp ../MATLAB.cpp ../mathp.cpp ../timer.cpp visualizer.cpp ../opengl.cpp ../Rotation3.cpp
OBJECTS=$(SOURCES:.cpp=.o)
COMPILE=-c -w -O3
#PRINT_TO_FILE_FLAG - this will print all useful information to a .OUT file
#PRINT_TIME - print time to stdout
#BOOST_FLAG - use separate threads
#OPENGL_FLAG - runs visualizer - uses boost_flag by default
#REALTIME_FLAG - runs simulation in realtime to test HIL
#DEBUG_MODE - prints some useful stuff and uses PAUSE statements to step though
#KALMAN_FLAG - runs the kalman filter instead of using perfect sensor measurement
DIRECTIVES= -DPRINT_TIME -DPRINT_TO_FILE_FLAG
INCLUDE = -IModel_Resources/Include #This points to the .hpp for the gravity and magnet models
LIB = -L/usr/local/lib 
BOOST = -lpthread -lboost_system -lboost_thread -lboost_date_time
OPENGL = -lGL -lGLU -lglut

all: $(SOURCES) $(EXECUTABLE)
$(EXECUTABLE): $(OBJECTS)
	$(CC) $(INCLUDE) $(OBJECTS) -o $@ $(LIB) $(BOOST) $(OPENGL)
.cpp.o:
	$(CC) $(INCLUDE) $(COMPILE) $(DIRECTIVES) $< -o $@ $(LIB)
clean:
	rm -rf $(EXECUTABLE) *.o
	cd Model_Resources/Resources; rm -rf *.o
rebuild:
	make clean
	make
