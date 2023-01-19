CC=g++
EXECUTABLE=SUAMCPP.exe
#MATLAB and mathp are symbolic links to BlackBox.git/C++
#In order for the gravity/magnetic model to compile you need to include
#all cpp files in Model_Resources/Resources
SOURCES=$(wildcard source/*.cpp) 
OBJECTS=$(SOURCES:.cpp=.o)
INCLUDES = -I ../../Navio
COMPILE=-c -w -O3

all: $(SOURCES) $(EXECUTABLE)
$(EXECUTABLE): $(OBJECTS)
	$(CC) $(INCLUDES) $(OBJECTS) -o $@ 
.cpp.o:
	$(CC) $(INCLUDES) $(COMPILE) $< -o $@ 
clean:
	rm -rf $(EXECUTABLE) *.o
	rm -rf source/*.o
rebuild:
	make clean
	make
