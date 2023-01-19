OBJ=	sixdof.f90
CC = gfortran #Define compiler (CC) to be gfortran compiler
FFLAGS= -w -fno-automatic -O3 -ffree-form -ffree-line-length-none -fdollar-ok -fcray-pointer -fdefault-real-8 -fdefault-double-8
MAINFLAGS = -c -w -ffree-form -ffree-line-length-none -fdollar-ok -fcray-pointer -fdefault-real-8 -fdefault-double-8 -O3 #-fcheck=all
INCLUDE = #-I./libs/Linux -I./libs/Linux/modfiles - Archaic from AREA-I days
LIBS = #-L./libs/Linux -lWingsX - Archaic from AREA-I days
EXECUTABLE=LinuxRun.exe

all: $(OBJ) $(EXECUTABLE)

$(EXECUTABLE): $(OBJ)
	$(CC) $(FFLAGS) -o $(EXECUTABLE) $(OBJ)
source/sixdof.o: sixdof.f90
#	mkdir Output_Files
	rm -rf *.o LinuxRun.exe
	$(CC) $(MAINFLAGS) sixdof.f90
clean:
	rm -rf *.o LinuxRun.exe *.mod
rebuild:
	make clean
	make
