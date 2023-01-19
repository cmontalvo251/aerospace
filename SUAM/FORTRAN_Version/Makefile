OBJ=source/suam.o
CC = gfortran #Define compiler (CC) to be gfortran compiler
FFLAGS= -w -fno-automatic -O3 -ffree-form -ffree-line-length-none -fdollar-ok -fcray-pointer -fdefault-real-8 -fdefault-double-8
MAINFLAGS = -c -w -ffree-form -ffree-line-length-none -fdollar-ok -fcray-pointer -fdefault-real-8 -fdefault-double-8 -O3 #-fcheck=all
INCLUDE = #-I./libs/Linux -I./libs/Linux/modfiles - Archaic from AREA-I days
LIBS = #-L./libs/Linux -lWingsX - Archaic from AREA-I days
EXECUTABLE=SUAM.exe

all: $(OBJ) $(EXECUTABLE)

$(EXECUTABLE): $(OBJ)
	$(CC) $(FFLAGS) -o $(EXECUTABLE) $(OBJ)
source/suam.o: source/suam.f90
	rm -rf *.o SUAM.exe
	$(CC) $(MAINFLAGS) source/suam.f90
	mv suam.o source/suam.o
	mv bodydatatypes.mod source/bodydatatypes.mod
clean:
	rm -rf *.o SUAM.exe *.mod
	rm -rf source/*.o source/*.mod
rebuild:
	make clean
	make
