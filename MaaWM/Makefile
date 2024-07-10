OBJ=	source/MAAWM.o \
        source/dgeco.o \
        source/idamax.o source/dscal.o source/dgefa.o source/dgesl.o \
        source/dasum.o source/daxpy.o source/ddot.o   \
	source/rbf1.o source/surrogen.o \
	source/vecpoly.o \
	source/quadcopter.o
FFLAGS= -c -fno-automatic -O3
#MAINFLAGS = -c -w -ffree-form -ffree-line-length-none -fdollar-ok -fcray-pointer -fcheck=all -O3
MAINFLAGS = -c -w -ffree-form -ffree-line-length-none -fdollar-ok -fcray-pointer -fdefault-real-8 -fdefault-double-8 -O3 #-fcheck=all
EXECUTABLE = Run.exe
CC = gfortran #Define compiler (CC) to be gfortran compiler

all: $(OBJ) $(EXECUTABLE)

$(EXECUTABLE): $(OBJ)
	gfortran -O3 -o  Run.exe $(OBJ)
source/MAAWM.o: source/MAAWM.f90
	rm -rf Run.exe
	gfortran $(MAINFLAGS) source/MAAWM.f90
	mv MAAWM.o source/MAAWM.o
	mv maawmdatatypes.mod source/maawmdatatypes.mod
source/quadcopter.o: source/quadcopter.f90
	$(CC) $(MAINFLAGS) source/quadcopter.f90 -o source/quadcopter.o
source/vecpoly.o: source/vecpoly.f
	gfortran $(FFLAGS) source/vecpoly.f -o source/vecpoly.o
source/dasum.o: source/dasum.f
	gfortran $(FFLAGS) source/dasum.f -o source/dasum.o
source/daxpy.o: source/daxpy.f
	gfortran $(FFLAGS) source/daxpy.f -o source/daxpy.o
source/ddot.o: source/ddot.f
	gfortran $(FFLAGS) source/ddot.f -o source/ddot.o
source/dgeco.o: source/dgeco.f
	gfortran $(FFLAGS) source/dgeco.f -o source/dgeco.o
source/dgesl.o: source/dgesl.f
	gfortran $(FFLAGS) source/dgesl.f -o source/dgesl.o
source/dscal.o: source/dscal.f
	gfortran $(FFLAGS) source/dscal.f -o source/dscal.o
source/dgefa.o: source/dgefa.f
	gfortran $(FFLAGS) source/dgefa.f -o source/dgefa.o
source/idamax.o: source/idamax.f
	gfortran $(FFLAGS) source/idamax.f -o source/idamax.o
source/rbf1.o: source/rbf1.f
	gfortran $(FFLAGS) source/rbf1.f -o source/rbf1.o
source/surrogen.o: source/surrogen.f
	gfortran $(FFLAGS) source/surrogen.f -o source/surrogen.o
clean:
	rm -rf source/*.o Run.exe source/*.mod
rebuild:
	make clean
	make
