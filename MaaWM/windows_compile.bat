REM %comment?

REM %%%Run a make clean every time
del Run.exe
cd source
del *.o *.mod
cd ..\
REM %%%Assume that BlackBox/c++ is in this directory

REM %%%Code should be clean now. Time to compile
REM Compile each submodule
FOR %%A IN (vecpoly dasum daxpy ddot dgeco dgesl dscal dgefa idamax rbf1 surrogen) DO gfortran -c -fno-automatic -O3 source\%%A.f -o source\%%A.o
REM Compile the Main routine
gfortran -c -w -ffree-form -ffree-line-length-none -fdollar-ok -fcray-pointer -fcheck=all -O3 source/MAAWM.f90
gfortran -c -w -ffree-form -ffree-line-length-none -fdollar-ok -fcray-pointer -fdefault-real-8 -fdefault-double-8 -O3 source/quadcopter.f90 -o source/quadcopter.o
move MAAWM.o source/MAAWM.o
move maawmdatatypes.mod source/maawmdatatypes.mod
REM Link all object files into one mega massive exe
gfortran -O3 -o Run.exe source/MAAWM.o source/quadcopter.o source/vecpoly.o source/dasum.o source/daxpy.o source/ddot.o source/dgeco.o source/dgesl.o source/dscal.o source/dgefa.o source/idamax.o source/rbf1.o source/surrogen.o
REM Run.exe 
