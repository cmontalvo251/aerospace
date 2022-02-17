!!!!!!!!!!!!!!!!Module!!!!!!!!!!!!!!!!

module SIMDATATYPES

integer, parameter :: NOSTATE = 6						!Number of States (Default = 2)
real*8, parameter :: PI = 3.14159265358979323846		!Pi

type SATELLITESTRUCTURE
  real*8 :: m = 80.0					!mass of diver(kg) 
  real*8 :: A = 0.7						!Cross-Sectional area of parachute (m^2)
  real*8 :: cd = 0.2					!drag coefficient
  real*8 :: psi = 0.785398163
  real*8 :: theta = 1.2
  real*8 :: phi = 0
  real*8 :: STATE(12)					!Position 
  real*8 :: STATEDOT(12)				!Velocity 
  real*8 :: NEWSTATE(12)				!Array to hold updated values
  real*8 :: TIB(3,3) = 0.0				!Transformation from Inertial to Body
  real*8 :: TBI(3,3) = 0.0				!Body to Inertial Transformation Matrix
  real*8 :: Mult(3,1) = 0.0				!Array to use in force calculations
  real*8 :: FORCEX = 0.0
  real*8 :: FORCEY = 0.0
  real*8 :: FORCEZ = 0.0
end type SATELLITESTRUCTURE

type ENVIRONMENTSTRUCTURE
  real*8 :: GRAVITY = 9.81				!gravity(kg/m*s^2)
  real*8 :: RHO = 1.225					!density of air
end type ENVIRONMENTSTRUCTURE

type SIMSTRUCTURE
  real*8 :: T							!Time
  real*8 :: TFINAL
  real*8 :: TIMESTEP					!Timestep
  type(SATELLITESTRUCTURE) :: SAT       !This calls SATELLITESTRUCTURE to be used inside SIMSTRUCTURE
  type(ENVIRONMENTSTRUCTURE) :: ENV		!This call ENVIRONMENTSTRUCTURE to be used inside SIMSTRUCTURE
end type SIMSTRUCTURE

end module

!!!!!!!!!!!!!!!!PROGRAM!!!!!!!!!!!!!!!!

PROGRAM Satellite
Use SIMDATATYPES
IMPLICIT NONE

integer openflag	 
real*8 i				!# of iterations for do loop
real*8 n				!ending iteration
type(SIMSTRUCTURE) EVE  !EVE stands for all variables in module

openflag = 0

!Initial Conditions
EVE%SAT%STATE(1) = 0.0		!X													
EVE%SAT%STATE(2) = 0.0      !Y
EVE%SAT%STATE(3) = 0.0		!Z
EVE%SAT%STATE(4) = 100.0		!U
EVE%SAT%STATE(5) = 0.0		!V 
EVE%SAT%STATE(6) = 0.0		!W

EVE%T = 0                   !Time   
EVE%TIMESTEP = 0.001         !Timestep
EVE%TFINAL = 10             !Final Time

call Derivatives(EVE)		!Calls subroutine to calculate variables

open(unit=25,file ='SatelliteResults.txt',iostat=openflag) !Create/Replace file for log purposes

n = EVE%TFINAL/EVE%TIMESTEP !Define ending iteration for do loop as TFINAL/TIMESTEP

write(*,*) 'npoints = ', n

do i = 1,n

   write(*,*) 'Simulation ', i/n*100,' % Complete'                      	!Print to screen program completion percentage

   write(25,fmt = '(1000F30.10)'), EVE%T,EVE%SAT%STATE(1:NOSTATE)  !Log results to previously opened file
	
   call Derivatives(EVE)  !Call the derivative values calculated in sub

   !Euler's Intergration
   EVE%SAT%NEWSTATE = EVE%SAT%STATE + EVE%TIMESTEP*EVE%SAT%STATEDOT 
  
   EVE%SAT%STATE = EVE%SAT%NEWSTATE

   EVE%T = EVE%T + EVE%TIMESTEP

end do

END PROGRAM Satellite                   

Subroutine Derivatives(EVE)                                              !This subroutine calculates the derivatives
use SIMDATATYPES
implicit none
real*8 x, y, z, u, v, w, xdot, ydot, zdot, udot, vdot, wdot, ctheta, stheta, cphi, sphi, cpsi, spsi, velocity 
type(SIMSTRUCTURE) EVE

!Helpful variables
  ctheta = cos(EVE%SAT%theta)
  stheta = sin(EVE%SAT%theta)
  cphi = cos(EVE%SAT%phi)
  sphi = sin(EVE%SAT%phi)
  cpsi = cos(EVE%SAT%psi)
  spsi = sin(EVE%SAT%psi)

!Inertial to Body Transformation matrix
EVE%SAT%TIB(1,1) = ctheta*cpsi
EVE%SAT%TIB(1,2) = (sphi*stheta*spsi) - (cphi*spsi)
EVE%SAT%TIB(1,3) = (cphi*stheta*cpsi) + (sphi*spsi)
EVE%SAT%TIB(2,1) = (ctheta*spsi) 
EVE%SAT%TIB(2,2) = (sphi*stheta*spsi) + (ctheta*cpsi)
EVE%SAT%TIB(2,3) = (cphi*stheta*spsi) - (sphi*cpsi)
EVE%SAT%TIB(3,1) = -stheta
EVE%SAT%TIB(3,2) = sphi*ctheta
EVE%SAT%TIB(3,3) = cphi*ctheta

!Transpose Transformation matrix
EVE%SAT%TBI = transpose(EVE%SAT%TIB)


!Unwrap State Vector
x = EVE%SAT%STATE(1)
y = EVE%SAT%STATE(2)
z = EVE%SAT%STATE(3)
u = EVE%SAT%STATE(4)
v = EVE%SAT%STATE(5)
w = EVE%SAT%STATE(6)

!Call Gravity and Density Model
call GRAVITYMODEL(EVE)
call DENSITYSOS(EVE)

!Kinematic Relationship
!write(*,*) EVE%SAT%STATE(4:6)
EVE%SAT%STATEDOT(1:3) = matmul(EVE%SAT%TIB,EVE%SAT%STATE(4:6))
!write(*,*) EVE%SAT%STATEDOT(1:3)
!PAUSE
!Calculate velocity
velocity = sqrt(EVE%SAT%STATE(4)**2 + EVE%SAT%STATE(5)**2 + EVE%SAT%STATE(6)**2)

!Calculate Forces
EVE%SAT%FORCEX = EVE%SAT%m*EVE%ENV%GRAVITY*EVE%SAT%TBI(1,3) - 0.5*EVE%ENV%RHO*velocity*EVE%SAT%A*EVE%SAT%cd*EVE%SAT%STATE(4)
EVE%SAT%FORCEY = EVE%SAT%m*EVE%ENV%GRAVITY*EVE%SAT%TBI(2,3) - 0.5*EVE%ENV%RHO*velocity*EVE%SAT%A*EVE%SAT%cd*EVE%SAT%STATE(5)
EVE%SAT%FORCEZ = EVE%SAT%m*EVE%ENV%GRAVITY*EVE%SAT%TBI(3,3) - 0.5*EVE%ENV%RHO*velocity*EVE%SAT%A*EVE%SAT%cd*EVE%SAT%STATE(6)

!Dynamic Relationship
EVE%SAT%STATEDOT(4) = EVE%SAT%FORCEX/EVE%SAT%m
EVE%SAT%STATEDOT(5) = EVE%SAT%FORCEY/EVE%SAT%m
EVE%SAT%STATEDOT(6) = EVE%SAT%FORCEZ/EVE%SAT%m

end Subroutine Derivatives

Subroutine GRAVITYMODEL(EVE)
use SIMDATATYPES
implicit none
type(SIMSTRUCTURE) EVE

!Gravity Model
EVE%ENV%GRAVITY = 9.81

end Subroutine GRAVITYMODEL

Subroutine DENSITYSOS(EVE)
use SIMDATATYPES
implicit none
type(SIMSTRUCTURE) EVE

!Density Model
EVE%ENV%RHO = 1.225

end Subroutine DENSITYSOS
