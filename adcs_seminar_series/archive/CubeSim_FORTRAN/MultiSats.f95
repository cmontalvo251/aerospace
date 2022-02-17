!This program contains a 6DOF Simulator with attitude control and magnet model for multiple CubeSats
!Authors: Nicholas Carroll and Carlos Montalvo
!Year: 2015

!!!!!!!!!!!!!!!!!!!!!!!!!!!Module!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module SIMDATATYPES

integer, parameter :: NOSTATE = 12					!Number of States (Default = 12)
integer, parameter :: MAXSATS = 4
real*8, parameter :: PI = 3.14159265358979323846			!Pi

type SATELLITESTRUCTURE
  real*8 :: m = 0.0							!mass of cubesat(kg) 
  real*8 :: A = 0.0							!Cross-Sectional area of parachute (m^2)
  real*8 :: cd = 0.0							!drag coefficient
  real*8 :: chord = 0.0							!Chord (m)
  real*8 :: b = 0.0							!Base (m)
  real*8 :: I(3,3) = 0.0 						!Moment of Inertia
  real*8 :: clp = 0.0
  real*8 :: cmq = 0.0
  real*8 :: cnr = 0.0
  real*8 :: STATE(12) = 0.0						!Original 
  real*8 :: STATEDOT(12) = 0.0						!Derivative 
  real*8 :: XDOT1(12) = 0.0
  real*8 :: XDOT2(12) = 0.0
  real*8 :: XDOT3(12) = 0.0
  real*8 :: XDOT4(12) = 0.0
  real*8 :: STATENOMINAL(12)
  real*8 :: NEWSTATE(12)						!Array to hold updated values
  real*8 :: TIB(3,3) = 0.0						!Transformation from Inertial to Body
  real*8 :: TBI(3,3) = 0.0						!Body to Inertial Transformation Matrix
  real*8 :: FORCE(3) = 0.0
  real*8 :: MOMENT(3) = 0.0
  real*8 :: H(3,3) = 0.0
  real*8 :: wtilde(3,3) = 0.0
  real*8 :: THETACOMMAND = 0.0
  real*8 :: PSICOMMAND = 0.0
  real*8 :: Rm = 0.0
 
end type SATELLITESTRUCTURE

type ENVIRONMENTSTRUCTURE
  real*8 :: GRAVITY = 0.0			        		!gravity(kg/m*s^2)
  real*8 :: RHO = 0.0							!density of air
end type ENVIRONMENTSTRUCTURE

type SIMSTRUCTURE
  real*8 :: T								!Time
  real*8 :: TFINAL
  real*8 :: TIMESTEP							!Timestep
  integer :: NUMSATS							!Number of Satellites
  type(SATELLITESTRUCTURE) :: SAT(MAXSATS)       			!This calls SATELLITESTRUCTURE to be used inside SIMSTRUCTURE
  type(ENVIRONMENTSTRUCTURE) :: ENV					!This call ENVIRONMENTSTRUCTURE to be used inside SIMSTRUCTURE
end type SIMSTRUCTURE

end module
!!!!!!!!!!!!!!!!!!!!!!!!!!!End Module!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!PROGRAM!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM Satellite
Use SIMDATATYPES
IMPLICIT NONE

integer openflag	 
integer i								!# of iterations for do loop
integer j								!# of iterations for aircraft
real*8 n							
type(SIMSTRUCTURE) EVE                  				!EVE stands for all variables in module

openflag = 0
EVE%T = 0                   						!Time   
EVE%TIMESTEP = 0.001      						!Timestep
EVE%TFINAL = 50             						!Final Time
EVE%NUMSATS = 2								!Number of cubesats 

call ReadInputs(EVE)


do j = 1,EVE%NUMSATS
	
	call Derivatives(EVE,j)						!Calls subroutine to calculate variables
	
end do

open(unit=25,file ='SatelliteResults.txt',iostat=openflag) 		!Create/Replace file for log purposes

n = EVE%TFINAL/EVE%TIMESTEP 						!Define ending iteration for do loop as TFINAL/TIMESTEP

do i = 1,n

	write(25,fmt = '(1000F30.10)',advance = 'no'), EVE%T 		!Write time for entire SatelliteResults.txt file
	
	!Write Satellite 1 results
	do j = 1, EVE%NUMSATS - 1 

		write(25,fmt = '(1000F30.10)', advance = 'no') ,EVE%SAT(j)%STATE(1:NOSTATE)  
	
	end do

	!Write Satellite 2 results
	write(25,fmt = '(1000F30.10)'),EVE%SAT(EVE%NUMSATS)%STATE(1:NOSTATE)  

	do j = 1,EVE%NUMSATS

   		EVE%SAT(j)%STATENOMINAL = EVE%SAT(j)%STATE

   		!Runge-Kutta 4
   		call Derivatives(EVE,j)  				!Call the derivative values calculated in sub
   		EVE%SAT(j)%XDOT1 = EVE%SAT(j)%STATEDOT
   		EVE%SAT(j)%STATE = EVE%SAT(j)%STATENOMINAL + EVE%TIMESTEP*0.5*EVE%SAT(j)%XDOT1
   		call Derivatives(EVE,j)
   		EVE%SAT(j)%XDOT2 = EVE%SAT(j)%STATEDOT
   		EVE%SAT(j)%STATE = EVE%SAT(j)%STATENOMINAL + EVE%TIMESTEP*0.5*EVE%SAT(j)%XDOT2
   		call Derivatives(EVE,j)
   		EVE%SAT(j)%XDOT3 = EVE%SAT(j)%STATEDOT
   		EVE%SAT(j)%STATE = EVE%SAT(j)%STATENOMINAL + EVE%TIMESTEP*EVE%SAT(j)%XDOT3
   		call Derivatives(EVE,j)
   		EVE%SAT(j)%XDOT4 = EVE%SAT(j)%STATEDOT
   
   		EVE%SAT(j)%STATEDOT = (EVE%SAT(j)%XDOT1 + 2*EVE%SAT(j)%XDOT2 + 2*EVE%SAT(j)%XDOT3 + EVE%SAT(j)%XDOT4)/6.0

   		!Euler's Intergration
   		EVE%SAT(j)%NEWSTATE = EVE%SAT(j)%STATENOMINAL + EVE%TIMESTEP*EVE%SAT(j)%STATEDOT 
  
   		EVE%SAT(j)%STATE = EVE%SAT(j)%NEWSTATE

	end do
    
    !Update time
    EVE%T = EVE%T + EVE%TIMESTEP

end do

END PROGRAM Satellite 
!!!!!!!!!!!!!!!!!!!!!!!!!!!End Program!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!Subroutines!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                  

!!!!!!!!!!!!!!!!!!!!!!!!!!!Derivatives!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Subroutine Derivatives(EVE,aflag)                                              
use SIMDATATYPES
implicit none
integer aflag
real*8 x, y, z, phi, theta, psi, u, v, w, p, q, r 
real*8 xdot, ydot, zdot, udot, vdot, wdot, ctheta, stheta, ttheta
real*8 cphi, sphi, cpsi, spsi, velocity, uvw(3)
real*8 chord, b, clp, cmq, cnr, wtilde, phat, qhat, rhat, Iinv(3,3), sub1, sub2, sub3
real*8 Kp, Kd, Km, phidot, thetadot, psidot, pc, qc, rc, Rm, FX, xh, yh, zh
type(SIMSTRUCTURE) EVE

!Unwrap State Vector
x = EVE%SAT(aflag)%STATE(1)
y = EVE%SAT(aflag)%STATE(2)
z = EVE%SAT(aflag)%STATE(3)
phi = EVE%SAT(aflag)%STATE(4)
theta = EVE%SAT(aflag)%STATE(5)
psi = EVE%SAT(aflag)%STATE(6)
u = EVE%SAT(aflag)%STATE(7)
v = EVE%SAT(aflag)%STATE(8)
w = EVE%SAT(aflag)%STATE(9)
p = EVE%SAT(aflag)%STATE(10)
q = EVE%SAT(aflag)%STATE(11)
r = EVE%SAT(aflag)%STATE(12)


!Helpful variables
  ctheta = cos(theta)
  stheta = sin(theta)
  ttheta = tan(theta)
  cphi = cos(phi)
  sphi = sin(phi)
  cpsi = cos(psi)
  spsi = sin(psi)

!Inertial to Body Transformation matrix
EVE%SAT(aflag)%TIB(1,1) = ctheta*cpsi
EVE%SAT(aflag)%TIB(1,2) = (sphi*stheta*cpsi) - (cphi*spsi)
EVE%SAT(aflag)%TIB(1,3) = (cphi*stheta*cpsi) + (sphi*spsi)
EVE%SAT(aflag)%TIB(2,1) = (ctheta*spsi) 
EVE%SAT(aflag)%TIB(2,2) = (sphi*stheta*spsi) + (cphi*cpsi)
EVE%SAT(aflag)%TIB(2,3) = (cphi*stheta*spsi) - (sphi*cpsi)
EVE%SAT(aflag)%TIB(3,1) = -stheta
EVE%SAT(aflag)%TIB(3,2) = sphi*ctheta
EVE%SAT(aflag)%TIB(3,3) = cphi*ctheta

!Transpose Transformation matrix
EVE%SAT(aflag)%TBI = transpose(EVE%SAT(aflag)%TIB)

!wtilde matrix
EVE%SAT(aflag)%wtilde(1,1) = 0
EVE%SAT(aflag)%wtilde(1,2) = -r
EVE%SAT(aflag)%wtilde(1,3) = q
EVE%SAT(aflag)%wtilde(2,1) = r
EVE%SAT(aflag)%wtilde(2,2) = 0
EVE%SAT(aflag)%wtilde(2,3) = -p
EVE%SAT(aflag)%wtilde(3,1) = -q
EVE%SAT(aflag)%wtilde(3,2) = p
EVE%SAT(aflag)%wtilde(3,3) = 0

!H Matrix
EVE%SAT(aflag)%H(1,1) = 1
EVE%SAT(aflag)%H(1,2) = sphi*ttheta
EVE%SAT(aflag)%H(1,3) = cphi*ttheta
EVE%SAT(aflag)%H(2,1) = 0
EVE%SAT(aflag)%H(2,2) = cphi
EVE%SAT(aflag)%H(2,3) = -sphi
EVE%SAT(aflag)%H(3,1) = 0
EVE%SAT(aflag)%H(3,2) = sphi/ctheta
EVE%SAT(aflag)%H(3,3) = cphi/ctheta 

!Call Gravity and Density Model
call GRAVITYMODEL(EVE)
call DENSITYSOS(EVE)

!Kinematic Relationship
EVE%SAT(aflag)%STATEDOT(1:3) = matmul(EVE%SAT(aflag)%TIB,EVE%SAT(aflag)%STATE(7:9))

!Calculate velocity
velocity = sqrt(EVE%SAT(aflag)%STATE(7)**2 + EVE%SAT(aflag)%STATE(8)**2 + EVE%SAT(aflag)%STATE(9)**2)

sub1 = 0.5*EVE%ENV%RHO*velocity*EVE%SAT(aflag)%A*EVE%SAT(aflag)%cd*u
sub2 = 0.5*EVE%ENV%RHO*velocity*EVE%SAT(aflag)%A*EVE%SAT(aflag)%cd*v
sub3 = 0.5*EVE%ENV%RHO*velocity*EVE%SAT(aflag)%A*EVE%SAT(aflag)%cd*w

!Calculate Forces
EVE%SAT(aflag)%FORCE(1) = EVE%SAT(aflag)%m*EVE%ENV%GRAVITY*EVE%SAT(aflag)%TBI(1,3) - sub1 !X Force
EVE%SAT(aflag)%FORCE(2) = EVE%SAT(aflag)%m*EVE%ENV%GRAVITY*EVE%SAT(aflag)%TBI(2,3) - sub2 !Y Force
EVE%SAT(aflag)%FORCE(3) = EVE%SAT(aflag)%m*EVE%ENV%GRAVITY*EVE%SAT(aflag)%TBI(3,3) - sub3 !Z Force

!Magnet control variables
Km = 1
xh = EVE%SAT(2)%STATE(1) - EVE%SAT(1)%STATE(1)
yh = EVE%SAT(2)%STATE(2) - EVE%SAT(1)%STATE(2)
zh = EVE%SAT(2)%STATE(3) - EVE%SAT(1)%STATE(3)
EVE%SAT(1)%Rm = sqrt(xh**2 + yh**2 + zh**2)

!If t > 8 seconds then turn on magnet control
if (EVE%T .gt. 8) then
	FX = Km/EVE%SAT(1)%Rm**2
	EVE%SAT(aflag)%FORCE(1) = EVE%SAT(1)%FORCE(1) + FX
	if (FX .gt. 1000) then
		STOP
	end if
	
	else 
	
	FX = 0

end if

!Dynamic Relationship (udot, vdot, wdot)
EVE%SAT(aflag)%STATEDOT(7:9) = EVE%SAT(aflag)%FORCE/EVE%SAT(aflag)%m - matmul(EVE%SAT(aflag)%wtilde,EVE%SAT(aflag)%STATE(7:9))

!Phidot, Thetadot, Psidot
EVE%SAT(aflag)%STATEDOT(4:6) = matmul(EVE%SAT(aflag)%H,EVE%SAT(aflag)%STATE(10:12))
phidot = EVE%SAT(aflag)%STATEDOT(4)
thetadot = EVE%SAT(aflag)%STATEDOT(5)
psidot = EVE%SAT(aflag)%STATEDOT(6) 

!Phat, Qhat, Rhat
if (abs(velocity) .gt. 0) then
	phat = (p*b)/(2*velocity)
	qhat = (q*chord)/(2*velocity)
	rhat = (r*b)/(2*velocity)

	!Calculate Moments
	EVE%SAT(aflag)%MOMENT(1) = -0.5*EVE%ENV%RHO*(v**2)*EVE%SAT(aflag)%A*chord*clp*phat !L 
	EVE%SAT(aflag)%MOMENT(2) = -0.5*EVE%ENV%RHO*(v**2)*EVE%SAT(aflag)%A*chord*cmq*qhat !M
	EVE%SAT(aflag)%MOMENT(3) = -0.5*EVE%ENV%RHO*(v**2)*EVE%SAT(aflag)%A*chord*cnr*rhat !N

	else

	EVE%SAT(aflag)%MOMENT(1:3) = 0.0

end if

!Compute Attitude Control
call CONTROL(EVE,aflag)
Kp = 1
Kd = 0.5
pc = Kp*-phi + Kd*-phidot
qc = Kp*(EVE%SAT(aflag)%THETACOMMAND - theta) + Kd*-thetadot
rc = Kp*(EVE%SAT(aflag)%PSICOMMAND - psi) + Kd*-psidot
EVE%SAT(aflag)%MOMENT(1) = EVE%SAT(aflag)%MOMENT(1) - Km*(p-pc)
EVE%SAT(aflag)%MOMENT(2) = EVE%SAT(aflag)%MOMENT(2) - Km*(q-qc)
EVE%SAT(aflag)%MOMENT(3) = EVE%SAT(aflag)%MOMENT(3) - Km*(r-rc)


!Angular Velocity derivative (pdot, qdot, rdot)
Iinv = 0.0
Iinv(1,1) = 1/EVE%SAT(aflag)%I(1,1)
Iinv(2,2) = 1/EVE%SAT(aflag)%I(2,2)
Iinv(3,3) = 1/EVE%SAT(aflag)%I(3,3)
EVE%SAT(aflag)%STATEDOT(10:12) = matmul(Iinv,EVE%SAT(aflag)%MOMENT(1:3)-matmul(EVE%SAT(aflag)%wtilde,matmul(EVE%SAT(aflag)%I,EVE%SAT(aflag)%STATE(10:12))))

end Subroutine Derivatives
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!Gravity!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
Subroutine GRAVITYMODEL(EVE)
use SIMDATATYPES
implicit none
type(SIMSTRUCTURE) EVE

!Gravity Model
EVE%ENV%GRAVITY = EVE%ENV%GRAVITY

end Subroutine GRAVITYMODEL
!!!!!!!!!!!!!!!!!!!!!!!!!!!End Gravity!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!Density!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Subroutine DENSITYSOS(EVE)
use SIMDATATYPES
implicit none
type(SIMSTRUCTURE) EVE

!Density Model
EVE%ENV%RHO = 1.225

end Subroutine DENSITYSOS
!!!!!!!!!!!!!!!!!!!!!!!!!!!End Density!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!Control!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Subroutine CONTROL(EVE, aflag)
use SIMDATATYPES
implicit none
type(SIMSTRUCTURE) EVE

integer aflag
real*8 x, y, z, phi, theta, psi

!Determine R vector from one satellite to the other
if (aflag .eq. 2) then
	!From Sat 2 to Sat 1
	x = EVE%SAT(1)%STATE(1) - EVE%SAT(2)%STATE(1)
	y = EVE%SAT(1)%STATE(2) - EVE%SAT(2)%STATE(2)
	z = EVE%SAT(1)%STATE(3) - EVE%SAT(2)%STATE(3)
	EVE%SAT(aflag)%Rm = sqrt(x**2 + y**2 + z**2) 
else
	!From Sat 1 to Sat 2
	x = EVE%SAT(2)%STATE(1) - EVE%SAT(1)%STATE(1)
	y = EVE%SAT(2)%STATE(2) - EVE%SAT(1)%STATE(2)
	z = EVE%SAT(2)%STATE(3) - EVE%SAT(1)%STATE(3)
	EVE%SAT(aflag)%Rm = sqrt(x**2 + y**2 + z**2)
end if

theta = atan(-z, sqrt(x**2 + y**2))
psi = atan(y, x)
EVE%SAT(aflag)%THETACOMMAND = theta
EVE%SAT(aflag)%PSICOMMAND = psi
	
end Subroutine Control
!!!!!!!!!!!!!!!!!!!!!!!!!!!End Control!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!Read Inputs!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Subroutine ReadInputs(EVE)
use SIMDATATYPES
implicit none
type(SIMSTRUCTURE) EVE
integer i

open(97,file='Satellite.in',status='old')
read(97,*) EVE%SAT(1)%m		!Mass (kg)
read(97,*) EVE%SAT(1)%A		!Area (m^2)
read(97,*) EVE%SAT(1)%cd	!Drag Coefficient
read(97,*) EVE%ENV%RHO		!Density
read(97,*) EVE%ENV%GRAVITY	!Gravity (m/s^2)
read(97,*) EVE%SAT(1)%STATE(1)  !x
read(97,*) EVE%SAT(1)%STATE(2)  !y
read(97,*) EVE%SAT(1)%STATE(3)  !z
read(97,*) EVE%SAT(1)%STATE(4)  !psi
read(97,*) EVE%SAT(1)%STATE(5)  !theta
read(97,*) EVE%SAT(1)%STATE(6)  !phi
read(97,*) EVE%SAT(1)%STATE(7)  !u
read(97,*) EVE%SAT(1)%STATE(8)  !v
read(97,*) EVE%SAT(1)%STATE(9)  !w
read(97,*) EVE%SAT(1)%STATE(10) !p
read(97,*) EVE%SAT(1)%STATE(11) !q
read(97,*) EVE%SAT(1)%STATE(12) !r
read(97,*) EVE%SAT(1)%chord	!Height (m)
read(97,*) EVE%SAT(1)%b		!Width (m)
read(97,*) EVE%SAT(1)%I(1,1)	!Ixx
read(97,*) EVE%SAT(1)%I(2,2)	!Iyy	
read(97,*) EVE%SAT(1)%I(3,3)	!Izz
read(97,*) EVE%SAT(1)%clp 	!clp
read(97,*) EVE%SAT(1)%cmq 	!cmq
read(97,*) EVE%SAT(1)%cnr 	!cnr

!Set all of Cubesat2 variables equal to Cubesat1
do i = 2, EVE%NUMSATS
   EVE%SAT(i) = EVE%SAT(1)
end do

!Specifically set x,y,z,phi,theta,psi,u,v,and w of Cubesat2
EVE%SAT(2)%STATE(1) = EVE%SAT(1)%STATE(1) + 5
EVE%SAT(2)%STATE(2) = EVE%SAT(1)%STATE(2) + 4
EVE%SAT(2)%STATE(3) = EVE%SAT(1)%STATE(3) + 2
EVE%SAT(2)%STATE(4) = EVE%SAT(1)%STATE(4) + 0.2
EVE%SAT(2)%STATE(5) = EVE%SAT(1)%STATE(5) + 0.5 
EVE%SAT(2)%STATE(6) = EVE%SAT(1)%STATE(6) + 0.7
EVE%SAT(2)%STATE(10) = EVE%SAT(1)%STATE(10) + 0.3
EVE%SAT(2)%STATE(11) = EVE%SAT(1)%STATE(11) + 0.4
EVE%SAT(2)%STATE(12) = EVE%SAT(1)%STATE(12) + 0.1

end Subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!End Read Inputs!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!