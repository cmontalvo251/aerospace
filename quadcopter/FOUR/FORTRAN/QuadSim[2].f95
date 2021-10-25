PROGRAM SMD
  implicit none
  include 'Globals.f'
  real*8 k1(12,1),k2(12,1),k3(12,1),k4(12,1),Rkphi(12,1),State(12,1),TIME
  integer TIMENEXT

  write(*,*) 'Quadcopter Simulation'
  !Read Inputs
  call READINPUTS
  State = IState
  call WRITELOG
  !Intialize State Derivatives
  !write(*,*) State
  call Derivatives(0,State)
  !Open the output file
  open(unit=95, file='SimulationResults.txt')
  
  omegavec(1,1) = 0.0
  omegavec(2,1) = 0.0
  omegavec(3,1) = 0.0
  omegavec(4,1) = 0.0
  
  
  WAYPOINT=1
  
  !Loop until TIME > tfinal
  TIMENEXT = 0
  do TIME= tinit,tfinal+timestep,timestep
     if (int(TIME/tfinal*100) .eq. TIMENEXT) then
		write(*,*) 'Simulation ', int(TIME/tfinal*100),'Percent Complete'
		TIMENEXT = TIMENEXT + 1
	 end if
     !Print State to file
	  !write(*,*) State(12,1)
     write(95, fmt='(1000F30.10)') TIME,State(1:12,1),omegavec
   
     call Derivatives(TIME,State)
	 !Runge-Kutta Integration
	 k1 = StateDot
	 !write(*,*) StateDot
	 !PAUSE
	 call Derivatives(TIME+timestep/2,State+k1*timestep/2)
	 k2=StateDot
	 call Derivatives(TIME+timestep/2,State+k2*timestep/2)
	 k3=StateDot
	 call Derivatives(TIME+timestep,State+k3*timestep)
	 k4=StateDot
	 Rkphi=1.0D0/6.0D0*(k1+2.0D0*k2+2.0D0*k3+k4)
	 State = State + Rkphi*timestep
   
  end do
  close(95)
END PROGRAM SMD

SUBROUTINE WRITELOG
  implicit none
  include 'Globals.f'
  open(unit=92,file='Logfile.txt')
  write(92,*) 'Input File Read'
  write(92,*) 'Tinitial = ',tinit
  write(92,*) 'Tfinal = ',tfinal
  write(92,*) 'Timestep = ',timestep
  close(92)
END SUBROUTINE

SUBROUTINE READINPUTS
  implicit none 
  include 'Globals.f'
  open(unit=96,file ='inputs.txt',status='old')
  read(96,*) tinit
  read(96,*) tfinal
  read(96,*) timestep
  read(96,*) IState(1,1)
  read(96,*) IState(2,1)
  read(96,*) IState(3,1)
  read(96,*) IState(4,1)
  read(96,*) IState(5,1)
  read(96,*) IState(6,1)
  read(96,*) IState(7,1)
  read(96,*) IState(8,1)
  read(96,*) IState(9,1)
  read(96,*) IState(10,1)
  read(96,*) IState(11,1)
  read(96,*) IState(12,1)
  read(96,*) m
  read(96,*) g
  read(96,*) TIMEON
  read(96,*) xcommand
  read(96,*) ycommand
  read(96,*) zcommand
  read(96,*) phicommand
  read(96,*) thetacommand
  read(96,*) psicommand
  read(96,*) kpphi
  read(96,*) kdphi
  read(96,*) kptheta
  read(96,*) kdtheta
  read(96,*) kppsi
  read(96,*) kdpsi
  read(96,*) kpx
  read(96,*) kdx
  read(96,*) kpy
  read(96,*) kdy
  read(96,*) kpz
  read(96,*) kdz
  close(96)
  
END SUBROUTINE READINPUTS

SUBROUTINE Derivatives(TIME,State)
  implicit none 
  include 'Globals.f'
  real*8 phi,theta,psi,uvw(3,1),p,q,r,pqr(3,1),u,v,w
  real*8 ctheta,stheta,ttheta,cphi,sphi,spsi,cpsi,Tib(3,3),xyzdot(3,1),H(3,3),ptpdot(3,1)
  real*8 skew(3,3),XYZA(3,1), XYZW(3,1),XYZWb(3,1),XYZ(3,1),uvwdot(3,1),Ixx,Iyy,Izz,Imat(3,3),Iinv(3,3),LMN(3,1),pqrdot(3,1)
  real*8 Qnew(3,3),MATT(3,1), TIME,State(12,1),forcevec(4,1)
  real*8 kt, b, Rnew, rho, pi, Ct, Ctau, lphi, ltheta, omegaHover
  real*8 thrust, fa, fb, fc, fd, tauphi, tautheta, taupsi, tauaerovec(3,4), Tauaeromat(3,1)
  real*8 sumomega, Ir, omegar, Als, Alc, dxD, dyD, Gamma(3,1), Tauw(3,1)
  real*8 omegaNot, domegaLeft, domegaRight, omegaLeft, omegaRight, add, domegaFront, omegaFront, omegaBack
  real*8 domegaDiag, omegaDiag, omegaOpp, addroll, addpitch, addyaw, kp, kd
  real*8 xc(4,1), yc(4,1), zc(4,1), dist
  
	  phi = State(4,1)
	  theta = State(5,1)
	  psi= State(6,1)
	  uvw(1,1)=State(7,1)
	  uvw(2,1)=State(8,1)
	  uvw(3,1)=State(9,1)
	  u = State(7,1)
	  v = State(8,1)
	  w = State(9,1)
	  p = State(10,1)
	  q = State(11,1)
	  r = State(12,1)
	  pqr(1,1)=State(10,1)
	  pqr(2,1)=State(11,1)
	  pqr(3,1)=State(12,1)

	  !Notation for trigonometric funtions 
	  ctheta = cos(theta); stheta = sin(theta); ttheta = stheta / ctheta;
	  cphi = cos(phi); sphi = sin(phi); spsi = sin(psi); cpsi = cos(psi);
	  
	  !Tranformation matrix (Tib inertial to body)
	  Tib(1,1) = ctheta * cpsi
	  Tib(2,1) = ctheta * spsi
	  Tib(3,1) = -stheta
	  Tib(1,2) = sphi * stheta * cpsi - cphi * spsi
	  Tib(2,2) = sphi * stheta * spsi + cphi * cpsi
	  Tib(3,2) = sphi * ctheta
	  Tib(1,3) = cphi * stheta * cpsi + sphi * spsi
	  Tib(2,3) = cphi * stheta * spsi - sphi * cpsi
	  Tib(3,3) = cphi * ctheta
	 
	 !Kinematics
	  xyzdot = matmul(Tib,uvw)
	  
	!Rotational Kinematics
	  !H matrix
	  H(1,1) = 1 ; H(1,2) = sphi*ttheta ; H(1,3) = cphi*ttheta
	  H(2,1) = 0 ; H(2,2) = cphi        ; H(2,3) = -sphi
	  H(3,1) = 0 ; H(3,2) = sphi/ctheta ; H(3,3) = cphi/ctheta

	  ptpdot = matmul(H,pqr)
	  
	  !Skew Symmetric Matrix
	  skew(1,1) = 0  ; skew(1,2) = -r ;  skew(1,3) = q
	  skew(2,1) = r  ; skew(2,2) = 0  ;  skew(2,3) = -p
	  skew(3,1) = -q ; skew(3,2) = p  ;  skew(3,3) = 0
	  
	  !Body frame forces
	  !Aerodynamic forces
	  Alc = 0.5D0 ; Als = 0.5D0 ; dxD = 0.5D0 ; dyD = 0.5D0
	  Rnew = 0.12065D0; rho = 1.2D0; pi = 3.141592654D0; Ct = 0.05D0; Ctau = 0.1D0
	  kt = Ct*((4.0D0*rho*(Rnew**4))/(pi**2)) ; b = Ctau*((4.0D0*rho*(Rnew**5))/(pi**3))
	  lphi = 0.3302D0 ; ltheta = 0.2413D0
	  
	! ! !Roll controller
	    ! ! !For a pure roll manuever, we want F_T=mg
	     ! domegaRight = -kpphi*(State(4,1)-phicommand) - kdphi*State(10,1)
		   ! add = sqrt((m*g - 2.0D0*kt*(omegaNot+domegaRight)**2) /(2.0D0*kt))
	    ! ! !!!!!!!!!!!!!!domegaLeft = omegaNot - add
		   ! omegaRight = omegaNot + domegaRight
		  ! omegaLeft = add !!!omegaLeft= omegaNot - domegaLeft
	     ! omegavec(1,1) = omegaLeft
	     ! omegavec(2,1) = omegaRight
	    ! omegavec(3,1) = omegaRight
	    ! omegavec(4,1) = omegaLeft
		   
		   ! ! ! ! !Pitch Controller
	       ! domegaFront = -kptheta*(State(5,1)-thetacommand) - kdtheta*State(11,1)
	       ! add = sqrt((m*g - 2.0D0*kt*(omegaNot+domegaFront)**2) /(2.0D0*kt))
	       ! omegaFront = omegaNot + domegaFront
	       ! omegaBack = add
	       ! omegavec(1,1) = OmegaFront
	       ! omegavec(2,1) = OmegaFront
	       ! omegavec(3,1) = OmegaBack
	       ! omegavec(4,1) = OmegaBack
	  
	! ! !Yaw Controller
	   ! domegaDiag = -kppsi*(State(6,1)-psicommand) - kdpsi*State(12,1)
	   ! add = sqrt((m*g - 2.0D0*kt*(omegaNot+domegaDiag)**2) /(2.0D0*kt))
	   ! omegaDiag = omegaNot + domegaDiag
	   ! omegaOpp = add
	   ! omegavec(1,1) = OmegaDiag
	   ! omegavec(2,1) = OmegaOpp
	   ! omegavec(3,1) = OmegaDiag
	   ! omegavec(4,1) = OmegaOpp
	   
	    !Waypoint Control
	 xc(1,1) = 0.0D0
	 xc(2,1) = 10.0D0
	 xc(3,1) = 10.0D0
	 xc(4,1) = 0.0D0
	 yc(1,1) = 10.0D0
	 yc(2,1) = 10.0D0
	 yc(3,1) = 0.0D0
	 yc(4,1) = 0.0D0
	 zc(1,1) = -10.0D0
	 zc(2,1) = -20.0D0
	 zc(3,1) = -30.0D0
	 zc(4,1) = 0.0D0

	 xcommand = xc(WAYPOINT,1)
	 ycommand = yc(WAYPOINT,1)
	 zcommand = zc(WAYPOINT,1)
	 
	 dist= sqrt((xcommand - State(1,1))**2 + (ycommand - State(2,1))**2 + (zcommand - State(3,1))**2)
	 
	 if (dist .lt. 0.2D0) then
	    WAYPOINT= WAYPOINT+1
		if (WAYPOINT .gt. 4) then
			WAYPOINT=1
		end if
	end if

	  !Nominal thrust so that F_T=mg and altitude hold controller
	  omegaHover = sqrt(m*g/(4.0D0*kt))
	  omegaNot = 1.0D0*sqrt(m*g/(4.0D0*kt)) - kpz*(zcommand-State(3,1))-kdz*(-xyzdot(3,1))

	  !Maximum speed of the rotors is 1100 rad/s
	  if (omegaNot .gt. 1100D0) then
	  omegaNot=1100D0
	  end if 
	  
	  !This means that the rotors can't spin in reverse
	  if (omegaNot .lt. 0.0D0) then
	    omegaNot = 0.0D0
	  end if 
	  
	  phicommand = kpy*(State(2,1)-ycommand) + kdy*(xyzdot(2,1))
	  !phicommand = 0.01D0
	  if (abs(phicommand) .gt. 0.523599D0) then
			phicommand = sign(0.523599D0,phicommand)
	  end if
	  thetacommand = -kpx*(State(1,1)-xcommand)- kdx*(xyzdot(1,1))
	  !thetacommand = 0.0D0
	   if (abs(thetacommand) .gt. 0.523599D0) then
			thetacommand = sign(0.523599D0,thetacommand)
	  end if
	  psicommand = 0.0D0
	  
		  domegaLeft = 0.0D0
	  	  omegaRight = 0.0D0
		  domegaFront = 0.0D0
		  omegaBack= 0.0D0
		  domegaDiag= 0.0D0
		  omegaOpp = 0.0D0
	  
	 !Roll, Pitch and Heading Controller 
	    domegaLeft = kpphi*(State(4,1)-phicommand) + kdphi*State(10,1)
		  addroll = sqrt((m*g - 2.0D0*kt*(omegaHover+domegaLeft)**2) /(2.0D0*kt))
		  omegaRight = addroll - omegaHover
	   domegaFront = -kptheta*(State(5,1)-thetacommand) - kdtheta*State(11,1)
		  addpitch = sqrt((m*g - 2.0D0*kt*(omegaHover+domegaFront)**2) /(2.0D0*kt))
		  omegaFront = omegaHover + domegaFront
		  omegaBack = addpitch - omegaHover
		
	  domegaDiag = -kppsi*(State(6,1)-psicommand) - kdpsi*State(12,1)
		addyaw = sqrt((m*g - 2.0D0*kt*(omegaHover+domegaDiag)**2) /(2.0D0*kt))
	    omegaDiag = omegaHover + domegaDiag
		omegaOpp = addyaw - omegaHover
		
		!write(*,*) State
		!write(*,*) phicommand,kpphi,State(4,1),kdphi,State(10,1)
		!write(*,*) omegaNot,domegaDiag,omegaOpp,domegaFront,omegaBack,domegaLeft,omegaRight
		!write(*,*) addyaw
		!PAUSE
	
		omegavec(1,1) = omegaNot + domegaLeft + domegaFront + domegaDiag
		omegavec(2,1) = omegaNot + omegaRight + domegaFront + omegaOpp
		omegavec(3,1) = omegaNot + omegaRight + omegaBack + domegaDiag
		omegavec(4,1) = omegaNot + domegaLeft + omegaBack + omegaOpp
	  
	   if (omegavec(1,1) .gt. 1100D0) then
	   omegavec(1,1) = 1100D0
	   end if 
	   if (omegavec(2,1) .gt. 1100D0) then
	   omegavec(2,1) = 1100D0
	   end if 
	   if (omegavec(3,1) .gt. 1100D0)then 
	   omegavec(3,1) = 1100D0
	   end if 
	   if (omegavec(4,1) .gt. 1100D0) then
	   omegavec(4,1) = 1100D0
	   end if 
	  
	  !Testing Yaw Response
	  !omegavec(1,1) = 1.0*omegavec(1,1) -!1.048808848170152*omegavec(1,1)
	  !omegavec(2,1) = 1.0*omegavec(2,1) !0.948683298050514*omegavec(2,1)
	  !omegavec(3,1) = 1.0*omegavec(3,1) !1.048808848170152*omegavec(3,1)
	  !omegavec(4,1) = 1.0*omegavec(4,1) !0.948683298050514*omegavec(4,1)
	  
	  sumomega = sum(omegavec);
	  forcevec = kt*omegavec**2
	  thrust = kt*sum(omegavec**2)
	  
	  !Why is everything wrong?!
	  !write(*,*) 'omegavec = ',omegavec
	  !write(*,*) 'omegaleft = ',omegaLeft
	  !write(*,*) 'omegaright = ',omegaRight
	  !write(*,*) 'sumomega = ',sumomega
	  !write(*,*) 'forcevec = ',forcevec
	  !write(*,*) 'omega0 = ',omegaNot
	  !write(*,*) 'thrust = ',thrust
	  !write(*,*) 'gravity = ',m*g
	   !write(*,*) 'Omega should be ',sqrt(m*g/(4.0D0*kt)),thrust/4
	   
	  if (sumomega .lt. 1e-2) then
        XYZA(1,1) = 0
		XYZA(2,1) = 0
		XYZA(3,1) = 0
	  else		
		XYZA(1,1)= -thrust*(((Alc/(sumomega*Rnew))+dxD)*u - ((Als)/(sumomega*Rnew))*v)
		XYZA(2,1)= -thrust*(((Als)/(sumomega*Rnew))*u + (((Alc)/(sumomega*Rnew))+dyD)*v)
		XYZA(3,1)= -thrust
	  end if
		
	  !Weight
	  XYZW(1,1) = 0
	  XYZW(2,1) = 0
	  XYZW(3,1) = m*g
		
	  XYZWb=matmul(Tib,XYZW)

	  XYZ = XYZWb +XYZA
		
	  !Inertial properties
	  Ixx = 0.005D0
	  Iyy = 0.005D0
	  Izz = 0.009D0
	  
	  !Imatrix
	  Imat(1,1)= Ixx ; Imat(1,2) = 0   ; Imat(1,3) = 0
	  Imat(2,1) = 0  ; Imat(2,2) = Iyy ; Imat(2,3) = 0
	  Imat(3,1) = 0  ; Imat(3,2) = 0   ; Imat(3,3) = Izz
	 
	  !I inverse matrix
	  Iinv(1,1)= 1/Ixx ; Iinv(1,2) = 0     ; Iinv(1,3) = 0
	  Iinv(2,1) = 0    ; Iinv(2,2) = 1/Iyy ; Iinv(2,3) = 0
	  Iinv(3,1) = 0    ; Iinv(3,2) = 0     ; Iinv(3,3) = 1/Izz
	  
	  !!!!!!!!! Moments {LMN}^T = Gamma + Tau_w + Tau_A
	  
	  Ir = 0.005D0
	  
	  omegar = omegavec(1,1) - omegavec(2,1) + omegavec(3,1) - omegavec(4,1)
	  
	 
	  Gamma(1,1) = Ir * omegar * q
	  Gamma(2,1) = -Ir * omegar * p
	  Gamma(3,1) = 0
	 
	  !!!!!!!!! Aerodynamics
	  tauaerovec(1,1) = -kt*lphi  ; tauaerovec(1,2) = kt*lphi    ; tauaerovec(1,3) = kt*lphi  ; tauaerovec(1,4) = -kt*lphi
	  tauaerovec(2,1) = kt*ltheta ; tauaerovec(2,2) = -kt*ltheta ; tauaerovec(2,3) = -kt*ltheta ; tauaerovec(2,4) = kt*ltheta
	  tauaerovec(3,1) = b		; tauaerovec(3,2) = -b	 ; tauaerovec(3,3) = b	 ; tauaerovec(3,4) = -b
	  
	  Tauaeromat = matmul(tauaerovec,omegavec**2)
	  tauphi = tauaeromat(1,1) ; 
	  Tautheta = tauaeromat(2,1) ; 
	  Taupsi = tauaeromat(3,1)
	  
	  LMN(1,1) = Gamma(1,1)+ lphi*((forcevec(2,1)+forcevec(3,1))-(forcevec(1,1)+forcevec(4,1)))
	  LMN(2,1)= Gamma(2,1)+ ltheta*((forcevec(1,1)+forcevec(2,1)) -(forcevec(3,1)+forcevec(4,1)))
	  LMN(3,1)= Gamma(3,1)+ Tauaeromat(3,1)
	  !LMN(1,1)= Gamma(1,1)  + Tauaeromat(1,1)
	  !LMN(2,1)= Gamma(2,1)  + Tauaeromat(2,1)
	  !LMN(3,1)= Gamma(3,1)  + Tauaeromat(3,1)
	  
	  !write(*,*) 'L = ',LMN(1,1)
	  !write(*,*) 'M = ',LMN(2,1)
	  !write(*,*) 'N = ',LMN(3,1)
	  !write(*,*) 'calc = ',-kt*lphi*omega1**2 + kt*lphi*omega2**2 - kt*lphi*omega3**2 + kt*lphi*omega4**2
	  !write(*,*) 'forces = ',fa,fb,fc,fd
	  !write(*,*) 'My moment = '
	  !write(*,*) 'Tauaeromat = ',Tauaeromat
	  !write(*,*) 'Tauaerovec = ',Tauaerovec
	  !write(*,*) 'Omegavec = ',omegavec
	  
	  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  !Translational Dynamics
	  uvwdot=((1.0D0/m)*XYZ) - matmul(skew,uvw)
	  
	  !Rotational Dynamics
	  Qnew=matmul(skew,Imat)
	  MATT=LMN- matmul(Qnew,pqr)
	  pqrdot=matmul(Iinv,MATT)
	  
	  StateDot(1,1)= xyzdot(1,1)
	  StateDot(2,1)= xyzdot(2,1)
	  StateDot(3,1)= xyzdot(3,1)
	  StateDot(4,1)= ptpdot(1,1)
	  StateDot(5,1)= ptpdot(2,1)
	  StateDot(6,1)= ptpdot(3,1)
	  StateDot(7,1)= uvwdot(1,1)
	  StateDot(8,1)= uvwdot(2,1)
	  StateDot(9,1)= uvwdot(3,1)
	  StateDot(10,1)= pqrdot(1,1)
	  StateDot(11,1)= pqrdot(2,1)
	  StateDot(12,1)= pqrdot(3,1)
	 
	  
END SUBROUTINE Derivatives







