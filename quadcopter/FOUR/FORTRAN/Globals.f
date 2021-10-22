	  common xcommand,ycommand,zcommand,phicommand,thetacommand,psicommand,StateDot,m,g
      common tinit,tfinal,timestep,TIMEON,IState
	  common kpphi,kdphi,kptheta,kdtheta,kppsi,kdpsi,kpx,kdx,kpy,kdy,kpz,kdz,omegavec,omegavecS
	  common WAYPOINT
	  
      real*8 xcommand,ycommand,zcommand,phicommand,thetacommand,psicommand,IState(12,1),StateDot(12,1),m,g
      real*8 tinit,tfinal,timestep,TIMEON,omegavec(4,1),omegavecS(4,1)
	  real*8 kpphi,kdphi,kptheta,kdtheta,kppsi,kdpsi,kpx,kdx,kpy,kdy,kpz,kdz 
	  real*8 WAYPOINT