        subroutine SURROGEN(ikind)
C	program surrogen
        include 'preamble'
        real*8 xknot(NDIM,3),alf(NDIM,3)  !Three on alf for u,v,w
        real*8 x(999999,4),u(999999,3)
        real*8 A(NDIM,NDIM),zzz(ndim),rcond,rhs(ndim),chek(ndim) 
        real*8 xmax,xmin,ymax,ymin,zmax,zmin
        real*8 delx,dely,delz
	!!!!CHECK THESE VARIABLES AND ADD
	!!! idx,ax,ctrnext,ctr,rctr
        real*8 fac1,fac2,idx,ac,ctrnext,ctr,rctr
        real*8 epsilon
        real*8 rbf1,pi(3),qk(3),qj(3),val,val2
        real*8 vec(35)  ! Vector of quartic polynomial terms in 3-D.
        real*8 veeze(35)
        integer i,j,k,L,M,ndata,nknot,nrank,JOB,ipvt(ndim)
        integer jay,kay
        integer ikind
        integer ierr !!!!!ADD ierr as a variable
c$$$        open(unit=47,
c$$$     A   file='Constant_Time_UVW_Data_Aircraft_Ideal.dat',
c$$$     B   status='old')

	!!!!!!!!!!!!COMMENT OPEN COMMAND ABOVE AND REPLACE WITH THIS
        open(unit=47,
     A   file='Output_Files/Winds.OUT',
     B   status='old')

        open(unit=50,file='Output_Files/Xcen.dat',status='unknown')
        open(unit=51,file='Output_Files/Ualf.dat',status='unknown')
        open(unit=52,file='Output_Files/Valf.dat',status='unknown')
        open(unit=53,file='Output_Files/Walf.dat',status='unknown')

        nrank=1000
        epsilon=0.D0
c        ikind=1 !!!!!!!COMMMENT OUT IKIND
        JOB=0
c       read(47,*) ndata !!!Delete this line of code and make ndata an input
c	ndata=4000 !!!!COMMENT OUT NDATA AS WELL
        xmax=-77.7777D280
        xmin=77.77777D280
        ymax=-77.7777D280
        ymin=77.77777D280
        zmax=-77.7777D280
        zmin=77.77777D280 
        k = 1
        do while (1 .eq. 1) !!!Change this for loop to a while loop
	   !!!!!!!!!!EDIT THIS LINE OF CODE AC is the AIRCRAFT NUMBER
           !!!!!Add iostat and end
         read(unit=47,fmt=*,iostat=ierr,end=555) x(k,1),x(k,2),x(k,3),
     A          x(k,4),u(k,1),u(k,2),u(k,3),ac
         if(x(k,1).gt.xmax) xmax=x(k,1)
         if(x(k,1).lt.xmin) xmin=x(k,1)
         if(x(k,2).gt.ymax) ymax=x(k,2)
         if(x(k,2).lt.ymin) ymin=x(k,2)
         if(x(k,3).gt.zmax) zmax=x(k,3)
         if(x(k,3).lt.zmin) zmin=x(k,3)
         k = k + 1 !!!!ADD THIS k = k + 1 down here
        enddo  

        !!!!!!!!ADD THESE LINES OF CODE
 555    if (ierr .ne. 0) then
           ndata = k-1
           ierr = 0
        endif
        write(*,*) 'NDATA = ',ndata
         
c$$$         write(*,*) xmin, xmax !!!!!!COMMENT OUT XMIN WRITE STATEMENTS
c$$$         write(*,*) ymin,ymax
c$$$         write(*,*) zmin,zmax

         delx=xmax-xmin
         dely=ymax-ymin
         delz=zmax-zmin

         nrank=0
         do jay=1,32
          fac1=dfloat(jay)/dfloat(33+1)
         do kay=1,32
          nrank=nrank+1
          fac2=dfloat(kay)/dfloat(33+1)
          xknot(nrank,1)=xmin+delx*fac1
          xknot(nrank,2)=ymin+dely*fac2
          xknot(nrank,3)=zmin+delz*rand()
         enddo
         enddo
         write(50,*) nrank
         write(51,*) nrank
         write(52,*) nrank
         write(53,*) nrank
         do i=1,nrank
          write(50,*) xknot(i,1),xknot(i,2),
     A     xknot(i,3)
         enddo


C Build system matrix....
         do k=1,nrank
         rhs(k)=0.d0
          ipvt(k)=0
          zzz(k)=0.d0
         do j=1,nrank
          a(j,k)=0.d0
         enddo
         enddo

	 !!!!ADDD THIS HERE

	 write(*,*) 'SURROGATE MODEL ROUTINE'
	 write(*,*) ' '

	 ctrnext = 10  !!!NEW VARIABLE

          do i=1,nrank

	     !!!!!ADD THESE FOR HELP WITH PROGRESS
	     rctr = i
	     ctr = (rctr*100)/nrank
	     if (ctr .ge. ctrnext) then 
		write(*,*) 'Surrogate Model: ',int(ctrnext), ' % Complete'
		ctrnext = ctrnext + 10
	     end if


           pi(1)=xknot(i,1)
           pi(2)=xknot(i,2)
           pi(3)=xknot(i,3) 
           do k=1,ndata
            qk(1)=x(k,1)
            qk(2)=x(k,2)
            qk(3)=x(k,3)
            val=rbf1(pi,qk,ikind,epsilon)
            rhs(i)=rhs(i)+u(k,1)*val
            do j=1,nrank
            qj(1)=xknot(j,1)
            qj(2)=xknot(j,2)
            qj(3)=xknot(j,3)
            val2=rbf1(qj,qk,ikind,epsilon)
            A(i,j)=A(i,j)+val*val2
            enddo
           enddo
          enddo
C FACTOR....
        CALL DGECO(A,Ndim,Nrank,IPVT,RCOND,ZZZ)
        JOB=0
        write(*,*) 'Here is rcond: ',rcond
        CALL DGESL(A,Ndim,Nrank,IPVT,RHS,JOB)
        do k=1,nrank
         alf(k,1)=rhs(k)
         write(50+1,*) Alf(K,1)
        enddo
C GET THE OTHER TWO RHSs AND SOLVE .....
        DO Jay=2,3 
        do k=1,nrank
         rhs(k)=0.d0
        enddo
        do i=1,nrank
           pi(1)=xknot(i,1)
           pi(2)=xknot(i,2)
           pi(3)=xknot(i,3)
           do k=1,ndata
            qk(1)=x(k,1)
            qk(2)=x(k,2)
            qk(3)=x(k,3)
            val=rbf1(pi,qk,ikind,epsilon)
            rhs(i)=rhs(i)+u(k,Jay)*val
           enddo
          enddo
        JOB=0
        CALL DGESL(A,Ndim,Nrank,IPVT,RHS,JOB)
        do k=1,nrank
         alf(k,Jay)=rhs(k)
         write(50+Jay,*) Alf(K,Jay)
        enddo

        ENDDO ! ON JAY

	!!!!!!!!!!!!!COMMENT ALL THIS AS WELL
c$$$C Check a few:
c$$$
c$$$         qk(1)=x(514,1)
c$$$         qk(2)=x(514,2)
c$$$         qk(3)=x(514,3)
c$$$
c$$$         val=0.d0
c$$$         do i=1,nrank
c$$$          pi(1)=xknot(i,1)
c$$$          pi(2)=xknot(i,2)
c$$$          pi(3)=xknot(i,3) 
c$$$           val=val+rbf1(pi,qk,ikind,epsilon)*alf(i,2)
c$$$         enddo 
c$$$
c$$$         write(*,*) 'Check: ',val, 'against: ',u(514,2)
 
        stop
        end
