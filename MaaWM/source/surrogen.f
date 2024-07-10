        subroutine surrogen(ikind,k35,time,ncenin)
c	program surrogen
        include 'preamble'
        real*8 xknot(NDIM,3),alf(NDIM,3)  !Three on alf for u,v,w
        real*8 x(999999,4),u(999999,3)
        real*8 A(NDIM,NDIM),zzz(ndim),rcond,rhs(ndim),chek(ndim) 
        real*8 xmax,xmin,ymax,ymin,zmax,zmin
        real*8 delx,dely,delz,time
	!!!!CHECK THESE VARIABLES AND ADD
	!!! idx,ax,ctrnext,ctr,rctr
        real*8 fac1,fac2,idx,ac,ctrnext,ctr,rctr
        real*8 epsilon
        real*8 rbf1,pi(3),qk(3),qj(3),val,val2
        real*8 vec(35)  ! Vector of quartic polynomial terms in 3-D.
        real*8 veeze(35,3)     
        integer i,j,k,L,M,ndata,nknot,nrank,JOB,ipvt(ndim)
        integer ncen,ncenx,nceny,ncenin
        integer jay,kay
        integer ikind
        integer ierr !!!!!ADD ierr as a variable
        integer k35

        write(*,*) 'TIME = ',time
        write(*,*) 'IKIND = ',ikind
        write(*,*) 'K35 = ',k35
        write(*,*) 'NCENINx = ',ncenin

c$$$	k35 = 35
	
c$$$        open(unit=47,
c$$$     A   file='Output_Files/Constant_Time_UVW_Data_Aircraft_Ideal.dat',
c$$$     B   status='old')

	!!!!!!!!!!!!COMMENT OPEN COMMAND ABOVE AND REPLACE WITH THIS
        open(unit=47,
     A   file='Output_Files/Winds.OUT',
     B   status='old')

	!!!!ADD OUTPUT FILES TO THIS and change the 4000 to blank
        open(unit=50,file='Output_Files/Xcen.dat',status='unknown')
        open(unit=51,file='Output_Files/Ualf.dat',status='unknown')
        open(unit=52,file='Output_Files/Valf.dat',status='unknown')
        open(unit=53,file='Output_Files/Walf.dat',status='unknown')
        open(unit=54,file='Output_Files/Uveeze.dat',status='unknown')
        open(unit=55,file='Output_Files/Vveeze.dat',status='unknown')
        open(unit=56,file='Output_Files/Wveeze.dat',status='unknown')

        epsilon=0.D0  ! Not used ATM.
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

         ncenx = ncenin
         nceny = ncenin
         !ncenx = 32
         !nceny = 32

         ncen=0
         do jay=1,ncenx
          fac1=dfloat(jay)/dfloat(ncenx+2)
         do kay=1,nceny
          ncen=ncen+1
          fac2=dfloat(kay)/dfloat(nceny+2)
          xknot(ncen,1)=xmin+delx*fac1
          xknot(ncen,2)=ymin+dely*fac2
          xknot(ncen,3)=zmin+delz*rand()
         enddo
         enddo
         write(50,*) ncen
         write(51,*) ncen
         write(52,*) ncen
         write(53,*) ncen
         write(*,*) 'NCEN = ',ncen
         do i=1,ncen
          write(50,*) xknot(i,1),xknot(i,2),
     A     xknot(i,3)
         enddo

         nrank=ncen+k35
 
c        nrank=ncen

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

C The idea is we need an extra 35 rows/columns.
C The first step will be to add in the 35 columns
C for the ncen rows...
           do k=1,ndata

	     !!!!!ADD THESE FOR HELP WITH PROGRESS
	     rctr = k
	     ctr = (rctr*100)/ndata
	     if (ctr .ge. ctrnext) then 
		write(*,*) 'Surrogate Model: ',int(ctrnext), ' % Complete'
		ctrnext = ctrnext + 10
	     end if

            qk(1)=x(k,1)
            qk(2)=x(k,2)
            qk(3)=x(k,3)
            call vecpoly(qk,vec)
           do i=1,ncen
C =========================
            pi(1)=xknot(i,1)
            pi(2)=xknot(i,2)
            pi(3)=xknot(i,3)
            val=rbf1(pi,qk,ikind,epsilon)
            rhs(i)=rhs(i)+u(k,1)*val
            do j=1,ncen
            qj(1)=xknot(j,1)
            qj(2)=xknot(j,2)
            qj(3)=xknot(j,3)
            val2=rbf1(qj,qk,ikind,epsilon)
            A(i,j)=A(i,j)+val*val2
            enddo
            do L=1,k35
             j=ncen+L
             A(i,j)=A(i,j)+val*vec(L)
            enddo
C========================
           enddo
C Now, we need the final 35 rows.....
            do L=1,k35
             i=ncen+L
             val=vec(L)
             rhs(i)=rhs(i)+u(k,1)*val
             do j=1,ncen
            qj(1)=xknot(j,1)
            qj(2)=xknot(j,2)
            qj(3)=xknot(j,3)
            val2=rbf1(qj,qk,ikind,epsilon)
            A(i,j)=A(i,j)+val*val2
             enddo
            do M=1,k35
             j=ncen+M
             A(i,j)=A(i,j)+val*vec(M)
            enddo 
            enddo
          enddo
C FACTOR....
        CALL DGECO(A,Ndim,Nrank,IPVT,RCOND,ZZZ)
        JOB=0
        write(*,*) 'Here is rcond: ',rcond
        open(unit=57,file='Output_Files/RCOND.txt',position='append',
     &       status='old')
        write(57,*) ndata,rcond
        close(57)
        CALL DGESL(A,Ndim,Nrank,IPVT,RHS,JOB)
        do k=1,ncen
         alf(k,1)=rhs(k)
         write(50+1,*) Alf(K,1)
        enddo
        do k=1,k35
         veeze(k,1)=rhs(k+ncen)
         write(54,*) veeze(k,1)
        enddo
         
          
C GET THE OTHER TWO RHSs AND SOLVE .....
        DO Jay=2,3 
        do k=1,nrank
         rhs(k)=0.d0
        enddo
           do k=1,ndata
           do i=1,ncen
           pi(1)=xknot(i,1)
           pi(2)=xknot(i,2)
           pi(3)=xknot(i,3)
            qk(1)=x(k,1)
            qk(2)=x(k,2)
            qk(3)=x(k,3)
            call vecpoly(qk,vec)
            val=rbf1(pi,qk,ikind,epsilon)
            rhs(i)=rhs(i)+u(k,Jay)*val
           enddo
           do L=1,k35
            i=ncen+L
            val=vec(L) 
            rhs(i)=rhs(i)+u(k,Jay)*val
           enddo
          enddo
        JOB=0
        CALL DGESL(A,Ndim,Nrank,IPVT,RHS,JOB)
        do k=1,ncen
         alf(k,Jay)=rhs(k)
         write(50+Jay,*) Alf(K,Jay)
        enddo
        do k=1,k35
         veeze(k,Jay)=rhs(k+ncen)
         write(53+Jay,*) veeze(k,Jay)
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
