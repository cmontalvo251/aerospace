	subroutine SURROGEN(ikind)
        Implicit none
	INTEGER NDIM
	PARAMETER (NDIM=9999)
        real*8 xknot(NDIM,3),alf(NDIM,3)  !Three on alf for u,v,w
        real*8 x(999999,4),u(999999,3)
        real*8 A(NDIM,NDIM),zzz(ndim),rcond,rhs(ndim),chek(ndim) 
        real*8 xmax,xmin,ymax,ymin,zmax,zmin
        real*8 delx,dely,delz
        real*8 fac1,fac2,idx,ac,ctrnext,ctr,rctr
        real*8 rbf1,pi(3),qk(3),qj(3),val,val2
        integer i,j,k,L,M,ndata,nknot,nrank,JOB,ipvt(ndim)
        integer jay,kay,intdx,ikind
c$$$        open(unit=47,
c$$$     A   file='Output_Files/Constant_Time_UVW_Data_Aircraft_Ideal.txt',
c$$$     B   status='old')

        open(unit=47,
     A   file='Output_Files/Winds.OUT',
     B   status='old')

        open(unit=50,file='Output_Files/Xcen.dat',status='unknown')
        open(unit=51,file='Output_Files/Ualf.dat',status='unknown')
        open(unit=52,file='Output_Files/Valf.dat',status='unknown')
        open(unit=53,file='Output_Files/Walf.dat',status='unknown')

        nrank=1000
        JOB=0
        read(47,*) ndata,ndata,ndata,ndata,ndata,ndata,ndata
c         ndata=4000
        xmax=-77.7777D280
        xmin=77.77777D280
        ymax=-77.7777D280
        ymin=77.77777D280
        zmax=-77.7777D280
        zmin=77.77777D280
        do k=1,ndata
         read(47,*) x(k,1),x(k,2),x(k,3),x(k,4),u(k,1),u(k,2),u(k,3),ac
         if(x(k,1).gt.xmax) xmax=x(k,1)
         if(x(k,1).lt.xmin) xmin=x(k,1)
         if(x(k,2).gt.ymax) ymax=x(k,2)
         if(x(k,2).lt.ymin) ymin=x(k,2)
         if(x(k,3).gt.zmax) zmax=x(k,3)
         if(x(k,3).lt.zmin) zmin=x(k,3)
        enddo  
         
c$$$         write(*,*) xmin, xmax
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

	 write(*,*) 'SURROGATE MODEL ROUTINE'
	 write(*,*) ' '

	 ctrnext = 10
          do i=1,nrank
	     rctr = i
	     ctr = (rctr*100)/nrank
	     if (ctr .ge. ctrnext) then 
		write(*,*) 'Surrogate Model: ',int(ctrnext), ' % Complete'
		ctrnext = ctrnext + 10
	     end if
	     
           pi(1)=xknot(i,1)
           pi(2)=xknot(i,2)
           pi(3)=xknot(i,3) 
           do k=1,NdAta
            qk(1)=x(k,1)
            qk(2)=x(k,2)
            qk(3)=x(k,3)
            val=rbf1(pi,qk)
            rhs(i)=rhs(i)+u(k,1)*val
            do j=1,nrank
            qj(1)=xknot(j,1)
            qj(2)=xknot(j,2)
            qj(3)=xknot(j,3)
            val2=rbf1(qj,qk)
            A(i,j)=A(i,j)+val*val2
            enddo
           enddo
          enddo
C FACTOR....
        CALL DGECO(A,Ndim,Nrank,IPVT,RCOND,ZZZ)
        JOB=0
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
           do k=1,NdAta
            qk(1)=x(k,1)
            qk(2)=x(k,2)
            qk(3)=x(k,3)
            val=rbf1(pi,qk)
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

	write(*,*) ' '
	write(*,*) ' SURROGATE MODEL COMPLETE'

	close(50)
	close(51)
	close(52)
	close(53)
	close(47)

        stop
        end
