       function rbf1(P,Q,IK,EPS)
C ========================
C  IK=1  Linear spline
C  IK=2   
C  IK=3  Cubic spline
C  IK=4
C  IK=5  Quintic spline
C  IK=6
C  IK=7  Septic spline
C  IK=8  TBD
C  IK=9  TBD
C  IK=10 TBD
C  IK=11 TBD
C  IK=12 TBD
C ========================
       implicit none
       real*8 p(3),q(3),rbf1,EPS,rval
       integer k,IK
       rval=0.d0
       do k=1,3
       rval=rval+(p(k)-q(k))**2
       enddo
       rval=dsqrt(rval)
C  ODD ORDERS...
C  ===========================
       IF(IK.EQ.1) THEN
       rbf1=rval
       ENDIF
       IF(IK.EQ.3) THEN
       rbf1=rval**3
       ENDIF
       IF(IK.EQ.5) THEN
       rbf1=rval**5
       ENDIF
       IF(IK.EQ.7) THEN
       rbf1=rval**7
       ENDIF
C  =================================
C  EVEN ORDERS....
       IF(RVAL.GT.0.01D0) THEN
       IF(IK.EQ.2) THEN
        RBF1=RVAL**2*DLOG(RVAL)
       ENDIF
       IF(IK.EQ.4) THEN
        RBF1=RVAL**4*DLOG(RVAL)
       ENDIF
       IF(IK.EQ.6) THEN
        RBF1=RVAL**6*DLOG(RVAL)
       ENDIF
       ELSE !==============================
       IF(IK.EQ.2) THEN
        RBF1=RVAL**1*DLOG(RVAL**RVAL)
       ENDIF
       IF(IK.EQ.4) THEN
        RBF1=RVAL**3*DLOG(RVAL**RVAL)
       ENDIF
       IF(IK.EQ.6) THEN
        RBF1=RVAL**5*DLOG(RVAL**RVAL)
       ENDIF
       ENDIF
C =================================



       return
       end

