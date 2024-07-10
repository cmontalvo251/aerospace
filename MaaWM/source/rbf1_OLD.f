       function rbf1(P,Q)
       implicit none
       real*8 p(3),q(3),rbf1
       integer k
       rbf1=0.d0
       do k=1,3
       rbf1=rbf1+(p(k)-q(k))**2
       enddo
       rbf1=dsqrt(rbf1)
       return
       end

