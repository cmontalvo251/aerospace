	subroutine vecpoly(Q,vec)
        implicit none
        real*8 Q(3),x,y,z,vec(35)
        x=Q(1)
        Y=Q(2)
        Z=Q(3)
        vec(1)=1.d0
        vec(2)=x
        vec(3)=y
        vec(4)=z        ! Linear
        vec(5)=x**2
        vec(6)=y**2
        vec(7)=z**2
        vec(8)=x*y
        vec(9)=x*z
        vec(10)=y*z     ! Quadratic
        vec(11)=x*x*y
        vec(12)=x*x*z 
        vec(13)=y*y*x
        vec(14)=y*y*z
        vec(15)=z*z*x
        vec(16)=z*z*y
        vec(17)=x**3
        vec(18)=y**3
        vec(19)=z**3 
        vec(20)=x*y*z   ! Cubic
        vec(21)=x**4
        vec(22)=y**4
        vec(23)=z**4
        vec(24)=y*x**3
        vec(25)=z*x**3
        vec(26)=x*y**3
        vec(27)=z*y**3
        vec(28)=x*z**3
        vec(29)=y*z**3
        vec(30)=(x**2)*(y**2)
        vec(31)=(x**2)*(z**2)
        vec(32)=(z**2)*(y**2)
        vec(33)=x**2*y*z
        vec(34)=y**2*x*z
        vec(35)=z**2*x*y      ! Quartic
        return
        end 
