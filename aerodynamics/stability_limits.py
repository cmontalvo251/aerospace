import numpy as np
import matplotlib.pyplot as plt


def Cm(alpha,de,xcg):
    cbar = 0.064028838569138949
    xacw = cbar/4.0
    CMAC = -0.035
    CL0w = 0.21
    lt = (18.0/12.0)/3.28
    AR = (58**2)/(5.08*9.0)
    S = 0.30164839140151456
    b = np.sqrt(AR*S)
    bt = 0.5*b
    ARt = AR
    St = bt**2/ARt
    xsm = xcg/cbar - xacw/cbar
    Cm0 = CMAC + xsm*CL0w
    print Cm0
    VH = lt*St/(cbar*S)
    eta = 1.0
    Claw = (1.2-0.21)/(8.74-0) * 180.0/np.pi
    e = 1.0
    CLaw = Claw/(1+Claw/(np.pi*e*AR))
    Clat = Claw
    CLat = Clat/(1+Clat/(np.pi*e*ARt))
    Cma = xsm*(CLaw + eta*St/S*CLat) - VH*CLat
    CLtde = 0.355*St/S
    Cmde = (xsm*St/S - VH)*eta*CLtde
    
    return Cm0 + Cma*alpha + Cmde*de

xcg = np.linspace(0.0,0.2,10)
alpha = np.linspace(-20*np.pi/180.0,20*np.pi/180.0,100)
plt.close("all")

for cg in xcg:
    de = 0
    #xcg = 0.015
    Cm_curve = Cm(alpha,de,cg)
    plt.plot(alpha,Cm_curve,label=str(cg))

plt.legend()
plt.grid()
plt.show()