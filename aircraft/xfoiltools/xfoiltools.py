import numpy as np
import matplotlib.pyplot as plt

def computeclcd(AoA,saveFlnmAF,saveFlnmCp,pdfhandle,PLOTAIRFOIL):
    # %% READ DATA FILE: AIRFOIL
    dataBuffer = np.loadtxt(saveFlnmAF, skiprows=0)

    # Extract data from the loaded dataBuffer array
    XB = dataBuffer[:,0]
    YB = dataBuffer[:,1]
    
    # %% READ DATA FILE: PRESSURE COEFFICIENT
    dataBuffer = np.loadtxt(saveFlnmCp, skiprows=1)
    # Extract data from the loaded dataBuffer array
    X_0  = dataBuffer[:,0]
    Cp_0 = dataBuffer[:,1]

    # %% EXTRACT UPPER AND LOWER AIRFOIL DATA
    # Split airfoil into (U)pper and (L)ower
    XB_U = XB[YB >= 0]
    XB_L = XB[YB < 0]
    YB_U = YB[YB >= 0]
    YB_L = YB[YB < 0]
    # Split XFoil results into (U)pper and (L)ower
    Cp_U = Cp_0[YB >= 0]
    Cp_L = Cp_0[YB < 0]
    X_U  = X_0[YB >= 0]
    X_L  = X_0[YB < 0]

    ##The upper surface is backwards which messes with lift calcs
    XB_U = XB_U[-1::-1]
    YB_U = YB_U[-1::-1]
    Cp_U = Cp_U[-1::-1]
    X_U = X_U[-1::-1]

    # %% PLOT DATA
    # Plot airfoil
    if PLOTAIRFOIL == 1:
        plt.figure()
        plt.plot(XB_U,YB_U,'b.-',label='Upper')
        plt.plot(XB_L,YB_L,'r.-',label='Lower')
        plt.xlabel('X-Coordinate')
        plt.ylabel('Y-Coordinate')
        plt.title('Airfoil')
        plt.axis('equal')
        plt.grid()
        plt.legend()
        pdfhandle.savefig()

    # Plot pressure coefficient
    plt.figure()
    plt.plot(X_U,Cp_U,'b.-',label='Upper')
    plt.plot(X_L,Cp_L,'r.-',label='Lower')
    plt.xlim(0,1)
    plt.xlabel('X-Axis')
    plt.ylabel('Y-Axis')
    plt.title('Pressure Coefficient = ' + AoA + ' deg')
    plt.legend()
    plt.grid()
    plt.gca().invert_yaxis()
    pdfhandle.savefig()

    # plt.figure()
    # plt.plot(XB_U)
    # plt.plot(XB_L)
    # pdfhandle.savefig()

    ##Compute Lift
    cn = 0
    ca = 0
    for ctr in range(0,len(XB_L)-1):
        #print(XB_L[ctr])
        delX = XB_L[ctr+1]-XB_L[ctr]
        cn+=(Cp_L[ctr]-Cp_U[ctr])*delX
        delYu = YB_U[ctr+1]-YB_U[ctr]
        delYl = YB_L[ctr+1]-YB_L[ctr]
        dyudx = delYu/delX
        dyldx = delYl/delX
        ca+=(Cp_U[ctr]*dyudx - Cp_L[ctr]*dyldx)*delX
    print('AoA = ',AoA)
    print('Normal Coefficient = ',cn)
    print('Axial Coefficient = ',ca)
    alfa = np.float(AoA)*np.pi/180.0
    cl = cn * np.cos(alfa) - ca * np.sin(alfa)
    cd = cn * np.sin(alfa) + ca * np.cos(alfa)
    print('Lift Coefficient = ',cl)
    print('Drag Coefficient = ',cd)

    return cl,cd
