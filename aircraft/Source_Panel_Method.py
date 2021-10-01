import numpy as np
import matplotlib.pyplot as plt

##User Inputs
R = 3.0 ##Radius of the cylinder
Npanels = 10 #Must be an even number and an integer!!!
Vinf = 20.0 ##Speed of the flow

##This is the equation of the surface. You must change this if you want something other than a cylinder
def surface(x):
    return np.sqrt(R**2-x**2)

##You must also define the limits of the body. In this case the LE is -R and the TE is R
xsmooth = np.linspace(-R,R,1000)
yu_smooth = surface(xsmooth)
yl_smooth = -surface(xsmooth)

##You also need to figure out a way to discretize the panels
##Get the Panel Coordinates
dx = 4*R/(Npanels)
xpanels = np.arange(-R+dx/2,R,dx)
yu_panels = surface(xpanels)
yl_panels = -surface(xpanels)

#########################EVERYTHING BELOW HERE SHOULD RUN JUST FINE##############################

#Combine upper and lower surfaces
xpanels = np.append(xpanels,xpanels[-1::-1])
ypanels = np.append(yu_panels,yl_panels[-1::-1])
xpanels = np.append(xpanels,xpanels[0])
ypanels = np.append(ypanels,yu_panels[0])
#Compute the control points
xcontrols = []
ycontrols = []
##Compute normal vectors as well
normalxs = []
normalys = []
##Compute PHI_i as well and beta
phis = []
betas = []
theta = []
for i in range(0,len(xpanels)-1):
    xi = (xpanels[i+1]+xpanels[i])/2.0
    yi = (ypanels[i+1]+ypanels[i])/2.0
    dx = (xpanels[i+1]-xpanels[i])
    dy = (ypanels[i+1]-ypanels[i])
    xcontrols.append(xi)
    ycontrols.append(yi)
    #Compute the normal vector component by taking the unit vector and adding pi/2
    r = np.asarray([dx,dy])
    rhat = r/np.linalg.norm(r)
    A = np.asarray([[np.cos(np.pi/2),-np.sin(np.pi/2)],[np.sin(np.pi/2),np.cos(np.pi/2)]])
    nhat = np.matmul(A,rhat)
    beta = np.arctan2(nhat[1],nhat[0])
    phi = beta - np.pi/2.0
    betas.append(beta)
    theta.append(beta)
    phis.append(phi)
    normalxs.append(nhat[0])
    normalys.append(nhat[1])

##Just for plotting purposes
xcontrols = np.asarray(xcontrols)
ycontrols = np.asarray(ycontrols)

#Once we have all our matrices we now need a double for loop to compute each control point and each panel
I = np.zeros((Npanels,Npanels))
Is = np.zeros((Npanels,Npanels))
Vs = np.zeros(Npanels)
RHS = np.zeros(Npanels)
for idx in range(0,Npanels):
    ##The control point is
    xi = xcontrols[idx]
    yi = ycontrols[idx]
    phii = phis[idx]
    ##Populate the right hand side
    betai = betas[idx]
    RHS[idx] = -Vinf*np.cos(betai)
    Vs[idx] = Vinf*np.sin(betai)
    for jdx in range(0,Npanels):
        if idx == jdx:
            #Remember the contribution to itself is just 0.5
            Iij = np.pi ##Normal velocity
            Iijs = 0.0 ##Tangential velocity
        else:
            #Get the other variables
            Xj = xpanels[jdx]
            Yj = ypanels[jdx]
            #Get the next panel
            if jdx == Npanels-1:
                Xj1 = xpanels[0]
                Yj1 = ypanels[0]
            else:
                Xj1 = xpanels[jdx+1]
                Yj1 = ypanels[jdx+1]
            phij = phis[jdx]
            #Compute the wild Integral crap
            A = -(xi - Xj)*np.cos(phij) - (yi - Yj)*np.sin(phij)
            B = (xi-Xj)**2 + (yi-Yj)**2
            C = np.sin(phii-phij)
            D = (yi - Yj)*np.cos(phii) - (xi-Xj)*np.sin(phii)
            Sj = np.sqrt((Xj1-Xj)**2 + (Yj1-Yj)**2)
            E = np.sqrt(B-A**2)
            Iij = C/2.0 * np.log((Sj**2 + 2*A*Sj + B)/B) + (D-A*C)/E * (np.arctan((Sj+A)/E) - np.arctan(A/E))
            Iijs = (D-A*C)/(2*E) * np.log((Sj**2 + 2*A*Sj + B)/B) - C * (np.arctan((Sj+A)/E) - np.arctan(A/E))
        ##Populate the big NxN matrix
        I[idx][jdx] = Iij/(2.0*np.pi)
        Is[idx][jdx] = Iijs/(2.0*np.pi)
        #and print
        #print(idx,jdx,Iij)
        
##Now we use the linear algebra toolbox to solve for the circulations
circs = np.linalg.solve(I,RHS)

#Once we have the circulations we can compute the Velocity at each panel
Vi = np.matmul(Is,circs) + Vs

##Then of course the pressure coefficient
Cp = 1 - (Vi/Vinf)**2
    
#Let's print the big matrix,RHS and circs but only if Npanels < 10:
if Npanels < 10:
    print('I=',I)
    print('Is=',Is)
    print('RHS=',RHS)
    print('Circs=',circs)
    print('Vs = ',Vs)
    print('Vi = ',Vi)
    print('Cp = ',Cp)
    
##Let's make sure our solution is valid
print('This number (',np.sum(circs),') should be really close to zero')

##Plotting
plt.figure()
##Plot the body
plt.plot(xsmooth,yu_smooth,'b-',label='Upper Surface')
plt.plot(xsmooth,yl_smooth,'r-',label='Lower Surface')
##Plot the panel points
plt.plot(xpanels,ypanels,'kx')
plt.plot(xpanels,ypanels,'k-',label='Panels')
##Plot the control points
plt.plot(xcontrols,ycontrols,'ro',label='Control Points')
##Plot the normal vectors
for i in range(0,len(normalxs)):
    xstart = xcontrols[i]
    xend = xcontrols[i] + 0.5*R*normalxs[i]
    ystart = ycontrols[i]
    yend = ycontrols[i] + 0.5*R*normalys[i]
    nvecx = np.asarray([xstart,xend])
    nvecy = np.asarray([ystart,yend])
    plt.plot(nvecx,nvecy,'g-')
plt.grid()
plt.legend()
plt.axis('equal')


##Let's plot pressure as a function of beta? No beta
plt.figure()
plt.plot(theta,Cp,'r*',label='Numeric')
##Get analytic solution as well?
theta_A = np.linspace(-np.pi,np.pi,100)
Cp_A = 1.0-4.0*np.sin(theta_A)**2
plt.plot(theta_A,Cp_A,'b-',label='Analytic')
plt.grid()
plt.xlabel('Theta (rad)')
plt.ylabel('Cp')

plt.show()