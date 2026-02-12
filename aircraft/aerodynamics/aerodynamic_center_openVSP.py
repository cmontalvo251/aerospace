import numpy as np

##Rectangular Wing Example
S = 500.0 #in^2
cbar = 10.0 #in
xcg = 4.155 #in

##Theoretical xacW
xacW = cbar/4
print('Using c/4 the aero center of the main wing = ',xacW)

#Using OpenVSP you need Alfa, Cm and Cl
#Equation below is sensitive so be sure to export data to CSV
alfa = 10.0
Cl = 0.74156
Cm = 0.13410

#Using Statics
xacW_VSP = xcg - cbar*Cm/Cl
print('Using openVSP aero center of the main wing = ',xacW_VSP)

##Adding in a Tail moves some things around
xcg = 6.433 #in
St = 30.0 #in^2
xac_t = 25 - 2.34 + 2.34/2.0
#Nose to nose the tails are 25 inches apart.
#The xac_w is 2.34" and the root chord of the tail is 5"
#so I'm assuming the aero center is half of the main wing aero center
print('Distance between aero centers of wings = ',xac_t)

#Theoretical aerodynamic center is
xMAC = xacW + (St/S)*xac_t  ###<---this equation is wrong
print('Theoretical location of xMAC = ',xMAC)

#Running VSP you get the following data for Cl and drag
Cl = 0.8323
Cm = 0.15779
xMAC_VSP = xcg - cbar*Cm/Cl ###I think this equation is correct.
print('xMAC using VSP = ',xMAC_VSP)

##At a cg location 6.433 the system is unstable which means the xMAC is less than 6.433
##At a cg location 4.155 the system is stable which means the xMAC is greater than 4.155

##So it's more than likely 4.53 which is in between (4.155, 6.433)
##If xMAC was at 3.9 the syste would still be unstable at a cg of 4.155
