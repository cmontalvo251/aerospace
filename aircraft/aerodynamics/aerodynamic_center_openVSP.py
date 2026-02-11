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
Cl = 0.71655
Cm = 0.1419366

#Using Statics
xacW_VSP = xcg - cbar*Cm/Cl
print('Using openVSP aero center of the main wing = ',xacW_VSP)