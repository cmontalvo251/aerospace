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
CLw = 0.74156
CMw = 0.13410

#Using Statics
xacW_VSP = xcg - cbar*CMw/CLw
print('Using openVSP aero center of the main wing = ',xacW_VSP)

##Adding in a Tail moves the cg
xcg = 6.816 #in

##Tail parameters
ST = 100.0 #in^2
cT = 5.0 #in
xacT = 50 + xacW_VSP*(cT/cbar-1)
#Nose to nose the tails are 50 inches apart.
#The xac_w is 2.34" and the root chord of the tail is 5"
#so I'm assuming the aero center is half of the main wing aero center
print('Distance between aero centers of wings = ',xacT)

#Theoretical aerodynamic center is Complicated because it's highly sensitive to the lift on the tail
#Thus it's recommended you run a stability sweep on the tail independently of the main wing
CLT = 0.66571*0.8 #add in "wake" effects
CL = (CLw + ST/S*CLT) #This is essentially the lift of the entire airplane
print('Theoretical Total Lift = ',CL)
VH = ST*xacT/(S*cbar) #This is the volume coefficient of the horizontal stabilizer
xMAC = xacW_VSP + CLT/CL*VH*cbar
print('Theoretical location of xMAC ignoring wake effects = ',xMAC)

#Running VSP you get the following data for Cl and drag
CL = 0.85513 #Notice this is lower tan the theoretical lift due to wake effects
CM = -0.06617 #moving the tail back far enough actually creates stability
xMAC_VSP = xcg - cbar*CM/CL
print('xMAC using VSP = ',xMAC_VSP)

##At a cg location 6.433 the system is unstable which means the xMAC is less than 6.433
##At a cg location 4.155 the system is stable which means the xMAC is greater than 4.155

##So it's more than likely 4.53 which is in between (4.155, 6.433)
