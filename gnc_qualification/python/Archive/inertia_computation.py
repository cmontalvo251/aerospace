##THis code will compute the inertia of ta CubeSAT with Solar Panels
import numpy as np
import matplotlib.pyplot as plt

#Mass - 25 kg - This includes the mass of solar panels
total_mass = 25.0
#Mass of Solar panel = 1.3 kg (each panel is 0.65 kg)
mass_solar_panel = 0.65
number_of_panels = 2.0
mass_cuboid = total_mass - mass_solar_panel*number_of_panels

##Size of Cuboid
Length = 366.0/1000.0 #meters
Width = 226.0/1000.0 #meters
Depth = 226.0/1000.0 #meters

##Size of Solar Panels
Length_sp = 3.0/1000.0
Width_sp = 296.90/1000.0
Depth_sp = 543.90/1000.0

###Distance to Solar Panel Centroid from Satellite Centroid
#Center of mass of the solar panel = lx = 100 mm, 113 mm
lx = 0.0
ly = 261.45/1000.0 #meters
lz = -183.0/1000.0 #meters

##Inertias about centroid
Inertia_cuboid = mass_cuboid/12.0 * np.array([[Length**2 + Width**2,0,0],[0,Length**2+Depth**2,0],[0,0,Width**2+Depth**2]])
Inertia_sp_centroid = mass_solar_panel/12.0 * np.array([[Length_sp**2 + Width_sp**2,0,0],[0,Length_sp**2+Depth_sp**2,0],[0,0,Width_sp**2+Depth_sp**2]])
##Use the parallel axis theorem and move the sp to the centroid of the cuboid
rsp1 = np.array([[0,-lz,ly],[lz,0,-lx],[-ly,lx,0]])
rspT1 = np.transpose(rsp1)
Inertia_sp1 = Inertia_sp_centroid + mass_solar_panel*np.matmul(rsp1,rspT1)
rsp2 = np.array([[0,-lz,-ly],[lz,0,-lx],[ly,lx,0]])
rspT2 = np.transpose(rsp2)
Inertia_sp2 = Inertia_sp_centroid + mass_solar_panel*np.matmul(rsp2,rspT2)

print('Inertia of Cuboid (kg-m^2) = ',Inertia_cuboid)
print('Inertia of Solar Panel about Centroid (kg-m^2) = ',Inertia_sp_centroid)
print('Inertia of Solar Panel about Cuboid Center (1) (kg-m^2) = ',Inertia_sp1)
print('Inertia of Solar Panel about Cuboid Center (2) (kg-m^2) = ',Inertia_sp2)

##Compute the inertia of total Satellite
Inertia = Inertia_cuboid + Inertia_sp1 + Inertia_sp2
print('Total Inertia (kg-m^2) = ',Inertia)

##
Ixx = Inertia[0][0]
Iyy = Inertia[1][1]
Izz = Inertia[2][2]
tip_off_rate = np.linspace(0,10,100)
FS = 2.0
Hxx = Ixx*tip_off_rate*np.pi/180.0*FS
Hyy = Iyy*tip_off_rate*np.pi/180.0*FS
Hzz = Izz*tip_off_rate*np.pi/180.0*FS

## kg-m^2 * rad/s -> N-m-s -> kg m/s^2 -m -s -> kg-m^2/s

plt.plot(tip_off_rate,Hxx*1000.0,label='X Axis')
plt.plot(tip_off_rate,Hyy*1000.0,label='Y Axis')
plt.plot(tip_off_rate,Hzz*1000.0,label='Z Axis')
plt.xlabel('Tip Off Rate (deg/s)')
plt.ylabel('Inertia Storage (mN-m-s)')
plt.title('Safety Factor = '+str(FS))
plt.legend()
plt.grid()
plt.show()
