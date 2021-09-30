#!/usr/bin/python

##In order to import this toolbox into a python script you need to 
##do the following. Copy the following lines of code below
# import sys
# sys.path.append('/home/carlos/Dropbox/BlackBox/pdf')
# from pdf import *

# or

# In order to get python to search for all of your lovely blackbox 
# python routines. Add this to your .bashrc file

# for d in /home/carlos/Dropbox/BlackBox/*/; do
# 	PYTHONPATH+=:$d
# done
# export PYTHONPATH

# In order to get this to work in Thonny you need to navigate to
# /home/user/.thonny/BundledPython36/lib/python3.6/site-packages and place
# a symbolic link here

# In Enthough you need to make symbolic links here
# /home/carlos/.local/share/canopy/edm/envs/User/lib/python2.7/site-packages

import matplotlib.pyplot as plt
import datetime
import csv
from matplotlib.dates import DateFormatter
import numpy as np
from pdf import *
import plotting as P

def create_text_file(quad_data,fileoutname):
    #Now we need to output all the data to a text file for MATLAB use
    # GPS+Arduino time in hours - 1
    # Latitude - 2
    # Longitude - 3
    # Altitude (for quad output the autopilot altitude, for pitot output the the pressure altitude offset by GPS?) - 4 - Ya. That will have to do
    # Phi - for pitot just output all zeros for this - 5
    # Theta - 6
    # Psi - 7
    # Lateral speed - for pitot just output GPS.speed - 8
    # Lateral speed- 9 - crap. GPS.speed got commented out because of space issues
    # Vertical speed - fp4 all -99 - 10
    # P - pitot all -99 - 11
    # Q - 12
    # R - 13
    # Ax - pitot all -99 - 14
    # Ay - 15
    # Az - 16
    # A0 - for quad all -99 - 17
    # A1 - 18
    # A2 - raw bits for pitot probe so we can use Matlab to calibrate - 19
    # A3 - 20
    # Temperature - all -99 for quad - 21
    # Pressure - 22 - except for this one
    # Humidity - 23

    baroRead = quad_data[0]
    gpsRead = quad_data[1]
    IMURead = quad_data[3]
    ATTRead = quad_data[4]

    if gpsRead.GPSINIT == 0:
        print('Logfile not written because GPS not initialized')
        return

    #open file for writing
    csvwriter = csv.writer(open(fileoutname,"wb"),delimiter=",")

    #Ok so here's the problem. All of these sensors update at different rates.
    #GPS at 4 Hz or so so for a 1600 msond flight the GPS has 7530 data points or about 4 Hz
    #The Attitude sensor on the other hand is operating at like 10 Hz so the Attittude sensor has 14800 data points
    #Finally the IMU is absurdly fast and is operating at like 40 Hz so it has like 74000 data points. So we need to interpolate
    #these things like so
    #ALmost forgot about the barometer. It looks like the barometer updates around 10 Hz since it has 14800 data points
    print("Length of Data Streams = ",len(gpsRead.time_quad),len(ATTRead.ROLL),len(IMURead.AccX),len(baroRead.pressurekPA))
    print("Durations of Data Streams = ",gpsRead.timeSec[-1]-gpsRead.timeSec[0],baroRead.timeSec[-1]-baroRead.timeSec[0],ATTRead.timeMS_ATT[-1]-ATTRead.timeMS_ATT[0],IMURead.timeMS_IMU[-1]-IMURead.timeMS_IMU[0],baroRead.timeSec[-1]-baroRead.timeSec[0])

    #Interpolate ROLL, PITCH and YAW !!!!!!!!
    #TIME IS ALL FUCKED on QUAD NEED TO INTERPOLATE ON ALL OF THESE HERE AN EXAMPLE
    #Major problem - This makes the code take FORRRRRRRRRRRRRREEEEEVVVVVVEERRRR to write to a .out file
    #Not sure what to do here. Fuck. Ok it looks like if you interpolate all at once it's way faster. Crisis averted?
    # roll = np.interp(gpsRead.timeSec[idx],ATTRead.timeMS_ATT,ATTRead.ROLL)*np.pi/180.0
    # pitch = np.interp(gpsRead.timeSec[idx],ATTRead.timeMS_ATT,ATTRead.PITCH)*np.pi/180.0
    # yaw = np.interp(gpsRead.timeSec[idx],ATTRead.timeMS_ATT,ATTRead.YAW)*np.pi/180.0
    print("Interpolating Roll,pitch,yaw")
    print(gpsRead.timeSec[0],gpsRead.timeSec[-1]) #This is 91.928 and 183.758???
    roll_interp = np.interp(gpsRead.timeSec,ATTRead.timeMS_ATT,ATTRead.ROLL)*np.pi/180.0
    pitch_interp = np.interp(gpsRead.timeSec,ATTRead.timeMS_ATT,ATTRead.PITCH)*np.pi/180.0
    yaw_interp = np.interp(gpsRead.timeSec,ATTRead.timeMS_ATT,ATTRead.YAW)*np.pi/180.0
    print("Interpolating IMU")
    gyrX_interp = np.interp(gpsRead.timeSec,IMURead.timeMS_IMU,IMURead.GyrX)*np.pi/180.0
    gyrY_interp = np.interp(gpsRead.timeSec,IMURead.timeMS_IMU,IMURead.GyrY)*np.pi/180.0
    gyrZ_interp = np.interp(gpsRead.timeSec,IMURead.timeMS_IMU,IMURead.GyrZ)*np.pi/180.0
    accX_interp = np.interp(gpsRead.timeSec,IMURead.timeMS_IMU,IMURead.AccX)
    accY_interp = np.interp(gpsRead.timeSec,IMURead.timeMS_IMU,IMURead.AccY)
    accZ_interp = np.interp(gpsRead.timeSec,IMURead.timeMS_IMU,IMURead.AccZ)
    print("Interpolating Barometer")
    pressure_interp = np.interp(gpsRead.timeSec,baroRead.timeSec,baroRead.pressurekPA)
    print("Length of New data Streams",len(gpsRead.time_quad),len(roll_interp),len(gyrX_interp),len(pressure_interp))
    print("interpolation done....")
    
    for idx in range(0,len(gpsRead.timeSec)):
        # GPS Data
        # tot_time_hr - 2
        #print "Writing Row = ",idx,' out of ',len(gpsRead.timeSec)
                
        row = [#baroRead.timeHr[idx],
            #I think the barometer timer is kind of off so I'm going to use the gps timer
            gpsRead.time_quad[idx],
        # Latitude - 0
               gpsRead.latitude[idx],
        # Longitude - 1
               gpsRead.longitude[idx],
        # x_vec_np - 3 
        # y_vec_np - 4 
        # alt_vec_np - 5
               gpsRead.altitude_CTUN[idx]+gpsRead.altitude[0], ##Altitude from CTUN (AUTOPILOT) offset by GPS altitude so we get MSL and not AGL
        # Note that GPS Start = data[2][0] and GPS End = data[2][-1]
        # Phi,theta,psi
            roll_interp[idx],
            pitch_interp[idx],
            yaw_interp[idx],
        # ,xdot,ydot,zdot
               gpsRead.speed[idx], #shit. Is this in knots? No this is in m/s 
               gpsRead.speed[idx], ##http://ardupilot.org/copter/docs/common-downloading-and-analyzing-data-logs-in-mission-planner.html
               gpsRead.CLIMB_RATE[idx], ##According to this website they are in m/s
        # p,q,r
               gyrX_interp[idx], ##http://ardupilot.org/copter/docs/common-downloading-and-analyzing-data-logs-in-mission-planner.html
               gyrY_interp[idx], ###This website says they are in deg/sec
               gyrZ_interp[idx], #So this converts them from deg/s to rad/s
        # Ax,Ay,Az
               accX_interp[idx], #m/s^2
               accY_interp[idx],
               accZ_interp[idx],
        # Pitot Data 
        # tot_time_sec_zero - 0 
        # airspeed_ms_filtered_all - 1 
        # airspeed_ms_all - 2 
        # raw_data - 3
               -99,
               -99,
               -99,
               -99,
        # print row
        # numpitot = len(airspeed_ms_all) 
        # Pitot Orientation = West,North,South,East (Qaud Ref) 
    
        # Temperature/Pressure/Humidity 
        # tot_time_sec_zero - 0
        # Temperature (C) - 1
               -99,
        # Pressure (kPa) - 2
               pressure_interp[idx],
        # Humidity (%) - 3
               -99]
        # Pressure Altitude (m) - 4
        csvwriter.writerow(row)

    print('Text File written to '+fileoutname)

def get_quad_data(fileName):

    try:
        with open(fileName) as file:
            fileExists = True
    except IOError as e:
        print("Unable to open qaudcopter file") #Does not exist OR no read permissions)
        fileExists = False

    baroRead = barometer()
    gpsRead  = gps()
    rcinoutRead = rc()
    IMURead = IMU()
    ATTRead = ATTITUDE()
    BATTRead = battery()    

    if (fileExists):
        fileObject = open(fileName, "r")
        gpsRead.GPSINIT = 0 #This means we haven't acquired GPS yet
        for line in fileObject:
            splitLine = line.split(',')
            #print splitLine
            if (splitLine[0] == "FMT") and (splitLine[3].replace(' ','') == "GPS"):
                formatGPS = splitLine[4].replace(' ','')
                if formatGPS[0] == "B":
                    #APM ver 2.0.18
                    #Status,TimeMS,Week,NSats,HDop,Lat,Lng,RelAlt,Alt,Spd,GCrs,VZ,T
                    apmVersion = '2.0.18'
                    rad2deg = 100.0 #angles are in centi-degrees???
                    ##http://ardupilot.org/copter/docs/common-downloading-and-analyzing-data-logs-in-mission-planner.html
                    ##According to this website the angles are in centi-degrees
                else:
                    #APM ver 2.0.21
                    #TimeMS,Status,GPSTimeMS,Week,NSats,HDop,Lat,Lng,RelAlt,Alt,Spd,GCrs,VZ,T
                    apmVersion = '2.0.21'
                    rad2deg = 1.0 #angles are already in degrees - But according to just my basic understanding and an experiment I did
                    #at the Irvington Airfield it looks as though the angles here are already in degrees. Odd.
                apmVersion = '2.0.26'
                print("APM VERSION = ",apmVersion,formatGPS[0])
                print("Resetting GPSINIT")
                gpsRead.GPSINIT = 0
                
            if (splitLine[0] == "BARO"):
                timeMS = int(splitLine[1])
                # print('timeMSBARO = ',timeMS)
                if gpsRead.GPSON_Sec == -99:
                    gpsRead.GPSON_Sec = timeMS/1000.0
                    print('GPS Time = ',gpsRead.GPSON_Sec)
                #print timeMS
                if len(baroRead.timeMS) > 0:
                    if timeMS + baroRead.timeOffset < baroRead.timeMS[-1]:
                        baroRead.timeOffset = baroRead.timeMS[-1]
                baroRead.timeMS.append(timeMS + baroRead.timeOffset)
                if len(gpsRead.timeMS_GPS) == 0:
                    baroRead.timeMS_GPS.append(-99)
                else:
                    baroRead.timeMS_GPS.append(gpsRead.timeMS_GPS[-1])
                baroRead.altitude.append(float(splitLine[2]))
                baroRead.pressure.append(float(splitLine[3]))
                baroRead.temperature.append(float(splitLine[4]))
            elif (splitLine[0] == "CTUN"):
                timeMS = int(splitLine[1])/1000.0
                # print('timeMSCTUN = ',timeMS)
                if len(gpsRead.timeMS_CTUN) > 0:
                    if timeMS + gpsRead.CTUN_OFFSET < gpsRead.timeMS_CTUN[-1]:
                        gpsRead.CTUN_OFFSET = gpsRead.timeMS_CTUN[-1]
                gpsRead.timeMS_CTUN.append(timeMS + gpsRead.CTUN_OFFSET)
                gpsRead.D_altitude.append(np.float(splitLine[5]))
                gpsRead.altitude_CTUN.append(np.float(splitLine[6]))
                gpsRead.baro_altitude.append(np.float(splitLine[7]))
                gpsRead.DS_altitude.append(np.float(splitLine[8]))
                gpsRead.SON_altitude.append(np.float(splitLine[9]))
                gpsRead.CLIMB_RATE_DES.append(np.float(splitLine[10])/100.0)
                gpsRead.CLIMB_RATE.append(np.float(splitLine[11])/100.0)
            elif (splitLine[0] == "GPS"):
                #Parse File
                #APM ver 2.0.18
                #GPS,Status,TimeMS,Week,NSats,HDop,Lat,Lng,RelAlt,Alt,Spd,GCrs,VZ,T
                #APM ver 2.0.21
                #GPS,TimeMS,Status,GPSTimeMS,Week,NSats,HDop,Lat,Lng,RelAlt,Alt,Spd,GCrs,VZ,T
                #APM ver 2.0.26
                #GPS,QBIHBcLLefffB,TimeUS,Status,GMS,GWk,NSats,HDop,Lat,Lng,Alt,Spd,GCrs,VZ,U
                #GPS,TimeUS,Status,GMS,GWk,NSats,HDop,Lat,Lng,Alt,Spd,GCrs,VZ,U
                # 0,  1      2       3 4    5     6    7  8   9   10  11   12
                if apmVersion == '2.0.18':
                    timeMS_internal = int(splitLine[1])
                    timeMS = int(splitLine[2])
                    week = int(splitLine[3])
                    numSats = int(splitLine[4])
                    #dummy = int(splitLine[5])
                    lat = float(splitLine[6])
                    lon = float(splitLine[7])
                    altitude = float(splitLine[9])*100
                    rel_alt = float(splitLine[8])*100
                    speed = float(splitLine[10])
                    VZ = -np.float(splitLine[12])
                elif apmVersion == '2.0.21':
                    timeMS_internal = int(splitLine[1])
                    timeMS = int(splitLine[3]) #GPSTimeMS is shifted here. ugh.
                    week = int(splitLine[4])
                    numSats = int(splitLine[5])
                    #dummy = int(splitLine[6])
                    lat = float(splitLine[7])
                    lon = float(splitLine[8])
                    rel_alt = float(splitLine[9])
                    altitude = float(splitLine[10])
                    speed = float(splitLine[11])
                    VZ = -np.float(splitLine[13])
                elif apmVersion == '2.0.26':
                    timeMS_internal = float(splitLine[1])
                    status = splitLine[2]
                    timeMS = int(float(splitLine[3])) #GPSTimeUS is shifted here. ugh.
                    week = int(splitLine[4])
                    numSats = int(splitLine[5])
                    #dummy = int(splitLine[6])
                    lat = float(splitLine[7])
                    lon = float(splitLine[8])
                    altitude = float(splitLine[9])
                    rel_alt = altitude
                    speed = float(splitLine[10])
                    gcrs = splitLine[11]
                    VZ = -np.float(splitLine[12])

                #It's possible that we never get GPS. If that's the case we need to
                #make a note of it.
                    
                #Need to run a few checks and make sure that we have a valid GPS lock                
                if gpsRead.GPSINIT == 1: #This means we have initialized GPS
                    if len(gpsRead.timeMS_GPS) > 0:
                        if timeMS + gpsRead.GPS_OFFSET < gpsRead.timeMS_GPS[-1]:
                            print("TimeMS = ",timeMS)
                            print("Old OFFSET = ",gpsRead.GPS_OFFSET)
                            gpsRead.GPS_OFFSET = gpsRead.timeMS_GPS[-1]
                            print("gpsRead.timeMS_GPS[-1] = ",gpsRead.timeMS_GPS[-1])
                            print("New Offset = ",gpsRead.GPS_OFFSET)
                    if len(gpsRead.timeMS_internal) > 0:
                        if (timeMS_internal + gpsRead.internal_OFFSET < gpsRead.timeMS_internal[-1]):
                            gpsRead.internal_OFFSET = gpsRead.timeMS_internal[-1]
                    gpsRead.timeMS_internal.append(timeMS_internal+gpsRead.internal_OFFSET)
                    gpsRead.timeMS_GPS.append(timeMS + gpsRead.GPS_OFFSET)
                    gpsRead.week_GPS.append(week)
                    if week != 0:
                        gpsRead.WEEK = week
                    gpsRead.numSats.append(numSats)
                    gpsRead.latitude.append(lat)
                    gpsRead.longitude.append(lon)
                    #Ok here we have a problem if altitude goes negative
                    if altitude < 0:
                        altitude = 0
                    gpsRead.altitude.append(altitude)
                    gpsRead.rel_altitude.append(rel_alt)
                    gpsRead.speed.append(speed)
                    gpsRead.VZ.append(VZ)
                else:
                    #print lat
                    if abs(lat) > 28 and abs(lat) < 32:
                        #print lon
                        if abs(lon) > 86 and abs(lon) < 90:
                            #print altitude
                            if abs(altitude) < 1000:
                                #Make sure MS time is initialized as well
                                if timeMS > 0:
                                    gpsRead.GPSINIT = 1
                                    #Grab when GPS turned on in MS
                                    #It's possible that we have GPS as soon as the system turned on.
                                    #Insane
                                    if gpsRead.GPSON_Sec == 0:
                                        if len(baroRead.timeMS) > 0:
                                            gpsRead.GPSON_Sec = np.float(baroRead.timeMS[-1])/1000.0
                                            print('GPS Time = ',gpsRead.GPSON_Sec)
                                        else:
                                            gpsRead.GPSON_Sec = -99 #Is this ok? or should we wait? yea we need to wait
                                    print('GPS Acquired = ',gpsRead.GPSINIT)
            elif (splitLine[0] == "RCIN"):
                # Channel 1 - Roll
                # Channel 2 - Pitch
                # Channel 3 - Throttle
                # Chanlle 4 - Yaw
                timeMS = int(splitLine[1])
                # print('timeMSRCIN = ',timeMS)
                if len(rcinoutRead.timeMSIN) > 0:
                    if timeMS + rcinoutRead.IN_OFFSET < rcinoutRead.timeMSIN[-1]:
                        rcinoutRead.IN_OFFSET = rcinoutRead.timeMSIN[-1]
                val1 = int(splitLine[1])+rcinoutRead.IN_OFFSET
                rcinoutRead.timeMSIN.append(val1)
                val2 = int(splitLine[2])
                rcinoutRead.RollChannel.append(val2)
                val3 = int(splitLine[3])
                rcinoutRead.PitchChannel.append(val3)
                val4 = int(splitLine[4])
                rcinoutRead.ThrottleChannel.append(val4)
                val5 = int(splitLine[5])
                rcinoutRead.YawChannel.append(val5)
            elif (splitLine[0] == "RCOU"):
                timeMS = int(splitLine[1])
                if len(rcinoutRead.timeMSOUT) > 0:
                    if timeMS + rcinoutRead.OUT_OFFSET < rcinoutRead.timeMSOUT[-1]:
                        rcinoutRead.OUT_OFFSET = rcinoutRead.timeMSOUT[-1]
                rcinoutRead.timeMSOUT.append(int(splitLine[1])+rcinoutRead.OUT_OFFSET)
                rcinoutRead.TopRight.append(int(splitLine[2])) #Motor 1 - Top Right - This configuration comes from 
                rcinoutRead.BottomLeft.append(int(splitLine[3])) #Motor 2 - Bottom Left - the pixhawk Quad-X profile
                rcinoutRead.TopLeft.append(int(splitLine[4])) #Motor 3 - Top Left - But This could be wrong
                rcinoutRead.BottomRight.append(int(splitLine[5])) #Motor 4 - Bottom Right
            elif (splitLine[0] == "IMU"):
                timeMS = np.float(splitLine[1])/1000.0
                if len(IMURead.timeMS_IMU) > 0:
                    if timeMS + IMURead.IMU_OFFSET < IMURead.timeMS_IMU[-1]:
                        IMURead.IMU_OFFSET = IMURead.timeMS_IMU[-1]
                timeMS_IMU = int(splitLine[1])/1000.0+IMURead.IMU_OFFSET
                # print('timeMS_IMU = ',timeMS_IMU)
                IMURead.timeMS_IMU.append(timeMS_IMU)
                ##http://ardupilot.org/copter/docs/common-downloading-and-analyzing-data-logs-in-mission-planner.html
                ##According to this website the values here are in deg/sec
                IMURead.GyrX.append(np.float(splitLine[2]))
                IMURead.GyrY.append(np.float(splitLine[3])) 
                IMURead.GyrZ.append(np.float(splitLine[4]))
                IMURead.AccX.append(np.float(splitLine[5]))
                IMURead.AccY.append(np.float(splitLine[6]))
                IMURead.AccZ.append(np.float(splitLine[7])+9.81) #This is how you know Acc is in m/s^2
            elif (splitLine[0] == "IMU2"):
                timeMS = np.float(splitLine[1])/1000.0
                if len(IMURead.timeMS_IMU2) > 0:
                    if timeMS + IMURead.IMU2_OFFSET < IMURead.timeMS_IMU2[-1]:
                        IMURead.IMU2_OFFSET = IMURead.timeMS_IMU2[-1]
                IMURead.timeMS_IMU2.append(np.float(splitLine[1])/1000.0+IMURead.IMU2_OFFSET)
                IMURead.GyrX2.append(np.float(splitLine[2]))
                IMURead.GyrY2.append(np.float(splitLine[3]))
                IMURead.GyrZ2.append(np.float(splitLine[4]))
                IMURead.AccX2.append(np.float(splitLine[5]))
                IMURead.AccY2.append(np.float(splitLine[6]))
                IMURead.AccZ2.append(np.float(splitLine[7])+9.81)
            elif (splitLine[0] == "EKF2"):
                timeMS = np.float(splitLine[1])/1000.0
                if len(IMURead.timeMS_EKF2) > 0:
                    if timeMS + IMURead.EKF2_OFFSET < IMURead.timeMS_EKF2[-1]:
                        IMURead.EKF2_OFFSET = IMURead.timeMS_EKF2[-1]
                IMURead.timeMS_EKF2.append(np.float(splitLine[1])/1000.0+IMURead.EKF2_OFFSET)
                IMURead.AX.append(np.float(splitLine[2]))
                IMURead.AY.append(np.float(splitLine[3]))
                IMURead.AZ.append(np.float(splitLine[4]))
                IMURead.MX.append(np.float(splitLine[5]))
                IMURead.MY.append(np.float(splitLine[6]))
                IMURead.MZ.append(np.float(splitLine[7]))
            elif (splitLine[0] == "AHR2"):
                timeMS = np.float(splitLine[1])/1000.0
                if len(ATTRead.timeMS_AHR2) > 0:
                    if timeMS + ATTRead.AHR2_OFFSET < ATTRead.timeMS_AHR2[-1]:
                        ATTRead.AHR2_OFFSET = ATTRead.timeMS_AHR2[-1]
                ATTRead.timeMS_AHR2.append(np.float(splitLine[1])/1000.0+ATTRead.AHR2_OFFSET)
                ATTRead.ROLL2.append(np.float(splitLine[2])*rad2deg)
                ATTRead.PITCH2.append(np.float(splitLine[3])*rad2deg)
                ATTRead.YAW2.append(np.float(splitLine[4])*rad2deg)
                ATTRead.altitude.append(np.float(splitLine[5]))
            elif (splitLine[0] == "EKF1"):
                timeMS = np.float(splitLine[1])/1000.0
                if len(ATTRead.timeMS_EKF1) > 0:
                    if timeMS + ATTRead.EKF1_OFFSET < ATTRead.timeMS_EKF1[-1]:
                        ATTRead.EKF1_OFFSET = ATTRead.timeMS_EKF1[-1]
                ATTRead.timeMS_EKF1.append(np.float(splitLine[1])/1000.0+ATTRead.EKF1_OFFSET)
                ATTRead.ROLL1.append(np.float(splitLine[2])*rad2deg)
                ATTRead.PITCH1.append(np.float(splitLine[3])*rad2deg)
                ATTRead.YAW1.append(np.float(splitLine[4])*rad2deg)
            elif (splitLine[0] == "ATT"):
                timeMS = np.float(splitLine[1])/1000.0
                if len(ATTRead.timeMS_ATT) > 0:
                    if timeMS + ATTRead.ATT_OFFSET < ATTRead.timeMS_ATT[-1]:
                        ATTRead.ATT_OFFSET = ATTRead.timeMS_ATT[-1]
                ATTRead.timeMS_ATT.append(np.float(splitLine[1])/1000.0+ATTRead.ATT_OFFSET)
                roll_command = np.float(splitLine[2])*rad2deg
                roll = np.float(splitLine[3])*rad2deg
                pitch_command = np.float(splitLine[4])*rad2deg
                pitch = np.float(splitLine[5])*rad2deg
                yaw_command = np.float(splitLine[6])*rad2deg
                yaw = np.float(splitLine[7])*rad2deg
                    
                ATTRead.ROLL_DES.append(roll_command)
                ATTRead.ROLL.append(roll)
                ATTRead.PITCH_DES.append(pitch_command)
                ATTRead.PITCH.append(pitch)
                ATTRead.YAW_DES.append(yaw_command)
                ATTRead.YAW.append(yaw)
            elif (splitLine[0] == "CURR"):
                timeMS = np.float(splitLine[1])/1000.0
                if len(BATTRead.timeSec) > 0:
                    if timeMS + BATTRead.OFFSET < BATTRead.timeSec[-1]:
                        BATTRead.OFFSET = BATTRead.timeSec[-1]
                BATTRead.timeSec.append(np.float(splitLine[1])/1000.0+BATTRead.OFFSET)
                BATTRead.Voltage.append(np.float(splitLine[4])/100.0)
                BATTRead.Current.append(np.float(splitLine[5])/100.0) 

    else:
        print("Sorry, Quadcopter Log File: " + fileName + " Does not exist")
        return 0

    if gpsRead.GPSINIT == 1:
        gpsRead.calibrate_time()
        gpsRead.fix_gps()
    baroRead.establish_units(gpsRead)
    rcinoutRead.establish_units()
    IMURead.convert2np()
    ATTRead.convert2np()   

    #Do a check on Yaw. If you encounter a yaw angle greater than 360 degrees something is wrong.
    yaw_max = np.max(ATTRead.YAW)
    if yaw_max > 360:
        print("Something looks goofed in the yaw axis. Check the rad2deg value.")
        sys.exit()
    
    directory = "Quad Data Directory \n" + "baro Class - 0 \n" + "gps Class - 1 \n" + "RCIN/OUT - 2 \n" + "IMU - 3 \n" + "ATTITUDE - 4 \n" + "Battery - 5 \n"

    return [baroRead,gpsRead,rcinoutRead,IMURead,ATTRead,BATTRead,directory]


#********************************************************************************
# Functions for time conversions...
#********************************************************************************
#-----------------------------------------------------
# * Return Modified Julian Day given calendar year,
# * month (1-12), and day (1-31).
# * - Valid for Gregorian dates from 17-Nov-1858.
# * - Adapted from sci.astro FAQ.
#-----------------------------------------------------
def DateToMJD(year, month, day):
    Y = 367 * year
    YM = int(year + (month + 9) / 12)
    YM2 = int(year + (month - 9) / 7)
    C  = int(YM2 / 100 + 1)
    J = int(7 * YM / 4)
    A = int(3 * C / 4)
    B = int(275 * month / 9)
    val = Y - J - A + B + day + 1721028 - 2400000
    # print('Y = ',Y)
    # print('YM = ',YM)
    # print('YM2 = ',YM2)
    # print('C = ',C)
    # print('J = ',J)
    # print('A = ',A)
    # print('B = ',B)
    # print('DateToMJD = ',val)
    return val


#-----------------------------------------------------
# * Convert Modified Julian Day to calendar date.
# * - Assumes Gregorian calendar.
# * - Adapted from Fliegel/van Flandern ACM 11/#10 p 657 Oct 1968.
#-----------------------------------------------------
def MJDToDate(MJD):
    J = int(MJD + 2400001 + 68569)
    # print('J=',J)
    C = int(4 * J / 146097)
    # print('C=',C)
    A = int((146097 * C + 3) / 4)
    # print('A = ',A)
    J = J - A
    # print('J2 = ',J)
    Y = int(4000 * (J + 1) / 1461001)
    # print('Y=',Y)
    B = int(31 - 1461 * Y / 4)
    # print('B = ',B)
    J = J + B
    # print('J3=',J)
    M = int(80 * J / 2447)
    # print('M=',M)
    D = int(2447 * M / 80)
    day = J - D
    # print('Day=',day)
    J = int(M / 11)
    month = M + 2 - (12 * J)
    year = 100 * (C - 49) + Y + J
    
    return year, month, day

#********************************************************************************
# Classes to hold each set of data.
#********************************************************************************

#/////////////////////////////////////////////////////////
#Class for RC Channels

#////////////////////////////////////////////////////////

class rc(object):
    def __init__(self):
        self.timeMSIN = []
        self.timeMSOUT = []
        self.RollChannel = []
        self.PitchChannel = []
        self.ThrottleChannel = []
        self.YawChannel = []
        # self.RollChannelOut = []
        # self.PitchChannelOut = []
        # self.ThrottleChannelOut = []
        # self.YawChannelOut = []
        self.TopRight = []
        self.TopLeft = []
        self.BottomRight = []
        self.BottomLeft = []
        self.IN_OFFSET = 0
        self.OUT_OFFSET = 0

    def establish_units(self):
        self.timeSecIN = []
        self.timeSecOUT = []
        for i in self.timeMSIN:
            self.timeSecIN.append(i / 1000.0)
        for i in self.timeMSOUT:
            self.timeSecOUT.append(i / 1000.0)
        #After you read in channels you will run self.channels_np = np.array(self.channels)
        #print self.channels_np[0,:] #- this is the entire first data point of all channels
        #print raw_np[:,0] #- this is the first channel

    def plot_rc(self,x0,xf,pp):
        plt.figure()
        plt.plot(self.timeSecIN,self.RollChannel)
        #plt.plot(self.timeSecOUT,self.RollChannelOut,label='Roll SERVO')
        plt.xlabel('Time (ms)')
        plt.ylabel('Roll Channel (ms)')
        #plt.legend(loc=2)
        plt.grid(); plt.xlim([x0,xf])
        pp.savefig()

        plt.figure()
        plt.plot(self.timeSecIN,self.PitchChannel)
        #plt.plot(self.timeSecOUT,self.PitchChannelOut,label='Pitch SERVO')
        plt.xlabel('Time (ms)')
        plt.ylabel('Pitch Channel (ms)')
        #plt.legend(loc=2)
        plt.grid(); plt.xlim([x0,xf])
        pp.savefig()

        plt.figure()
        plt.plot(self.timeSecIN,self.YawChannel)
        #plt.plot(self.timeSecOUT,self.YawChannelOut,label='Yaw SERVO')
        plt.xlabel('Time (ms)')
        plt.ylabel('Yaw Channel (ms)')
        #plt.legend(loc=2)
        plt.grid(); plt.xlim([x0,xf])
        pp.savefig()
        
        plt.figure()
        plt.plot(self.timeSecIN,self.ThrottleChannel)
        #plt.plot(self.timeSecOUT,self.ThrottleChannelOut,label='Throttle SERVO')
        plt.xlabel('Time (ms)')
        plt.ylabel('Throttle Channel (ms)')
        #plt.legend(loc=2)
        plt.grid(); plt.xlim([x0,xf])
        pp.savefig()

        plt.figure()
        plt.plot(self.timeSecOUT,self.TopRight,label='Top Right')
        plt.plot(self.timeSecOUT,self.TopLeft,label='Top Left')
        plt.plot(self.timeSecOUT,self.BottomRight,label='Bottom Right')
        plt.plot(self.timeSecOUT,self.BottomLeft,label='Bottom Left')
        plt.xlabel('Time (ms)')
        plt.ylabel('Motor Outputs (ms)')
        plt.legend(loc=2)
        plt.grid(); plt.xlim([x0,xf])
        pp.savefig()

        plt.figure()
        # plt.plot(self.timeSecOUT,self.TopRight,label='Top Right')
        plt.plot(self.timeSecOUT,self.TopLeft,label='Top')
        # plt.plot(self.timeSecOUT,self.BottomRight,label='Bottom Right')
        plt.plot(self.timeSecOUT,self.BottomLeft,label='Bottom')
        plt.xlabel('Time (ms)')
        plt.ylabel('Left Motor Outputs (ms)')
        plt.legend(loc=2)
        plt.grid(); plt.xlim([x0,xf])
        pp.savefig()

        plt.figure()
        plt.plot(self.timeSecOUT,self.TopRight,label='Top')
        #plt.plot(self.timeSecOUT,self.TopLeft,label='Top')
        plt.plot(self.timeSecOUT,self.BottomRight,label='Bottom')
        #plt.plot(self.timeSecOUT,self.BottomLeft,label='Bottom')
        plt.xlabel('Time (ms)')
        plt.ylabel('Right Motor Outputs (ms)')
        plt.legend(loc=2)
        plt.grid(); plt.xlim([x0,xf])
        pp.savefig()

        plt.figure()
        plt.plot(self.timeSecOUT,self.TopRight,label='Right')
        plt.plot(self.timeSecOUT,self.TopLeft,label='Left')
        #plt.plot(self.timeSecOUT,self.BottomRight,label='Bottom Right')
        #plt.plot(self.timeSecOUT,self.BottomLeft,label='Bottom Left')
        plt.xlabel('Time (ms)')
        plt.ylabel('Top Motor Outputs (ms)')
        plt.legend(loc=2)
        plt.grid(); plt.xlim([x0,xf])
        pp.savefig()

        plt.figure()
        #plt.plot(self.timeSecOUT,self.TopRight,label='Top Right')
        #plt.plot(self.timeSecOUT,self.TopLeft,label='Top Left')
        plt.plot(self.timeSecOUT,self.BottomRight,label='Right')
        plt.plot(self.timeSecOUT,self.BottomLeft,label='Left')
        plt.xlabel('Time (ms)')
        plt.ylabel('Bottom Motor Outputs (ms)')
        plt.legend(loc=2)
        plt.grid(); plt.xlim([x0,xf])
        pp.savefig()
        

#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# Class for ATTITUDE(self):
#
#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

class ATTITUDE(object):
    def __init__(self):
        #AHRS2
        self.timeMS_AHR2 = []
        self.ROLL2 = []
        self.PITCH2 = []
        self.YAW2 = []
        self.altitude = []
        self.AHR2_OFFSET = 0
        #EKF1
        self.timeMS_EKF1 = []
        self.ROLL1 = []
        self.PITCH1 = []
        self.YAW1 = []
        self.EKF1_OFFSET = 0
        #ATT
        self.timeMS_ATT = []
        self.ROLL = []
        self.PITCH = []
        self.YAW = []
        self.ROLL_DES = []
        self.PITCH_DES = []
        self.YAW_DES = []
        self.ATT_OFFSET = 0
        
    def convert2np(self):
        self.timeMS_ATT = np.array(self.timeMS_ATT)
        self.ROLL = np.array(self.ROLL)
        self.PITCH = np.array(self.PITCH)
        self.YAW = np.array(self.YAW)

    def plot_ATTITUDE(self,x0,xf,pp):
        plt.figure()
        print("Roll and Roll Desired = ",self.ROLL[0],self.ROLL_DES[0])
        plt.plot(self.timeMS_AHR2,self.ROLL2,label='AHRS2')
        plt.plot(self.timeMS_EKF1,self.ROLL1,label='EKF1')
        plt.plot(self.timeMS_ATT,self.ROLL,label='AUTOPILOT')
        plt.plot(self.timeMS_ATT,self.ROLL_DES,label='Command')
        plt.grid(); plt.xlim([x0,xf]); #plt.ylim([-50,50])
        plt.xlabel('Time (ms)')
        plt.ylabel('Roll Angle (deg)')
        plt.legend(loc=2)
        pp.savefig()

        plt.figure()
        print("Pitch and Pitch Desired = ",self.PITCH[0],self.PITCH_DES[0])
        plt.plot(self.timeMS_AHR2,self.PITCH2,label='AHRS2')
        plt.plot(self.timeMS_EKF1,self.PITCH1,label='EKF1')
        plt.plot(self.timeMS_ATT,self.PITCH,label='AUTOPILOT')
        plt.plot(self.timeMS_ATT,self.PITCH_DES,label='Command')
        plt.grid(); plt.xlim([x0,xf]) ; #plt.ylim([-50,50])
        plt.xlabel('Time (ms)')
        plt.ylabel('Pitch Angle (deg)')
        plt.legend(loc=2)
        pp.savefig()

        plt.figure()
        plt.plot(self.timeMS_AHR2,self.YAW2,label='AHRS2')
        plt.plot(self.timeMS_EKF1,self.YAW1,label='EKF1')
        plt.plot(self.timeMS_ATT,self.YAW,label='AUTOPILOT')
        plt.plot(self.timeMS_ATT,self.YAW_DES,label='Command')
        plt.grid(); plt.xlim([x0,xf]) ; 
        plt.xlabel('Time (ms)')
        plt.ylabel('Yaw Angle (deg)')
        plt.legend(loc=2)
        pp.savefig()
        

#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# Class for IMU
#
#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

class IMU(object):
    def __init__(self):
        #EKF2
        self.timeMS_EKF2 = []
        self.AX = []
        self.AY = []
        self.AZ = []
        self.MX = []
        self.MY = []
        self.MZ = []
        self.EKF2_OFFSET = 0
        #IMU
        self.timeMS_IMU = []
        self.GyrX = []
        self.GyrY = []
        self.GyrZ = []
        self.AccX = []
        self.AccY = []
        self.AccZ = []
        self.IMU_OFFSET = 0
        #IMU2
        self.timeMS_IMU2 = []
        self.GyrX2 = []
        self.GyrY2 = []
        self.GyrZ2 = []
        self.AccX2 = []
        self.AccY2 = []
        self.AccZ2 = []
        self.IMU2_OFFSET = 0

    def convert2np(self):
        self.timeMS_IMU = np.asarray(self.timeMS_IMU)
        self.GyrX = np.asarray(self.GyrX)
        self.GyrY = np.asarray(self.GyrY)
        self.GyrZ = np.asarray(self.GyrZ)
        self.AccX = np.asarray(self.AccX)
        self.AccY = np.asarray(self.AccY)
        self.AccZ = np.asarray(self.AccZ)

    def plot_IMU(self,x0,xf,pp,baroRead=None):
        if x0 == -99:
            x0 = baroRead.timeSec[0]
            xf = baroRead.timeSec[-1]
            
        plt.figure()
        plt.plot(self.timeMS_IMU,self.AccX,label='IMU')
        plt.plot(self.timeMS_IMU2,self.AccX2,label='IMU2')
        #plt.plot(self.timeMS_EKF2,self.AX,label='EKF2')
        P.xlim_auto(x0,xf,self.timeMS_IMU,self.AccX)
        plt.grid(); 
        plt.xlabel('Time (ms)')
        plt.ylabel('Acceleration X (m/s^2)')
        plt.legend(loc=2)
        pp.savefig()

        plt.figure()
        plt.plot(self.timeMS_IMU,self.AccY,label='IMU')
        plt.plot(self.timeMS_IMU2,self.AccY2,label='IMU2')
        #plt.plot(self.timeMS_EKF2,self.AY,label='EKF2')
        P.xlim_auto(x0,xf,self.timeMS_IMU,self.AccY)
        plt.grid(); 
        plt.xlabel('Time (ms)')
        plt.ylabel('Acceleration Y (m/s^2)')
        plt.legend(loc=2)
        pp.savefig()

        plt.figure()
        plt.plot(self.timeMS_IMU,self.AccZ,label='IMU')
        plt.plot(self.timeMS_IMU2,self.AccZ2,label='IMU2')
        #plt.plot(self.timeMS_EKF2,self.AZ,label='EKF2')
        P.xlim_auto(x0,xf,self.timeMS_IMU,self.AccZ)
        plt.grid()
        plt.xlabel('Time (ms)')
        plt.ylabel('Acceleration Z (m/s^2)')
        plt.legend(loc=2)
        pp.savefig()
        
        plt.figure()
        #plt.plot(self.timeMS_EKF2,self.MX,label='EKF2')
        plt.plot(self.timeMS_IMU,self.GyrX,label='IMU')
        plt.plot(self.timeMS_IMU2,self.GyrX2,label='IMU2')
        P.xlim_auto(x0,xf,self.timeMS_IMU,self.GyrX)
        plt.grid()
        plt.xlabel('Time (ms)')
        plt.ylabel('Rate Gyro X (deg/s)')
        plt.legend(loc=2)
        pp.savefig()

        plt.figure()
        #plt.plot(self.timeMS_EKF2,self.MY,label='EKF2')
        plt.plot(self.timeMS_IMU,self.GyrY,label='IMU')
        plt.plot(self.timeMS_IMU2,self.GyrY2,label='IMU2')
        P.xlim_auto(x0,xf,self.timeMS_IMU,self.GyrY)
        plt.grid()
        plt.xlabel('Time (ms)')
        plt.ylabel('Rate Gyro Y (deg/s)')
        plt.legend(loc=2)
        pp.savefig()

        plt.figure()
        #plt.plot(self.timeMS_EKF2,self.MZ,label='EKF2')
        plt.plot(self.timeMS_IMU,self.GyrZ,label='IMU')
        plt.plot(self.timeMS_IMU2,self.GyrZ2,label='IMU2')
        P.xlim_auto(x0,xf,self.timeMS_IMU,self.GyrZ)
        plt.grid()
        plt.xlabel('Time (ms)')
        plt.ylabel('Rate Gyro Z (deg/s)')
        plt.legend(loc=2)
        pp.savefig()

#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# Class for barometer readings
#
# CRt = climb rate
#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

class barometer(object):
    def __init__(self):
        self.timeMS      = []
        self.altitude    = []
        self.pressure    = []
        self.temperature = []
        self.timeHr      = []
        self.timeHr_GPS  = []
        self.timeMS_GPS  = []
        self.timeOffset  = 0
        
    def establish_units(self,gpsRead):
        self.pressureMB = []
        self.pressurekPA = []
        for i in self.pressure:
            self.pressureMB.append(i / 100.0)
            self.pressurekPA.append(i/100.0 * 0.1)
        self.timeSec = []
        #Find gps hour but this only works if GPSINIT == 1
        GPSON = 0
        if gpsRead.GPSINIT == 1:
            for i in gpsRead.time_quad:
                if abs(i) > 0.01 and GPSON == 0:
                    GPSON = i
        else:
            GPSON = 0
        #print GPSON
        for i in self.timeMS:
            self.timeSec.append(i / 1000.0)
            #Offset timeHr by gpstime
            self.timeHr.append(i/1000.0 * (1.0/3600.0) + GPSON)

        for i in range(len(self.timeMS_GPS)):
            MJD_DateSpace = np.floor(gpsRead.WEEK * 7 + self.timeMS_GPS[i] / 86400000)
            # print('Week = ',gpsRead.WEEK)
            # print('timeMS_GPS = ',self.timeMS_GPS[i])
            # print('MJD_DateSpace = ',MJD_DateSpace)
            MJD_Start = DateToMJD(1980, 1, 6) + MJD_DateSpace
            # print('MJD_Start = ',MJD_Start)
            year, month, day = MJDToDate(MJD_Start)
            # print('ymd = ',year,month,day)
            self.month = month
            self.day = day
            self.year = year
            
            if self.timeMS_GPS[i] == -99:
                val = GPSON
            else:
                val = self.timeMS_GPS[i]
            secondsLeft = (val / 1000) % 86400
            # print('SecondsLeft = ',secondsLeft)
            hours = secondsLeft / 3600
            secondsLeft = secondsLeft - hours * 3600
            minutes = secondsLeft / 60
            seconds = secondsLeft - minutes * 60
            milli = self.timeMS_GPS[i] % 1000

            currentDT = datetime.datetime(int(year), int(month), int(day), int(hours), int(minutes), int(seconds), int(milli))
            #print hours,minutes,seconds,milli
            mytime = hours+minutes/60.0+seconds/3600.0+(milli/1000.0)/3600.0
            #print mytime
            #print self.timeMS_GPS[i],hours,minutes,seconds,milli,mytime

            self.timeHr_GPS.append(mytime)
        self.timeSec_GPS_np = np.array(self.timeHr_GPS)
        self.timeSec_GPS_np = (self.timeSec_GPS_np-self.timeSec_GPS_np[0])*3600.0

        #THis is for the barometer
        self.altitude_np = np.array(self.altitude)
        self.altitude_MSL_np = self.altitude_np
        if gpsRead.GPSINIT == 1:
            self.altitude_MSL_np += gpsRead.altitude[0]

        self.timeSec = np.array(self.timeSec)
            
    def plot_barometer(self, figObj,x0,xf,pp):
        axis1 = figObj.add_subplot(1, 1, 1)
        axis1.plot(self.timeSec, self.pressureMB, color = "blue")
        # axis1.set_ylim([1010.0, 1020.0])
        axis1.get_yaxis().get_major_formatter().set_useOffset(False)
        axis1.set_ylabel("Pressure (mb)")
        axis1.set_xlabel("Time (ms)")
        axis1.set_title("Drone flight ")
        #axis1.set_xlim([x0,xf])
        axis1.grid()
        # for t1 in axis1.get_yticklabels():
        #     t1.set_color("blue")
        
        # axis2 = axis1.twinx()
        # axis1.plot(self.timeSec, self.altitude, color = "red")
        #axis1.set_ylim([-2.0, 30.0])
        # axis1.set_ylabel("Altitude [m]")
        #axis2.set_xlim([500,600])
        
        # for t1 in axis2.get_yticklabels():
        #     t1.set_color("red")

        pp.savefig()

#/////////////////////////
#Clas for battery health
#///////////////////////

class battery(object):
    def __init__(self):
        self.timeSec = []
        self.Voltage = []
        self.Current = []
        self.OFFSET = 0

    def plot_battery(self,pp):
        plt.figure()
        plt.plot(self.timeSec,self.Voltage)
        plt.xlabel('Time (ms)')
        plt.ylabel('Voltage (V)')
        #plt.xlim
        plt.grid()
        pp.savefig()

        plt.figure()
        plt.plot(self.timeSec,self.Current)
        plt.xlabel('Time (ms)')
        plt.ylabel('Current (A)')
        #plt.xlim
        plt.grid()
        pp.savefig()
        
#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# Class for GPS readings
#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

class gps(object):
    def __init__(self):
        self.timeMS_GPS = []
        self.timeMS_internal = []
        #self.rawTIME_GPS = []
        self.timeMS_CTUN = []
        self.week_GPS   = []
        self.WEEK = -99
        self.numSats    = []
        self.latitude   = []
        self.longitude  = []
        self.altitude   = []
        self.rel_altitude = []
        self.speed      = []
        self.VZ         = []
        self.month      = []
        self.day        = []
        self.D_altitude = []
        self.baro_altitude = []
        self.altitude_CTUN = []
        self.altitude_MSL = []
        self.DS_altitude = []
        self.SON_altitude = []
        self.CLIMB_RATE_DES = []
        self.CLIMB_RATE = []
        self.internal_OFFSET = 0
        self.CTUN_OFFSET = 0
        self.GPS_OFFSET = 0
        self.GPSINIT = 0
        self.GPSON_Sec = 0

    def fix_gps(self):
        #altitude is MSL altitude but uses GPS and is inaccurate
        #altitude_CTUN is AGL and uses the barometer to get more accurate
        #data but is at a higher frequency. 
        
        for i in range(len(self.altitude_CTUN)):
            self.altitude_MSL.append(self.altitude_CTUN[i] + self.altitude[0])

        #Ok so now altitude_MSL is a combination of GPS and barometer so it's highly
        #accurate and reports MSL instead of AGL but it's at a higher frequency
        #so we need to interpolate
        # timeSec and altitude have the same amount of data points
        # timeMS_CTUN and altitude_CTUN have the same amount of data points
        # so altitude_MSL has the same number of data point as timeMS_CTUN
        #EXAMPLE = np.interp(gpsRead.timeSec,ATTRead.timeMS_ATT,ATTRead.ROLL)
        #ROLL and timeMS_ATT have the same number of data points
        self.altitude_MSL_interp = np.interp(self.timeSec,self.timeMS_CTUN,self.altitude_MSL)

        ##Fix Lat and Lon
        m_lat = np.mean(self.latitude)
        m_lon = np.mean(self.longitude)
        tol = 0.00002
        for i in range(0,len(self.latitude)):
            #self.latitude[i] = m_lat
            #self.longitude[i] = m_lon
            lat = self.latitude[i]
            if lat > (1.0+tol)*m_lat or lat < (1.0-tol)*m_lat:
                self.latitude[i] = m_lat
            lon = self.longitude[i]
            if abs(lon) > (1+tol)*abs(m_lon) or abs(lon) < (1-tol)*abs(m_lon):
                self.longitude[i] = m_lon
            
    def calibrate_time(self):
        self.timeUTC = []
        self.time_quad = []

        #print self.timeMS_GPS

        for i in range(len(self.timeMS_GPS)):
            MJD_DateSpace = self.week_GPS[i] * 7 + (self.timeMS_GPS[i] / 86400000)
            # print('MJD_DateSpace2 = ',MJD_DateSpace)
            MJD_Start = DateToMJD(1980, 1, 6) + MJD_DateSpace
            year, month, day = MJDToDate(MJD_Start)
            self.month = month
            self.day = day
            self.year = year

            secondsLeft = (self.timeMS_GPS[i] / 1000) % 86400
            # print('SecondsLeft2 = ',secondsLeft)
            hours = int(secondsLeft / 3600)
            secondsLeft = secondsLeft - hours * 3600
            minutes = int(secondsLeft / 60)
            seconds = int(secondsLeft - minutes * 60)
            milli = self.timeMS_GPS[i] % 1000

            #print('YMDHMSM = ',year,month,day,hours,minutes,seconds,milli)
            currentDT = datetime.datetime(year, month, day, hours, minutes, seconds, milli)
            #print hours,minutes,seconds,milli
            mytime = hours+minutes/60.0+seconds/3600.0+(milli/1000.0)/3600.0
            #print mytime
            #print self.timeMS_GPS[i],hours,minutes,seconds,milli,mytime

            if len(self.time_quad) > 0:
                if mytime == self.time_quad[-1]:
                    mytime+=1e-8
            self.time_quad.append(mytime)
            self.timeUTC.append(currentDT)

        #print self.time_quad[0]
        #This convert time_quad to seconds
        #self.timeSec = np.array(self.time_quad)
        #self.timeSec = (self.timeSec-self.timeSec[0])*3600.0 + self.GPSON_Sec
        self.timeSec = np.array(self.timeMS_internal)/1000.0

        print(self.timeSec[-1],self.timeSec[0])

    def plotGPS_Speed(self,x0,xf,pp):
        plt.figure()
        plt.plot(self.timeSec,self.speed)
        plt.grid(); plt.xlim([x0,xf])
        plt.xlabel('Time (ms)')
        plt.ylabel('Lateral Speed (m/s)') #It says m/s here but is this right? Yes this is right. 
        pp.savefig() ##According to the following website http://ardupilot.org/copter/docs/common-downloading-and-analyzing-data-logs-in-mission-planner.html

        plt.figure()
        plt.plot(self.timeSec,self.VZ,label="GPS")
        plt.plot(self.timeMS_CTUN,self.CLIMB_RATE_DES,label="Desired")
        plt.plot(self.timeMS_CTUN,self.CLIMB_RATE,label="AUTOPILOT")
        plt.legend(loc="best")
        plt.grid(); plt.xlim([x0,xf])
        plt.xlabel('Time (ms)')
        plt.ylabel('Vertical Speed (m/s)')
        pp.savefig()

    def plot_Time(self,pp):
        plt.figure()
        plt.plot(self.timeSec,self.time_quad,marker='s')
        plt.grid(); 
        plt.xlabel('Time (ms)')
        plt.ylabel('Time GPS (HH.MM/100)')
        pp.savefig()
    
    def plot_latlon(self, figObj,pp):
        axis1 = figObj.add_subplot(1, 1, 1)
        axis1.plot(self.longitude, self.latitude, color = "green")
        rangeLon = max(self.longitude) - min(self.longitude)
        rangeLat = max(self.latitude) - min(self.latitude)

        plt.ticklabel_format(style = "plain", useOffset = False)
        
        axis1.set_title("Lat/Lon during " + str(self.month) + "/" + str(self.day))
        axis1.set_ylabel("Latitude")
        axis1.set_xlabel("Longitude")
        axis1.set_xlim([min(self.longitude) - 0.25 * rangeLon, max(self.longitude) + 0.25 * rangeLon])
        axis1.set_ylim([min(self.latitude) - 0.25 * rangeLat, max(self.latitude) + 0.25 * rangeLat])
        axis1.grid()
        
        #axis2 = plt.axes([0.2, 0.7, 0.15, 0.15])
        #axis2.plot(self.timeUTC, self.speed, color = "purple")
        #axis2.set_title("Speed", fontsize = 6)
        # #axis2.set_xlabel("UTC Time", fontsize = 6)
        # axis2.set_ylabel("Speed [??]", fontsize = 6)
        # axis2.xaxis.set_major_formatter(DateFormatter('%H:%M:%S'))
        # axis2.grid(which = "both", axis = "both", color = "gray", linestyle = "dotted", linewidth = 0.5)
        # for t1 in axis2.get_yticklabels():
        #     t1.set_fontsize(6)
        # for t1 in axis2.get_xticklabels():
        #     t1.set_fontsize(6)
        #     t1.set_rotation(60)

        pp.savefig()
        
#Cross Class Plotting///////////////////////////////////////////////

def plot_altitude(quad_data,x0,xf,pp):
    baroRead = quad_data[0]
    gpsRead = quad_data[1]
    rcinoutRead = quad_data[2]
    IMURead = quad_data[3]
    ATTRead = quad_data[4]

    plt.figure()
    plt.plot(baroRead.timeSec,baroRead.altitude,label="Barometer")
    #plt.plot(ATTRead.timeMS_AHR2,ATTRead.altitude,label="AHR2")
    if gpsRead.GPSINIT == 1:
        #plt.plot(gpsRead.timeSec,gpsRead.altitude,label="GPS (MSL)")
        plt.plot(gpsRead.timeSec,gpsRead.rel_altitude,label="GPS (AGL)")
        plt.plot(gpsRead.timeMS_CTUN,gpsRead.altitude_CTUN,label="AUTOPILOT") ##Yea so altitude_CTUN is AUTOPILOT which stiches barometer and GPS
        plt.plot(gpsRead.timeMS_CTUN,gpsRead.D_altitude,label="Desired")
        #plt.plot(gpsRead.timeMS_CTUN,gpsRead.SON_altitude,label="Sonar")
        #plt.plot(gpsRead.timeMS_CTUN,gpsRead.DS_altitude,label="Desired S")
        #plt.plot(gpsRead.timeMS_CTUN,gpsRead.baro_altitude,label="Barometer 2")
    
    plt.xlabel("Time (ms)")
    plt.xlim([x0,xf])
    plt.ylabel("Altitude (m)")
    plt.grid()
    plt.legend(loc='best')
    pp.savefig()

##########################################################################################
##########################################################################################
# Main program
#--------------------------------------------------------------------------------

#------------------------------------------------------------
# Hardcode file name
#------------------------------------------------------------

if __name__ == "__main__":

    SHOWPLOTS = 0
    pp = PDF(SHOWPLOTS,plt)

    if len(sys.argv) > 1:
        fileName = sys.argv[1]
    else:
        print(sys.argv)
        print('No Filename found or given')
        sys.exit()

    #fileName = "Data_Files/10_25_2016/FlightALL.log"
    #fileName = "Data_Files/11_29_2016/Flight2.log"

    quad_data = get_quad_data(fileName)

    if not quad_data:
        sys.exit()

    print("Quad Data Directory = ",quad_data[-1])
    baroRead = quad_data[0]
    gpsRead = quad_data[1]
    rcinoutRead = quad_data[2]
    IMURead = quad_data[3]
    ATTRead = quad_data[4]
    BATTRead = quad_data[5]

    #Defaults. DO NOT TOUCH
    x0 = baroRead.timeSec[0]
    xf = baroRead.timeSec[-1]

    # x0 = 530
    # xf = 550
    # x0 = 40
    # xf = 50

    figure1 = plt.figure()#figsize = (11, 8))
    print('Plotting Barometer')
    baroRead.plot_barometer(figure1,x0,xf)

    if gpsRead.GPSINIT == 1:
        print('Plotting GPS')
        figure2 = plt.figure()
        gpsRead.plot_latlon(figure2)
        gpsRead.plotGPS_Speed(x0,xf)
        gpsRead.plot_Time()
    else:
        print("GPSINIT NOT FOUND")

    print('Plotting RC')
    rcinoutRead.plot_rc(x0,xf)

    print('Plotting IMU')
    IMURead.plot_IMU(x0,xf,pp)

    print('Plotting Attitude')
    ATTRead.plot_ATTITUDE(x0,xf)

    print('Plotting Altitude')
    plot_altitude(quad_data,x0,xf)

    print('Plot battery')
    BATTRead.plot_battery()

    ##Output stuff to textfile
    print('Creating text file...')
    create_text_file(quad_data,'Quad_Data.out')
    
    pp.close()

# http://ardupilot.org/copter/docs/common-downloading-and-analyzing-data-logs-in-mission-planner.html

#Currently PLotting - Remember this changed in 2.0.18 to 2.0.21
#FMT, 138, 25, AHR2, IccCfLL, TimeMS,Roll,Pitch,Yaw,Alt,Lat,Lng
#FMT, 1, 23, ATT, IccccCCCC, TimeMS,DesRoll,Roll,DesPitch,Pitch,DesYaw,Yaw,ErrRP,ErrYaw
#FMT, 136, 21, BARO, Iffcf, TimeMS,Alt,Press,Temp,CRt
#FMT, 140, 43, EKF1, IccCffffffccc, TimeMS,Roll,Pitch,Yaw,VN,VE,VD,PN,PE,PD,GX,GY,GZ
#FMT, 130, 45, GPS, BIHBcLLeeEefI, Status,TimeMS,Week,NSats,HDop,Lat,Lng,RelAlt,Alt,Spd,GCrs,VZ,T
#FMT, 4, 33, CTUN, Ihhhffecchh, TimeMS,ThrIn,AngBst,ThrOut,DAlt,Alt,BarAlt,DSAlt,SAlt,DCRt,CRt
#FMT, 131, 31, IMU, Iffffff, TimeMS,GyrX,GyrY,GyrZ,AccX,AccY,AccZ
#FMT, 135, 31, IMU2, Iffffff, TimeMS,GyrX,GyrY,GyrZ,AccX,AccY,AccZ
#FMT, 133, 35, RCIN, Ihhhhhhhhhhhhhh, TimeMS,C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,C13,C14
#FMT, 134, 23, RCOU, Ihhhhhhhh, TimeMS,Chan1,Chan2,Chan3,Chan4,Chan5,Chan6,Chan7,Chan8
#FMT, 9, 23, CURR, IhIhhhf, TimeMS,ThrOut,ThrInt,Volt,Curr,Vcc,CurrTot

#############################Still need#####################################


#FMT, 141, 26, EKF2, Ibbbcchhhhhh, TimeMS,AX,AY,AZ,VWN,VWE,MN,ME,MD,MX,MY,MZ
#FMT, 142, 27, EKF3, Icccccchhhc, TimeMS,IVN,IVE,IVD,IPN,IPE,IPD,IMX,IMY,IMZ,IVT
#FMT, 143, 25, EKF4, IcccccccbbBB, TimeMS,SV,SP,SH,SMX,SMY,SMZ,SVT,OFN,EFE,FS,DS

#FMT, 15, 25, MAG, Ihhhhhhhhh, TimeMS,MagX,MagY,MagZ,OfsX,OfsY,OfsZ,MOfsX,MOfsY,MOfsZ
#FMT, 27, 25, MAG2, Ihhhhhhhhh, TimeMS,MagX,MagY,MagZ,OfsX,OfsY,OfsZ,MOfsX,MOfsY,MOfsZ

##########################NEED TO FIND DOCUMENTATION ON THESE###################

#FMT, 137, 13, POWR, ICCH, TimeMS,Vcc,VServo,Flags
#FMT, 146, 16, RAD, IBBBBBHH, TimeMS,RSSI,RemRSSI,TxBuf,Noise,RemNoise,RxErrors,Fixed

#######################Don't thinki I'll ever need these################

#FMT, 129, 23, PARM, Nf, Name,Value
#FMT, 23, 8, DU32, BI, Id,Value
#FMT, 19, 5, ERR, BB, Subsys,ECode
#FMT, 3, 6, MODE, Mh, Mode,ThrCrs - Maybe could use this one day
#FMT, 0, 129, 23, Nf, Name,Value
#FMT, 6, 17, PM, HHIhBHB, NLon,NLoop,MaxT,PMT,I2CErr,INSErr,INAVErr
#FMT, 151, 12, UBX1, IBHBB, TimeMS,Instance,noisePerMS,jamInd,aPower
#FMT, 152, 12, UBX2, IBbBbB, TimeMS,Instance,ofsI,magI,ofsQ,magQ

################Only needed for Waypoint control    ####################33

#FMT, 145, 41, CMD, IHHHfffffff, TimeMS,CTot,CNum,CId,Prm1,Prm2,Prm3,Prm4,Lat,Lng,Alt
#FMT, 5, 47, NTUN, Iffffffffff, TimeMS,DPosX,DPosY,PosX,PosY,DVelX,DVelY,VelX,VelY,DAccX,DAccY
