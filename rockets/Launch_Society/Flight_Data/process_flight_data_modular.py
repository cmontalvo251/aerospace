#!/usr/bin/python

##Import toolboxes
import numpy as np
import sys
import os

##Grab filename
if len(sys.argv) > 1:
    filename = sys.argv[1]
    raw_filename, file_extension = os.path.splitext(filename)
    print "Filename = ",filename
else:
    print "Filename not found"
    sys.exit();

##Open file name a parse through file
file_in = open(raw_filename+'.csv','rb')
file_out = open(raw_filename+'_processed.csv','wb')

for line in file_in:
    #Check for comment
    out_row = ''
    c = ','
    COMMENT = '%'
    DELIMITER = ' '
    if line[0] != COMMENT:
        #Split row
        row = line.split(DELIMITER)
        
        #Taken from Beecher Faust's output file
        # % Time (s) Altitude (ft) Vertical velocity (ft/s) Vertical acceleration (ft/s2) Total velocity (ft/s)
        # Total acceleration (ft/s2) Position East of launch (ft) Position North of launch (ft) Lateral distance (ft)
        # Lateral velocity (ft/s) Lateral acceleration (ft/s2) Latitude (deg) Longitude (deg) Angle of attack (deg) Roll rate (r/s)
        # Pitch rate (r/s) Yaw rate (r/s) Mass (oz) Longitudinal moment of inertia (lb-ft2) Rotational moment of inertia (lb-ft2)
        # % Thrust (N) Pitch moment coefficient (?) Yaw moment coefficient (?)
        for val in row:
            if len(val) > 0:
                out_row += str(val) + c
        out_row = out_row[:-1]
        print out_row
        file_out.write(out_row)

file_in.close()
file_out.close()
        
        
    
    
