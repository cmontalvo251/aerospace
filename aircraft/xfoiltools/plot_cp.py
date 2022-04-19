#!/usr/bin/python

##Import modules
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from xfoiltools import computeclcd

fAoA = ['-15.7894736842','15.7894736842']

### Get ready for straight to pdf
os.system('rm plots.pdf')
pdfhandle = PdfPages('plots.pdf')

for fctr in range(0,len(fAoA)):
    AoA = fAoA[fctr]
    saveFlnmAF = 'Airfoil' + AoA + '.txt'
    saveFlnmCp = 'Cp' + AoA + '.txt'
    cl,cd = computeclcd(AoA,saveFlnmAF,saveFlnmCp,pdfhandle,0)

pdfhandle.close()
os.system('evince plots.pdf')
