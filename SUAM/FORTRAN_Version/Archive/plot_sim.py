#!/usr/bin/python

import numpy as np 
import matplotlib.pyplot as pyplot
import sys
from pdf import *
import mymath as M 


def get_state_data(inputfilename):

	print('Loading File...')

	statedata = open(inputfilename)

	time = []
	x = []
	y = []
	z = []
	phi = []
	theta = []
	psi = []
	u = []
	v = []
	w = []
	p = []
	q = []
	r = []
	T1 = []
	T1dot = []
	T2 = []
	T2dot =[]
	T3 = []
	T3dot = []
	T4 = []
	T4dot = []

	lenfile = 0

        datastate = []

	for line in statedata:
		lenfile+=1

		if len(line) > 2:
			row = line.split(',')
                else:
                        print 'whoops'

		if len(row) > 1:
                        row_np = []
                        for x in row:
                                val = np.float(x)
                                row_np.append(val)
                        datastate.append(np.asarray(row_np))

	print('datastate made')
	return np.asarray(datastate)

def create_state_plots(datastate,pp):
	#create sense hat plots
	time = datastate[:,0]
	x = datastate[:,1]
	y = datastate[:,2]
	z = datastate[:,3]
	u = datastate[:,4]
	v = datastate[:,5] 
	w = datastate[:,6]
	phi= datastate[:,7]
	theta = datastate[:,8]
	psi = datastate[:,9]
	p = datastate[:,10]
	q = datastate[:,11]
	r = datastate[:,12]
	T1 = datastate[:,13]
	T1dot = datastate[:,14]
	T2 = datastate[:,15]
	T2dot = datastate[:,16]
	T3 = datastate[:,17]
	T3dot = datastate[:,18]
	T4 = datastate[:,19]
	T4dot = datastate[:,20]

	#plot x
	plt.figure()
	plt.plot(time,x,label = 'x')
	plt.xlabel('Time (s)')
	plt.ylabel('x (m)')
	plt.grid()
	plt.legend()
	plt.title('x')
	pp.savefig()
	

	#plot y

	plt.figure()
	plt.plot(time,y,label = 'y')
	plt.xlabel('Time (s)')
	plt.ylabel('y (m)')
	plt.grid()
	plt.legend()
	plt.title('y vs. t')
	pp.savefig()
	#plot z
	plt.figure()
	plt.plot(time,z,label = 'z')
	plt.xlabel('Time (s)')
	plt.ylabel('z (m)')
	plt.grid()
	plt.legend()
	plt.title('z vs. t')
	pp.savefig()
	#plot phi
	plt.figure()
	plt.plot(time,phi,label = r'\phi')
	plt.xlabel('Time (s)')
	plt.ylabel('x (m)')
	plt.grid()
	plt.legend()
	plt.title('x')
	pp.savefig()

		#plot phi
	plt.figure()
	plt.plot(time,theta,label = r'\theta')
	plt.xlabel('Time (s)')
	plt.ylabel(r'\theta (rad)')
	plt.grid()
	plt.legend()
	plt.title('x')
	pp.savefig()

		#plot phi
	plt.figure()
	plt.plot(time,psi,label = r'\phi')
	plt.xlabel('Time (s)')
	plt.ylabel('x (m)')
	plt.grid()
	plt.legend()
	plt.title('x')
	pp.savefig()

		#plot phi
	plt.figure()
	plt.plot(time,u,label = r'\phi')
	plt.xlabel('Time (s)')
	plt.ylabel('x (m)')
	plt.grid()
	plt.legend()
	plt.title('x')
	pp.savefig()

		#plot phi
	plt.figure()
	plt.plot(time,v,label = r'\phi')
	plt.xlabel('Time (s)')
	plt.ylabel('x (m)')
	plt.grid()
	plt.legend()
	plt.title('x')
	pp.savefig()

		#plot phi
	plt.figure()
	plt.plot(time,w,label = r'\phi')
	plt.xlabel('Time (s)')
	plt.ylabel('x (m)')
	plt.grid()
	plt.legend()
	plt.title('x')
	pp.savefig()

		#plot phi
	plt.figure()
	plt.plot(time,p,label = r'\phi')
	plt.xlabel('Time (s)')
	plt.ylabel('x (m)')
	plt.grid()
	plt.legend()
	plt.title('x')
	pp.savefig()

		#plot phi
	plt.figure()
	plt.plot(time,q,label = r'\phi')
	plt.xlabel('Time (s)')
	plt.ylabel('x (m)')
	plt.grid()
	plt.legend()
	plt.title('x')
	pp.savefig()

		#plot phi
	plt.figure()
	plt.plot(time,r,label = r'\phi')
	plt.xlabel('Time (s)')
	plt.ylabel('x (m)')
	plt.grid()
	plt.legend()
	plt.title('x')
	pp.savefig()

		#plot phi
	plt.figure()
	plt.plot(time,T1,label = r'\phi')
	plt.xlabel('Time (s)')
	plt.ylabel('x (m)')
	plt.grid()
	plt.legend()
	plt.title('x')
	pp.savefig()

		#plot phi
	plt.figure()
	plt.plot(time,T1dot,label = r'\phi')
	plt.xlabel('Time (s)')
	plt.ylabel('x (m)')
	plt.grid()
	plt.legend()
	plt.title('x')
	pp.savefig()

		#plot phi
	plt.figure()
	plt.plot(time,T2,label = r'\phi')
	plt.xlabel('Time (s)')
	plt.ylabel('x (m)')
	plt.grid()
	plt.legend()
	plt.title('x')
	pp.savefig()

	#plot phi
	plt.figure()
	plt.plot(time,T2dot,label = r'\phi')
	plt.xlabel('Time (s)')
	plt.ylabel('x (m)')
	plt.grid()
	plt.legend()
	plt.title('x')
	pp.savefig()
	#plot phi
	plt.figure()
	plt.plot(time,T3,label = r'\phi')
	plt.xlabel('Time (s)')
	plt.ylabel('x (m)')
	plt.grid()
	plt.legend()
	plt.title('x')
	pp.savefig()
		#plot phi
	plt.figure()
	plt.plot(time,T3dot,label = r'\phi')
	plt.xlabel('Time (s)')
	plt.ylabel('x (m)')
	plt.grid()
	plt.legend()
	plt.title('x')
	pp.savefig()
		#plot phi
	plt.figure()
	plt.plot(time,T4,label = r'\phi')
	plt.xlabel('Time (s)')
	plt.ylabel('x (m)')
	plt.grid()
	plt.legend()
	plt.title('x')
	pp.savefig()	

		#plot phi
	plt.figure()
	plt.plot(time,T4dot,label = r'\phi')
	plt.xlabel('Time (s)')
	plt.ylabel('x (m)')
	plt.grid()
	plt.legend()
	plt.title('x')
	pp.savefig()
        
if __name__ == '__main__':

	print('Processing State output File')

	if len(sys.argv) < 2:
		print('Not enough input arguments. Need File location')
		sys.exit()
	else:
		print(sys.argv)
		inputfilenames = [sys.argv[1]]
	
	#inputfilename = 'State.OUTPUT'

	SHOWPLOTS = 0
	pp = PDF(SHOWPLOTS,plt)

	data_all = []
	inputfilename = inputfilenames[0]
	data = get_state_data(inputfilename)
	#data_all.append(data)
	
	#datastate = data_all
	#print len(datastate[0])

	create_state_plots(data,pp)

	pp.close()
