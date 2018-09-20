# parallel_radius.py
# James Widdicombe
# Last Updated 20/09/2018
# Script to calculate radius of an axion star

# Load the modules
import yt
import os
import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import fsolve
from yt import derived_field
import time
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

# Enable Parallelism
yt.enable_parallelism()

# Timings
start_time = time.time()

# CHANGE ME
# Loading dataset (Load from one folder up)
data_location = '../../plt*.hdf5' # Data file location

#Loading dataset
ts = yt.load(data_location)

# Plot parameters
line = 2.5 # Line width
alp = 0.6 # alpha

# Other factors
total_box_size = float(ts[0].domain_right_edge[0])
center = total_box_size/2.
Quality = 'cubic'

# Matplotlib Settings
rcParams.update({'figure.autolayout': True})
rcParams['axes.formatter.limits'] = [-3,3]

# Define an empty storage dictionary for collecting information
# in parallel through processing
storage = {}

for sto, i in ts.piter(storage=storage):
	#Timings
	L_start = time.time()
	
	# Look at initial Stars
	# Line out from 0 to center (x), at the center of (y and z)
	c = i.ray((0,center,center),(center,center,center))

	# Name x coordinate and rho from lineout
	x = c["x"]
	rho = c["rho"]

	# Find maximum rho value within lineout
	rho_max = float(rho.max())

	# Find the coordinate of where the maximum value is
	x_max = np.where(rho == rho_max)[0][0]
	x_max_val = x[x_max]
	
	 # Create a function that goes to negative when we are below 95% of rho max
	rhofun = interp1d(x, rho-0.05*rho_max, kind=Quality, fill_value="extrapolate")

	# Solve the function for where it goes to zero using a best guess past and
	# before the the max value of x
	x_1 = fsolve(rhofun,float(x[x_max])-1)[0]
	x_2 = fsolve(rhofun,float(x[x_max])+1)[0]

	# Size of the 95% rho in center
	size_x = (x_2-x_1)/2.
	
	# Store the frames information
	array = [i.current_time,time.time()-L_start,x_max_val,size_x]
	sto.result = array
	sto.result_id = str(i)

if yt.is_root():
	looptime = []
	Ftime = []
	rhomaxpos = []
	rhosize = []
	for L in sorted(storage.items()):
		looptime.append(float(L[1][1]))
		Ftime.append(float(L[1][0]))
		rhomaxpos.append(float(L[1][2]))
		rhosize.append(float(L[1][3]))
		#print L[1]
		
	# Max rho pos
	plt.figure(1)
	plt.plot(Ftime,rhomaxpos)
	plt.xlabel('Time $[1/m]$')
	plt.ylabel('Position x axis $[1/m]$')
	plt.grid()
	plt.savefig('max_rho_pos.png')
	plt.close()
	
	# Star Radius
	plt.figure(2)
	plt.plot(Ftime,rhosize)
	plt.xlabel('Time $[1/m]$')
	plt.ylabel('Star Radius $[1/m]$')
	plt.ylim(0,12)
	plt.grid()
	plt.savefig('star_radius.png')
	plt.close()
	
	total_compute = sum(looptime)
	speedup =  total_compute/(time.time()-start_time)
	np.savetxt('time.out',Ftime)
	np.savetxt('rho_max_pos.out',rhomaxpos)
	np.savetxt('rho_size.out',rhosize)
	np.savetxt('loop_time.out',looptime)
	np.savetxt('speed_up.out',[speedup])
	