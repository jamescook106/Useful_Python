# parallel_convergence.py
# James Widdicombe
# Last Updated 08/10/2018
# Calculate Convergence
# Currently only calculates rho

# Load the modules
import yt
import os
import numpy as np
from yt import derived_field
import time
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams
from scipy.interpolate import interp1d
from scipy.optimize import fsolve
from scipy import ndimage
from scipy.interpolate import RegularGridInterpolator

# Enable Parallelism
yt.enable_parallelism()

# Plot parameters
line = 2.5 # Line width
alp = 0.6 # alpha

# Matplotlib Settings
rcParams.update({'figure.autolayout': True})
rcParams['axes.formatter.limits'] = [-3,3]

# CHANGE ME
# Loading dataset (Load from one folder up)
data_location = '../'
low_path = 'low/pltstate*'
medium_path = 'med/pltstate*'
high_path = 'high/pltstate*'

#Loading dataset
low = yt.load(data_location+low_path)
medium = yt.load(data_location+medium_path)
high = yt.load(data_location+high_path)

# Resolution
low_layers = 2
med_layers = 3
high_layers = 4
low_res_base = float(low[0].domain_right_edge[0])/float(low[0].domain_dimensions[0])
low_res_highest = low_res_base/(2**low_layers)
med_res_highest = low_res_base/(2**med_layers)
high_res_highest = low_res_base/(2**high_layers)

# Program Parameters 
center =  (low[0].domain_right_edge/2.)[0]
adjusted_right = int(center)*2 - 8
convergence_test_x_point = float(int(center)-8. )
convergence_test_y_point = float(int(center))
convergence_test_z_point = float(int(center))
scheme = 1

c = [convergence_test_x_point- (med_res_highest/2.),convergence_test_y_point- (med_res_highest/2.),convergence_test_z_point- (med_res_highest/2.)]

# Dimension of interpolation array
dim = np.array([10,10,10])

# Define an empty storage dictionary for collecting information
# in parallel through processing
storage_low = {}
storage_medium={}
storage_high = {}

for sto, i in low.piter(storage=storage_low):

	# All Data
	ptn = i.point(c)
	rho_basic = float(ptn["rho"][0])

	le_x = convergence_test_x_point - 5.*low_res_highest
	le_yz =  convergence_test_y_point - 5.*low_res_highest
	low_covering_grid = i.covering_grid(level=low_layers, left_edge=[le_x,le_yz,le_yz],dims=dim)

	xdata = np.array(low_covering_grid['x'][:,0,0])
	ydata = np.array(low_covering_grid['y'][0,:,0])
	zdata = np.array(low_covering_grid['z'][0,0,:])

	x0 = xdata[0]
	y0 = ydata[0]
	z0 = zdata[0]

	xcoord = (convergence_test_x_point-x0-(med_res_highest/2.))/low_res_highest
	ycoord = (convergence_test_y_point-y0-(med_res_highest/2.))/low_res_highest
	zcoord = (convergence_test_z_point-z0-(med_res_highest/2.))/low_res_highest

	#print xdata
	#print ydata
	#print zdata
	#print c
	#print xcoord
	#print ycoord
	#print xcoord

	data = np.array(low_covering_grid['rho'])
	var_low = ndimage.map_coordinates(data, [[xcoord], [ycoord], [zcoord]], order=scheme)[0]

	array = [i.current_time,rho_basic,var_low]

	
	sto.result = array
	sto.result_id = str(i)

for sto, i in medium.piter(storage=storage_medium):

	ptn = i.point(c)
	rho_basic = float(ptn["rho"][0])

	le_x = convergence_test_x_point - 5.*med_res_highest
	le_yz =  convergence_test_y_point - 5.*med_res_highest
	med_covering_grid = i.covering_grid(level=med_layers, left_edge=[le_x,le_yz,le_yz],dims=dim)

	xdata = np.array(med_covering_grid['x'][:,0,0])
	ydata = np.array(med_covering_grid['y'][0,:,0])
	zdata = np.array(med_covering_grid['z'][0,0,:])

	x0 = xdata[0]
	y0 = ydata[0]
	z0 = zdata[0]

	xcoord = (convergence_test_x_point-x0-(med_res_highest/2.))/med_res_highest
	ycoord = (convergence_test_y_point-y0-(med_res_highest/2.))/med_res_highest
	zcoord = (convergence_test_z_point-z0-(med_res_highest/2.))/med_res_highest

	#print xdata
	#print ydata
	#print zdata
	#print c
	#print xcoord
	#print ycoord
	#print xcoord

	data = np.array(med_covering_grid['rho'])
	var_med = ndimage.map_coordinates(data, [[xcoord], [ycoord], [zcoord]], order=scheme)[0]

	array = [i.current_time,rho_basic,var_med]

	sto.result = array
	sto.result_id = str(i)
	
for sto, i in high.piter(storage=storage_high):

	ptn = i.point(c)
	rho_basic = float(ptn["rho"][0])

	le_x = convergence_test_x_point - 5.*high_res_highest
	le_yz =  convergence_test_y_point - 5.*high_res_highest
	high_covering_grid = i.covering_grid(level=high_layers, left_edge=[le_x,le_yz,le_yz],dims=dim)

	xdata = np.array(high_covering_grid['x'][:,0,0])
	ydata = np.array(high_covering_grid['y'][0,:,0])
	zdata = np.array(high_covering_grid['z'][0,0,:])

	x0 = xdata[0]
	y0 = ydata[0]
	z0 = zdata[0]

	xcoord = (convergence_test_x_point-x0-(med_res_highest/2.))/high_res_highest
	ycoord = (convergence_test_y_point-y0-(med_res_highest/2.))/high_res_highest
	zcoord = (convergence_test_z_point-z0-(med_res_highest/2.))/high_res_highest

	#print xdata
	#print ydata
	#print zdata
	#print c
	#print xcoord
	#print ycoord
	#print xcoord

	data = np.array(high_covering_grid['rho'])
	var_med = ndimage.map_coordinates(data, [[xcoord], [ycoord], [zcoord]], order=scheme)[0]

	array = [i.current_time,rho_basic,var_med]

	sto.result = array
	sto.result_id = str(i)

# Sort data and plot graphs on a single rank
if yt.is_root():
	timedata_low = []
	low_rho_basic = []
	low_rho_ndimage=[]
	timedata_med = []
	med_rho_basic = []
	med_rho_ndimage=[]
	timedata_high = []
	high_rho_basic = []
	high_rho_ndimage = []
	convergence = []
	convergence_time = []

	# Low Results
	for L in sorted(storage_low.items()):
		timedata_low.append(L[1][0])
		low_rho_basic.append(L[1][1])
		low_rho_ndimage.append(L[1][2])

	# Medium Results
	for L in sorted(storage_medium.items()):
		timedata_med.append(L[1][0])
		med_rho_basic.append(L[1][1])
		med_rho_ndimage.append(L[1][2])

	# High Results
	for L in sorted(storage_high.items()):
		timedata_high.append(L[1][0])
		high_rho_basic.append(L[1][1])
		high_rho_ndimage.append(L[1][2])
		
	# Calculate Convergence
	for i in range (0, min(len(timedata_low),len(timedata_med), len(timedata_high))):
		convergence_time.append(timedata_high[i])
		temp = abs(low_rho_ndimage[i]-med_rho_ndimage[i])/abs(med_rho_ndimage[i]-high_rho_ndimage[i])
		convergence.append(np.log2(temp))
		
	convergence_mean = []
	convergence_plus = []
	convergence_minus = []
	for i in range (0, min(len(timedata_low),len(timedata_med), len(timedata_high))):
		convergence_mean.append(np.mean(convergence))
		convergence_plus.append(np.mean(convergence)+np.std(convergence))
		convergence_minus.append(np.mean(convergence)-np.std(convergence))

	# rho_basic
	plt.figure(1)
	plt.plot(timedata_low,low_rho_basic,label='Basic Low')
	plt.plot(timedata_med,med_rho_basic,label='Basic Medium')
	plt.plot(timedata_high,high_rho_basic,label='Basic High')
	plt.xlabel('Time $[1/m]$')
	#plt.xlim(0,220)
	#plt.ylim(0,0.000001)
	plt.grid()
	plt.legend()
	plt.savefig('rho_basic.png')
	plt.close()

	# rho_ndimage
	plt.figure(2)
	plt.plot(timedata_low,low_rho_ndimage,label='ndimage Low')
	plt.plot(timedata_med,med_rho_ndimage,label='ndimage Medium')
	plt.plot(timedata_high,high_rho_ndimage,label='ndimage High')
	plt.xlabel('Time $[1/m]$')
	#plt.xlim(0,200)
	#plt.ylim(0,0.000001)
	plt.grid()
	plt.legend()
	plt.savefig('rho_ndimage_'+str(scheme)+'.png')
	plt.close()
	
	# Convergence
	plt.figure(2)
	plt.plot(convergence_time,convergence)
	plt.xlabel('Time $[1/m]$')
	plt.ylim(0,4.3)
	plt.grid()
	#plt.legend()
	plt.savefig('convergence_'+str(scheme)+'.png')
	plt.close()
	
	# Convergence
	plt.figure(2)
	plt.plot(convergence_time,convergence)
	plt.plot(convergence_time,convergence_mean,color='orange')
	plt.plot(convergence_time,convergence_plus,linestyle='--',alpha=0.5,color='orange')
	plt.plot(convergence_time,convergence_minus,linestyle='--',alpha=0.5,color='orange')
	plt.xlabel('Time $[1/m]$')
	plt.ylim(0,4.3)
	plt.grid()
	#plt.legend()
	plt.savefig('convergence_'+str(scheme)+'_average.png')
	plt.close()

#np.savetxt('time.out',timedata)_low
