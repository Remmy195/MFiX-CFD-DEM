import os, glob
import numpy as np
from glob import glob

basepath = os.getcwd()
j_slices = 5
dirs = ['Phi0p1','Phi0p15','Phi0p2']
VFcut_min = [0.65, 0.625, 0.60] #Average of (0.4,0.9), (0.4,0.85), (0.4,0.8)
VFcut_max = [0.95, 0.925, 0.90] #Average of (0.9,1.0), (0.85,1.0), (0.8,1.0)
i = 0

for entry in dirs:
	fname=str(entry)+str('/VF.dat')
	dir_name = os.path.basename(os.path.dirname(entry))

	with open(fname) as f1:
		lines1 = f1.readlines()[5:]
	
	f1.close()
	
	start=0
	end = len(lines1)
	
	x = []
	y = []
	z = []
	vf = []
	
	for line in lines1[start:end]:
		values = [float(s) for s in line.split()]
		x.append(values[0])
		y.append(values[1])
		z.append(values[2])
		vf.append(values[3])
	
	startk = 0
	endk = int((end + 1)/j_slices - 1)
	
	ycut_settling = []
	ycut_filling = []
	vf_filling = []
		
	for j in range(j_slices):
		for k in range(endk-2,startk,-1):
			if (vf[k] < VFcut_max[i]):	# average of ep_g = ep_g0 and 1.0
				break
		ycut_settling.append(y[k])
		startk = endk + 1
		endk = startk + int((end + 1)/j_slices - 1)

	ycut_settling_mean = np.mean(np.asarray(ycut_settling))

	startk = 0
	endk = int((end + 1)/j_slices - 1)

	for j in range(j_slices):
		for k in range(startk+2,endk):
			if (vf[k] > VFcut_min[i]):	# average of ep_g = ep* and 0.85
				break
		ycut_filling.append(y[k])
		vf_filling.append(vf[startk])
		startk = endk + 1
		endk = startk + int((end + 1)/j_slices - 1)
	
	ycut_filling_mean = np.mean(np.asarray(ycut_filling))
	vf_filling_mean = np.mean(np.asarray(vf_filling))
	i = i + 1				# to update min. anad max. cut-off
	
	print (dir_name, ycut_settling_mean,ycut_filling_mean,vf_filling_mean)

