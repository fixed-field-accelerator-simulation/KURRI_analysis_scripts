#!/usr/bin/python
from __future__ import division
import pylab as plt
import numpy as np
import sys
import analysedata_150619 as adata
import math

plt.rcParams.update({'axes.labelsize':14})
plt.rcParams.update({'xtick.labelsize':14})
plt.rcParams.update({'ytick.labelsize':14})
plt.rcParams['legend.fontsize']=14

PROTON_MASS = 938.27203e6
SPEED_OF_LIGHT = 299792458

read_svk_file = True
if read_svk_file:
	svkfile = "../svk20.dat"
	ff = open(svkfile, 'r')
	
	time_limit = 17e-3
	
	time_svk = []
	ke_svk = []
	rf_freq_svk = []
	for indf ,line in enumerate(ff):
		
		lspl = line.split()
		
		time_svk.append(float(lspl[1]))
		ke_svk.append(float(lspl[2]))
		rf_freq_svk.append(float(lspl[3]))
		
		if time_svk[-1] > time_limit:
			break
			
def ke_lookup(time):
	""" Find kinetic energy at specified time in svk file """
	i = 0 
	ke = None
	for t in time_svk:
		if t > time:
			ke = ke_svk[i]	
			break	
		i = i + 1
	
	return ke


probenames=["F01","F02","F03","F05","F12"] #2015

basedir_l = ["../..//KURRI_analysis_scripts/qinbin_proc/", ""]

probe_pos_l = [[1,2,3,5,12],[1,5]]

dir_data_l = ["2015-06-24/700A/", "180921/"]#900A data (June 2015)

probe_color =  {'F01': 'b', 'F02': 'k', 'F03':'m','F05':'g','F12':'c'}

style_l = ['-','o--']
offset_l = [0, 0.0]

timeaxis = False

index_d = 0
for basedir, probe_pos, dir_data, style, offset in zip(basedir_l, probe_pos_l, dir_data_l, style_l, offset_l):
	probe_names = []
	for p in probe_pos:
		if p < 10:
			probe_names.append("F0"+str(p))
		else:
			probe_names.append("F"+str(p))

	files = ["QinBin_"+nam+".txt" for nam in probe_names]	
	data_all = adata.readset(basedir+dir_data, files)

	pline = []
	for data_probe, probe in zip(data_all, probe_names):
	
		time = data_probe[1]
		
		ke_a = np.array([ke_lookup(t) for t in 1e-6*data_probe[1]])
		te_a = ke_a + PROTON_MASS
		p_mev = 1e-6*np.array([np.sqrt(te**2 - PROTON_MASS**2) for te in te_a])
	
		if timeaxis:
			xax = 1e-3*time
		else:
			xax = p_mev
			
		pl = plt.plot(xax - offset, data_probe[0],probe_color[probe[:3]]+style, linewidth=index_d+1, label=probe)
		pline.append(pl)
		
	index_d = index_d + 1

plt.xlim(xmax=350)
plt.legend(loc = "upper left")

if timeaxis:
	plt.xlabel('loss time [ms]')
else:
	plt.xlabel('p [MeV/c]')
	
plt.ylabel('closed orbit [m]')
plt.savefig('co_compare_paxis')
plt.show()



