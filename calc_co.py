from __future__ import division
import pylab as plt
import numpy as np
import sys
import os

probe_offset_map =  {'F1': 4705.5, 'F5': 4708}
PROTON_MASS = 938.27203e6
SPEED_OF_LIGHT = 299792458


dir_data = "180921/"
files = os.listdir(dir_data)

qinbin_files = []
probe_id_l = []
for fl in files:
	if "QinBin" in fl:
		qinbin_files.append(fl)
		if "F1" in fl:
			probe_id_l.append("F1")
		else:
			probe_id_l.append("F5")

print "found following files in directory ",dir_data
print "qinbin_files ",qinbin_files
print "probe_id ",probe_id_l


read_svk_file = True
if read_svk_file:
	svkfile = "../svk20.dat"
	ff = open(svkfile, 'r')
	
	time_limit = 17e-3
	
	time_svk = []
	ke_svk = []
	for indf ,line in enumerate(ff):
		
		lspl = line.split()
		
		time_svk.append(float(lspl[1]))
		ke_svk.append(float(lspl[2]))
		
		if time_svk[-1] > time_limit:
			break
	
	
	#plt.plot(time_svk, ke_svk, 'k-', linestyle='-', linewidth=2)
	#plt.xlim(time_svk[0], time_svk[-1])
	#plt.show()
	
	
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

F1_relpos_mm = np.array([100, 150, 200])
F1_times = np.array([0.00828123995, 0.00994051994, 0.0118635999])


F5_relpos_mm = np.array([200, 250, 300])
F5_times = np.array([0.00833383995, 0.0101358799, 0.0119313199])

col_l = ['b','g']
ip = 0
for fname, probe in zip(qinbin_files, probe_id_l):
	
	print dir_data+fname
	ff = open(dir_data+fname,'r')
	
	pos = []
	times_ms = []
	il = 0
	for line in ff:
		if il > 0:
			lspl = line.split()
		
			pos.append(float(lspl[0]))
			times_ms.append(1e-3*float(lspl[1]))
		
		il = il + 1

	
	ke = np.array([ke_lookup(1e-3*t) for t in times_ms])
	print ke
	te = ke + PROTON_MASS
	p = 1e-6*np.sqrt(te**2 - PROTON_MASS**2)

	plt.plot(p, pos, col_l[ip]+'o-', linewidth=2, label=probe)
	
	plt.xlabel('p [MeV/c]')
	plt.ylabel('Closed orbit [m]')
	ip = ip + 1
	
plt.legend(loc = 'upper left')
plt.savefig('closedorb')
plt.show()