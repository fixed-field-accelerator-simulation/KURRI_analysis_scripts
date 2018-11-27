from __future__ import division
import pylab as plt
import numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import quadrature
import sys
import os
import matplotlib.cm as cm
import analysedata_180919 as adata
import tomo_defs_181101 as tdefs
import math
from subprocess import call
import time 
import matplotlib.animation as animation
from subprocess import call
from scipy import ndimage
from scipy.interpolate import griddata




#set up formatting for movie files
#Writer = animation.writers['ffmpeg']
#writer = Wriiter(fps=15, metadata=dict(artist='Me'), bitrate=1800)

read_rf_waveform = True
read_phis_file = True #read file containnig phis timings


add_time_offset = False

#calculate offset in time w.r.t RF lookup table. 
calc_rf_time_offset = False #special run to calculate offset of timebase from RF table time
if calc_rf_time_offset:
	tomo_times = np.array([-0.05e-3])
	nturns_sort = 10000 #10000 for full intensity file

	write_offset_file = False
	show_time_plot = True
	
	intensity_calc = True
	write_intensity_file = False

	add_time_offset = False
	use_file_offset = False
else:
	tomo_times = np.array([0.8e-3]) #this is with respect to offset if add_time_offset is True
	use_file_offset = True
	nturns_sort = 200

	intensity_calc = False

	write_offset_file = False
	write_intensity_file = False

#intensity_calc = False
#if intensity_calc:
#	tomo_times = np.array([-0.05e-3])	
#	nturns_sort = 10000


#tomo_times = np.array([0.08e-3])

dirname = "../data/2018/9/20/machida-san_pattern/proc/"
#dirname = "../data/2018/9/20/proc/"
#dirname = "../data/2018/9/21/machida_san_foil_escape/proc/"
#dirname = "../data/2018/9/25/machida-san_2018_09_25/proc/"
#dirname = "../data/2018/9/27/JB_2018_09_27/proc/"

#dirname = "../data/2018/9/21/machida_san_foil_escape/" # "../data/2018/9/20/"

filter_string = '' #look for csv files with this string in file name
phi_s_deg = 20 #set synchronous phase
turns_syncfrac = 0.5 #fraction of synchrotron oscillation to include
recon_start = 1 #reconstruct phase space at this turn first
recon_step_frac = 0.1 #0.1 #reconstruction steps in terms of fraction of sychrotron period

#recon_step = 25 #step size terms of turn number between recon_start and recon_stop
#recon_stop = 100 #reconstruct phase space at this turn last
animate = False

ncheck = 0 #number of synch point settings to check on either side of that set

idelay =  -13 #-13 #-23 for flattop 21/9	 #-22#-22 #-21	 #-20, -26 #apply this delay to rf zerocrossing time
#if idelay=None, use of the following methods to establish it
idelay_corr_phi_s = False #if True, correct above idelay for nonzero phi_s


#...................................
#these methods are applied if idelay = None
shift_data_sym = True #shift data by looking for symmetry point
shift_data_winmin = False #shift minimum of moving window average where window is out-of-bucket duration
shift_data_max = False #shifts data so that maximum is at centre 
shift_data_min = False #shifts data so that minimum is at edge

if idelay !=None:
	print "applying the following shift to rf waveform to find synchronous time ",idelay

else:
	shift_names = ["window min","Shift max", "Shift min", "symmetry search"]
	shift_switches = [shift_data_winmin, shift_data_max, shift_data_min, shift_data_sym]
	true_count = shift_switches.count(True)
	if true_count > 1:
		print "ensure just one shift method is selected"
		sys.exit()
	elif true_count == 0:
		print "no data shifting"

	if shift_data_winmin and phi_s_deg == 0:
		print "window method requires non-zero phi_s"
		sys.exit()
	print "using the following method to find synchronous time: ", shift_names[shift_switches.index(True)]


show_data_switch = False #show raw data
show_time_fit = False #show fit used to calculate ToF
mountain_range_plot = False
plot_data_shift = False

print "phi_s set to ",phi_s_deg, " degrees"

npfac = 1.0 #interpolate onto time axis with npfac times the number of points as the data
interpolate_data = True #if False interpolate data onto time axis which is determined by rf waveform

subtract_baseline = False #subtract no beam baseline signal

tomo_outdir = "tomo_output"
files = os.listdir(os.curdir)
if tomo_outdir not in files:
	print "create output directory "+tomo_outdir
	os.mkdir(tomo_outdir)

print "reconstruct at following times [ms] ", 1e3*tomo_times
print "phi_s set to [deg] ",phi_s_deg
if phi_s_deg == 0:
	flattop = True
	flattop_time = 10e-3#0.3e-3 #estimate
	print "flattop time is [ms] ", 1e3*flattop_time
else:
	flattop = False

def sorted_ls(path):
    mtime = lambda f: os.stat(os.path.join(path, f)).st_mtime
    return list(sorted(os.listdir(path), key=mtime))

def first_chars(x):
	return(x[:40])

sort_by_name = True
if sort_by_name:
	files = sorted(os.listdir(dirname), key = first_chars)
else:
	files = sorted_ls(dirname)

print "the following csv files were found"
#files =  os.listdir(dirname)
csv_files = []
icsv = 0 
for i1, f1 in enumerate(files):
	if 'csv' in f1 and filter_string in f1:
		print "index ",icsv," file ", f1
		csv_files.append(f1)
		icsv = icsv + 1

rawin1 =  raw_input("select bunch monitor data file ")
try:
	i_file_sel = int(rawin1)
except:
	print "input is not integer"
	sys.exit()

if read_rf_waveform:
	
	rf_file_guess = csv_files[i_file_sel].replace("wbeam","rf")

	if rf_file_guess in csv_files:
		i_rf_file_sel = csv_files.index(rf_file_guess)
		print "found following matching rf file ",rf_file_guess
		
	else:
		rawin2 = raw_input("select RF data file (enter to select "+ str(i_file_sel + 1) + ")  ")
		
		try:
			i_rf_file_sel = int(rawin2)
		except:
			if rawin2 == '':
				i_rf_file_sel = i_file_sel + 1
			else:
				i_rf_file_sel = None
				print "no valid RF file selected"
				

if i_rf_file_sel != None:
	print "selected rf file index ",i_rf_file_sel


if subtract_baseline:
	rawin3 = raw_input("select baseline (no beam) file  ")

	try:
		i_baseline_file_sel = int(rawin3)
	except:
		print "no valid baseline file selected"



if phi_s_deg == 0:
	binslowerignore = 0 #add zeroes before profile (earlier time)
	binsupperignore = 0 #add zeroes after profile (later time)
else:
	pass
	#binslowerignore = int(npfac*35) #add zeroes before profile (earlier time)
	#binsupperignore = int(npfac*19) #add zeroes after profile (later time)

#longitudinal tomography based on bunch monitor data.
#dirname = "../data/2015/6/24/700A/"



#svk file from Tom. rf burst starts -0.5ms, i.e. before injection.
#format is index, time[s], kinetic energy [eV], rf freq [Hz], mean radius [m]
read_svk_file = True
if read_svk_file:
	svkfile = "../sep18_exp/svk20.dat"
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
	
	
	#plt.plot(time_svk, ke_svk, 'k-', linestyle='-', linewidth=2)
	#plt.xlim(time_svk[0], time_svk[-1])
	#plt.show()
	
	
def ke_lookup(time, toffset = 0):
	""" Find kinetic energy at specified time in svk file """
	i = 0 
	ke = None
	
	if time_svk[0] > time - toffset:
		print "time less than minimum in table : ",time-toffset
		print "use initial energy in table"
		return ke_svk[0]
		
	for t in time_svk:
		if t > time - toffset:
			ke = ke_svk[i]	
			break	
		i = i + 1
	
	return ke
	
def f_rf_lookup(time, toffset = 0):
	""" Find rf frequency at specified time in svk file """
	i = 0 
	rf_freq = None
	for t in time_svk:
		if t > time - toffset:
			rf_freq = rf_freq_svk[i]	
			break	
		i = i + 1
	
	return rf_freq

def phis_file_read(fname):
	"""Read phis file"""
	ff = open(dirname+fname,'r')
	t_phis_file = []
	phis_file = []
	for lraw in ff:
		line = lraw.split()
		t_phis_file.append(float(line[0]))
		phis_file.append(float(line[1]))
	ff.close()

	return t_phis_file, phis_file
	
		
pre2018_data = False 
if pre2018_data:
	channel_id = ["F01"] #identifier for probes to be analysed 

	#indices_select used to select particular data files associated with each probe
	indices_select = [15] #process all files
	#indices_select = range(1,5) #select files with index 1 to 4 for each probe
	#indices_selec
	t = [0,1,2] #select first three files for each probe

	#select channel files
	ch_files_sel = adata.select_filenames(dirname, indices_select, channel_id)
	print "selected data file ",ch_files_sel
	tdat_chf, data_chf = tdefs.read_scope(dirname, ch_files_sel[0][0])

else:
	fname_sig = csv_files[i_file_sel]
	#fname = "20180920-004.csv"#"20180920-004.csv"
	#read scope data
	tdat_chf, data_chf = tdefs.read_scope(dirname, fname_sig)

	if read_rf_waveform and i_rf_file_sel != None:
		fname_rf = csv_files[i_rf_file_sel]
		tdat_rf, sig_rf = tdefs.read_scope(dirname, fname_rf)
		tdat_rf_ms = 1e3*np.array(tdat_rf)
		sig_rf_sc = 0.01*np.array(sig_rf)

#plt.plot(tdat_chf, data_chf)
#plt.show()
#sys.exit()

t_phis_file = None
if read_phis_file:
	if fname_sig[:4] == 'case':
		
		fname_phis = 'case'+fname_sig[4]+'_phis.txt'
		print "look for phis file named ",fname_phis
		if fname_phis in files:
			print "read ",fname_phis


			t_phis_file, phis_file = phis_file_read(fname_phis)

			print "phis_file ",phis_file

			#plt.plot(t_phis_file, phis_file, 'ko-')
			#plt.show()
	

		else:
			print "phis_file not found"
			sys.exit()



time_interval = tdat_chf[2] - tdat_chf[1]
time_interval_ns = 1e9*time_interval
tdat_chf_ms = 1e3*np.array(tdat_chf)



if subtract_baseline:
	fname_base = csv_files[i_baseline_file_sel]
	#fname = "20180920-004.csv"#"20180920-004.csv"
	#read scope data
	tdat_base, data_base = tdefs.read_scope(dirname, fname_base)

	dirproc = dirname + 'proc/'


	fdiff = open(dirproc+'data_less_baseline.txt','w')
	for tb, data_sig, data_base in zip(tdat_base, data_chf, data_base):
		data_diff = data_sig - data_base
		print >>fdiff, repr(tb),",", repr(data_diff)
		


i1 = 0
for i1, t1 in enumerate(tdat_chf):
	if t1 == 0:
		index_zero = i1
		break
		
tof_guess = 600e-9


def show_data(time_indices, nturns, shift=0):

	print "tdat_rf_ms step ",tdat_rf_ms[1] - tdat_rf_ms[0]
	tpad = 2e-4

	nt = 150

	if nt > len(time_indices):
		nt = len(time_indices) - 1
	tfirst = tdat_rf_ms[time_indices[0]] -tpad
	tlast =  tdat_rf_ms[time_indices[nt]] + tpad

	plt.subplot(211)
	plt.plot(tdat_rf_ms, sig_rf_sc, 'r')
	plt.plot(tdat_chf_ms, data_chf,'k')
	
	#plt.axvline(x=t_hminus, color='r')
	
	for lower_lim in tomo_times_data_l:
		plt.axvline(x=tfirst, color='r')
	
	#plt.xlim(0.5,1.2)
	#plt.xlim(0, tdat_chf_ms[-1])

	#plt.ylim(-0.2, 0.2)
	plt.ylim(-0.05, 0.05)
	plt.ylabel('bunch monitor (all)')
	plt.axvline(x=tdat_rf_ms[time_indices[0]])
	

	plt.subplot(212)
	plt.plot(tdat_chf_ms + 1e3*time_interval*shift, data_chf,'k')
	if read_rf_waveform:
		plt.plot(tdat_rf_ms , sig_rf_sc, 'r')
	plt.xlim(tfirst, tlast)
	plt.ylim(-0.1, 0.1)

	plt.axvline(x=tdat_rf_ms[time_indices[0]])
	plt.ylabel('bunch monitor (zoom)')
	plt.xlabel('time [ms]')
	plt.savefig('bunchmonitor_data')
	plt.show()

	sys.exit()
	


def sort_data(istart, tof_fit_a, time_indices, nturns, n_slices, interpolate_data = True, nshift = 0, show_mrange=True, isym=0):
	"""Divide bunch monitor data into individual turns. There is the option to interpolate onto regular time points. 
	istart: the index of the first time point
	tof_fit_a: times at the start of each turn
	time_indices: indices associated with tof_fit_a
	nturns: limit on the number of turns to sort
	n_slices: number of time points per turn
	interpolate_data: interpolate onto a regular grid delimited by tof_fit_a if True
	nshift: shift time by this number of points to deal with offset of rf signal (if shift_time is True)
	
	
"""
	
	#iend = istart + nturns*points_per_turn
	 
	
	#tstart = tof_fit_a[0] #tdat_chf[istart] #t_dat_all[0]
	
	#plt.plot(t_dat_all, sig_dat_all)
	#plt.show()
	#sys.exit()
	#use ToF to define time.
	#t_dat_adj = tstart + np.linspace(0, nturns*tof, nturns*n_slices)

	if npfac > 1 and not interpolate_data:
		print "must interpolate data if npfac != 1"
		sys.exit()
	

	#interpolate onto time base Trev/n_slices
	t_dat_adj = []
	for it in range(1, nturns+1):
		time1 = tof_fit_a[it-1]
		time2 = tof_fit_a[it]
		t_oneturn_data = np.linspace(time1, time2 ,n_slices) #fitdiff_a[it]
		t_dat_adj  = np.concatenate((t_dat_adj, t_oneturn_data))


	t_dat_adj = t_dat_adj + time_interval*nshift

	t_adj_last = t_dat_adj[-1]
	i1 = 0
	for i1, t1 in enumerate(tdat_chf):
		if t1 >= t_adj_last:
			index_t_last = i1
			break	

	
	t_dat_orig = np.array(tdat_chf[istart:index_t_last])
	sig_dat_orig = np.array(data_chf[istart:index_t_last])
		
	sig_interp = np.interp(t_dat_adj, t_dat_orig, sig_dat_orig)
	
	sig_dat_split = np.split(sig_interp, nturns)
	sig_interp = np.interp(t_dat_adj, t_dat_orig, sig_dat_orig)


	if interpolate_data:
		t_turn_data = np.linspace(0, tof, n_slices)
		sig_dat_split = np.split(sig_interp, nturns)
	else:
		print "interpolate_data = False option needs work"
		sys.exit()
		sig_dat_split = []
		for it in range(1, nturns+1):
			sig_turn = data_chf[time_indices[it-1]:time_indices[it]]
			print "len sig_turn ",len(sig_turn)
			sig_dat_split.append(sig_turn)
		#sys.exit()
		#sig_dat_split = np.array(sig_dat_split)
	
		t_turn_data = t_dat_orig[:n_slices]

		#sig_dat = np.array(data_chf[istart:istart + nturns*n_slices])
		#sig_dat_split = np.split(sig_dat, nturns)
		

		#plt.plot(t_dat_adj, sig_interp,'r-')
		#plt.plot(t_dat_orig, sig_dat_orig,'ko')
		#plt.axvline(x=t_dat_all[-1])
		#plt.show()
		#sys.exit()
	
	
	if mountain_range_plot and show_mrange:
		fig = plt.figure(figsize=(8, 8), facecolor='w')
		ax = plt.subplot(111, frameon=False)

		
	it = 1
	data_win_sig_sel = []
	for data_win in sig_dat_split[:nturns]:
		#print "data_win ",data_win
		data_win_sig =  max(data_win) - np.array(data_win)
		data_win_sig_sel.append(data_win_sig)
		
		if mountain_range_plot and show_mrange:
			lw = 0.5 - (it / 3000)
			ax.plot(range(1,n_slices+1), 700*data_win_sig + it,'k', lw=lw)
		it = it + 1*1
	
	if mountain_range_plot and show_mrange:		
		ax.set_ylabel('turn number')
		plt.axvline(x=binslowerignore+0.5)
		plt.axvline(x=n_slices - binsupperignore+0.5)
		if nshift == 0:
			plt.axvline(x=isym, color='r', linestyle='-')
		plt.axvline(x=sync_bin, color='r', linestyle='--')
		plt.savefig('mountainrange')
		plt.show()



	prof_data_a = np.array(data_win_sig_sel)	
	
	return prof_data_a, t_turn_data


def phase_extrema(phi_s):
	"""For a given phi_s find max and min phases along separatrix"""
	#calculate phase extrema along separtrix
	phi_l1 = math.pi - phi_s

	if phi_s == 0.0:
		phi_l2 = -math.pi
	else:
		#guess turning point 
		#may break down if phi_s is close to 90 degrees
		if phi_s < 50*math.pi/180:
			phi_g = 0
		elif phi_s >= 50*math.pi/180 and phi_s < 0.5*math.pi:
			phi_g = 80*math.pi/180
		else:
			phi_g = math.pi

		phi_l2 =  tdefs.turning_point(phi_s, phi_guess = phi_g)
		
	return phi_l1, phi_l2
	
def find_symmetry_point_shape(data):
	"""Look for point where data on either side is most similar"""

	ld = len(data)
	ldh = int(0.5*ld)

	if ld%2 == 1:
		oddcorr = 1
	else:
		oddcorr = 0

	
	chisq_l = []
	for d in range(ldh):
		#indh1 = range(d, ldh+d)
		#indh2 = range(ldh+d, ld+d)
		indh1 = range(ldh+d, ld+d - oddcorr)
		indh2 = range(d, ldh+d)
		half1 = data.take(indh1, mode='wrap')
		half2 = data.take(indh2, mode='wrap')

		half2rev = half2[::-1]		
		halfdiff = half1 - half2rev

		chisq = np.sum(halfdiff**2)
		chisq_l.append(chisq)
					
	chisq_min = min(chisq_l)
	isym = chisq_l.index(chisq_min)
	
	return isym, chisq_min


def find_symmetry_point_integral(data):
	"""Look for point where data on either side is most similar"""

	ld = len(data)
	ldh = int(0.5*ld)

	if ld%2 == 1:
		oddcorr = 1
	else:
		oddcorr = 0

	
	chisq_l = []
	for d in range(ldh):
		#indh1 = range(d, ldh+d)
		#indh2 = range(ldh+d, ld+d)
		indh1 = range(ldh+d, ld+d - oddcorr)
		indh2 = range(d, ldh+d)
		half1 = data.take(indh1, mode='wrap')
		half2 = data.take(indh2, mode='wrap')

		half1_sum = np.sum(half1)
		half2_sum = np.sum(half2)
		halfdiff = abs(half1_sum - half2_sum)

		chisq = halfdiff
		chisq_l.append(chisq)
					
	chisq_min = min(chisq_l)
	isym = chisq_l.index(chisq_min)
	#print isym, chisq_minke_lookup(0.03e-3)

	#plt.plot(chisq_l)
	#plt.show()
	#sys.exit()
	
	return isym, chisq_min


def read_data_file(filen):
	
	f1 = open(filen,'r')
	data_l = []
	for data in f1:
		data_l.append(float(data))
	data_a = np.array(data_l)

	return data_a


def filter_isolated_cells(array, struct):
	""" Return array with completely isolated single cells removed
	:param array: Array with completely isolated single cells
	:param struct: Structure array for generating unique regions
	:return: Array with minimum region size > 1
	"""

	nonzero_indices = np.nonzero(image_data_tr)
	
	nd = len(image_data_tr[0])
	array_ones = np.zeros((nd, nd))
	array_ones[nonzero_indices] = 1

	filtered_array = np.copy(array)
	
	id_regions, num_ids = ndimage.label(array_ones, structure = struct)

	
	id_sizes = np.array(ndimage.sum(array_ones, id_regions, range(num_ids + 1)))
	
	#area_mask = (id_sizes == 1)
	area_mask = (id_sizes <= 3)
	filtered_array[area_mask[id_regions]] = 0
	
	#indices_filtered = np.nonzero(filtered_array)
	
	return filtered_array

#Hamiltonian
H_fn = lambda phi, dp, phi_s, a1, a2: a1*(dp**2) + a2*(np.cos(phi) - np.cos(phi_s) + (phi - phi_s)*np.sin(phi_s))

#set parameters
phi_s = phi_s_deg*math.pi/180 
kindex = 7.6 #
alpha_c = 1/(kindex+1)
gamma_tr = (1/alpha_c)**0.5
V_rf = 4e3 # in [V]
harmonic = 1
bending_radius=0.489396551697 # in m
Rnom = 4.41
circumference = Rnom*2*np.pi
phi_offset = 0 #np.arcsin(bending_radius*circumference*Bdot/V_rf) # measured from aligned focussing phase (0 or pi)

print "alpha_c ",alpha_c

if use_file_offset:
	case_spl = fname_sig.split("_")
	case_num = case_spl[0]
	print "case_num ",case_num
	print "file ",dirname+str(case_num)+"_toffset.txt"
	try:
		
		fo = open(dirname+str(case_num)+"_toffset.txt",'r')
		fo.readline()
		offset_file_ms = 1e3*float(fo.readline())
		
	except:
		print "offset file not found"
		pass
	
	print "offset_file_ms ",offset_file_ms

#calculate kinetic energy, beta, synchrotron period at each time

if add_time_offset:
	toffset_rf = 1e-3*offset_file_ms
else:
	toffset_rf = 0

if not flattop:
	Ekin_a = np.array([ke_lookup(t, toffset_rf) for t in tomo_times])
else:
	ke_flattop = ke_lookup(flattop_time, toffset_rf)

	Ekin_l = []
	for t in tomo_times:
		if t < flattop_time:
			Ekin_l.append(ke_lookup(t, toffset_rf))
		else:
			Ekin_l.append(ke_flattop)
	
	Ekin_a = np.array(Ekin_l)

PROTON_MASS = 938.27203e6
SPEED_OF_LIGHT = 299792458

print "Ekin_a ",Ekin_a 
if Ekin_a[0] != None:
	Etot_a = Ekin_a + PROTON_MASS
	beta_a = adata.ke_to_relativistic_beta(Ekin_a, PROTON_MASS)
	betagam_a = adata.ke_to_relativistic_beta_gamma(Ekin_a, PROTON_MASS)
	gamma_a = betagam_a/beta_a
	eta_a = alpha_c - gamma_a**(-2)
	Qs_a = np.sqrt(np.abs(eta_a) * V_rf / (2 * np.pi * beta_a**2 * Etot_a))
	Ts_a = 1/Qs_a

	print "times ",[t - toffset_rf for t in tomo_times] 
	print "kinetic energies at selected times ",Ekin_a
	print "rel. beta at selected times ",beta_a
	print "eta_a ",eta_a 
	print "Qs ",Qs_a
	print "Ts ",Ts_a

#print "f_rf_lookup at t selected ",f_rf_lookup(0)


#plt.plot(tomo_times, Ts_a,'ko-')
#plt.show()

tomo_times_data_l = []
istart_guess_l = []
for tt in tomo_times:

	print "tt ",tt
	if add_time_offset:
		tt = tt + toffset_rf
	
	print "requested time ",tt

	i = 0
	icrude = None
	for t in tdat_chf:
		if t > tt:
			icrude = i
			break
		i = i + 1
	
	print "icrude ",icrude
	#icrude = index_zero + int(round(tt/time_interval))
	t_crude = tdat_chf[icrude]
	print "t_crude ",t_crude

	istart_guess_l.append(icrude)
	tomo_times_data_l.append(t_crude)



points_per_turn = int(round(tof_guess/time_interval)) 


print "specified initial times ",tomo_times




if not read_rf_waveform:
	#calculate TOF
	tof, time_indices, ttime_a = tdefs.calc_tof_bm_data(lower_lim, points_per_turn, nturns_sort+1)
	istart = time_indices[0]
else:
	#tof_m1, indices_m1, ttime_a_m1 = tdefs.calc_tof_bm_data(lower_lim, points_per_turn, nturns_sort+1)
	tof, time_indices, ttime_a = tdefs.calc_tof_rf_waveform(tdat_rf, sig_rf, nturns_sort+1, istart_guess_l[0])
	istart = time_indices[0]

	if len(ttime_a) < nturns_sort+1:					
		print "nturns_fit limited by data to ",nturns_fit

	print "tof rf ",tof


	tof_a = [ttime_a[i] - ttime_a[i-1] for i in range(1, len(ttime_a))]
	freq_a = np.array([1/t for t in tof_a])
	freq_MHz_a = 1e-6*freq_a

	
	if calc_rf_time_offset:

		print "calc_rf_time_offset section"
		
		if t_phis_file != None:
			t1 = 1e-3*t_phis_file[0]
			if len(t_phis_file) > 1:
				t2 = 1e-3*t_phis_file[1]
			else:
				t2 = None
		 
	 

		fzero_ideal = f_rf_lookup(0) #frequency at zero time in table
		f_zero_ideal_MHz = 1e-6*fzero_ideal
		print "ideal t=0 freq ",fzero_ideal

		
		npfit = 3
		nturns_fit = 2000
		ffit_coef = np.polyfit(ttime_a[1:nturns_fit], freq_a[0:nturns_fit-1], npfit)
		ffit = np.poly1d(ffit_coef)

		offset =  (fzero_ideal - ffit_coef[-1])/ffit_coef[npfit-1]
		print "offset [ms]",1e3*offset

		t_rf_table = np.linspace(0, ttime_a[-1], 100)

		if write_offset_file:
			#construct offset file name
			fspl = csv_files[i_file_sel].split("_")
			foffset = ''
			for str1 in fspl[:-1]:
				foffset = foffset + str1+"_"
			offset_file = foffset + "toffset.txt"
			print "offset_file ",offset_file
			
			fo = file(dirname+offset_file,'w')
			print >>fo, offset
		
		#for t in t_phis_file:
		#	t_s = 1e-3*t
		#	if t_s < ttime_a[-1]:
		#		plt.axvline(x=offset + t_s, color='gray')	 


		if show_time_plot:
			plt.plot(t_rf_table, [1e-6*f_rf_lookup(t) for t in t_rf_table],'k-',label='lookup table')
			
			plt.plot(ttime_a[1:], freq_MHz_a, 'b.',label='data (1/ToF)')
			plt.plot(ttime_a[:nturns_fit], 1e-6*ffit(ttime_a[:nturns_fit]),'r', label='data fit')

			
			
			plt.xlabel('time [s]')
			plt.ylabel('frequency [MHz]')
			plt.axvline(x= offset, color='r',label='calc. offset')
			
			plt.legend()
			plt.savefig('frequency_offset')
			plt.show()
			




					
#loop over times, do tomography on each loop
for index_time in range(len(tomo_times)):

	
	#calculate points per turn based on average tof
	points_per_turn = int(round(tof/time_interval)) 
	n_slices = int(npfac*points_per_turn)
	
	
	
	if Ekin_a[0] != None:
		T_rf = tof/harmonic#2*np.pi/(harmonic*omegarev) #rf period
		T_rf_ns = 1e9*T_rf
		omega_rev = 2*math.pi/T_rf
		omega_rf = omega_rev*harmonic


		#n_slices_active = n_slices-(binslowerignore + bip0_mev = 1e-6*np.sqrt(Etot_a[index_time]**2 - PROTON_MASS**2)nsupperignore)
		if idelay != None and idelay_corr_phi_s:
			#calculate further shift in delay to account for phi_s !=0
			idelay_phi_s = int(n_slices*phi_s_deg/360)

			idelay = idelay - idelay_phi_s
			print idelay

		#bucket phase extrema calculation
		if phi_s_deg == 0:
			phi_l1 = math.pi
			phi_l2 = -math.pi
		else:
			phi_l1, phi_l2 = phase_extrema(phi_s)

		print "phi_l1, phi_l2 ",phi_l1, phi_l2
		if phi_l1 < phi_l2:
			phi_lower = phi_l1
			phi_upper = phi_l2label='data'
		else:
			phi_lower = phi_l2
			phi_upper = phi_l1

		print "phi_lower ",phi_lower

		bucket_phase_width = abs(phi_l1 - phi_l2)
		bucket_dur_ns = T_rf_ns*(bucket_phase_width/(2*math.pi))
		outbucket_dur_ns = T_rf_ns - bucket_dur_ns

		tadj = T_rf_ns*(phi_l2/(2*math.pi)) - outbucket_dur_ns 

		synch_index = int(round(n_slices*(phi_s - phi_l2)/(2*math.pi)))
		print "synch_index ",synch_index

		n_slices_active = int(n_slices*bucket_dur_ns/T_rf_ns) #number of slices in bucket
		n_slices_outb = n_slices - n_slices_active
		print "bucket points, out-of-bucket points, bucket points/all points ",n_slices_active, n_slices_outb, n_slices_active/n_slices

		if shift_data_winmin and phi_s_deg == 0:
			print "shift_data_winmin must be False if phi_s_deg is zero"

		if phi_s_deg != 0:
			binslowerignore = n_slices_outb#abs(int(n_slices*phi_lower/(math.pi))) #n_slices_outb
			binsupperignore =   n_slices - (binslowerignore + n_slices_active)

		#omega_rf = 2*math.pi*f_rf_lookup(tomo_times[index_time])
		npts = 1000
		a1 = 0.5*harmonic*omega_rf*eta_a[index_time]
		a2 = omega_rf*V_rf/(2*math.pi*(beta_a[index_time]**2)*Etot_a[index_time])
		dp_sep, t_sep_ns = tdefs.separatrix(phi_s, eta_a[index_time], T_rf, harmonic, npts, a1, a2)

		phase_sep = 2*math.pi*(t_sep_ns)/T_rf_ns	
		t_bucket = np.linspace(t_sep_ns[0], t_sep_ns[-1], n_slices)

		#t_data = np.linspace(-0.5*T_rf_ns, 0.5*T_rf_ns, n_slices) #data fills rf period

		t_phi_sync_ns = phi_s*T_rf_ns/(2*math.pi) #time at synchronous phase
		t_phi_l2 = T_rf_ns*(phi_l2/(2*math.pi)) #time at phi_l2

		t_sep_ns = np.insert(t_sep_ns, 0, t_phi_l2)
		dp_sep = np.insert(dp_sep,0,0)
		phase_sep = np.insert(phase_sep, 0, phi_l2) 

		phase_sep_deg = (180/math.pi)*phase_sep
		sync_frac =  (t_phi_sync_ns- t_sep_ns[0])/bucket_dur_ns
		
		print "phis ",phi_s_deg
		print "sync_frac ",sync_frac

		sync_bin = binslowerignore + int(n_slices_active*sync_frac)
		
		print "int(n_slices_active*sync_frac) ",int(n_slices_active*sync_frac)
		print "sync bin ",sync_bin
		
		
		#Hamiltonian on separatrix
		Hsep = H_fn(math.pi - phi_s, 0, phi_s, a1, a2)


	
	prof_data_a, t_turn_data = sort_data(istart, ttime_a, time_indices, nturns_sort, n_slices, interpolate_data, show_mrange = True)

	#tphase1_ms = 0.5
	#tphase2_ms = 0.3
	#print "ke at specified times ",ke_lookup(1e-3*tphase1_ms, toffset_rf), ke_lookup(1e-3*(tphase1_ms + tphase2_ms), toffset_rf)
	#print "ke at 0.18 ms ",ke_lookup(1e-3*0.18, toffset_rf)
	#sys.exit()
	show_imshow = False
	yaxis_time = True
	if show_imshow:
		if yaxis_time:
			tlim1 = 1e3*(ttime_a[0] - toffset_rf)
			tlim2 = 1e3*(ttime_a[-2] - toffset_rf)
			plt.imshow(prof_data_a, origin='lower', aspect='auto', extent = [0, n_slices, tlim1, tlim2])
			plt.ylabel('time')
		else:
			plt.imshow(prof_data_a, origin='lower', aspect='auto')
			plt.ylabel('turn')
		plt.xlabel('time bin')

		for t in t_phis_file:
			t = t
			if t > tlim1 and t < tlim2:
				plt.axhline(y=t, color='gray')

		#if 1e3*ttime_a[-1] > tphase1_ms + tphase2_ms:
		#	plt.axhline(y=tphase1_ms)
		#	plt.axhline(y=tphase1_ms + tphase2_ms)
		
		plt.savefig('imshow_prof0')
		plt.show()
		
	
	
	
	if intensity_calc:
		#calculate intensity and write to file
		
		
		if write_intensity_file:
			#construct offset file name
			fspl = csv_files[i_file_sel].split("_")
			foffset = ''
			for str1 in fspl[:-1]:
				foffset = foffset + str1+"_"
			offset_file = foffset + "intensity.txt"
			print "intensity ",offset_file
			
			fi = file(dirname+offset_file,'w')
			#fi = open('intensity.txt','w')
		
		intensity_raw = []
		fwhm_raw = []
		i1 = 0
		for prof_data, t1 in zip(prof_data_a, ttime_a[:-1]):
			int1 = time_interval_ns*np.trapz(prof_data)
			
			intensity_raw.append(int1)
		
			#calculate FWHM
			diff = max(prof_data) - min(prof_data)
			HM = diff/2
			pos_extremum = prof_data.argmax()
			try:
				nearest_above = (np.abs(prof_data[pos_extremum:-1] - HM)).argmin()
				nearest_below = (np.abs(prof_data[0:pos_extremum] - HM)).argmin()
				
				FWHM = pos_extremum+nearest_above - nearest_below
				fwhm_raw.append(FWHM)
			except:
				FWHM = 0
				fwhm_raw.append(0)
			
			
			#print "nearest below,above, FWHM ",nearest_below, nearest_above, FWHM
			show_intensity_fwhm_example = False
			if show_intensity_fwhm_example:
				if t1 > 0.5e-3 and FWHM*time_interval_ns > 200: 

					#plt.plot(np.abs(prof_data[0:pos_extremum] - HM))
					#plt.show()
					print "t1 ",t1
					print "FWHM points, time [ns] ",FWHM, FWHM*time_interval_ns

					time_ax = time_interval_ns*np.arange(len(prof_data))
					plt.plot(time_ax, prof_data,'k')
					plt.axvline(x=nearest_below*time_interval_ns, color='r')
					plt.axvline(x=(pos_extremum+nearest_above)*time_interval_ns, color='r')
					plt.ylabel('signal [V]')
					plt.xlabel('time [ns]')
					plt.title('integral [nVs] '+str(int1)[:5]+ ' FWHM [ns] '+str(time_interval_ns*FWHM))

					plt.savefig('integral_fwhm.png')								
					plt.show()
					#sys.exit()


			i1 = i1 + 1
			if write_intensity_file:
				print >> fi, repr(t1), repr(int1), repr(FWHM)

		if write_intensity_file:
			fi.close()	

		
		print "time_interval_ns ",time_interval
		
		max_intensity = max(intensity_raw)
		imax = intensity_raw.index(max_intensity)

		t_imax = ttime_a[imax]
		print "intensity maximum at ",t_imax
		print "offset ",t_imax - 5e-5
		#print "frequency at max ",1e-6*ffit(t_imax)
		#print "t_phis_file ", t_phis_file
		if write_offset_file:
			print >>fo, t_imax - 5e-5		
			fo.close()
			

		offset_intensity_ms = 1e3*(t_imax - 5e-5)	

		fig = plt.figure()
		
		#plt.plot(np.arange(1,nturns_sort+1),intensity_raw,'k-')
		#plt.subplot(211)
		#plt.plot(intensity_raw,'k.')
		
		#plt.plot(1e3*np.array(tdat_chf), data_chf,'k-')
		#plt.xlim(1e3*ttime_a[0],1e3*ttime_a[-1])
		#plt.ylabel('bm signal [V]')
		#plt.xlabel('turn')
		fwhm_raw_a = np.array(fwhm_raw)

		print "minimum FWHM ",min(fwhm_raw)*time_interval_ns
		plt.subplot(211)
		plt.plot(1e3*ttime_a[:-1],intensity_raw,'k.')
		plt.ylabel('intensity [nVs]')
		
		plt.ylim(ymin=0)
		plt.xlim(offset_intensity_ms,1e3*ttime_a[-1])
		
		plt.axvline(x=offset_intensity_ms, color='m')

		if t_phis_file != None:
			for t in t_phis_file:
				if t != 0.0:
					plt.axvline(x=t+offset_intensity_ms, color='r')

		plt.subplot(212)
		plt.plot(1e3*ttime_a[:-1],time_interval_ns*fwhm_raw_a,'k.')
		plt.xlabel('time [ms]')
		plt.ylabel('FWHM [nVs]')
		plt.xlim(offset_intensity_ms,1e3*ttime_a[-1])
		plt.ylim(ymin=0, ymax = 350)

		plt.axvline(x=offset_intensity_ms, color='m')

		if t_phis_file != None:
			for t in t_phis_file:
				if t != 0.0:
					plt.axvline(x=t+offset_intensity_ms, color='r')

		#plt.axvline(x=0.18)
		plt.savefig('intensity')
		plt.show()

		#sys.exit()
		
		
		
	#number of turns specified in terms of sync period
	nturns = int(turns_syncfrac*Ts_a[index_time])
	recon_stop = nturns

	
	p0_mev = 1e-6*np.sqrt(Etot_a[index_time]**2 - PROTON_MASS**2)
	bh = ((2*Qs_a)/(harmonic*abs(eta_a[index_time]))) #bucket height
	
	time_str = str(1e3*tomo_times[index_time])
	lower_lim = tomo_times_data_l[index_time]
	
	print "bh ",bh
	print "lower_lim ",lower_lim

	#tof_fit_a = t_zero_cross

	#this is not a real Bdot - it is a parameter used to specify phis in the tomo code
	Bdot = (V_rf*np.sin(phi_s))/(bending_radius*2*math.pi*Rnom)

	show_sep = False
	if show_sep:
		plt.plot(phase_sep_deg, dp_sep, 'ko-')
		plt.axvline(x=180*phi_l1/math.pi)
		plt.axvline(x=180*phi_l2/math.pi)
		plt.show()

	#sort data into turns covering roughly a synchrotron period
	nturns_sync = int(round(Ts_a[index_time]))
	
	

	

	if shift_data_min or shift_data_max or shift_data_winmin or shift_data_sym or idelay != None:
		#prof_data_tr = np.transpose(prof_data_a)
		#data_integ_0 = [sum(d) for d in prof_data_tr]
		
		data_integ_0 = np.sum(prof_data_a, axis=0)
		imax_col = np.argmax(data_integ_0)
		imin_col = np.argmin(data_integ_0)
		#imax_col = data_integ_0.index(max(data_integ_0))		
		#imin_col = data_integ_0.index(min(data_integ_0))

		
		if idelay != None:
			print "use hardcoded idelay ",idelay
		elif shift_data_min:
			if imin_col != 0:
				idelay = int(imin_col/npfac)
			print "imin_col ",imin_col
		elif shift_data_max:
			if imax_col - 0.5*n_slices!= 0:
				idelay = int((imax_col - 0.5*n_slices)/npfac)
			print "imax_col ", imax_col
		elif shift_data_winmin:

			plt.plot(data_integ_0)
			plt.show()
			runmean = tdefs.running_mean_wrap(data_integ_0, n_slices_outb)
			imax_wincol = np.argmax(runmean)
			imin_wincol = np.argmin(runmean)
			idelay = int(imin_wincol/npfac) - int(0.5*n_slices_outb)
			print "calculated delay is ",idelay
		elif shift_data_sym:

			
			ld = len(data_integ_0)

			chisq_sym_turns = []
			
			isym_shape_list = []
			isym_integ_list = []
			
			turnhalfwin = 50
			turn_list = np.arange(turnhalfwin, nturns_sort)

			intensity = []
			for i1 in turn_list:
				data_int = np.sum(prof_data_a[i1-turnhalfwin: i1+turnhalfwin], axis=0)

				isym_shape, _ = find_symmetry_point_shape(data_int*(nturns_sync/i1))
				isym_integ, chisq_min = find_symmetry_point_integral(data_int)
				
				
				binl = range(ld)
				
				chisq_sym_turns.append(chisq_min)

				isym_shape_list.append(isym_shape)
				isym_integ_list.append(isym_integ)

				indsh = np.arange(ld) + isym_integ - sync_bin
				data_integ_shift = data_int.take(indsh, mode='wrap')

				intensity.append(np.trapz(data_integ_shift))

				if i1 in [turn_list[0], turn_list[500], turn_list[1000], turn_list[2000],turn_list[2700]]:
					plt.plot(data_integ_shift)
					#plt.plot(data_integ_shift*(220/i1))

			
			plt.show()

			
			chisq_sym_turns = np.array(chisq_sym_turns)

			nwin = int(nturns_sync)
			range_chisq_win_l = []
			turn_chisq_check = nturns_sync + np.arange(len(turn_list) - nwin)
			for d in range(len(turn_list) - nwin):
				chisq_win = chisq_sym_turns.take(np.arange(d, nwin+d))
				range_chisq_win = max(chisq_win) - min(chisq_win)
				range_chisq_win_l.append(range_chisq_win)
			

			#plt.plot(range_chisq_win_l)
			#plt.show()

			threshold = 10
			ithres = None
			for i, r in enumerate(range_chisq_win_l):
				if r < threshold:
					ithres = i
					break

			if ithres == None:
				print "threshold not reached, used final value of symmetry point"
				ithres = len(range_chisq_win_l) - 1
			else:
				print "converges at turn ",turn_chisq_check[ithres]

			isym_integ_av = sum(isym_integ_list)/len(isym_integ_list)
			isym_final =  int(isym_integ_av)
			idelay =  isym_final - sync_bin #- int(0.5*n_slices))	
		
						
			print "symmetry point av",isym_integ_av
			print "midpoint bin",int(0.5*n_slices)
			print "sync_bin ",sync_bin
			print "synch phase w.r.t midpoint bin ",(sync_bin - int(0.5*n_slices))

			
			print "calculated delay is ",idelay
					
			show_sym_plots = True
			if show_sym_plots:


				indsh = np.arange(ld) + idelay
				data_integ_shift = data_integ_0.take(indsh, mode='wrap')

				
				indsh1 = np.arange(ld) -26
				data_integ_shift_fixed = data_integ_0.take(indsh1, mode='wrap')


				plt.plot(data_integ_0,'k', label='original')
				plt.plot(data_integ_shift,'b', label='shift sym')
				#plt.plot(data_integ_shift_fixed,'b--', label='shift 26')
				plt.axvline(x = isym_final)
				plt.axvline(x = sync_bin,color='gray',linestyle='--',label='sync phase')
				plt.legend()
				plt.savefig('shift_integ')
				plt.show()


				#plt.subplot(311)
				#plt.plot(turn_list, chisq_sym_turns)
				#plt.ylabel('asymmetry')
				#plt.ylim(ymin=0)
				#plt.subplot(211)
				#plt.plot(turn_chisq_check, range_chisq_win_l,'ro')
				#plt.axvline(x=turn_chisq_check[ithres])	
				#plt.ylabel('window range')
				#plt.ylim(ymin=0)
				plt.subplot(211)
				plt.plot(turn_list, isym_integ_list, 'k',label='int')
				plt.plot(turn_list, isym_shape_list, 'r',label='shp')
				plt.axhline(y=isym_integ_av)
				plt.xlabel('turn number')	
				plt.ylabel('sym point')
				#plt.plot(data_int)
				#plt.plot(data_integ_shift)
				plt.legend()
				plt.subplot(212)
				plt.plot(turn_list, intensity, 'ko-')
				plt.xlabel('turn')
				plt.ylabel('intensity')
				plt.savefig('symcalc')
				plt.show()



				#sys.exit()


		print "idelay used ",idelay
		if show_data_switch and index_time == 0:
			nturns_show = nturns_sync
			tind = np.array(time_indices)
			show_data(tind, nturns_show, shift=idelay)

		
		if idelay != 0:

			istart = istart + idelay
			prof_data_a, t_turn_data = sort_data(istart, ttime_a, time_indices, nturns_sort, n_slices, interpolate_data, nshift= idelay)

			

			print "len prof_data_a ",len(prof_data_a)
			row_max = np.max(prof_data_a, axis=1)
			prof_data_norm = prof_data_a/row_max	[:,np.newaxis]

			#prof_data_norm_a = n
			#prof_data_tr = np.transpose(prof_data_a)
			#data_integ_1o = [sum(d) for d in prof_data_tr]
		
			integral_by_turn = [np.trapz(prof) for prof in prof_data_a]

			data_integ_1 = np.sum(prof_data_a, axis=0)
			data_integ_1_norm = np.sum(prof_data_norm, axis = 0)

			#data_integ_0 = data_integ_0/np.max(data_integ_0)
			#data_integ_1 = data_integ_1/np.max(data_integ_1)
			#data_integ_1_norm = data_integ_1_norm/max(data_integ_1_norm)

			show_imshow_2 = False
			if show_imshow_2:
				plt.imshow(prof_data_a, origin='lower', aspect='auto')
				plt.xlabel('time bin')
				plt.ylabel('turn')
				plt.savefig('imshow_prof1')
				plt.show()

			
			prof_data_flat = prof_data_a.flatten()
			fft_prof_data = np.fft.fft(prof_data_flat)
			fftabs_prof_data = np.abs(fft_prof_data)
			
			from scipy import signal
			widths = np.arange(1,1000)
			cwtmatr = signal.cwt(prof_data_flat, signal.ricker, widths)
			#plt.plot(prof_data_flat)
			plt.imshow(cwtmatr, extent=[-1, 1, widths[-1], 1], origin='lower', aspect='auto')
			plt.show()
			
			sys.exit()
			#plt.matshow(prof_data_a)
			#plt.imshow(prof_data_a, origin='lower')
			#plt.show()
			
			#plt.plot(prof_data_a[0], 'k')
			#plt.plot(prof_data_a[100], 'g')
			#plt.plot(prof_data_a[-1], 'r')
			#plt.show()
			#plt.plot(integral_by_turn)
			#plt.ylim(ymin=0)
			#plt.show()

			xd = np.arange(1,n_slices+1)
			
			t_turn_data_ns = 1e9*t_turn_data + tadj
			
			if plot_data_shift:

				if shift_data_winmin:
					plt.subplot(211)

				plt.plot(xd, data_integ_0,'ko-',label='before shift')
				plt.plot(xd, data_integ_1,'ro-',label='after shift')
				#plt.plot(xd, data_integ_1_norm, 'bo-')
				plt.axvline(x = synch_index, color='gray')
				plt.xlim(xmax=len(data_integ_0))
				
				plt.ylabel('integrated signal')
				plt.title('shift = '+str(idelay))

				if shift_data_winmin:
					plt.subplot(212)
					plt.plot(runmean, 'ko-')
				plt.xlabel('slice number')
				
				plt.xlim(xmax=len(data_integ_0))

				

				plt.legend()
				plt.savefig('sig_integ')
				plt.show()
				#sys.exit()
				
			
		else:
			t_turn_data_ns = 1e9*t_turn_data

	else:
		t_turn_data_ns = 1e9*t_turn_data + tadj

	
	#normalise each turn of data by sum at that turn
	#prof_row_sum =  np.sum(prof_data_a, axis=1)
	#prof_data_norm_a =  np.transpose(np.transpose(prof_data_a)/prof_row_sum)

	print "write ", nturns, " turns to file"
	#f1 = open('kurridata.txt', 'w')
	#for it in range(nturns):
	#	yprof = prof_data_a[it]
	#	yprof_norm = yprof/max(yprof)
	#	#plt.plot(yprof,'ko-')
	#	for y in yprof:
	#		print >>f1,y
	#f1.close()


	chisq_sum_l = []

	if ncheck == 0:
		idelay_a = np.array([idelay])
	else:
		#idelay_a = np.arange(idelay - ncheck, idelay+ncheck+1)

		#idelay_a = np.arange(idelay - ncheck, idelay+3)
		#idelay_a = np.arange(-35, -18)

		if phi_s_deg == 0:
			idelay_a = np.arange(-30, -18)
		else:
			idelay_a = np.arange(-35, -18)

	print "idelaya ", idelay_a


	for idelay1 in idelay_a:

		istart = istart + idelay1
		prof_data_a, t_turn_data = sort_data(istart, ttime_a, time_indices, nturns_sort, n_slices, interpolate_data, nshift= idelay1, show_mrange=False)
		

		prof_row_sum =  np.sum(prof_data_a, axis=1)
		prof_data_norm_a =  np.transpose(np.transpose(prof_data_a)/prof_row_sum)

		print "write ", nturns, " turns to file"
		print "first, last time in ttime_a [ms] ",1e3*ttime_a[0] - offset_file_ms, 1e3*ttime_a[nturns] - offset_file_ms
		#sys.exit()

		f1 = open('kurridata.txt', 'w')
			
		profmaxind = [np.argmax(prof_data_a[it]) for it in range(nturns)]
		intensity_totomo = []
		for it in range(nturns):
			int1 = time_interval_ns*np.trapz(prof_data_a[it])
			intensity_totomo.append(int1)
			#plt.plot(prof_data_a[it]+0.01*it)
		#plt.show()
	
		
		#plt.plot(profmaxind)
		#plt.ylabel("index of maximum signal")
		#plt.show()

		
		#plt.plot(intensity_totomo)
		#plt.ylabel("intensity")
		#plt.show()
		
		
		for it in range(nturns):
			yprof = prof_data_a[it]
			yprof_norm = yprof/max(yprof)
			#plt.plot(yprof,'ko-')
			for y in yprof:
				print >>f1,y
		f1.close()



		run_tomo_code = True
		if run_tomo_code:
		
			recon_step = int(recon_step_frac*Ts_a[index_time])
			input_settings = tdefs.get_input("input_default.dat")
		
			binwidth_time = time_interval

			fname = 'input_v2.dat'

			descrip = 'KURRI'
			datafile = 'kurridata.txt'
			input_settings["numframes"] = nturns
			input_settings["numbins"] = n_slices
			input_settings["synchbin"] = synch_index #int(0.5*n_slices)
			input_settings["framebinwidth"] = binwidth_time/npfac	
			input_settings["binslowerignore"] = binslowerignore
			input_settings["binsupperignore"] = binsupperignore
			input_settings["filmstart"] = recon_start
			input_settings["filmstop"] = recon_stop
			input_settings["filmstep"] = recon_step
			input_settings['Bdot'] = Bdot
			input_settings['VRF1ref'] = V_rf
			input_settings['Rnom'] = Rnom
			input_settings['rhonom'] = bending_radius
			input_settings['gammatnom'] = gamma_tr
			input_settings['pickup'] = 1.0
			input_settings['selffields'] = 0
			tdefs.write_inputfile(fname,descrip, datafile, input_settings)

			temp_files = ['plotinfo.data','profiles.data','d001.data','d100.data','jmax.data','image001.data','image100.data']
			#for file1 in temp_files:
			#	if os.path.isfile(file1):
			#		os.remove(file1)
			#		print "delete ",file1


			time_cerncode0 = time.time()
			call(["tomo"])
			time_cerncode1 = time.time()

			print "CERN code execution time ",time_cerncode1-time_cerncode0

			files = os.listdir(os.curdir)
			for fn in files:
				if fn[-4:] == 'data':
					os.rename(fn,tomo_outdir+"/"+fn)

					 
	
		show_tomo_result = True
		if show_tomo_result:

			#read output parameters from plotinfo.data
			plotinfo = tdefs.get_plotinfo(input_file=tomo_outdir+"/plotinfo.data")
		
			x0 = plotinfo['xat0'][0]
			y0 = plotinfo['yat0'][0]
			dEbin =  plotinfo['dEbin'][0]
			dtbin =  plotinfo['dtbin'][0]
			phi_s_out = plotinfo['phi_s'][0] #synch phase

			fname_prof = 'profiles.data'
			profiles_norm_raw = read_data_file(tomo_outdir+'/'+fname_prof)
			profiles_norm_a = np.split(profiles_norm_raw, nturns)


			show_profile_out = False
			if show_profile_out:

				fname_raw = 'kurridata.txt'

				bunchdata= read_data_file(fname_raw)
				raw_spl = np.split(bunchdata, nturns)

				#plt.plot(profiles_a)
				for raw, prof in zip(raw_spl, profiles_norm_a):
					plt.plot(raw,'k')
					plt.plot(prof,'ro')
				plt.show()

				sys.exit()

		
			inputp = tdefs.get_input()
			E0, gamma0, eta0, beta0, omega_rev0, Vpeak = tdefs.set_param(inputp)
			T = 2*np.pi/omega_rev0 #rf period

			max_dp_sep = max(dp_sep)

	
		#len_img_data = len(image_data_a)
		ndim = n_slices_active  #int(len_img_data**0.5)

		#construct time, dp arrays taking synchronous point x0, y0 into account.
		dE_a_eV = dEbin*(np.array(range(ndim)) - y0)
		dE_a_mev = 1e-6*dE_a_eV
		E_tot_eV_a = Etot_a[index_time] + dE_a_eV
		p_a = [1e-6*(Et**2 - PROTON_MASS**2)**0.5 for Et in E_tot_eV_a]
		dp_a = (p_a - p0_mev)/p0_mev
		t_synch = 1e9*phi_s_out*T/(2*math.pi)
		t_tomo_ns = 1e9*dtbin*(np.array(range(ndim)) - x0) + t_synch

		phase_tomo = t_tomo_ns*(2*math.pi)/T_rf_ns
		phase_tomo_deg = phase_tomo*180/math.pi


		phase_mesh, dp_mesh = np.meshgrid(phase_tomo, dp_a)

		H_mesh = -H_fn(phase_mesh, dp_mesh, phi_s, a1, a2)

		show_H_mesh = False
		if show_H_mesh:
			plt.contour(phase_tomo_deg, dp_a, H_mesh)
			plt.plot(phase_sep_deg, dp_sep,'r') #bucket
			plt.plot(phase_sep_deg, -dp_sep,'r') #bucket
			plt.show()
			sys.exit()
		
		t_shift_ns = -1e9*dtbin*(x0 - 0.5*ndim) #+ t_synch
		#print "time shift applied ",t_shift_ns
		
		t_data_step = t_turn_data[1] - t_turn_data[0]
		#t_tomo_ns = np.linspace(t_sep_ns[0], t_sep_ns[-1], n_slices)
		
		recon_id = np.arange(recon_start, recon_stop+1, recon_step)
		print "recon_id ",recon_id

		#recon_id = recon_id[:-1]
		#plt.ion()

		image_all_a = []
		#image_std = []
		#tomo_time_proj_l = []

		image_scaled = True #scale image data by measured profile integral

		image_std_l = []
		image_max_l = []

		iloop = 0 
		chisq_l = []
		intensity_scale = []

		#create regular grid in Hamiltonian
		ngrid = 20
		H_grid = np.linspace(0, -Hsep, ngrid)
		H_grid_delta = H_grid[1] - H_grid[0]

		for irc in recon_id:

			#if run_pyheadtail:

			#counts, xbins, ybins,_ = plt.hist2d(bunch_t_ini_ns, bunch_dp_ini, bins=60, norm=LogNorm())

			print "recon at turn ",irc
			#process main tomography output file
			#fname_img = 'image001.data'

			fname_img = tomo_outdir+'/image000'[:-len(str(irc))]+ str(irc) + '.data'
			print "fname_img ",fname_img
		
			time_fread0 = time.time()
			image_data_a = read_data_file(fname_img) #1D data

			time_fread1 = time.time()
			#print "time to read file ",time_fread1 - time_fread0

			#image_trimzeroes = np.trim_zeros(image_data_a)
			#mean_nonzero = np.mean(image_trimzeroes)
			#print image_data_a.shape, image_trimzeroes.shape
			#print "mean, std w/o zeroes ",mean_nonzero, np.std(image_trimzeroes), np.std(image_data_a)
			#hist, bin_edges = np.histogram(1e3*image_data_a)
			
	
			#pts_above = np.where(image_data_a > 4*np.std(image_data_a))[0]
			#print "numpy of points above limit ",pts_above.shape
			
			#split into array
			image_data_spl = np.split(image_data_a, ndim)
			#transpose to get phase on horizontal axis, energy on vertical
			image_data_tr = np.transpose(image_data_spl)

		

			#flatten, sort, interpolate and rebin image and Hamiltonian
			image_flat = image_data_tr.flatten()
			H_flat = H_mesh.flatten()
			argHsort = H_flat.argsort()
			image_flat_sort = image_flat[argHsort]
			H_flat_sort = H_flat[argHsort]
			#map onto regular grid			
			image_Hgrid = np.zeros(ngrid)
			for H1, image1 in zip(H_flat_sort, image_flat_sort):
				ig =  int(H1/H_grid_delta)
				if ig < ngrid:
					image_Hgrid[ig] = image_Hgrid[ig] + image1
			#normalise by maximum value
			image_flat_maxn = image_flat/np.max(image_flat)
			image_Hgrid_maxn = image_Hgrid/np.max(image_Hgrid)
			#sum points
			image_Hgrid_integ = np.array([sum(image_Hgrid[:i]) for i in range(ngrid)])
			image_Hgrid_integ_rev = 1 - image_Hgrid_integ

			plot_image_H = False
			if plot_image_H:
				plt.plot(H_flat, image_flat_maxn,'ko')
				plt.plot(H_grid, image_Hgrid_maxn,'ro-')
				plt.plot(H_grid, image_Hgrid_integ_rev,'mo-')
				plt.xlim(0, -Hsep)
				plt.show()
					
			#print "image_data_ones shape ",image_data_ones.shape
			#image_data_f = filter_isolated_cells(image_data_ones, struct=np.ones((3,3)))

			filter_isolates = True
			if filter_isolates:
				#image_data_fil  = filter_isolated_cells(image_data_tr, image_data_ones, struct=np.array([[0,1,0],[1,1,1],[0,1,0]]))
				image_data_fil  = filter_isolated_cells(image_data_tr, struct=np.array([[0,0,0],[1,1,1],[0,0,0]]))
				image_data_fil  = filter_isolated_cells(image_data_fil, struct=np.array([[0,1,0],[0,1,0],[0,1,0]]))

			#print "indices_f shape ",indices_f
			#image_data_fil = image_data_tr[indices_f]

			#print "filtered image ",image_data_fil.shape
			#plt.subplot(121)
			#plt.spy(image_data_tr)
			#plt.subplot(122)
			#plt.spy(image_data_fil)
			#plt.savefig('filter_isolates')
			#plt.show()
			#sys.exit()

			if filter_isolates:
				image_data_tr = image_data_fil
				image_data_spl = np.transpose(image_data_tr)

			image_sum_all = np.sum(image_data_a)

			if image_scaled:
				intensity_scale.append(intensity_totomo[irc-1])
				image_std_l.append(np.std(image_data_tr)*intensity_totomo[irc-1])
				image_max_l.append(np.amax(image_data_tr)*intensity_totomo[irc-1])
			else:
				image_std_l.append(np.std(image_data_tr))
				image_max_l.append(np.amax(image_data_tr))

			#print np.std(image_data_tr), np.std(image_data_tr.flatten)
			#sys.exit()

			
			#tomo_profile = [sum(d) for d in image_data_tr]
			tomo_time_proj = [sum(d) for d in image_data_spl]
			tomo_mom_proj = [sum(d) for d in image_data_tr]

			image_all_a.append(image_data_tr)		
			#tomo_time_proj_l.append(tomo_time_proj)

			

			#compare with data
			prof_diff = profiles_norm_a[irc-1] - tomo_time_proj
			#prof_diff = prof_data_norm_a[irc-1][binslowerignore: n_slices - binsupperignore] - tomo_time_proj
			chisq = sum(prof_diff**2)

			#print "chisq ",chisq
			chisq_l.append(chisq)
		
			if iloop == 0:
				tomo_time_proj_a = np.array([tomo_time_proj])
				image_Hgrid_a = image_Hgrid
			else:
				tomo_time_proj_a = np.vstack((tomo_time_proj_a, tomo_time_proj))
				image_Hgrid_a = np.vstack((image_Hgrid_a, image_Hgrid))


			iloop = iloop + 1

				
		chisq_sum = sum(chisq_l)
		chisq_sum_l.append(chisq_sum)

		print "idelay, chisq summed over turns ",idelay1, chisq_sum
	

	if len(idelay_a) > 10:

		plt.plot(idelay_a, chisq_sum_l,'ko-')
		plt.xlabel('time point shift')
		plt.ylabel('chisq sum')
		plt.ylim(ymin=0)
		plt.axvline(x=idelay, color='gray')
		plt.savefig('timeshiftscan')
		plt.show()


compare_image_Hgrid = True
if compare_image_Hgrid:
	for image_H in image_Hgrid_a:
		plt.plot(H_grid, image_H)
	plt.show()
	sys.exit()

compare_sel_profiles = False
if compare_sel_profiles:
	col_l = ['k','r','b']
	ip = 0
	
	for proj, irc in zip(tomo_time_proj_a, recon_id):
		plt.subplot(211)
		plt.plot(proj,col_l[ip])
		plt.subplot(212)
		plt.plot(prof_data_a[irc-1], col_l[ip]) 
		ip = ip + 1
	plt.show()


#cont_limit = 4.2*max(image_std_l)
cont_limit = max(image_max_l)

if animate:



	#plt.subplot(211)
	levels = np.linspace(0, cont_limit, 30)
	#plt.title('reconstructed turn '+sttomo_time_proj_lr(irc)+ ' time '+time_str+' ms')
	fig, ax = plt.subplots()
	cont = ax.contour(t_tomo_ns, dp_a, image_all_a[0], 40)
	def animate(i):	
		ax.clear()	
		if image_scaled:
			cont = ax.contour(t_tomo_ns, dp_a, intensity_scale[i]*image_all_a[i], levels = levels)
		else:
			cont = ax.contour(t_tomo_ns, dp_a, image_all_a[i], levels = levels)
		ax.plot(t_sep_ns, dp_sep,'r') #bucket
		ax.plot(t_sep_ns, -dp_sep,'r') 
		ax.set_xlabel('time [ns]')
		ax.set_ylabel(r"$\Delta$p/p")
		ax.set_xlim(1.01*t_sep_ns[0],1.01*t_sep_ns[-1])
		ax.set_ylim(1.01*bh, -1.01*bh)
		#cont.set_zdata(image_all_a[1])
		return cont


	Nt = len(image_all_a)
	anim = animation.FuncAnimation(fig, animate, frames=Nt, interval=500, blit=False)
	anim.save('tomo.gif', dpi=80, writer='imagemagick')
	#anim.save('tomo.mp4', fps=30, extra_args=['-vcodec', 'libx264'])
	#anim.save("tomo.mp4", writer=writer)
	plt.show()

else:
	
	plot_surface = False
	plot_contour = True
	
	if plot_surface:
		from mpl_toolkits.mplot3d import Axes3D
		t_mesh, dp_mesh = np.meshgrid(t_tomo_ns, dp_a)

		print "len intensity_totomo ",len(intensity_totomo)
		print "image max over all profiles ",max(image_max_l)

		ipl = 0
		for irc in recon_id:
		
			#print "image sum ",np.sum(image_all_a[ipl])

			
			#print "measured profile integral ",intensity_totomo[irc-1]
			
			if image_scaled:
				image = image_all_a[ipl]*intensity_totomo[irc-1]
			else:
				image = image_all_a[ipl]

			fig = plt.figure()
			ax = fig.gca(projection='3d')			
			plt.title('reconstructed turn '+str(irc)+ ' time '+time_str+' ms')
			ax.plot_surface(t_mesh, dp_mesh, image, cmap="seismic") #or coolwarm

			ax.set_xlabel('time [ns]')
			ax.set_ylabel(r"$\Delta$p/p")
			
			
			ax.set_zlim(0, cont_limit)
			ipl = ipl + 1	

			plt.show()

		


	print "recon_id ",recon_id
	ipl = 0
	for irc in recon_id:

		versus_time = False
		if versus_time:
			xaxis = t_tomo_ns
		else:
			xaxis = phase_tomo_deg

		if image_scaled:
			image = image_all_a[ipl]*intensity_totomo[irc-1]
		else:
			image = image_all_a[ipl]
		
		show_imshow_spy = False
		if show_imshow_spy:
			plt.subplot(211)
			plt.imshow(image, vmax = cont_limit)
			plt.subplot(212)
			plt.spy(image)
			plt.show()
			sys.exit()
		

		#levels = np.linspace(0, 0.008, 20)
		levels = np.linspace(0, cont_limit, 30)
		print "plot phasespace at turn ",irc
		
		plt.subplot(211)

		#plt.subplot(211)
		#plt.title('reconstructed turn '+str(irc)+ ' time '+time_str+' ms')
		plt.title('reconstructed turn '+str(irc))
		if plot_contour:
			ctom = plt.contour(xaxis, dp_a, image, levels = levels)
			#ctom = plt.contour(t_tomo_ns, dp_a, image)
			#plt.colorbar(ctom)
		
		#print "ctom ",ctom.levels	
		#sys.exit()

		#plt.contourf(t_gen_ns, dp_bin_ini, Hgrid_inv, ctom.levels)
		#plt.contourf(t_gen_ns, dp_bin_ini, Hgrid_inv)

		if versus_time:
			plt.plot(t_sep_ns, dp_sep,'r') #bucket
			plt.plot(t_sep_ns, -dp_sep,'r') #bucket
			plt.axvline(x=t_tomo_ns[int(x0)]) #synch time
			plt.xlim(1.01*t_sep_ns[0],1.01*t_sep_ns[-1])
			plt.xlabel('time [ns]')
		else:
			plt.plot(phase_sep_deg, dp_sep,'r') #bucket
			plt.plot(phase_sep_deg, -dp_sep,'r') #bucket
			plt.axvline(x=phi_s_deg) #synch time
			plt.xlabel(r'$\phi$ [deg]')

		#plt.plot(xaxis, np.zeros(len(xaxis)),'mo') #check xaxis points
		
		plt.ylim(1.01*bh, -1.01*bh)

		x1,x2,y1,y2 = plt.axis()
		
		plt.ylabel(r"$\Delta$p/p")

		plt.subplot(212)

		#t_turn_data_ns_shift = t_turn_data_ns #- 20# t_shift_ns

		#plt.plot(t_turn_data_ns, prof_data_norm,'ko-',linewidth=2,label='data')
		#plt.plot(t_turn_data_ns, prof_data_norm_a[irc-1],'k.-',linewidth=2,label='data')

		plt.plot(xaxis, profiles_norm_a[irc-1],'ko-',linewidth=2,label='data')

		#plt.plot(t_tomo_ns, prof_data_norm,'ko-',linewidth=2,label='data')

		#print "len t_tomo_ns ",len(t_tomo_ns)
		#print "tomo_time_proj_a shape ",tomo_time_proj_a.shape
		
		plt.plot(xaxis, tomo_time_proj_a[ipl],'r-',linewidth=2,label='tomography')
		
		#plt.axvline(x=t_bucket[0], color='gray')
		#plt.axvline(x=t_bucket[-1], color='gray')

		if versus_time:
			if binslowerignore > 0:
				plt.axvline(x=t_turn_data_ns[binslowerignore] + 0.1*t_data_step, color='gray', linestyle='--')
			if binsupperignore > 0:
				plt.axvline(x=t_turn_data_ns[n_slices - binsupperignore]-0.1*t_data_step, color='gray', linestyle='--')		
			plt.axvline(x=t_synch, color='gray', linestyle = '--')
			plt.xlim(1.01*t_sep_ns[0],1.01*t_sep_ns[-1])
			plt.xlabel('time [ns]')		

		else:
			plt.axvline(x=phi_s_deg, color='gray', linestyle = '--')
			plt.xlabel(r'$\phi$ [deg]')		
		#plt.ylim(-1.1*max_dp_sep, 1.1*max_dp_sep) 


		#plt.xlim(t_tomo_ns[0],t_tomo_ns[-1])
		#plt.xlim(t_sep_ns[0], t_sep_ns[-1])
		#plt.ylim(-1.1*max_dp_sep, 1.1*max_dp_sep) 
		
		plt.ylabel('normalised intensity')
		plt.legend(loc = 'upper right')
		plt.savefig('phasespace_'+str(irc))
		plt.tight_layout()


		plt.show()

		ipl = ipl + 1








