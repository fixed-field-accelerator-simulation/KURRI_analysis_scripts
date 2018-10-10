from __future__ import division
import pylab as plt
import numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import quadrature
import sys
import os

import matplotlib.cm as cm
import analysedata_180919 as adata
import tomo_defs_181008 as tdefs
import math
from subprocess import call
import time 
import matplotlib.animation as animation
from subprocess import call

#set up formatting for movie files
#Writer = animation.writers['ffmpeg']
#writer = Wriiter(fps=15, metadata=dict(artist='Me'), bitrate=1800)

read_rf_waveform = True

tomo_times = np.array([12e-4]) #start time at which do do tomography #8.0e-4
dirname = "../data/2018/9/20/"
#dirname = "../data/2018/9/21/machida_san_foil_escape/" # "../data/2018/9/20/"
filter_string = '' #look for csv files with this string in file name
phi_s_deg = 20 #set synchronous phase
turns_syncfrac = 0.5 #fraction of synchrotron oscillation to include
recon_start = 1 #reconstruct phase space at this turn first
recon_step = 50 #step size terms of turn number between recon_start and recon_stop
recon_stop = 101 #reconstruct phase space at this turn last
animate = False

ncheck = 0 #number of synch point settings to check on either side of that set

idelay =  -27 #-21 #apply this delay to rf zerocrossing time
#if idelay=None, use of the following methods to establish it

#...................................
#these methods are applied if idelay = None
shift_data_sym = False #shift data by looking for symmetry point
shift_data_winmin = True #shift minimum of moving window average where window is out-of-bucket duration
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
mountain_range_plot = True
plot_data_shift = False

print "phi_s set to ",phi_s_deg, " degrees"

npfac = 1.0 #interpolate onto time axis with npfac times the number of points as the data
interpolate_data = True #if False interpolate data onto time axis which is determined by rf waveform

tomo_outdir = "tomo_output"
files = os.listdir(os.curdir)
if tomo_outdir not in files:
	print "create output directory "+tomo_outdir
	os.mkdir(tomo_outdir)
0.000334929534919

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
	print "selected rf file ",i_rf_file_sel


if phi_s_deg == 0:
	binslowerignore = 0 #add zeroes before profile (earlier time)
	binsupperignore = 0 #add zeroes after profile (later time)
else:
	pass
	#binslowerignore = int(npfac*35) #add zeroes before profile (earlier time)
	#binsupperignore = int(npfac*19) #add zeroes after profile (later time)

#longitudinal tomography based on bunch monitor data.
#dirname = "../data/2015/6/24/700A/"

PROTON_MASS = 938.27203e6
SPEED_OF_LIGHT = 299792458

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
	
def f_rf_lookup(time):
	""" Find kinetic energy at specified time in svk file """
	i = 0 
	rf_freq = None
	for t in time_svk:
		if t > time:
			rf_freq = rf_freq_svk[i]	
			break	
		i = i + 1
	
	return rf_freq
	
		
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
	tdat_chf, data_chf = adata.read_scope(dirname, ch_files_sel[0][0])

else:
	fname = csv_files[i_file_sel]
	#fname = "20180920-004.csv"#"20180920-004.csv"
	#read scope data
	tdat_chf, data_chf = adata.read_scope(dirname, fname)

	if read_rf_waveform and i_rf_file_sel != None:
		fname = csv_files[i_rf_file_sel]
		tdat_rf, sig_rf = adata.read_scope(dirname, fname)
		tdat_rf_ms = 1e3*np.array(tdat_rf)
		sig_rf_sc = 0.01*np.array(sig_rf)

		
time_interval = tdat_chf[2] - tdat_chf[1]
tdat_chf_ms = 1e3*np.array(tdat_chf)


i1 = 0
for i1, t1 in enumerate(tdat_chf):
	if t1 == 0:
		index_zero = i1
		break
		
tof_guess = 600e-9


def show_data(time_indices, nturns, shift=0):

	print "tdat_rf_ms step ",tdat_rf_ms[1] - tdat_rf_ms[0]
	tpad = 2e-4
	tfirst = tdat_rf_ms[time_indices[0]] -tpad
	t20 =  tdat_rf_ms[time_indices[20]] + tpad
	tm20 =  tdat_rf_ms[time_indices[nturns-20]] - tpad
	tlast = tdat_rf_ms[time_indices[nturns]] + tpad

	plt.subplot(311)
	#plt.plot(tdat_rf_ms, sig_rf_sc, 'r')
	plt.plot(tdat_chf_ms, data_chf,'k')
	
	#plt.axvline(x=t_hminus, color='r')
	
	for lower_lim in tomo_times_data_l:
		plt.axvline(x=tfirst, color='r')
		plt.axvline(x=tlast, color='r')
	
	plt.xlim(0, tdat_chf_ms[-1])
	plt.ylim(-0.2, 0.2)
	plt.ylabel('bunch monitor (all)')
	plt.axvline(x=tdat_rf_ms[time_indices[0]])
	
	plt.subplot(312)
	plt.plot(tdat_chf_ms + 1e3*time_interval*shift, data_chf,'k')
	if read_rf_waveform:
		plt.plot(tdat_rf_ms , sig_rf_sc, 'r')
	plt.xlim(tfirst, t20)
	plt.ylim(-0.1, 0.1)

	plt.axvline(x=tdat_rf_ms[time_indices[0]])
	plt.subplot(313)
	plt.plot(tdat_chf_ms + 1e3*time_interval*shift, data_chf,'k')
	if read_rf_waveform:
		plt.plot(tdat_rf_ms, sig_rf_sc, 'r')
	plt.xlim(tm20, tlast)
	plt.ylim(-0.1, 0.1)

	plt.axvline(x=tdat_rf_ms[time_indices[nturns-20]])
	plt.ylabel('bunch monitor (first 20 turns)')
	plt.xlabel('time [ms]')
	plt.savefig('bunchmonitor_data')
	plt.show()
	
	




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

	it = 0.1
	data_win_sig_sel = []
	for data_win in sig_dat_split[:nturns]:
		#print "data_win ",data_win
		data_win_sig =  max(data_win) - np.array(data_win)
		data_win_sig_sel.append(data_win_sig)
		
		if mountain_range_plot and show_mrange:
			lw = 0.5 - it / 3000.0
			ax.plot(range(1,n_slices+1), 150*data_win_sig + it,'k', lw=lw)
		it = it + 1
	
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
	
def find_symmetry_point(data):
	"""Look for point where data on either side is most similar"""

	ld = len(data)
	ldh = int(0.5*ld)
			
	chisq_l = []
	for d in range(ldh):
		#indh1 = range(d, ldh+d)
		#indh2 = range(ldh+d, ld+d)
		indh1 = range(ldh+d, ld+d)
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


def read_data_file(filen):
	
	f1 = open(filen,'r')
	data_l = []
	for data in f1:
		data_l.append(float(data))
	data_a = np.array(data_l)

	return data_a

	
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

#calculate kinetic energy, beta, synchrotron period at each time

if not flattop:
	Ekin_a = np.array([ke_lookup(t) for t in tomo_times])
else:
	ke_flattop = ke_lookup(flattop_time)

	Ekin_l = []
	for t in tomo_times:
		if t < flattop_time:
			Ekin_l.append(ke_lookup(t))
		else:
			Ekin_l.append(ke_flattop)
	
	Ekin_a = np.array(Ekin_l)

print "Ekin_a ",Ekin_a 
Etot_a = Ekin_a + PROTON_MASS
beta_a = adata.ke_to_relativistic_beta(Ekin_a, PROTON_MASS)
betagam_a = adata.ke_to_relativistic_beta_gamma(Ekin_a, PROTON_MASS)
gamma_a = betagam_a/beta_a
eta_a = alpha_c - gamma_a**(-2)
Qs_a = np.sqrt(np.abs(eta_a) * V_rf / (2 * np.pi * beta_a**2 * Etot_a))
Ts_a = 1/Qs_a

print "kinetic energies at selected times ",Ekin_a
print "beta at selected tims ",beta_a
print "eta_a ",eta_a 
print "Qs ",Qs_a
print "Ts ",Ts_a


#plt.plot(tomo_times, Ts_a,'ko-')
#plt.show()

tomo_times_data_l = []
istart_guess_l = []
for tt in tomo_times:
	
	icrude = index_zero + int(round(tt/time_interval))
	t_crude = tdat_chf[icrude]

	istart_guess_l.append(icrude)
	tomo_times_data_l.append(t_crude)



points_per_turn = int(round(tof_guess/time_interval)) 


print "specified initial times ",tomo_times

#loop over times, do tomography on each loop
for index_time in range(len(tomo_times)):

	#number of turns specified in terms of sync period
	nturns = int(turns_syncfrac*Ts_a[index_time])
	
	p0_mev = 1e-6*np.sqrt(Etot_a[index_time]**2 - PROTON_MASS**2)
	bh = ((2*Qs_a)/(harmonic*abs(eta_a[index_time]))) #bucket height
	
	time_str = str(1e3*tomo_times[index_time])
	lower_lim = tomo_times_data_l[index_time]
	
	nturns_fit = 1001
	if not read_rf_waveform:
		#calculate TOF
		tof, time_indices, ttime_a = tdefs.calc_tof_bm_data(lower_lim, points_per_turn, nturns_fit)
		istart = time_indices[0]
	else:
		#tof_m1, indices_m1, ttime_a_m1 = tdefs.calc_tof_bm_data(lower_lim, points_per_turn, nturns_fit)
		tof, time_indices, ttime_a = tdefs.calc_tof_rf_waveform(tdat_rf, sig_rf, nturns_fit, istart_guess_l[index_time])
		istart = time_indices[0]

		print "tof rf ",tof
		
		#tof_m1_a = [ttime_a_m1[i] - ttime_a_m1[i-1] for i in range(1, nturns_fit)]
		tof_a = [ttime_a[i] - ttime_a[i-1] for i in range(1, nturns_fit)]
		
		#plt.plot(tof_m1_a, 'ko-')
		#plt.plot(tof_a, 'ro-')
		#plt.show()
		#sys.exit()


	#tof_fit_a = t_zero_cross

	T = tof/harmonic#2*np.pi/(harmonic*omegarev) #rf period
	T_rf_ns = 1e9*T
	omega_rev = 2*math.pi/T
	omega_rf = omega_rev*harmonic

	#update points per turnt_turn_data
	points_per_turn = int(round(tof/time_interval)) 
	n_slices = int(npfac*points_per_turn)
	#n_slices_active = n_slices-(binslowerignore + binsupperignore)
	
	#bucket phase extrema calculation
	if phi_s_deg == 0:
		phi_l1 = math.pi
		phi_l2 = -math.pi
	else:
		phi_l1, phi_l2 = phase_extrema(phi_s)

	print "phi_l1, phi_l2 ",phi_l1, phi_l2
	if phi_l1 < phi_l2:
		phi_lower = phi_l1
		phi_upper = phi_l2
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
	print "bucket points, out-of-bucket points ",n_slices_active, n_slices_outb

	if shift_data_winmin and phi_s_deg == 0:
		print "shift_data_winmin must be False if phi_s_deg is zero"

	if phi_s_deg != 0:
		binslowerignore = n_slices_outb#abs(int(n_slices*phi_lower/(math.pi))) #n_slices_outb
		binsupperignore =   n_slices - (binslowerignore + n_slices_active)

	#omega_rf = 2*math.pi*f_rf_lookup(tomo_times[index_time])
	npts = 1000
	a1 = 0.5*harmonic*omega_rf*eta_a[index_time]
	a2 = omega_rf*V_rf/(2*math.pi*(beta_a[index_time]**2)*Etot_a[index_time])
	dp_sep, t_sep_ns = tdefs.separatrix(phi_s, eta_a[index_time], T, harmonic, npts, a1, a2)
	
	phase_sep = 2*math.pi*(t_sep_ns)/(1e9*T)
	t_bucket = np.linspace(t_sep_ns[0], t_sep_ns[-1], n_slices)

	#t_data = np.linspace(-0.5*T_rf_ns, 0.5*T_rf_ns, n_slices) #data fills rf period

	t_phi_sep_ns = phi_s*(1e9*T)/(2*math.pi) #time at synchronous phase
	t_phi_l2 = T_rf_ns*(phi_l2/(2*math.pi)) #time at phi_l2
		
	t_sep_ns = np.insert(t_sep_ns, 0, t_phi_l2)
	dp_sep = np.insert(dp_sep,0,0)

	phase_sep_deg = (180/math.pi)*phase_sep
	sync_frac =  (t_phi_sep_ns- t_sep_ns[0])/bucket_dur_ns
	print "synch_frac ",sync_frac

	sync_bin = binslowerignore + int(n_slices_active*sync_frac)
	print "sync bin ",sync_bin
	

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
	
	nturns_sort = 1000
	prof_data_a, t_turn_data = sort_data(istart, ttime_a, time_indices, nturns_sort, n_slices, interpolate_data, show_mrange = False)

	
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
			isym_list = []

			turn_list = np.arange(nturns_sync, nturns_sort)

			for i1 in turn_list:
				data_int = np.sum(prof_data_a[:i1], axis=0)
				isym, chisq_min = find_symmetry_point(data_int*(nturns_sync/i1))
				
				
				binl = range(ld)
				
				chisq_sym_turns.append(chisq_min)
				isym_list.append(isym)

				indsh = np.arange(ld) + isym
				data_integ_shift = data_int.take(indsh, mode='wrap')

				#if i1%10 == 0:
					#plt.plot(data_integ_shift*(220/i1))
			#plt.show()

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

			isym_final =  isym_list[ithres]
			idelay =  isym_final - sync_bin #- int(0.5*n_slices))
			
			print "midpoint bin",int(0.5*n_slices)
			print "synch phase w.r.t midpoint bin ",(sync_bin - int(0.5*n_slices))

			print "symmetry point ",isym_list[ithres]
			print "calculated delay is ",idelay
					
			show_sym_plots = True
			if show_sym_plots:
				plt.subplot(311)
				plt.plot(turn_list, chisq_sym_turns)
				plt.ylabel('asymmetry')
				plt.ylim(ymin=0)
				plt.subplot(312)
				plt.plot(turn_chisq_check, range_chisq_win_l,'ro')
				plt.axvline(x=turn_chisq_check[ithres])	
				plt.ylabel('window range')
				plt.ylim(ymin=0)
				plt.subplot(313)
				plt.plot(turn_list, isym_list, 'k')
				plt.ylabel('sym point')
				#plt.plot(data_int)
				#plt.plot(data_integ_shift)
				plt.savefig('symcalc')
				plt.show()

				indsh = np.arange(ld) + idelay
				data_integ_shift = data_integ_0.take(indsh, mode='wrap')
				plt.plot(data_integ_shift)
				plt.axvline(x = 0.5*ld)
				
				plt.show()

				sys.exit()



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

		idelay_a = np.arange(idelay - ncheck, idelay+3)
		#idelay_a = np.arange(-26, -15)
		
	print "idelaya ", idelay_a


	for idelay1 in idelay_a:

		istart = istart + idelay1
		prof_data_a, t_turn_data = sort_data(istart, ttime_a, time_indices, nturns_sort, n_slices, interpolate_data, nshift= idelay1, show_mrange=False)

		prof_row_sum =  np.sum(prof_data_a, axis=1)
		prof_data_norm_a =  np.transpose(np.transpose(prof_data_a)/prof_row_sum)

		print "write ", nturns, " turns to file"
		f1 = open('kurridata.txt', 'w')
		for it in range(nturns):
			yprof = prof_data_a[it]
			yprof_norm = yprof/max(yprof)
			#plt.plot(yprof,'ko-')
			for y in yprof:
				print >>f1,y
		f1.close()



		run_tomo_code = True
		if run_tomo_code:
		
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

			print "delete files before run "
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


			show_profile_out = False
			if show_profile_out:

				fname_prof = 'profiles.data'
				fname_raw = 'kurridata.txt'

				bunchdata= read_data_file(fname_raw)
			
				profiles_a = read_data_file(tomo_outdir+'/'+fname_prof)

				nprof = len(bunchdata)/n_slices
				print len(profiles_a), len(profiles_a)/64

				print "nprof ",nprof

				raw_spl = np.split(bunchdata, nprof)
				prof_spl = np.split(profiles_a, nprof)

				#plt.plot(profiles_a)
				for raw, prof in zip(raw_spl, prof_spl):
					plt.plot(raw,'k')
					plt.plot(prof,'ro')
				plt.show()

		
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
		
		t_shift_ns = -1e9*dtbin*(x0 - 0.5*ndim) #+ t_synch
		#print "time shift applied ",t_shift_ns
		
		t_data_step = t_turn_data[1] - t_turn_data[0]
		#t_tomo_ns = np.linspace(t_sep_ns[0], t_sep_ns[-1], n_slices)
		
		recon_id = np.arange(recon_start, recon_stop+1, recon_step)
		print "recon_id ",recon_id

		#plt.ion()

		image_all_a = []
		tomo_time_proj_l = []

		iloop = 0 
		chisq_l = []
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

			image_sum = np.sum(image_data_a)

			#split into array
			image_data_spl = np.split(image_data_a, ndim)
			#transpose to get phase on horizontal axis, energy on vertical
			image_data_tr = np.transpose(image_data_spl)

			#tomo_profile = [sum(d) for d in image_data_tr]
			tomo_time_proj = [sum(d) for d in image_data_spl]
			tomo_mom_proj = [sum(d) for d in image_data_tr]

			image_all_a.append(image_data_tr)		
			tomo_time_proj_l.append(tomo_time_proj)

			#compare with data
			prof_diff = prof_data_norm_a[irc-1][binslowerignore: n_slices - binsupperignore] - tomo_time_proj
			chisq = sum(prof_diff**2)

			#print "chisq ",chisq
			chisq_l.append(chisq)
		
			if iloop == 0:
				tomo_time_proj_a = tomo_time_proj
			else:
				tomo_time_proj_a = np.vstack((tomo_time_proj_a, tomo_time_proj))

			iloop = iloop + 1

				
		chisq_sum = sum(chisq_l)
		chisq_sum_l.append(chisq_sum)

		print "idelay, chisq summed over turns ",idelay1, chisq_sum
	

	if len(idelay_a) > 1:
		

		
		plt.plot(idelay_a, chisq_sum_l,'ko-')
		plt.xlabel('time point shift')
		plt.ylabel('chisq sum')
		plt.ylim(ymin=0)
		plt.axvline(x=idelay, color='gray')
		plt.savefig('timeshiftscan')
		plt.show()
		
	
		sys.exit()

		
	#plt.subplot(311)
	#plt.contour(t_tomo_ns, dp_a, image_all_a[0], 40)
	#plt.subplot(312)
	#plt.contour(t_tomo_ns, dp_a, image_all_a[1], 40)
	#plt.subplot(313)
	#plt.contour(t_tomo_ns, dp_a, image_all_a[2], 40)
	#plt.show()

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
		

	if animate:

		#plt.subplot(211)
		#plt.title('reconstructed turn '+sttomo_time_proj_lr(irc)+ ' time '+time_str+' ms')
		fig, ax = plt.subplots()
		cont = ax.contour(t_tomo_ns, dp_a, image_all_a[0], 40)
		def animate(i):	
			ax.clear()	
			cont = ax.contour(t__ns, dp_a, image_all_a[i], 40)
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

		ipl = 0
		for irc in recon_id:
	
			plt.subplot(211)
			
			plt.title('reconstructed turn '+str(irc)+ ' time '+time_str+' ms')
			ctom = plt.contour(t_tomo_ns, dp_a, image_all_a[ipl], 40)
			#print "ctom ",ctom.levels	
			#sys.exit()

			#plt.contourf(t_gen_ns, dp_bin_ini, Hgrid_inv, ctom.levels)
			#plt.contourf(t_gen_ns, dp_bin_ini, Hgrid_inv)

			plt.plot(t_sep_ns, dp_sep,'r') #bucket
			plt.plot(t_sep_ns, -dp_sep,'r') #bucket
			
			#if run_pyheadtail:
			#	plt.contour(-TT, DPP, HH, levels=[0], colors='magenta')
			#counts, xbins, ybins,_ = plt.hist2d(t_tomo_ns, dp_a, bins=60, norm=LogNorm())
			#plt.xlim(-300,300)
			plt.xlim(1.01*t_sep_ns[0],1.01*t_sep_ns[-1])
			plt.ylim(1.01*bh, -1.01*bh)
			
			plt.axvline(x=t_tomo_ns[int(x0)])
			
			x1,x2,y1,y2 = plt.axis()
			
			plt.xlabel('time [ns]')
			plt.ylabel(r"$\Delta$p/p")
			
			#plt.contour(counts, 20, extent=[xbins.min(), xbins.max(), ybins.min(), ybins.max()])
			#plt.ylabel(r"$\Delta$p [MeV]")

			plt.subplot(212)

			#t_turn_data_ns_shift = t_turn_data_ns #- 20# t_shift_ns

			#plt.plot(t_turn_data_ns, prof_data_norm,'ko-',linewidth=2,label='data')
			plt.plot(t_turn_data_ns, prof_data_norm_a[irc-1],'k.-',linewidth=2,label='data')
			#plt.plot(t_tomo_ns, prof_data_norm,'ko-',linewidth=2,label='data')
			plt.plot(t_tomo_ns, tomo_time_proj_a[ipl],'r-',linewidth=2,label='tomography')
			#plt.xlim(t_bucket[0], t_bucket[-1])

			plt.axvline(x=t_bucket[0], color='gray')
			plt.axvline(x=t_bucket[-1], color='gray')
			
			if binslowerignore:
				plt.axvline(x=t_turn_data_ns[binslowerignore] + 0.1*t_data_step, color='gray', linestyle='--')
			if binsupperignore > 0:
				plt.axvline(x=t_turn_data_ns[n_slices - binsupperignore]-0.1*t_data_step, color='gray', linestyle='--')		
			
			plt.axvline(x=t_synch, color='gray', linestyle = '--')
			#plt.ylim(-1.1*max_dp_sep, 1.1*max_dp_sep) 

			plt.xlim(1.01*t_sep_ns[0],1.01*t_sep_ns[-1])
			#plt.xlim(t_tomo_ns[0],t_tomo_ns[-1])
			#plt.xlim(t_sep_ns[0], t_sep_ns[-1])
			#plt.ylim(-1.1*max_dp_sep, 1.1*max_dp_sep) 
			plt.xlabel('time [ns]')
			plt.ylabel('normalised intensity')
			plt.legend(loc = 'upper left')
			plt.savefig('phasespace_'+str(index_time))
			plt.tight_layout()
			
			
			plt.show()

			ipl = ipl + 1








