from __future__ import division
import pylab as plt
import numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import quadrature
import sys
import os

import matplotlib.cm as cm
import analysedata_180919 as adata
import tomo_defs_180927 as tdefs
import math
from subprocess import call
import time 
import matplotlib.animation as animation

#set up formatting for movie files
#Writer = animation.writers['ffmpeg']
#writer = Wriiter(fps=15, metadata=dict(artist='Me'), bitrate=1800)

tomo_times = np.array([11e-4]) #start time at which do do tomography #8.0e-4
dirname = "../data/2018/9/20/"
#dirname = "../data/2018/9/21/machida_san_foil_escape/" # "../data/2018/9/20/"
filter_string = '' #look for csv files with this string in file name
phi_s_deg = 0 #set synchronous phase
turns_syncfrac = 0.5 #fraction of synchrotron oscillation to include
recon_start = 1 #reconstruct phase space at this turn first
recon_step = 50 #step size terms of turn number between recon_start and recon_stop
recon_stop = 101 #reconstruct phase space at this turn last
animate = True

shift_data_winmin = False #shift minimum of moving window average where window is out-of-bucket duration
shift_data_max = False #shifts data so that maximum is at centre 
shift_data_min = True #shifts data so that minimum is at edge

show_data_switch = False #show raw data
show_time_fit = True #show fit used to calculate ToF
mountain_range_plot = True
plot_data_shift = True

read_rf_waveform = True

npfac = 1.0 #interpolate onto time axis with npfac times the number of points as the data
interpolate_data = True #if False use interpolated data

tomo_outdir = "tomo_output"
files = os.listdir(os.curdir)
if tomo_outdir not in files:
	print "create output directory "+tomo_outdir
	os.mkdir(tomo_outdir)


print "reconstruct at following times [ms] ", 1e3*tomo_times
print "phi_s set to [deg] ",phi_s_deg
if phi_s_deg == 0:
	flattop = True
	flattop_time = 0.3e-3 #estimate
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

if read_rf_waveform:
	rawin2 = raw_input("select RF data file ")

	try:
		i_rf_file_sel = int(rawin2)
	except:
		i_rf_file_sel = None


if phi_s_deg == 0:
	binslowerignore = 0 #add zeroes before profile (earlier time)
	binsupperignore = 0 #add zeroes after profile (later time)
else:
	binslowerignore = int(npfac*35) #add zeroes before profile (earlier time)
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


def show_data(lower_lim_l):
	plt.subplot(211)
	#plt.plot(tdat_rf_ms, sig_rf_sc, 'r')
	plt.plot(tdat_chf_ms, data_chf,'k')
	
	#plt.axvline(x=t_hminus, color='r')
	
	for lower_lim in lower_lim_l:
		plt.axvline(x=1e3*lower_lim, color='r')
		plt.axvline(x=1e3*(lower_lim + nturns*tof_guess), color='r')
	
	plt.xlim(0, tdat_chf_ms[-1])
	plt.ylim(-0.2, 0.2)
	plt.ylabel('bunch monitor (all)')

	plt.subplot(212)
	plt.plot(tdat_chf_ms, data_chf,'k')
	if read_rf_waveform:
		plt.plot(tdat_rf_ms, sig_rf_sc, 'r')
	plt.xlim(1e3*lower_lim, 1e3*(lower_lim + 20*tof_guess))
	plt.ylabel('bunch monitor (first 20 turns)')
	plt.xlabel('time [ms]')
	plt.savefig('bunchmonitor_data')
	plt.show()
	
	
def calc_tof(lower_lim, points_per_turn):
	"""Calculate TOF by fitting the time of successive peaks"""
	
	upper_lim_fit = lower_lim + nturns_fit*tof_guess
	
	find_hminus = False
	_, pmx, _, i_peak_l, _, _ = adata.find_peaks(tdat_chf, data_chf, points_per_turn, find_hminus, x_low = lower_lim, x_high = upper_lim_fit, find_maxima = True)
	
	tp = np.arange(nturns_fit)
	npoly = 1
	fitpmx = np.polyfit(tp, pmx, npoly)
	tof_fit_mean = fitpmx[npoly - 1]
	tof_fit_ns = 1e9*tof_fit_mean
	polyfitpmx = np.poly1d(fitpmx)
	tof_fit_a = polyfitpmx(tp)

	tof_diff = pmx - polyfitpmx(tp)
	fit_tof_diff = np.polyfit(tp, tof_diff, 2)
	polyfitdiff = np.poly1d(fit_tof_diff)
	fitdiff_a = polyfitdiff(tp) - polyfitdiff[0]
	
	#ampt = max(tof_diff) - min(tof_diff)
	#timeosc = [ampt*np.sin(0.25*2*np.pi*Qs*i) - 0.5*ampt for i in np.arange(nturns_fit)]
	#timeosc = timeosc - timeosc[0]
	
	if show_time_fit:
		pmx_a = np.array(pmx)
		plt.subplot(211)
		plt.plot(1e6*pmx_a,'k-',linewidth=2,label='data')
		plt.plot(1e6*polyfitpmx(tp),'r',label='fit')
		plt.ylabel('turn-by-turn times [us]')
		plt.title('slope = '+str(tof_fit_ns)+ ' ns')
		plt.legend(loc = 'upper left')
		plt.subplot(212)
		plt.plot(1e6*(pmx_a - polyfitpmx(tp)))
		plt.plot(1e6*polyfitdiff(tp) ,'r')
		plt.xlabel('turn')
		plt.ylabel('difference [us]')
		plt.savefig('time_fit')
		plt.show()
	
	#print (tof_fit+1e9*fitdiff_a[0])/time_interval , (tof_fit+1e9*fitdiff_a[-1])/time_interval
			
	return tof_fit_mean, i_peak_l, fitdiff_a, tof_fit_a
		
	


def sort_data(istart, fitdiff_a, tof_fit_a, nturns,n_slices, interpolate_data = True, nshift = 0):
	"""Divide bunch monitor data into individual turns. There is the option to interpolate onto regular time points. """
	
	#iend = istart + nturns*points_per_turn
	 
	
	tstart = tof_fit_a[0] #tdat_chf[istart] #t_dat_all[0]
	
	#plt.plot(t_dat_all, sig_dat_all)
	#plt.show()
	#sys.exit()
	#use ToF to define time.
	#t_dat_adj = tstart + np.linspace(0, nturns*tof, nturns*n_slices)

	if npfac > 1 and not interpolate_data:
		print "must interpolate data if npfac != 1"
		sys.exit()
	
	if not interpolate_data:
		print "interpolate = False option deprecated. Code needs to be checked"
		sys.exit()
		t_dat_all = np.array(tdat_chf[istart:iend])
		sig_dat_all = np.array(data_chf[istart:iend])
		
		t_dat_split = np.split(t_dat_all, nturns)
		sig_dat_split = np.split(sig_dat_all, nturns)
		
		#turn data overlaps on final/first point
		index_turn_start_l= [istart]	
		for it in range(nturns):
			if it > 0:
				for j1, t1 in enumerate(tdat_chf[istart:iend]):
					if  abs(tstart + it*tof_fit + fitdiff_a[it] - t1) < 3e-8:
						index_turn_start = istart +j1
						index_turn_start_l.append(index_turn_start)
						break
		
		t_dat_split = []
		sig_dat_split = []
		for it, ist in enumerate(index_turn_start_l):
			t_dat_split.append(tdat_chf[ist:ist+points_per_turn])
			sig_dat_split.append(data_chf[ist:ist+points_per_turn])
						
		sig_dat_split = np.array(sig_dat_split)	
		
		t_turn_data = t_dat_split[0]
	else:
		#interpolate onto time base Trev/n_slices
		
		#nturns = 3

		t_turn_data = np.linspace(0, tof_fit_mean, n_slices)

		new_code = True
		if new_code:		
			t_dat_adj = []
			#tstart_update  = tstart

			for it in range(1, nturns+1):
				time1 = tof_fit_a[it-1]
				time2 = tof_fit_a[it]
				t_oneturn_data = np.linspace(time1, time2 ,n_slices) #fitdiff_a[it]
				t_dat_adj  = np.concatenate((t_dat_adj, t_oneturn_data))

			t_dat_adj = t_dat_adj + time_interval*nshift
				
		else:
			
			t_dat_adj = []
			for it in range(nturns):
				t_this_turn = tstart + it*tof_fit_mean + t_turn_data + fitdiff_a[it] 
				t_dat_adj  = np.concatenate((t_dat_adj, t_this_turn))
	

		


		t_adj_last = t_dat_adj[-1]
		i1 = 0
		for i1, t1 in enumerate(tdat_chf):
			if t1 >= t_adj_last:
				index_t_last = i1
				break	

	
		t_dat_orig = np.array(tdat_chf[istart:index_t_last])
		sig_dat_orig = np.array(data_chf[istart:index_t_last])
		
		sig_interp = np.interp(t_dat_adj, t_dat_orig, sig_dat_orig)
		t_dat_split = np.split(t_dat_adj, nturns)
		sig_dat_split = np.split(sig_interp, nturns)
		sig_interp = np.interp(t_dat_adj, t_dat_orig, sig_dat_orig)
		t_dat_split = np.split(t_dat_adj, nturns)
		sig_dat_split = np.split(sig_interp, nturns)
		

		#plt.plot(t_dat_adj, sig_interp,'r-')
		#plt.plot(t_dat_orig, sig_dat_orig,'ko')
		#plt.axvline(x=t_dat_all[-1])
		#plt.show()
		#sys.exit()
	
	it = 0
	data_win_sig_sel = []
	for data_win in sig_dat_split[:nturns]:
		data_win_sig =  max(data_win) - data_win
		data_win_sig_sel.append(data_win_sig)
		
		if mountain_range_plot:
			plt.plot(range(1,n_slices+1), 150*data_win_sig + it,'.-')
		it = it + 1
	
	if mountain_range_plot:		
		plt.ylabel('turn number')
		plt.axvline(x=binslowerignore+0.5)
		plt.axvline(x=n_slices - binsupperignore+0.5)
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
	
	
	
#set parameters
phi_s = phi_s_deg*math.pi/180 
kindex = 7.6
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

lower_lim_l = []
for tt in tomo_times:
	icrude = index_zero + int(round(tt/time_interval))
	t_crude = tdat_chf[icrude]
	lower_lim_l.append(t_crude)


#find zero crossing times in rf data
if read_rf_waveform:
	zero_cross_ind = []
	t_zero_cross = []
	plot_zero_calc = False
	for ind in range(icrude, len(tdat_rf)):
		if sig_rf[ind-1] < 0.0 and sig_rf[ind] > 0.0:
			slope = (sig_rf[ind] - sig_rf[ind-1])
			offset = sig_rf[ind] - slope
			zerot_fac = -offset/slope
			
			print tdat_rf[ind-1] + zerot_fac*time_interval
			if plot_zero_calc:
				xp = np.linspace(0,1,10)
				yp = slope*xp + offset
				plt.plot([0,1],[sig_rf[ind-1],sig_rf[ind]],'ko')
				plt.plot(xp, yp, 'r-')
				plt.axvline(x=zerot_fac)	
				plt.show()

			t_zero_cross.append(tdat_rf[ind-1] + zerot_fac*time_interval)

			zero_cross_ind.append(ind)
		if len(zero_cross_ind) > 400:
			break
	
	t_zero_cross = np.array(t_zero_cross)
	
	zero_cross_sep = [zero_cross_ind[i] - zero_cross_ind[i-1] for i in range(1, len(zero_cross_ind))]
	
	sig_zero_cross_sep = [sig_rf[zero_cross_ind[i]] - sig_rf[zero_cross_ind[i]-1] for i in range(1, len(zero_cross_ind))]
	t_zero_cross_sep = [t_zero_cross[i] - t_zero_cross[i-1] for i in range(1, len(t_zero_cross))]
	print "zero cross sep ",zero_cross_sep
 	print "time_zeros ",sig_zero_cross_sep
	print "t_zero_cross_sep ",t_zero_cross_sep

	#plt.plot(t_zero_cross_sep, 'ko')
	#plt.show()
	#sys.exit()
print "specified initial times ",tomo_times
print "chosen initial times ",lower_lim_l


#loop over times, do tomography on each loop
for index_time in range(len(tomo_times)):

	#number of turns specified in terms of sync period
	nturns = int(turns_syncfrac*Ts_a[index_time])
	
	if show_data_switch and index_time == 0:
		show_data(lower_lim_l)
		
	
	p0_mev = 1e-6*np.sqrt(Etot_a[index_time]**2 - PROTON_MASS**2)
	bh = ((2*Qs_a)/(harmonic*abs(eta_a[index_time]))) #bucket height
	
	time_str = str(1e3*tomo_times[index_time])

	lower_lim = lower_lim_l[index_time]
		
	points_per_turn = int(round(tof_guess/time_interval)) 
	
	
	#calculate TOF
	nturns_fit = 400
	tof_fit_mean, i_peak_l, fitdiff_a, tof_fit_a = calc_tof(lower_lim, points_per_turn)
	istart = i_peak_l[0]
	tof = tof_fit_mean

	tof_fit_a = t_zero_cross

	T = tof/harmonic#2*np.pi/(harmonic*omegarev) #rf period
	T_rf_ns = 1e9*T
	omega_rev = 2*math.pi/T
	omega_rf = omega_rev*harmonic

	#update points per turnt_turn_data
	points_per_turn = int(round(tof/time_interval)) 
	n_slices = int(npfac*points_per_turn)
	#n_slices_active = n_slices-(binslowerignore + binsupperignore)
	

	#bucket phase calculation
	if phi_s_deg == 0:
		phi_l1 = math.pi
		phi_l2 = -math.pi
	else:
		phi_l1, phi_l2 = phase_extrema(phi_s)

	bucket_phase_width = abs(phi_l1 - phi_l2)
	bucket_dur_ns = T_rf_ns*(bucket_phase_width/(2*math.pi))
	outbucket_dur_ns = T_rf_ns - bucket_dur_ns
	print "phase extrema phi_l1, phi_l2 ",180*phi_l1/math.pi, 180*phi_l2/math.pi
	print "bucket phase width [rad, deg]",bucket_phase_width,180*bucket_phase_width/math.pi
	print "bucket duration, out-of-bucket calc from phase [ns] ",bucket_dur_ns, outbucket_dur_ns
	print "time phi_l1, phi_l2 ", T_rf_ns*(phi_l1/(2*math.pi)), T_rf_ns*(phi_l2/(2*math.pi))
	
	tadj = T_rf_ns*(phi_l2/(2*math.pi)) - outbucket_dur_ns 


	synch_index = int(round(n_slices*(phi_s - phi_l2)/(2*math.pi)))
	print "synch_index ",synch_index
	
	
	n_slices_active = int(n_slices*bucket_dur_ns/T_rf_ns) #number of slices in bucket
	n_slices_outb = n_slices - n_slices_active
	print "bucket points, out-of-bucket points ",n_slices_active, n_slices_outb

	if shift_data_winmin and phi_s_deg == 0:
		print "shift_data_winmin must be False if phi_s_deg is zero"

	if phi_s_deg != 0:
		binslowerignore = n_slices_outb
		binsupperignore = n_slices - (binslowerignore + n_slices_active)
		
	#sort data into turns covering roughly a synchrotron period
	nturns_sync = int(round(Ts_a[index_time]))


	prof_data_a, t_turn_data = sort_data(istart, fitdiff_a, tof_fit_a, nturns_sync, n_slices, interpolate_data)	
	
	ishift = 0
	if shift_data_min or shift_data_max or shift_data_winmin or ishift != 0:
		prof_data_tr = np.transpose(prof_data_a)
		data_integ_0 = [sum(d) for d in prof_data_tr]


		imax_col = data_integ_0.index(max(data_integ_0))		
		imin_col = data_integ_0.index(min(data_integ_0))

		
		if shift_data_winmin:
			runmean = tdefs.running_mean_wrap(data_integ_0, n_slices_outb)
			imax_wincol = np.argmax(runmean)
			imin_wincol = np.argmin(runmean)

		
		if shift_data_min:
			if imin_col != 0:
				ishift = int(imin_col/npfac)
			print "imin_col ",imin_col
		elif shift_data_max:
			if imax_col - 0.5*n_slices!= 0:
				ishift = int((imax_col - 0.5*n_slices)/npfac)
			print "imax_col ", imax_col
		elif shift_data_winmin:
				ishift = int(imin_wincol/npfac) - int(0.5*n_slices_outb)
		
		print "ishift ",ishift
		if ishift != 0:

			istart = istart + ishift

			prof_data_a, t_turn_data = sort_data(istart, fitdiff_a, tof_fit_a, nturns_sync, n_slices, interpolate_data, nshift= ishift)
			prof_data_tr = np.transpose(prof_data_a)
			data_integ_1 = [sum(d) for d in prof_data_tr]
			
			xd = np.arange(1,n_slices+1)
			
			t_turn_data_ns = 1e9*t_turn_data + tadj

			if plot_data_shift:

				if shift_data_winmin:
					plt.subplot(211)

				plt.plot(xd, data_integ_0,'ko-',label='before shift')
				plt.plot(xd, data_integ_1,'r',label='after shift')
				plt.axvline(x = synch_index, color='gray')
				plt.xlim(xmax=len(data_integ_0))
				
				plt.ylabel('integrated signal')
				plt.title('signal integrated over turns specified')

				if shift_data_winmin:
					plt.subplot(212)
					plt.plot(runmean, 'ko-')
				plt.xlabel('slice number')
				
				plt.xlim(xmax=len(data_integ_0))

				

				plt.legend()
				plt.show()
			
		else:
			t_turn_data_ns = 1e9*t_turn_data

	else:
		t_turn_data_ns = 1e9*t_turn_data + tadj

	outbucket_dur_ns
	print "write ", nturns, " turns to file"
	f1 = open('kurridata.txt', 'w')
	for it in range(nturns):
		yprof = prof_data_a[it]
		yprof_norm = yprof/max(yprof)
		#plt.plot(yprof,'ko-')
		
		for y in yprof:
			print >>f1,y
	f1.close()
	
		
	#omega_rf = 2*math.pi*f_rf_lookup(tomo_times[index_time])
	npts = 1000
	a1 = 0.5*harmonic*omega_rf*eta_a[index_time]
	a2 = omega_rf*V_rf/(2*math.pi*(beta_a[index_time]**2)*Etot_a[index_time])
	dp_sep, t_sep_ns = tdefs.separatrix(phi_s, eta_a[index_time], T, harmonic, npts, a1, a2)
	
	phase_sep = 2*math.pi*(t_sep_ns)/(1e9*T)
	t_bucket = np.linspace(t_sep_ns[0], t_sep_ns[-1], n_slices)

	#t_data = np.linspace(-0.5*T_rf_ns, 0.5*T_rf_ns, n_slices) #data fills rf period

	t_phi_sep_ns = phi_s*(1e9*T)/(2*math.pi)
	t_phi_l2 = T_rf_ns*(phi_l2/(2*math.pi)) 
		
	t_sep_ns = np.insert(t_sep_ns,0,t_phi_l2)
	dp_sep = np.insert(dp_sep,0,0)


	phase_sep_deg = (180/math.pi)*phase_sep
	sync_frac =  (t_phi_sep_ns- t_sep_ns[0])/bucket_dur_ns
	print "synch_frac ",sync_frac
	sync_bin = n_slices_active*sync_frac
	print "sync bin ",sync_bin
	
	show_sep = False
	if show_sep:
		plt.plot(phase_sep_deg, dp_sep, 'ko-')
		plt.axvline(x=180*phi_l1/math.pi)
		plt.axvline(x=180*phi_l2/math.pi)
		plt.show()

	#this is not a real Bdot - it is a parameter used to specify phis in the tomo code
	Bdot = (V_rf*np.sin(phi_s))/(bending_radius*2*math.pi*Rnom)

	from subprocess import call
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


		call(["tomo"])

		files = os.listdir(os.curdir)
		for fn in files:
			if fn[-4:] == 'data':
				os.rename(fn,tomo_outdir+"/"+fn)
				 
		

	def read_data_file(filen):
	
		f1 = open(filen,'r')
		data_l = []
		for data in f1:
			data_l.append(float(data))
		data_a = np.array(data_l)

		return data_a

	
	show_tomo_result = True
	if show_tomo_result:

		#read output parameters from plotinfo.data
		plotinfo = tdefs.get_plotinfo(input_file=tomo_outdir+"/plotinfo.data")
	
		x0 = plotinfo['xat0'][0]
		y0 = plotinfo['yat0'][0]
		dEbin =  plotinfo['dEbin'][0]
		dtbin =  plotinfo['dtbin'][0]
		phi_s_out = plotinfo['phi_s'][0] #synch phase

		print "dtbin ",dtbin

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
	t_a_ns = 1e9*dtbin*(np.array(range(ndim)) - x0) + t_synch
	
	t_shift_ns = -1e9*dtbin*(x0 - 0.5*ndim) #+ t_synch
	
	print "time shift applied ",t_shift_ns
	
	t_data_step = t_turn_data[1] - t_turn_data[0]
	#t_a_ns = np.linspace(t_sep_ns[0], t_sep_ns[-1], n_slices)
	
	recon_id = np.arange(recon_start, recon_stop+1, recon_step)
	print "recon_id ",recon_id

	#plt.ion()

	image_all_a = []
	tomo_time_proj_l = []
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
		print "time to read file ",time_fread1 - time_fread0

		#split into array
		image_data_spl = np.split(image_data_a, ndim)
		#transpose to get phase on horizontal axis, energy on vertical
		image_data_tr = np.transpose(image_data_spl)

		#tomo_profile = [sum(d) for d in image_data_tr]
		tomo_time_proj = [sum(d) for d in image_data_spl]
		tomo_mom_proj = [sum(d) for d in image_data_tr]

		image_all_a.append(image_data_tr)		
		tomo_time_proj_l.append(tomo_time_proj)
		
		#t_turn_data_ns = 1e9*t_turn_data - 0.5*T_rf_ns 
	
	#plt.subplot(311)
	#plt.contour(t_a_ns, dp_a, image_all_a[0], 40)
	#plt.subplot(312)
	#plt.contour(t_a_ns, dp_a, image_all_a[1], 40)
	#plt.subplot(313)
	#plt.contour(t_a_ns, dp_a, image_all_a[2], 40)
	#plt.show()

	
	if animate:

		#plt.subplot(211)
		#plt.title('reconstructed turn '+sttomo_time_proj_lr(irc)+ ' time '+time_str+' ms')
		fig, ax = plt.subplots()
		cont = ax.contour(t_a_ns, dp_a, image_all_a[0], 40)
		def animate(i):	
			ax.clear()	
			cont = ax.contour(t_a_ns, dp_a, image_all_a[i], 40)
			ax.plot(t_sep_ns, dp_sep,'r') #bucket
			ax.plot(t_sep_ns, -dp_sep,'r') 
			#cont.set_zdata(image_all_a[1])
			return cont


		Nt = len(image_all_a)
		anim = animation.FuncAnimation(fig, animate, frames=Nt, interval=500, blit=False)
		#anim.save("tomo.mp4", writer=writer)
		plt.show()

	else:

		ipl = 0
		for irc in recon_id:
	
			plt.subplot(211)
			
			plt.title('reconstructed turn '+str(irc)+ ' time '+time_str+' ms')
			ctom = plt.contour(t_a_ns, dp_a, image_all_a[ipl], 40)
			#print "ctom ",ctom.levels	
			#sys.exit()

			#plt.contourf(t_gen_ns, dp_bin_ini, Hgrid_inv, ctom.levels)
			#plt.contourf(t_gen_ns, dp_bin_ini, Hgrid_inv)

			plt.plot(t_sep_ns, dp_sep,'r') #bucket
			plt.plot(t_sep_ns, -dp_sep,'r') #bucket
			
			#if run_pyheadtail:
			#	plt.contour(-TT, DPP, HH, levels=[0], colors='magenta')
			#counts, xbins, ybins,_ = plt.hist2d(t_a_ns, dp_a, bins=60, norm=LogNorm())
			#plt.xlim(-300,300)
			plt.xlim(1.01*t_sep_ns[0],1.01*t_sep_ns[-1])
			plt.ylim(1.01*bh, -1.01*bh)
			
			plt.axvline(x=t_a_ns[int(x0)])
			
			x1,x2,y1,y2 = plt.axis()
			
			plt.xlabel('time [ns]')
			plt.ylabel(r"$\Delta$p/p")
			
			#plt.contour(counts, 20, extent=[xbins.min(), xbins.max(), ybins.min(), ybins.max()])
			#plt.ylabel(r"$\Delta$p [MeV]")

			plt.subplot(212)

			prof_data_norm = prof_data_a[irc-1]/max(prof_data_a[0])
			tomo_t_proj_norm = tomo_time_proj_l[ipl]/max(tomo_time_proj_l[ipl])

			t_turn_data_ns_shift = t_turn_data_ns #+ t_shift_ns
			plt.plot(t_turn_data_ns, prof_data_norm,'ko-',linewidth=2,label='data')
			#plt.plot(t_a_ns, prof_data_norm,'ko-',linewidth=2,label='data')
			plt.plot(t_a_ns, tomo_t_proj_norm,'ro-',linewidth=2,label='tomography')
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
			#plt.xlim(t_a_ns[0],t_a_ns[-1])
			#plt.xlim(t_sep_ns[0], t_sep_ns[-1])
			#plt.ylim(-1.1*max_dp_sep, 1.1*max_dp_sep) 
			plt.xlabel('time [ns]')
			plt.ylabel('normalised intensity')
			plt.legend(loc = 'upper left')
			plt.savefig('phasespace_'+str(index_time))
			plt.tight_layout()
			
			
			plt.show()

			ipl = ipl + 1








