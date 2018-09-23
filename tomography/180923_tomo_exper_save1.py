from __future__ import division
import pylab as plt
import numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import quadrature
import sys
import os
import matplotlib.cm as cm
import analysedata_150619 as adata
import tomo_defs_180911 as tdefs
import math
from subprocess import call

tomo_times = np.array([13.6e-4]) #start time at which do do tomography
#dirname = "../data/2018/9/20/"
dirname = "../data/2018/9/21/machida_san_foil_escape/" # "../data/2018/9/20/"
filter = 'wb' #look for csv files with this string in file name
phi_s_deg = 0 #set synchronous phase
turns_syncfrac = 0.5 #fraction of synchrotron oscillation to include
recon_turn = 1 #reconstruct phase space at this turn
show_data_switch = True #show raw data
show_time_fit = False #show fit used to calculate ToF
mountain_range_plot = True

npfac = 1.0 #interpolate onto time axis with npfac times the number of points as the data
interpolate_data = True #if False use interpolated data

shift_data_max = False #shifts data so that maximum is at centre 
shift_data_min = False #shifts data so that minimum is at edge

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
	if 'csv' in f1 and filter in f1:
		print "index ",icsv," file ", f1
		csv_files.append(f1)
		icsv = icsv + 1

rawin1 =  raw_input("select file ")
try:
	i_file_sel = int(rawin1)
except:
	print "input is not integer"

if phi_s_deg == 0:
	binslowerignore = 0 #add zeroes before profile (earlier time)
	synch_index = int(npfac*75)
	binsupperignore = 0 #add zeroes after profile (later time)
else:
	binslowerignore = int(npfac*35) #add zeroes before profile (earlier time)
	synch_index = int(npfac*40)
	binsupperignore = int(npfac*19) #add zeroes after profile (later time)

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
	plt.xlim(1e3*lower_lim, 1e3*(lower_lim + 20*tof_guess))
	plt.ylabel('bunch monitor (first 20 turns)')
	plt.xlabel('time [ms]')
	plt.savefig('bunchmonitor_data')
	plt.show()
	
	
def calc_tof(lower_lim, nturn_points):
	"""Calculate TOF by fitting the time of successive peaks"""
	
	upper_lim_fit = lower_lim + nturns_fit*tof_guess
	
	find_hminus = False
	_, pmx, _, i_peak_l, _, _ = adata.find_peaks_new(tdat_chf, data_chf, nturn_points, find_hminus, x_low = lower_lim, x_high = upper_lim_fit, find_maxima = True)
	
	tp = np.arange(nturns_fit)
	fitpmx = np.polyfit(tp, pmx, 1)
	tof_fit = fitpmx[0]
	tof_fit_ns = 1e9*tof_fit
	polyfitpmx = np.poly1d(fitpmx)
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
			
	return tof_fit, i_peak_l, fitdiff_a
		
	


def sort_data(istart, fitdiff_a, interpolate_data = True):
	"""Divide bunch monitor data into individual turns. There is the option to interpolate onto regular time points. """
	
	
	#iend = istart + nturns*nturn_points
	 
	
	tstart = tdat_chf[istart] #t_dat_all[0]
	
	#plt.plot(t_dat_all, sig_dat_all)
	#plt.show()
	#sys.exit()
	#use ToF to define time.
	#t_dat_adj = tstart + np.linspace(0, nturns*tof, nturns*n_slices)

	if npfac > 1 and not interpolate_data:
		print "must interpolate data if npfac != 1"
		sys.exit()
	
	if not interpolate_data:
		t_dat_all = np.array(tdat_chf[istart:iend])
		sig_dat_all = np.array(data_chf[istart:iend])
		
		t_dat_split = np.split(t_dat_all, nturns_fit)
		sig_dat_split = np.split(sig_dat_all, nturns_fit)
		
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
			t_dat_split.append(tdat_chf[ist:ist+nturn_points])
			sig_dat_split.append(data_chf[ist:ist+nturn_points])
						
		sig_dat_split = np.array(sig_dat_split)	
		
		t_turn_data = t_dat_split[0]
	else:
		#interpolate onto time base Trev/n_slices
		
		t_turn_data = np.linspace(0, tof_fit, n_slices)
		t_dat_adj = []
		for it in range(nturns):
			t_this_turn = tstart + it*tof_fit + t_turn_data + fitdiff_a[it] 
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

		
			
	it = 0
	data_win_sig_sel = []
	for data_win in sig_dat_split[:nturns]:
		data_win_sig =  max(data_win) - data_win
		data_win_sig_sel.append(data_win_sig)
		
		if mountain_range_plot:
			plt.plot(range(1,n_slices+1), data_win_sig + it*0.008,'.-')
		it = it + 1
	
	if mountain_range_plot:		
		plt.axvline(x=binslowerignore+0.5)
		plt.axvline(x=n_slices - binsupperignore+0.5)
		plt.savefig('mountainrange')
		plt.show()
		
	
	prof_data_a = np.array(data_win_sig_sel)	


	
	return prof_data_a, t_turn_data




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
for time in tomo_times:
	icrude = index_zero + int(round(time/time_interval))
	t_crude = tdat_chf[icrude]
	lower_lim_l.append(t_crude)

print "specified initial times ",tomo_times
print "chosen initial times ",lower_lim_l


for index_time in range(len(tomo_times)):


	
	nturns = int(turns_syncfrac*Ts_a[index_time])
	
	if show_data_switch and index_time == 0:
		show_data(lower_lim_l)
	
	p0_mev = 1e-6*np.sqrt(Etot_a[index_time]**2 - PROTON_MASS**2)
	
	time_str = str(1e3*tomo_times[index_time])

	
	lower_lim = lower_lim_l[index_time]
	
	
	nturn_points = int(round(tof_guess/time_interval)) 
	nturns_fit = 400
	#calculate TOF
	tof_fit, i_peak_l, fitdiff_a = calc_tof(lower_lim, nturn_points)
	
	istart = i_peak_l[0]
	
	#update points per turn
	nturn_points = int(round(tof_fit/time_interval)) 
	n_slices = int(npfac*nturn_points)
	n_slices_active = n_slices-(binslowerignore + binsupperignore)
	print  "nturn_points ",nturn_points
	
		
	prof_data_a, t_turn_data = sort_data(istart, fitdiff_a, interpolate_data)	
	
	
	if shift_data_min or shift_data_max:
		prof_data_tr = np.transpose(prof_data_a)
		data_integ_0 = [sum(d) for d in prof_data_tr]
		
		imax_col = data_integ_0.index(max(data_integ_0))
		imin_col = data_integ_0.index(min(data_integ_0))
		
		ishift = 0
		if shift_data_min:
			if imin_col != 0:
				ishift = int(imin_col/npfac)
			print "imin_col ",imin_col, ""
		else:
			if imax_col - 0.5*n_slices!= 0:
				ishift = int((imax_col - 0.5*n_slices)/npfac)
		
			print "max_col - 0.5*n_slices) ", int(imax_col - 0.5*n_slices)
			
			
		if ishift != 0:
			istart = istart + ishift
			prof_data_a, t_turn_data = sort_data(istart, fitdiff_a, interpolate_data)
			prof_data_tr = np.transpose(prof_data_a)
			data_integ_1 = [sum(d) for d in prof_data_tr]
			
			xd = np.arange(1,n_slices+1)
			
			plt.plot(xd, data_integ_0,'k',label='before shift')
			plt.plot(xd, data_integ_1,'r',label='after shift')
			plt.axvline(x = 0.5*n_slices, color='gray')
			plt.xlim(xmax=len(data_integ_0))
			
			plt.xlabel('slice number')
			plt.ylabel('integrated signal')
			plt.title('signal integrated over turns specified')
			
			plt.legend()
			plt.show()

	
	print "write ", nturns, " turns to file"
	f1 = open('kurridata.txt', 'w')
	for it in range(nturns):
		yprof = prof_data_a[it]
		yprof_norm = yprof/max(yprof)
		#plt.plot(yprof,'ko-')
		
		for y in yprof:
			print >>f1,y
	f1.close()
	
	
	tof = tof_fit	
	omegarev = beta_a[index_time]*SPEED_OF_LIGHT/Rnom #angular rev. frequency
	T = tof/harmonic#2*np.pi/(harmonic*omegarev) #rf period
	T_rf_ns = 1e9*T
	bh = ((2*Qs_a)/(harmonic*abs(eta_a[index_time])))
	omega_rf = omegarev*harmonic
	a1 = 0.5*harmonic*omega_rf*eta_a[index_time]
	a2 = omega_rf*V_rf/(2*math.pi*(beta_a[index_time]**2)*Etot_a[index_time])

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
		#phi_l2 = tdefs.turning_point2(phi_s, a1, a2, 1000)
		
	print "phase extrema along separatrix, phi_l1, phi_l2 ",phi_l1, phi_l2
	print "in degrees ", 180*phi_l1/math.pi, 180*phi_l2/math.pi

	bucket_phase_width = abs(phi_l1 - phi_l2)

	npts = 1000
	dp_sep, t_phi_ns = tdefs.separatrix(phi_s, eta_a[index_time], T, harmonic, npts, a1, a2)
	phase_sep = 2*math.pi*(t_phi_ns)/(1e9*T)


	t_bucket = np.linspace(t_phi_ns[0], t_phi_ns[-1], n_slices)
	#t_data = np.linspace(-0.5*T_rf_ns, 0.5*T_rf_ns, n_slices) #data fills rf period

	print "bucket phase width [rad, deg]",bucket_phase_width,180*bucket_phase_width/math.pi
	print "bucket duration calc from phase [ns] ", 

	bucket_dur_ns = T_rf_ns*(bucket_phase_width/(2*math.pi))
	

	#bucket_dur_ns = t_phi_ns[-1] - t_phi_ns[0]
	outbucket_dur_ns = T_rf_ns - bucket_dur_ns
	print "rf period, bucket duration ",T_rf_ns, bucket_dur_ns, outbucket_dur_ns, 
	print "fraction time out of bucket ",outbucket_dur_ns/T_rf_ns


	bucket_points = int(n_slices*bucket_dur_ns/T_rf_ns)

	print "bucket points, out-of-bucket points ",bucket_points, n_slices - bucket_points


	t_phi_sep_ns = phi_s*(1e9*T)/(2*math.pi)

	phase_sep_deg = (180/math.pi)*phase_sep

	show_sep = False
	if show_sep:
		plt.plot(phase_sep_deg, dp_sep, 'ko-')
		plt.axvline(x=180*phi_l1/math.pi)
		plt.axvline(x=180*phi_l2/math.pi)
		plt.show()
	#sys.exit()


	print "bucket time limits ",t_phi_ns[0], t_phi_ns[-1]
	print "t_phi_sep_ns ",phi_s_deg, t_phi_sep_ns

	sync_frac =  (t_phi_sep_ns- t_phi_ns[0])/bucket_dur_ns
	print "synch_frac ",sync_frac
	sync_bin = bucket_points*sync_frac
	print "sync bin ",sync_bin


	#this is not a real Bdot - it is a parameter used to specify phis in the tomo code
	Bdot = (V_rf*np.sin(phi_s))/(bending_radius*2*math.pi*Rnom)


	from subprocess import call
	run_tomo_code = True
	if run_tomo_code:
	
		input_settings = tdefs.get_input("input_default.dat")
	
		#synch_index_ideal = 
		recon_start = recon_turn
		recon_stop = recon_turn
	
		binwidth_time = time_interval
	
		fname = 'input_v2.dat'

		descrip = 'KURRI'
		datafile = 'kurridata.txt'
		input_settings["numframes"] = nturns
		input_settings["numbins"] = n_slices
		input_settings["synchbin"] = synch_index#int(0.5*n_slices)
		input_settings["framebinwidth"] = binwidth_time/npfac	
		input_settings["binslowerignore"] = binslowerignore
		input_settings["binsupperignore"] = binsupperignore
		input_settings["filmstart"] = recon_start
		input_settings["filmstop"] = recon_stop
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
		for file1 in temp_files:
			if os.path.isfile(file1):
				os.remove(file1)
				print "delete ",file1

		call(["tomo"])
	

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
		plotinfo = tdefs.get_plotinfo()
	
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
		
			profiles_a = read_data_file(fname_prof)

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
		E0, gamma0, eta0, beta0, omegarev0, Vpeak = tdefs.set_param(inputp)
		T = 2*np.pi/omegarev0 #rf period

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
	#t_a_ns = np.linspace(t_phi_ns[0], t_phi_ns[-1], n_slices)
	
	recon_id = np.arange(recon_start, recon_stop+1)


	#plt.ion()

	for irc in recon_id:

		#if run_pyheadtail:

		#counts, xbins, ybins,_ = plt.hist2d(bunch_t_ini_ns, bunch_dp_ini, bins=60, norm=LogNorm())

		print "recon at turn ",irc
		#process main tomography output file
		#fname_img = 'image001.data'
		fname_img = 'image000'[:-len(str(irc))]+ str(irc) + '.data'
		print "fname_img ",fname_img
		image_data_a = read_data_file(fname_img) #1D data
		#split into array
		image_data_spl = np.split(image_data_a, ndim)
		#transpose to get phase on horizontal axis, energy on vertical
		image_data_tr = np.transpose(image_data_spl)

		#tomo_profile = [sum(d) for d in image_data_tr]
		tomo_time_proj = [sum(d) for d in image_data_spl]
		tomo_mom_proj = [sum(d) for d in image_data_tr]		
		
		t_turn_data_ns = 1e9*t_turn_data - 0.5*T_rf_ns 
		
		plt.subplot(211)
		
		plt.title('reconstructed turn '+str(recon_turn)+ ' time '+time_str+' ms')
		ctom = plt.contour(t_a_ns, dp_a, image_data_tr, 40)
		#print "ctom ",ctom.levels	
		#sys.exit()

		#plt.contourf(t_gen_ns, dp_bin_ini, Hgrid_inv, ctom.levels)
		#plt.contourf(t_gen_ns, dp_bin_ini, Hgrid_inv)

		plt.plot(t_phi_ns, dp_sep,'r') #bucket
		plt.plot(t_phi_ns, -dp_sep,'r') #bucket
		
		#if run_pyheadtail:
		#	plt.contour(-TT, DPP, HH, levels=[0], colors='magenta')
		#counts, xbins, ybins,_ = plt.hist2d(t_a_ns, dp_a, bins=60, norm=LogNorm())
		#plt.xlim(-300,300)
		plt.xlim(1.01*t_phi_ns[0],1.01*t_phi_ns[-1])
		plt.ylim(1.01*bh, -1.01*bh)
		
		plt.axvline(x=t_a_ns[int(x0)])
		
		x1,x2,y1,y2 = plt.axis()
		
		plt.xlabel('time [ns]')
		plt.ylabel(r"$\Delta$p/p")
		
		#plt.contour(counts, 20, extent=[xbins.min(), xbins.max(), ybins.min(), ybins.max()])
		#plt.ylabel(r"$\Delta$p [MeV]")

		plt.subplot(212)

		prof_data_norm = prof_data_a[irc-1]/max(prof_data_a[0])
		tomo_t_proj_norm = tomo_time_proj/max(tomo_time_proj)

		t_turn_data_ns_shift = t_turn_data_ns + t_shift_ns
		plt.plot(t_turn_data_ns_shift, prof_data_norm,'ko-',linewidth=2,label='data')
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

		plt.xlim(1.01*t_phi_ns[0],1.01*t_phi_ns[-1])
		#plt.xlim(t_a_ns[0],t_a_ns[-1])
		#plt.xlim(t_phi_ns[0], t_phi_ns[-1])
		#plt.ylim(-1.1*max_dp_sep, 1.1*max_dp_sep) 
		plt.xlabel('time [ns]')
		plt.ylabel('normalised intensity')
		plt.legend(loc = 'upper left')
		plt.savefig('phasespace_'+str(index_time))
		plt.tight_layout()
		
		
		plt.show()
		#plt.pause(1.0)
		#plt.cla()











