from __future__ import division
import pylab as plt
import numpy as np
import sys
import os
import matplotlib.cm as cm
import math 
import matplotlib.animation as animation
from subprocess import call
from scipy.interpolate import griddata
import tomo_defs_20191202 as tdefs


tomo_time = 0.00035 #start time data used in tomography (this scope time used in the data files, not w.r.t. injection)
#the injection time will be calculated by the script	

flattop_time = None #time at which energy stops increasing (Note: phi_s_deg is set independently of this)
animate = True #animate tomography results

nturns_sort = 200#Number of turns of data to read and sort (a subset is used in tomography reconstruction)
filter_data = False #remove high frequency noise from data before reconstruction

show_imshow = True #show turn-by-turn data using pylab.imshow
stop_after_imshow = False #True to show bunch monitor heat map (imshow). False to continue with tomography after imshow

time_offset_ms = None#If None the offset is calculated by finding the time at which the frequency of the RF waveform is equal to the ideal frequency at 11 MeV.
plot_time_offset = True #show how the injection time offset is calculated.

turn_recon_l = []

#specify the start and end turn and increment at which the profile is reconstructed
specify_absolute_turns = True #if False, specify turns in terms of fraction of a sync. period
if specify_absolute_turns:
	recon_start = 1 #recconstruct phase space at this turn first
	turns_tomo = 176#243#243#243#243
	recon_step = 50
	recon_stop = 176#120#turns_tomo
	
	turn = recon_start
	while turn < recon_stop:
		turn_recon_l.append(turn)
		turn = turn + recon_step
	turn_recon_l.append(turns_tomo)

else:
    turns_syncfrac = 0.5 #fraction of synchrotron oscillation to include
    recon_step_frac = 0.1 #0.1 #reconstruction steps in terms of fraction of sychrotron period
	

#set principal RF parameters
phi_s_deg = 20 #set synchronous phase
V_rf = 4e3#2.69e3#4.0e3 # in [V]#2.7e3, 20/9 
harmonic = 1
phi_s = phi_s_deg*math.pi/180 #synchronous phase

#second harmonic
V_rf2 = 0 #set to zero so no additional harmonic is added
Vratio = 0#V_rf2/V_rf
hratio = 3 #ratio of 2nd to 1st RF system
phi12 = 0.5*math.pi #phase between RF systems (in radians of the principal harmonic)

if flattop_time !=None:
	if tomo_time >= flattop_time and phi_s_deg !=0.0:
		print "WARNING: phi_s is ",phi_s_deg, "deg but tomo_time is later than flattop_time"


H_fn = lambda phi, dp, phi_s, a1, a2: a1*(dp**2) + a2*(np.cos(phi) - np.cos(phi_s) + (phi - phi_s)*np.sin(phi_s))
H_fn_gen_pot = lambda n, phi, phi_s, coef, ang_off: coef*(np.cos(n*phi + ang_off) - np.cos(n*phi_s + ang_off) + (phi - phi_s)*np.sin(n*phi_s + ang_off))

#March&April 2019
#dirname = "../data/2019/3/27/"
#dirname = "../data/2019/6/13/"
dirname = "../data/2019/6/20190620/"
multi_channel = True

col_bm_rf = [0,2]#first file on 20/6/2019
#col_bm_rf = [0,1]#other files on 20/6/2019


#time_offset_ms = 0.1143 #offset time in March&April 2019 was 0.1143


#Uesugi-san table containing constant phis=20 degree programme 
svkfilepath = "svk20.dat"

show_data_switch = False #show raw bunch monitor and rf data
show_mountain_range = False #show turn-by-turn data in mounntain range format


read_rf_waveform = True #read RF data

interpolate_data = True #if False interpolate data onto time axis which is determined by rf waveform

idelay_corr_phi_s = False #if True, correct above idelay for nonzero phi_s

#jump point 0.00221497s

nscan_param = 1 #number of synch point settings to check on either side of idelay setting

read_file_offset = False #no longer supported

#calculate offset in time w.r.t RF lookup table. 
calc_rf_time_offset = True #special run to calculate offset of timebase from RF table time


intensity_calc = True
write_intensity_file = False
write_offset_file = False

PROTON_MASS = 938.27231e6 
SPEED_OF_LIGHT = 299792458
mass = PROTON_MASS


#set lattice parameters
kindex = 7.6
alpha_c = 1/(kindex+1)
gamma_tr = (1/alpha_c)**0.5
rho=0.489396551697 #bending radius [m]
	
print "phi_s set to ",phi_s_deg, " degrees"

#set tomography output directory
tomo_outdir = "tomo_output"
files = os.listdir(os.curdir)
if tomo_outdir not in files:
	print "create output directory "+tomo_outdir
	os.mkdir(tomo_outdir)

print "reconstruct at following time [ms] ", 1e3*tomo_time


interactive_data_sel = True
if interactive_data_sel:
	#if "/data/2019/3" in dirname or "/data/2019/4" in dirname:
		 
	if "/data/2019/6/13" in dirname:
		print "June 2019 data"
		fname_sig, tdat_rf, sig_rf, time_data, signal_data, t_phis_file, phis_file, offset_file_ms = tdefs.read_data_files(dirname, read_rf_waveform = read_rf_waveform, read_file_offset = read_file_offset, multi_channel = multi_channel)

	else:
		fname_sig, time_data, signal_data, tdat_rf, sig_rf, t_phis_file, phis_file, offset_file_ms = tdefs.read_data_files(dirname, col_bm_rf, read_rf_waveform = read_rf_waveform, read_file_offset = read_file_offset, multi_channel = multi_channel)
else:
	if multi_channel:
		time_data, waveform_data = tdefs.read_scope(path, multi_channel = multi_channel)

		waveform_data = np.array(waveform_data)
		waveform_data_tr = np.transpose(waveform_data)
	
		tdat_rf = time_data
		sig_rf = waveform_data_tr[0]
		signal_data = waveform_data_tr[1]

#read svk file
time_svk, ke_svk, rf_freq_svk = tdefs.read_svk_file(svkfilepath)
svk_data = [time_svk, ke_svk, rf_freq_svk]

time_interval = time_data[2] - time_data[1]
time_interval_ns = 1e9*time_interval
time_data_ms = 1e3*np.array(time_data)

#idelay is delay in time bins between the real RF phase seen by the beam and that seen on the scope
if round(time_interval_ns) == 2:
	idelay = -47
else:
	idelay = int(round(-47*(2/time_interval_ns)))

print "time interval_ns ",time_interval_ns
print "idelay ",idelay
print "time data range ",time_data[0], time_data[-1]

#find time indices at tomo_times

i = 0
for t in time_data:
	if t > tomo_time:
		i_data = i
		break
	i = i + 1
istart_guess = i_data


#If injection time is not supplied, calculate it here using the RF frequency
if time_offset_ms == None:
	nturns_rf = 1000
	tof_mean, tof_l, time_indices, ttime_a = tdefs.calc_tof_rf_waveform(tdat_rf, sig_rf, nturns_rf+1, -0.05e-3)

	if len(ttime_a) < nturns_sort+1:					
		print "nturns limited by data to ",len(ttime_a)	
	if calc_rf_time_offset:
		time_offset_ms = tdefs.time_offset_calc(dirname, ttime_a, t_phis_file, svkfilepath, svk_data, fname_sig, write_offset_file, plot_time_offset)

	print "Calculated injection time (time_offset) [ms] ",time_offset_ms


time_offset = 1e-3*time_offset_ms
	
tof_mean, tof_l, time_indices, ttime_a = tdefs.calc_tof_rf_waveform(tdat_rf, sig_rf, nturns_sort+1, tomo_time)
istart = time_indices[0]

#calculate kinetic energy, beta, synchrotron period at each selected time

print "time_offset, tomo_time ",time_offset, tomo_time

if flattop_time != None:
	ke_flattop = tdefs.ke_lookup(flattop_time, svkfilepath, time_offset)

	if tomo_time < flattop_time:
		ke = tdefs.ke_lookup(tomo_time, svkfilepath, time_offset)
	else:
		ke = ke_flattop
else: 
	ke = tdefs.ke_lookup(tomo_time, svkfilepath, time_offset)



prepare_for_tomo = True			
if prepare_for_tomo:
	
	#calc harmonics of RF waveform
	calc_harmonics = False
	if calc_harmonics:
		tdefs.calc_waveform_harmonics(tdat_rf, sig_rf, ttime_a)
	
	#calculate points per turn based on average tof
	points_per_turn = int(round(tof_mean/time_interval)) 
	nslices = int(points_per_turn)
	
	#calculate mean tof
	t_turn_data = np.linspace(0, tof_mean, nslices)

	T_rf = tof_mean/harmonic
	T_rf_ns = 1e9*T_rf

	#work out longitudinal parameters of interest
	Ts, dtbin_phase,omega_rf, Rnom, Bnom, phase_sep_deg, dp_sep, bh, sync_bin, synch_index, binslz, binsuz, tadj, Hsep, phi_l1, phi_l2 = tdefs.longitudinal_param(V_rf, phi_s, harmonic, alpha_c, rho, ke, tof_mean, nslices, nturns_sort, mass)

	Etot_turn, dEbin_c_turn, p0_turn, eta_turn, beta_turn = tdefs.longitudinal_turn(V_rf, phi_s, harmonic, alpha_c, ke, Rnom, nturns_sort, dtbin_phase, mass)
	
	print "ke, Ts ",ke, Ts

	p0_mev_turn = 1e-6*p0_turn

	#calculate delay in rf signal
	if idelay == None:
		print "calculate delay in signal"
		prof_data_a = tdefs.sort_data(istart, time_data, signal_data, ttime_a, time_indices, nturns_sort, nslices, interpolate_data)
		#sort data into turns covering roughly a synchrotron period
		nt_Ts = int(round(Ts))
		idelay = tdefs.calc_sync_point(prof_data_a, nt_Ts, sync_bin, show_sym_plots = True)

	#introduce further delay for non-zero phi_s (assuming idelay applies to a stationary bucket)
	if idelay != None and idelay_corr_phi_s:
		#calculate further shift in delay to account for phi_s !=0
		idelay_phi_s = int(nslices*phi_s_deg/360)

		idelay = idelay - idelay_phi_s

		print "extra delay for non-zero phi_s [bins] ",idelay_phi_s
		print idelay

		#sys.exit()

	#sort data into turns, interpolate onto fixed time axis (phase)
	istart = istart + idelay

	if filter_data:
		signal_data = tdefs.filter_data(signal_data, time_interval_ns*1e-9)

	prof_data_a = tdefs.sort_data(istart, time_data, signal_data, ttime_a, time_indices, nturns_sort, nslices, interpolate_data, nshift= idelay) 
	
	phase_ax = np.linspace(-180,180,nslices) + 1.1544309246398268

	
	#find symmetry - look for bin where integral on each side is most equal
	sym_l = []
	for data in prof_data_a:
		sym, chisq = tdefs.find_symmetry_point_integral(data, plot_data = False)
		sym_l.append(sym)
	sym_a = np.array(sym_l)

	sym_phase_a = phase_ax[sym_a]


	if intensity_calc:
		intensity, phase_fmax, std_prof_max, fwhm = tdefs.intensity_calc_fn(prof_data_a, phase_ax, ttime_a, tof_l, time_interval_ns, dirname, fname_sig, sym_phase_a, filter_data, write_intensity_file = write_intensity_file)	
	
	if show_mountain_range:
		bcktlim1 = binslz 
		bcktlim2 = nslices - binsuz
		xaxis_deg = np.linspace(-180,180,nslices)
		tdefs.mountain_range(prof_data_a, nturns_sort, nslices, incl_xaxis = True, xaxis= xaxis_deg, xlabel='phase [deg]')			

	midpt_sym_phase = 0.5*(np.max(sym_phase_a) + np.min(sym_phase_a))
	midpt_phase_fmax = 0.5*(np.max(phase_fmax) + np.min(phase_fmax))
	print "calculate extrema in oscillation in terms of sym_phase and phase_fmax"
	print "*******************************************"
	print "sym_phase (black line) extrema"
	print "minimum: turn, time [ms] phase [deg] ",np.argmin(sym_phase_a)+1, 1e3*ttime_a[np.argmin(sym_phase_a)], np.min(sym_phase_a)
	print "maximum: turn, time [ms] phase [deg] ",np.argmax(sym_phase_a)+1, 1e3*ttime_a[np.argmax(sym_phase_a)], np.max(sym_phase_a)
	print "midpt sym_phase ",midpt_sym_phase
	print "suggested phase jump from minimum ",midpt_sym_phase-np.min(sym_phase_a), " at ", 1e3*ttime_a[np.argmin(sym_phase_a)], "ms"
	print "suggested phase jump from maximum ",midpt_sym_phase-np.max(sym_phase_a), " at ", 1e3*ttime_a[np.argmax(sym_phase_a)], "ms"

	print "*******************************************"
	print "phase_fmax (magenta line) extrema"
	print "minimum: turn, time [ms] phase [deg]",np.argmin(phase_fmax)+1, 1e3*ttime_a[np.argmin(phase_fmax)], np.min(phase_fmax)
	print "maximum: turn, time [ms] phase [deg] ",np.argmax(phase_fmax)+1, 1e3*ttime_a[np.argmax(phase_fmax)], np.max(phase_fmax)
	print "midpt phase_fmax ",midpt_phase_fmax

	print "suggested phase jump from minimum ",midpt_phase_fmax-np.min(phase_fmax), " at ", 1e3*ttime_a[np.argmin(phase_fmax)], "ms"
	print "suggested phase jump from minimum ",midpt_phase_fmax-np.max(phase_fmax), " at ", 1e3*ttime_a[np.argmax(phase_fmax)], "ms"

	if show_imshow:

		tlim1 = 1e3*(ttime_a[0])	
		tlim2 = 1e3*(ttime_a[-1])

		tomo_lim = [recon_start, turns_tomo]
		extent = [-180, 180, tlim1, tlim2]

		tdefs.imshow_plot(prof_data_a,extent,ylabel='time (ms)',xlabel='waveform phase [deg]', hline_l =turn_recon_l, hline_l_2 =tomo_lim, extra_data = [sym_phase_a, phase_fmax], interactive = animate)

	if stop_after_imshow:
		print "stop_after_imshow ",stop_after_imshow
		sys.exit()


	t_turn_data_ns = 1e9*t_turn_data + tadj

	chisq_sum_l = []

	if nscan_param == 1:
		idelay_a = np.array([idelay])
	else:

		if phi_s_deg == 0:
			idelay_a = np.arange(-30, -18)
		else:
			idelay_a = np.arange(-35, -18)


	#number of turns specified in terms of sync period
	if specify_absolute_turns:
		nturns = turns_tomo
	else:
		nturns = int(turns_syncfrac*Ts)

	#this is not a real Bdot - it is a parameter used to specify phis in the tomo code
	Bdot = (V_rf*np.sin(phi_s))/(rho*2*math.pi*Rnom)
	time_str = str(1e3*tomo_time)
	for iscan in range(nscan_param):


		if nscan_param > 1:
			idelay1 = idelay_a[iscan]
			istart = istart + idelay1
			prof_data_a, t_turn_data = sort_data(istart, ttime_a, time_indices, nturns_sort, nslices, interpolate_data, nshift= idelay1)
		

		prof_row_sum =  np.sum(prof_data_a, axis=1)
		prof_data_norm_a =  np.transpose(np.transpose(prof_data_a)/prof_row_sum)

		print "write ", nturns, " turns to file"
		print "first, last time in ttime_a [ms] ",1e3*ttime_a[0] - time_offset_ms, 1e3*ttime_a[nturns] -time_offset_ms
		#sys.exit()

		

		intensity_totomo = []
		for it in range(nturns):
			int1 = time_interval_ns*np.trapz(prof_data_a[it])
			intensity_totomo.append(int1)
			#plt.plot(prof_data_a[it]+0.01*it)
		#plt.show()
	
		f1 = open('rawdata.txt', 'w')
		for it in range(recon_start, recon_start+nturns):
			yprof = prof_data_a[it]
			yprof_norm = yprof/max(yprof)
			#plt.plot(yprof,'ko-')
			for y in yprof:
				print >>f1,y
		f1.close()

		if not specify_absolute_turns:
			recon_step = int(recon_step_frac*Ts)


		#***RUN TOMOGRAPHY RECONSTRUCTION***
		tomo_code_param = [nturns, nslices, synch_index, dtbin_phase, binslz, binsuz, recon_start, recon_stop, recon_step, Bnom, Bdot, V_rf, V_rf2, harmonic, hratio, phi12, Rnom, rho, gamma_tr]
		tdefs.run_tomo_code(tomo_outdir, tomo_code_param)
	
		show_tomo_result = True
		if show_tomo_result:

			#read output parameters from plotinfo.data
			plotinfo = tdefs.get_plotinfo(input_file=tomo_outdir+"/plotinfo.data")
		
			x0 = plotinfo['xat0'][0]
			y0 = plotinfo['yat0'][0]
			dEbin =  plotinfo['dEbin'][0]
			dtbin =  plotinfo['dtbin'][0]
			phi_s_out = plotinfo['phi_s'][0] #synch phase


			#dEbin_calc = dtbin*harmonic*omega_
			fname_prof = 'profiles.data'
			profiles_norm_raw = tdefs.file_read(tomo_outdir+'/'+fname_prof)
			profiles_norm_data_a = np.split(profiles_norm_raw, nturns)

			#fname_prof = 'profiles_recon.data'
			#profiles_norm_recon_raw = tdefs.file_read(tomo_outdir+'/'+fname_prof)
			#profiles_norm_recon_a = np.split(profiles_norm_recon_raw, nturns)

			

			show_profile_out = False
			if show_profile_out:

				fname_raw = 'rawdata.txt'

				bunchdata= tdefs.file_read(fname_raw)
				raw_spl = np.split(bunchdata, nturns)

				
				tdefs.mountain_range(profiles_norm_data_a, nturns, nslices, fmount='mountainrange_data')


				tdefs.mountain_range(profiles_norm_recon_a, nturns, nslices, fmount='mountainrange_recon')


				max_prof_index_data = np.argmax(profiles_norm_data_a, axis=1)
				max_prof_index_recon = np.argmax(profiles_norm_recon_a, axis=1)
				phase_profmax_data = [phase_ax[m] for m in max_prof_index_data]
				phase_profmax_recon = [phase_ax[m] for m in max_prof_index_recon]
				mean_phase_data_norm = [sum([p1*f for p1,f in zip(phase_ax, prof)])/sum(prof) for prof in profiles_norm_data_a]

				plt.plot(phase_profmax_data,'r-')
				plt.plot(phase_profmax_recon,'r--')
				plt.plot(mean_phase_data_norm, 'k-')
				plt.show()
				

				#plt.plot(profiles_a)
				#for raw, prof, prof_recon in zip(raw_spl, profiles_norm_data_a, profiles_norm_recon_a):
					#plt.plot(raw,'k')
					#plt.subplot(211)
				#	plt.plot(prof,'r')
					#plt.subplot(212)
				#	plt.plot(prof_recon,'b')
				#plt.show()


		
			inputp = tdefs.get_input()
			E0, gamma0, eta0, beta0, omega_rev0, Vpeak = tdefs.set_param(inputp)
			

			

		T = 2*np.pi/omega_rev0 #rf period
		max_dp_sep = max(dp_sep)


		#len_img_data = len(image_data_a)
		ndim = nslices - (binslz + binsuz)  #int(len_img_data**0.5)

		t_synch = 1e9*phi_s_out*T/(2*math.pi)

		t_tomo_ns = 1e9*dtbin*(np.array(range(ndim)) - x0) + t_synch

		phase_tomo = t_tomo_ns*(2*math.pi)/T_rf_ns
		phase_tomo_deg = phase_tomo*180/math.pi

		print "x0 from tomography code ",x0

		
		print "t_synch ",t_synch, " phase at t_synch ",t_synch*(2*math.pi)/T_rf_ns

		
		t_shift_ns = -1e9*dtbin*(x0 - 0.5*ndim) #+ t_synch
		#print "time shift applied ",t_shift_ns
		
		t_data_step = t_turn_data[1] - t_turn_data[0]
		#t_tomo_ns = np.linspace(t_sep_ns[0], t_sep_ns[-1], nslices)
		
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

		dp_a_recon = []
		dp_sep_recon = []
		discrep_total = 0
		for irc in recon_id:

			#if run_pyheadtail:

			#counts, xbins, ybins,_ = plt.hist2d(bunch_t_ini_ns, bunch_dp_ini, bins=60, norm=LogNorm())

			print "recon at turn ",irc
			#process main tomography output file
			#fname_img = 'image001.data'


			fname_discrep = tomo_outdir+'/d000'[:-len(str(irc))]+ str(irc) + '.data'
			discrep_data = tdefs.file_read(fname_discrep, sel_col=1)
			discrep_final = discrep_data[-1]
			discrep_total = discrep_total + discrep_final
			print "discrepancy final, total ",discrep_final, discrep_total  

			fname_img = tomo_outdir+'/image000'[:-len(str(irc))]+ str(irc) + '.data'
			print "fname_img ",fname_img
		
			image_data_a = tdefs.file_read(fname_img) #1D data
			
			#split into array
			image_data_spl = np.split(image_data_a, ndim)
			#transpose to get phase on horizontal axis, energy on vertical
			image_data_tr = np.transpose(image_data_spl)

			#construct dp pixel array taking synchronous point x0, y0 into account.
			delta_E_a = dEbin_c_turn[irc-1]*(np.array(range(ndim)) - y0)
			E_tot_eV_a = Etot_turn[irc-1] + delta_E_a
			p_a = [1e-6*(Et**2 - mass**2)**0.5 for Et in E_tot_eV_a]
			dp_a = (p_a - p0_mev_turn[irc-1])/p0_mev_turn[irc-1]

			#calculate separatrix at each energy

			npts = 1000
			a1 = 0.5*harmonic*omega_rf*eta_turn[irc-1]
			a2 = omega_rf*V_rf/(2*math.pi*(beta_turn[irc-1]**2)*Etot_turn[irc-1])
			dp_sep1, phase_sep1,_ = tdefs.separatrix(phi_s, eta_turn[irc-1], harmonic, npts, a1, a2)
			dp_sep1 = np.insert(dp_sep1,0,0)
			phase_sep1 = np.insert(phase_sep1, 0, phi_l2) 

			phase_mesh, dp_mesh = np.meshgrid(phase_tomo, dp_a)
			H_mesh = -H_fn(phase_mesh, dp_mesh, phi_s, a1, a2)
			#varphi_mesh = varphi_fn(phase_mesh, dp_mesh, phi_s)	

			show_H_mesh = False
			if show_H_mesh:
				plt.contour(phase_tomo_deg, dp_a, H_mesh)
				plt.plot(phase_sep_deg, dp_sep,'r') #bucket
				plt.plot(phase_sep_deg, -dp_sep,'r') #bucket
				plt.show()
				sys.exit()
	
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
				image_data_fil  = tdefs.filter_isolated_cells(image_data_tr, struct=np.array([[0,0,0],[1,1,1],[0,0,0]]))
				image_data_fil  = tdefs.filter_isolated_cells(image_data_fil, struct=np.array([[0,1,0],[0,1,0],[0,1,0]]))


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
			prof_diff = profiles_norm_data_a[irc-1] - tomo_time_proj
			chisq = sum(prof_diff**2)

			#print "chisq ",chisq
			chisq_l.append(chisq)
		
			if iloop == 0:
				tomo_time_proj_a = np.array([tomo_time_proj])
				image_Hgrid_a = image_Hgrid
				dp_a_recon = np.array([list(dp_a)])
				delta_E_a_recon = np.array([list(delta_E_a)])
				dp_sep_recon = np.array([list(dp_sep1)])
			else:
				tomo_time_proj_a = np.vstack((tomo_time_proj_a, tomo_time_proj))
				image_Hgrid_a = np.vstack((image_Hgrid_a, image_Hgrid))
				dp_a_recon = np.vstack((dp_a_recon, dp_a))
				delta_E_a_recon = np.vstack((delta_E_a_recon, delta_E_a))
				dp_sep_recon = np.vstack((dp_sep_recon, dp_sep1))

			iloop = iloop + 1
			

		if nscan_param > 1:
			chisq_sum = sum(chisq_l)
			chisq_sum_l.append(chisq_sum)

			print "idelay, chisq summed over turns ",idelay1, chisq_sum
		else:
			print "tomography vs input data chisq ",chisq
	

	if nscan_param > 1:
		plt.plot(idelay_a, chisq_sum_l,'ko-')
		plt.xlabel('time point shift')
		plt.ylabel('chisq sum')
		plt.ylim(ymin=0)
		plt.axvline(x=idelay, color='gray')
		plt.savefig('timeshiftscan')
		plt.show()


compare_image_Hgrid = False
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
#cont_limit = 0.7*max(image_max_l)

cont_limit = 5*0.0235
print "!!cont_limit fixed at",cont_limit

#select horizontal axis
versus_time = False
versus_phim = False #phim is phase - phis
versus_dp = True #if False vertical axis is delta_E

if versus_time:
	xaxis = t_tomo_ns
elif versus_phim:
	xaxis = phase_tomo_deg - phi_s_deg
	phi_off = phi_s_deg
else:
	xaxis = phase_tomo_deg
	phi_off = 0


if animate:
	#Display tomography reconstruction results in animated form

	animate_incl_subplot = True #include comparison of line integral with bunch monitor data if True

	#plt.subplot(211)
	levels = np.linspace(0, cont_limit, 30)
	#plt.title('reconstructed turn '+sttomo_time_proj_lr(irc)+ ' time '+time_str+' ms')
	
	if animate_incl_subplot:
		fig, (ax1, ax2) = plt.subplots(2,1)
	else:
		fig, ax1 = plt.subplots()

	def animate(ifrm):
		"""Function called by FuncAnimation
		ifrm: Frame index"""
	
		irc = recon_id[ifrm] #turn index

		ax1.clear()	
		if animate_incl_subplot:
			ax2.clear()

		if image_scaled:
			cont = ax1.contour(xaxis, dp_a_recon[ifrm], intensity_scale[ifrm]*image_all_a[ifrm], levels = levels)
		else:
			cont = ax1.contour(xaxis, dp_a_recon[ifrm], image_all_a[ifrm], levels = levels)
		
		npts = 1000
		dp_sep1, phase_sep1, t_sep_ns = tdefs.separatrix(phi_s, eta_turn[irc-1], harmonic, npts, a1, a2, T=T)

		#ax1.plot(t_sep_ns, dp_sep1,'r') #bucket
		#ax1.plot(t_sep_ns, -dp_sep1,'r') 

		ax1.plot(phase_sep_deg - phi_off, dp_sep_recon[ifrm],'r') #bucket
		ax1.plot(phase_sep_deg - phi_off, -dp_sep_recon[ifrm],'r') #bucket
		
		
		ax1.set_ylabel(r"$\Delta$p/p")
		ax1.set_xlim(phase_sep_deg[0], phase_sep_deg[-1])
		ax1.set_ylim(-1.01*bh, 1.01*bh)

		ax1.set_title('turn '+str(irc)+' (press SPACE to pause/play)')

		if animate_incl_subplot:
			ax2.set_xlabel(r'$\phi$ [deg]')
		else:
			ax1.set_xlabel(r'$\phi$ [deg]')			

		if animate_incl_subplot:
			ax2.plot(xaxis, profiles_norm_data_a[irc-1],'k-',linewidth=2,label='data')		
			ax2.plot(xaxis, tomo_time_proj_a[ifrm],'r-',linewidth=2,label='tomography')
			ax2.legend(loc = 'upper right')	
			
		#cont.set_zdata(image_all_a[1])
		return cont


	def on_press(event):
		"""Function that handles events. Pause animiation if space bar is pressed"""

		if event.key.isspace():
		    if anim.running:
		        anim.event_source.stop()
		    else:
		        anim.event_source.start()
		    anim.running ^= True
		elif event.key == 'left':
		    anim.direction = -1
		elif event.key == 'right':
		    anim.direction = +1

		# Manually update the plot
		if event.key in ['left','right']:
		    t = anim.frame_seq.next()
		    animate(t)
		    plt.draw()

	Nt = len(image_all_a)
	#Nt = recon_id

	fig.canvas.mpl_connect('key_press_event', on_press)
	anim = animation.FuncAnimation(fig, animate, frames=Nt, interval=500, blit=False)

	#save gif
	anim.save('tomo.gif', dpi=80, writer='imagemagick')

	anim.running = True
	anim.direction = +1

	plt.ioff()
	plt.show()

else:
	#Display tomography reconstruction results as a series of surface and/or contour plots
	plot_surface = False
	plot_contour = True
		
	dp_bin_ini = np.linspace(-bh, bh, nslices)
	ph_a = phase_tomo
	ph_mesh, dp_mesh = np.meshgrid(ph_a, dp_bin_ini, sparse = True)
	
	#Hgrid = H_fn(ph_mesh, dp_mesh, phi_s, a1, a2)
	Hgrid = H_fn(ph_mesh, dp_mesh, phi_s, a1, a2) + H_fn_gen_pot(hratio*harmonic, ph_mesh, phi_s, Vratio*a2, phi12)

	if plot_surface:
		from mpl_toolkits.mplot3d import Axes3D
		t_mesh, dp_mesh = np.meshgrid(t_tomo_ns, dp_a)

		ipl = 0
		for irc in recon_id:
			if image_scaled:
				image_surf = image_all_a[ipl]*intensity_totomo[irc-1]
			else:
				image_surf = image_all_a[ipl]

			fig = plt.figure()
			ax = fig.gca(projection='3d')			
			plt.title('reconstructed turn '+str(irc)+ ' time '+time_str+' ms')
			ax.plot_surface(t_mesh, dp_mesh, image_surf, cmap="seismic") #or coolwarm

			ax.set_xlabel('time [ns]')
			ax.set_ylabel(r"$\Delta$p/p")
			
			
			ax.set_zlim(0, cont_limit)
			ipl = ipl + 1	

			plt.show()

	if plot_contour:
		ipl = 0
		for irc in recon_id:

			if image_scaled:
				image = image_all_a[ipl]*intensity_totomo[irc-1]
			else:
				image = image_all_a[ipl]
			
			
			#levels = np.linspace(0, 0.008, 20)
			levels = np.linspace(0, cont_limit, 30)
			print "plot phasespace at turn & time ",irc, ttime_a[irc-1]
			#print "max dp/p along separatrix ",max(dp_sep_recon[ipl])

			show_H_contours = False

			plt.subplot(211)
			#plt.title('reconstructed turn '+str(irc)+' time [ms] '+str(1e3*ttime_a[irc-1]))
			#plt.title('reconstructed turn '+str(irc))
			if plot_contour:

				if show_H_contours:
					plt.contourf(xaxis, dp_bin_ini, Hgrid, 20,alpha=0.7,cmap = 'RdGy')

				if versus_dp:
					ctom = plt.contour(xaxis, dp_a_recon[ipl], image, levels = levels)
				else:
					ctom = plt.contour(xaxis, delta_E_a_recon[ipl], image, levels = levels)
				#ctom = plt.contour(t_tomo_ns, dp_a, image)
				#plt.colorbar(ctom)
			
			show_fund_separatrix = True
			if versus_time:
				plt.plot(t_sep_ns, dp_sep,'r') #bucket
				plt.plot(t_sep_ns, -dp_sep,'r') #bucket
				plt.axvline(x=t_tomo_ns[int(x0)]) #synch time
				plt.xlim(1.01*t_sep_ns[0],1.01*t_sep_ns[-1])
				plt.xlabel('time [ns]')
			else:
				
				if show_fund_separatrix:
					plt.plot(phase_sep_deg - phi_off, dp_sep_recon[ipl],'r') #bucket
					plt.plot(phase_sep_deg - phi_off, -dp_sep_recon[ipl],'r') #bucket

				#plt.axvline(x=phi_s_deg - phi_off) #synch time
				plt.xlim(phase_sep_deg[0], phase_sep_deg[-1])
				#plt.xlabel(r'$\phi$ [deg]')
			
			if versus_dp:
				plt.ylim(-1.3*bh, 1.3*bh)	
			else:
				pass
				#plt.ylim(1.01*bh, -1.01*p0_turn[irc-1]*bh/beta_turn[irc-1])

			x1,x2,y1,y2 = plt.axis()
			
			plt.ylabel(r"$\Delta$p/p")

			plt.subplot(212)

			plt.plot(xaxis, profiles_norm_data_a[irc-1],'k-',linewidth=2,label='data')		
			plt.plot(xaxis, tomo_time_proj_a[ipl],'r-',linewidth=2,label='tomography')
			
			if versus_time:
				if binslz > 0:
					plt.axvline(x=t_turn_data_ns[binslz] + 0.1*t_data_step, color='gray', linestyle='--')
				if binsuz > 0:
					plt.axvline(x=t_turn_data_ns[nslices - binsuz]-0.1*t_data_step, color='gray', linestyle='--')		
				plt.axvline(x=t_synch, color='gray', linestyle = '--')
				plt.xlim(1.01*t_sep_ns[0],1.01*t_sep_ns[-1])
				plt.xlabel('time [ns]')		

			else:
				#plt.axvline(x=phi_s_deg - phi_off, color='gray', linestyle = '--')
				plt.xlabel(r'$\phi$ [deg]')
				plt.xlim(phase_sep_deg[0], phase_sep_deg[-1])		

			
			plt.ylabel('normalised signal')
			#plt.legend(loc = 'upper left')
			plt.savefig('phasespace_'+str(irc))
			plt.tight_layout()

			plt.show()

			ipl = ipl + 1








