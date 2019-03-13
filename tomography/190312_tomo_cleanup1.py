from __future__ import division
import pylab as plt
import numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import quadrature
import sys
import os
import matplotlib.cm as cm
import tomo_defs_190311 as tdefs
import math
from subprocess import call
import time 
import matplotlib.animation as animation
from subprocess import call

from scipy.interpolate import griddata

H_fn = lambda phi, dp, phi_s, a1, a2: a1*(dp**2) + a2*(np.cos(phi) - np.cos(phi_s) + (phi - phi_s)*np.sin(phi_s))

phi_s_deg = 0 #set synchronous phase
phi_s = phi_s_deg*math.pi/180

#set data directory
#dirname = "../data/2018/9/20/proc/"
#dirname = "../data/2018/9/20/machida-san_pattern/proc/"
dirname = "../data/2018/9/21/machida_san_foil_escape/proc/"
#dirname = "../data/2018/9/25/machida-san_2018_09_25/proc/"
#dirname = "../data/2018/9/27/JB_2018_09_27/proc/"
#dirname = "../data/2018/9/21/machida_san_foil_escape/" # "../data/2018/9/20/"

#Uesugi-san table containing constant phis=20 degree programme 
svkfilepath = "../sep18_exp/svk20.dat"

show_data_switch = False #show raw bunch monitor and rf data
show_mountain_range = False #show turn-by-turn data in mounntain range format
show_imshow = False #show turn-by-turn data using pylab.imshow

read_rf_waveform = True #read RF data
read_phis_file = True #read file containnig phis timings

npfac = 1.0 #interpolate onto time axis with npfac times the number of points as the data
interpolate_data = True #if False interpolate data onto time axis which is determined by rf waveform
animate = False #animate phase space results

#idelay - offset time by this to get synchronous phase correct
idelay =  None#-24#-27#-24 #-13 #-23 for flattop 21/9	 #-38 if phi_s is 33 degrees?
idelay_corr_phi_s = False #if True, correct above idelay for nonzero phi_s


nscan_param = 1 #number of synch point settings to check on either side of idelay setting

#calculate offset in time w.r.t RF lookup table. 
calc_rf_time_offset = False #special run to calculate offset of timebase from RF table time
if calc_rf_time_offset:
	tomo_times = np.array([-0.05e-3])
	nturns_sort = 10000 #10000 for full intensity file

	write_offset_file = False
	plot_time_offset = True
	
	intensity_calc = True
	write_intensity_file = False

	add_time_offset = False
	read_file_offset = False
else:
	tomo_times = np.array([0.3e-3]) #this is with respect to offset if add_time_offset is True
	read_file_offset = True
	nturns_sort = 1000

	intensity_calc = False

	write_offset_file = False
	write_intensity_file = False
	add_time_offset = True

PROTON_MASS = 938.27231e6 #938.27203e6
SPEED_OF_LIGHT = 299792458
mass = PROTON_MASS

specify_absolute_turns = True #if False, specify turns in terms of fraction of a sync. period
if specify_absolute_turns:
	turns_absolute = 100
	recon_step = 50
	recon_stop = turns_absolute
else:
	turns_syncfrac = 0.5 #fraction of synchrotron oscillation to include
	recon_step_frac = 0.1 #0.1 #reconstruction steps in terms of fraction of sychrotron period

recon_start = 1 #reconstruct phase space at this turn first

#set parameters
kindex = 7.6
alpha_c = 1/(kindex+1)
gamma_tr = (1/alpha_c)**0.5
V_rf = 4.0e3 # in [V]#2.7e3, 20/9 
harmonic = 1
rho=0.489396551697 #bending radius [m]


#if V_rf != 4e3:
#	phi_s = np.arcsin((4e3/V_rf)*np.sin(phi_s))
#	phi_s_deg = phi_s*180/math.pi
#	print "modified phi_s (rad & deg) ",phi_s, phi_s_deg
	


if idelay !=None:
	print "applying the following shift to rf waveform to find synchronous time ",idelay

print "phi_s set to ",phi_s_deg, " degrees"

tomo_outdir = "tomo_output"
files = os.listdir(os.curdir)
if tomo_outdir not in files:
	print "create output directory "+tomo_outdir
	os.mkdir(tomo_outdir)

print "reconstruct at following times [ms] ", 1e3*tomo_times
print "phi_s set to [deg] ",phi_s_deg
if phi_s_deg == 0:
	flattop = True
	flattop_time = 0.3e-3#9.5e-4 #0.3e-3 #estimate
	print "flattop time is [ms] ", 1e3*flattop_time
else:
	flattop = False


#read the data
interactive_data_sel = True
if interactive_data_sel:
	fname_sig, time_data, signal_data, tdat_rf, sig_rf, t_phis_file, phis_file, offset_file_ms = tdefs.read_data_files(dirname, read_rf_waveform = read_rf_waveform, read_file_offset = read_file_offset)

#read svk file
time_svk, ke_svk, rf_freq_svk = tdefs.read_svk_file(svkfilepath)
svk_data = [time_svk, ke_svk, rf_freq_svk]

time_interval = time_data[2] - time_data[1]
time_interval_ns = 1e9*time_interval
time_data_ms = 1e3*np.array(time_data)

if npfac != 1:
	print "sorry, npfac != 1 option isn't supported"
	sys.exit()

			
def read_data_file_2D(filen):
	
	f1 = open(filen,'r')
	data_l = []
	for data in f1:
		dataspl = data.split()
		data_l.append(float(dataspl[1]))
	data_a = np.array(data_l)

	return data_a


if add_time_offset:
	toffset_rf = 1e-3*offset_file_ms
else:
	toffset_rf = 0

#find time indices at tomo_times
istart_guess_l = []
for tt in tomo_times:
	if add_time_offset:
		tt = tt + toffset_rf
	i = 0
	icrude = None
	for t in time_data:
		if t > tt:
			i_data = i
			break
		i = i + 1
	istart_guess_l.append(i_data)


if not read_rf_waveform:
	tof_guess = 600e-9
	points_per_turn = int(round(tof_guess/time_interval)) 
	#calculate TOF
	tof, time_indices, ttime_a = tdefs.calc_tof_bm_data(lower_lim, points_per_turn, nturns_sort+1)
	istart = time_indices[0]
else:
	#tof_m1, indices_m1, ttime_a_m1 = tdefs.calc_tof_bm_data(lower_lim, points_per_turn, nturns_sort+1)
	tof_mean, time_indices, ttime_a = tdefs.calc_tof_rf_waveform(tdat_rf, sig_rf, nturns_sort+1, istart_guess_l[0])
	istart = time_indices[0]

	if len(ttime_a) < nturns_sort+1:					
		print "nturns limited by data to ",len(ttime_a)
	
	print "tof mean ",tof_mean

	
	if calc_rf_time_offset:
		offset_file_ms = tdefs.time_offset_calc(dirname, ttime_a, t_phis_file, svk_data, fname_sig, write_offset_file, plot_time_offset)

	print "offset_file_ms ",offset_file_ms




#calculate kinetic energy, beta, synchrotron period at each selected time
if not flattop:
	Ekin_a = np.array([ke_lookup(t, svkfilepath, toffset_rf) for t in tomo_times])
else:
	ke_flattop = tdefs.ke_lookup(flattop_time, svkfilepath, toffset_rf)

	print "ke_flattop ",ke_flattop

	Ekin_l = []
	for t in tomo_times:
		if t < flattop_time:
			Ekin_l.append(ke_lookup(t, svkfilepath, toffset_rf))
		else:
			Ekin_l.append(ke_flattop)
	
	Ekin_a = np.array(Ekin_l)





				
#loop over times, do tomography on each loop
for index_time in range(len(tomo_times)):

	ke = Ekin_a[index_time] #kinetic energy
	
	#calculate points per turn based on average tof
	points_per_turn = int(round(tof_mean/time_interval)) 
	n_slices = int(npfac*points_per_turn)

	#calculate mean tof
	
	t_turn_data = np.linspace(0, tof_mean, n_slices)

	T_rf = tof_mean/harmonic
	T_rf_ns = 1e9*T_rf

	#work out longitudinal parameters of interest
	Ts, dtbin_phase,omega_rf, Rnom, Bnom, phase_sep_deg, dp_sep, bh, sync_bin, synch_index, binslz, binsuz, tadj, Hsep, phi_l1, phi_l2 = tdefs.longitudinal_param(V_rf, phi_s, harmonic, alpha_c, rho, ke, tof_mean, n_slices, nturns_sort, mass)

	Etot_turn, dEbin_c_turn, p0_turn, eta_turn, beta_turn = tdefs.longitudinal_turn(V_rf, phi_s, harmonic, alpha_c, ke, Rnom, nturns_sort, dtbin_phase, mass)

	p0_mev_turn = 1e-6*p0_turn

	#calculate delay in rf signal
	if idelay == None:
		prof_data_a = tdefs.sort_data(istart, time_data, signal_data, ttime_a, time_indices, nturns_sort, n_slices, interpolate_data)
		#sort data into turns covering roughly a synchrotron period
		nt_Ts = int(round(Ts))
		idelay = tdefs.calc_sync_point(prof_data_a, nt_Ts, sync_bin)

	#introduce further delay for non-zero phi_s (assuming idelay applies to a stationary bucket)
	if idelay != None and idelay_corr_phi_s:
		#calculate further shift in delay to account for phi_s !=0
		idelay_phi_s = int(n_slices*phi_s_deg/360)

		idelay = idelay - idelay_phi_s
		print idelay


	#sort data into turns, interpolate onto fixed time axis (phase)
	istart = istart + idelay
	prof_data_a = tdefs.sort_data(istart, time_data, signal_data, ttime_a, time_indices, nturns_sort, n_slices, interpolate_data, nshift= idelay) 

	if show_mountain_range:
		bcktlim1 = binslz 
		bcktlim2 = n_slices - binsuz
		tdefs.mountain_range(prof_data_a, nturns_sort, n_slices, bucketlim_low=bcktlim1, bucketlim_high=bcktlim2, sync_bin = sync_bin)			

	if show_imshow:
		tlim1 = 1e3*(ttime_a[0] - toffset_rf)
		tlim2 = 1e3*(ttime_a[-2] - toffset_rf)
		extent = [0, n_slices, tlim1, tlim2]
		tdefs.imshow_plot(prof_data_a,extent,ylabel='time')

	if intensity_calc:
		tdefs.intensity_calc_fn(prof_data_a, ttime_a, time_interval_ns, t_phis_file = t_phis_file)

	
	
	t_turn_data_ns = 1e9*t_turn_data + tadj

	chisq_sum_l = []

	if nscan_param == 1:
		idelay_a = np.array([idelay])
	else:

		if phi_s_deg == 0:
			idelay_a = np.arange(-30, -18)
		else:
			idelay_a = np.arange(-35, -18)

	print "idelaya ", idelay_a

	#number of turns specified in terms of sync period
	if specify_absolute_turns:
		nturns = turns_absolute
	else:
		nturns = int(turns_syncfrac*Ts_a[index_time])

	#this is not a real Bdot - it is a parameter used to specify phis in the tomo code
	Bdot = (V_rf*np.sin(phi_s))/(rho*2*math.pi*Rnom)
	time_str = str(1e3*tomo_times[index_time])
	for iscan in range(nscan_param):


		if nscan_param > 1:
			idelay1 = idelay_a[iscan]
			istart = istart + idelay1
			prof_data_a, t_turn_data = sort_data(istart, ttime_a, time_indices, nturns_sort, n_slices, interpolate_data, nshift= idelay1)
		

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
	
		for it in range(nturns):
			yprof = prof_data_a[it]
			yprof_norm = yprof/max(yprof)
			#plt.plot(yprof,'ko-')
			for y in yprof:
				print >>f1,y
		f1.close()

		run_tomo_code = True
		if run_tomo_code:
		
			if not specify_absolute_turns:
				recon_step = int(recon_step_frac*Ts_a[index_time])
			input_settings = tdefs.get_input("input_default.dat")

			fname = 'input_v2.dat'

			descrip = 'KURRI'
			datafile = 'kurridata.txt'
			input_settings["numframes"] = nturns
			input_settings["numbins"] = n_slices
			input_settings["synchbin"] = synch_index #int(0.5*n_slices)
			input_settings["framebinwidth"] = dtbin_phase/npfac	
			input_settings["binslowerignore"] = binslz
			input_settings["binsupperignore"] = binsuz
			input_settings["filmstart"] = recon_start
			input_settings["filmstop"] = recon_stop
			input_settings["filmstep"] = recon_step
			input_settings['Bref'] = Bnom
			input_settings['Bdot'] = Bdot
			input_settings['VRF1ref'] = V_rf
			input_settings['Rnom'] = Rnom
			input_settings['rhonom'] = rho
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

			
	

			#dEbin_calc = dtbin*harmonic*omega_
			fname_prof = 'profiles.data'
			profiles_norm_raw = tdefs.read_data_file(tomo_outdir+'/'+fname_prof)
			profiles_norm_data_a = np.split(profiles_norm_raw, nturns)

			fname_prof = 'profiles_recon.data'
			profiles_norm_recon_raw = tdefs.read_data_file(tomo_outdir+'/'+fname_prof)
			profiles_norm_recon_a = np.split(profiles_norm_recon_raw, nturns)


			

			show_profile_out = False
			if show_profile_out:

				fname_raw = 'kurridata.txt'

				bunchdata= tdefs.read_data_file(fname_raw)
				raw_spl = np.split(bunchdata, nturns)

				
				tdefs.mountain_range(profiles_norm_data_a, nturns, n_slices, fmount='mountainrange_data')


				tdefs.mountain_range(profiles_norm_recon_a, nturns, n_slices, fmount='mountainrange_recon')



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
		ndim = n_slices - (binslz + binsuz)  #int(len_img_data**0.5)

		t_synch = 1e9*phi_s_out*T/(2*math.pi)

		t_tomo_ns = 1e9*dtbin*(np.array(range(ndim)) - x0) + t_synch

		phase_tomo = t_tomo_ns*(2*math.pi)/T_rf_ns
		phase_tomo_deg = phase_tomo*180/math.pi

		
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
			discrep_data = read_data_file_2D(fname_discrep)
			discrep_final = discrep_data[-1]
			discrep_total = discrep_total + discrep_final
			print "discrepancy final, total ",discrep_final, discrep_total  
			
			
			
			fname_img = tomo_outdir+'/image000'[:-len(str(irc))]+ str(irc) + '.data'
			print "fname_img ",fname_img
		
			time_fread0 = time.time()
			image_data_a = tdefs.read_data_file(fname_img) #1D data

			time_fread1 = time.time()
			#print "time to read file ",time_fread1 - time_fread0


	
			#pts_above = np.where(image_data_a > 4*np.std(image_data_a))[0]
			#print "numpy of points above limit ",pts_above.shape
			
			#split into array
			image_data_spl = np.split(image_data_a, ndim)
			#transpose to get phase on horizontal axis, energy on vertical
			image_data_tr = np.transpose(image_data_spl)

			#construct dp pixel array taking synchronous point x0, y0 into account.
			E_tot_eV_a = Etot_turn[irc-1] + dEbin_c_turn[irc-1]*(np.array(range(ndim)) - y0)
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
				dp_sep_recon = np.array([list(dp_sep1)])
			else:
				tomo_time_proj_a = np.vstack((tomo_time_proj_a, tomo_time_proj))
				image_Hgrid_a = np.vstack((image_Hgrid_a, image_Hgrid))
				dp_a_recon = np.vstack((dp_a_recon, dp_a))
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
		versus_phim = True #phim is phase - phis
		if versus_time:
			xaxis = t_tomo_ns
		elif versus_phim:
			xaxis = phase_tomo_deg - phi_s_deg
			phi_off = phi_s_deg
		else:
			xaxis = phase_tomo_deg
			phi_off = 0


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
		#print "max dp/p along separatrix ",max(dp_sep_recon[ipl])


		plt.subplot(211)
		plt.title('reconstructed turn '+str(irc))
		if plot_contour:
			ctom = plt.contour(xaxis, dp_a_recon[ipl], image, levels = levels)
			#ctom = plt.contour(t_tomo_ns, dp_a, image)
			#plt.colorbar(ctom)
		
		if versus_time:
			plt.plot(t_sep_ns, dp_sep,'r') #bucket
			plt.plot(t_sep_ns, -dp_sep,'r') #bucket
			plt.axvline(x=t_tomo_ns[int(x0)]) #synch time
			plt.xlim(1.01*t_sep_ns[0],1.01*t_sep_ns[-1])
			plt.xlabel('time [ns]')
		else:
			plt.plot(phase_sep_deg - phi_off, dp_sep_recon[ipl],'r') #bucket
			plt.plot(phase_sep_deg - phi_off, -dp_sep_recon[ipl],'r') #bucket
			plt.axvline(x=phi_s_deg - phi_off) #synch time
			plt.xlabel(r'$\phi$ [deg]')
		
		plt.ylim(1.01*bh, -1.01*bh)

		x1,x2,y1,y2 = plt.axis()
		
		plt.ylabel(r"$\Delta$p/p")

		plt.subplot(212)

		plt.plot(xaxis, profiles_norm_data_a[irc-1],'ko-',linewidth=2,label='data')		
		plt.plot(xaxis, tomo_time_proj_a[ipl],'r-',linewidth=2,label='tomography')
		
		if versus_time:
			if binslz > 0:
				plt.axvline(x=t_turn_data_ns[binslz] + 0.1*t_data_step, color='gray', linestyle='--')
			if binsuz > 0:
				plt.axvline(x=t_turn_data_ns[n_slices - binsuz]-0.1*t_data_step, color='gray', linestyle='--')		
			plt.axvline(x=t_synch, color='gray', linestyle = '--')
			plt.xlim(1.01*t_sep_ns[0],1.01*t_sep_ns[-1])
			plt.xlabel('time [ns]')		

		else:
			#plt.axvline(x=phi_s_deg - phi_off, color='gray', linestyle = '--')
			plt.xlabel(r'$\phi$ [deg]')		

		
		plt.ylabel('normalised intensity')
		plt.legend(loc = 'upper right')
		plt.savefig('phasespace_'+str(irc))
		plt.tight_layout()

		plt.show()

		ipl = ipl + 1








