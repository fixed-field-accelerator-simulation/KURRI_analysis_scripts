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

#Calculate time to loss at various probe positions from bunch monitor data.
#30/9/2015 This version makes use of an updated algorithm signal_loss_time_new to find the loss time.
#7/8/2018 Option to process data from multiple probes.


dirname = "../data/2015/6/24/700A/"
channel_id = ["F01"] #identifier for probes to be analysed 

#indices_select used to select particular data files associated with each probe
indices_select = [15] #process all files
#indices_select = range(1,5) #select files with index 1 to 4 for each probe
#indices_selec
t = [0,1,2] #select first three files for each probe

#select channel files
ch_files_sel = adata.select_filenames(dirname, indices_select, channel_id)
print "selected data file ",ch_files_sel

clist = ['r','b','g','c']


	
def find_peaks(xdat, ydat, interval, find_hminus, x_low = None, x_high = None, find_maxima=False):
    
	"""find_maxima = True: Find maxima. Otherwise find minima"""

	pmxi = []
	pmx = []
	pmy = []
	pmi = []
	pmy_av = []
	pmy_std = []
    
	if x_low != None:
		for x in xdat:
			if x > x_low:
				index_low = xdat.index(x)
				break
	else:
		index_low = 0
        
	if x_high != None:
		for x in xdat:
			if x > x_high:
				index_high = xdat.index(x)
				break
	else:
		index_high = 0    

	npts = index_high - index_low + 1
    
	winsize = int(0.5*interval)

	for i in range(int(npts/interval)):

		my = 0.
		my_min = 100
		
            
		if i == 0 :
			istart = index_low 
			xwindow = xdat[index_low:index_low+2*winsize]
			ywindow = ydat[index_low:index_low+2*winsize]

			icen = index_low
			xwindow = xdat[icen-winsize: icen+winsize]
			ywindow = ydat[icen-winsize: icen+winsize]
            

			istart = icen - winsize      
		else:
			icen = mi + interval

			#test if in range
			if icen + winsize > len(xdat):
				print "reach end of data ",icen+winsize, len(xdat)
				break
            
			xwindow = xdat[icen-winsize: icen+winsize]
			ywindow = ydat[icen-winsize: icen+winsize]
            
 
			istart = icen - winsize
       
		if find_maxima:
			ywindow = [-y for y in ywindow]

		for xd,yd in zip(xwindow, ywindow):
			#look for positive H- peak
			if i == 0 and find_hminus:
				if yd > my:
					mx = xd
					my = yd
 					mi = ywindow.index(yd) + istart
			else:
			#look for negative proton peak
				if yd < my:
					mx = xd
					my = yd
					mi = ywindow.index(yd) + istart
		#print "iteration, mi ",i, mi   
		my = my - max(ywindow)
        
		my_av = sum(ywindow)/len(ywindow)
		my_std = np.std(ywindow)
                        
		if i == -1:            
			plt.plot(xwindow,ywindow,'k.-')
			plt.plot([mx],[my],'ro')
			plt.show()
			sys.exit()  

		pmxi.append(i)
		pmx.append(mx)
		pmy.append(my)
		pmi.append(mi)
		pmy_av.append(my_av)
		pmy_std.append(my_std)
        
	return pmxi, pmx, pmy, pmi, pmy_av, pmy_std

def profile_data(xdat, ydat, indices, interval):
	"""prepare profile data and write to file"""

	winsize = int(0.5*interval)

	xax = range(1,interval+1)
	ip = 0
	col_l = ['k','r','b','m','g']
	plot_proc = True
	
	f1 = open('kurridata.txt', 'w')
	y_prof_l = []
	for i1 in range(len(indices)):

		icen = indices[i1] + 7

		xwindow = xdat[icen-winsize: icen+winsize]
		ywindow = ydat[icen-winsize: icen+winsize]	

		ybase  = min(-ywindow)
		#ybase  = max(ywindow)
		#yprof = -(ywindow-ybase)
		
		yprof = -ywindow-ybase
		#yprof = -ywindow
		y_prof_l.append(yprof)
		
		
		for y in yprof:
			print >>f1, y
			
		if plot_proc:
			
			if (i1+1)%1 == 0 or i1 == 0:
			
				print "index plot ",i1
				plt.plot(xax, yprof, col_l[ip%len(col_l)]+'o-',label=str(i1+1), linewidth=2)
				ip = ip + 1

	if plot_proc:
		plt.xlim(1, interval)
		plt.legend(loc = 'upper left', title="turn")
		plt.ylabel('-signal - min(-signal)')
		plt.axvline(x=8, linestyle='--', color='gray')
		plt.savefig('bunchmonitor_proc')
		plt.show()
		#sys.exit()

	f1.close()	
	
	return y_prof_l

#read scope data
tdat_chf, data_chf = adata.read_scope(dirname, ch_files_sel[0][0])

index_hminus = data_chf.index(max(data_chf))
t_hminus = tdat_chf[index_hminus]

tof = 600e-9 #640e-9 #640e-9 #1.001*600e-9
time_interval = tdat_chf[2] - tdat_chf[1]
nturn_points = int(round(tof/time_interval)) # points = intervals + 1

#set parameters
phi_s_deg = 20#20
n_slices = nturn_points
phi_s = phi_s_deg*math.pi/180 
kindex = 7.6
alpha_c = 1/(kindex+1)
gamma_tr = (1/alpha_c)**0.5

#gamma_tr = 2.932
#alpha_c = gamma_tr**(-2)
V_rf = 4e3 # in [V]
harmonic = 1
#pipe_radius = 5e-2

#Bdot=0 # in T/s
bending_radius=0.489396551697 # in m
Rnom = 4.41
circumference = Rnom*2*np.pi
phi_offset = 0 #np.arcsin(bending_radius*circumference*Bdot/V_rf) # measured from aligned focussing phase (0 or pi)

Ekin = 11e6
PROTON_MASS = 938.27203e6
Etot = Ekin + PROTON_MASS
SPEED_OF_LIGHT = 299792458
	
beta = adata.ke_to_relativistic_beta(Ekin, PROTON_MASS)
betagam = adata.ke_to_relativistic_beta_gamma(Ekin, PROTON_MASS)
gamma = betagam/beta
eta = alpha_c - gamma**(-2)
p0_mev = 1e-6*np.sqrt(Etot**2 - PROTON_MASS**2)
Qs = np.sqrt(np.abs(eta) * V_rf / (2 * np.pi * beta**2 * Etot))

print "beta, gamma ",beta, gamma
print "eta ",eta
print "p0_mev ",p0_mev
	
omegarev = beta*SPEED_OF_LIGHT/Rnom #angular rev. frequency
T = tof#2*np.pi/(harmonic*omegarev) #rf period
T_ns = 1e9*T



bh = ((2*Qs)/(harmonic*abs(eta)))

omega_rf = omegarev*harmonic
a1 = 0.5*harmonic*omega_rf*eta
a2 = omega_rf*V_rf/(2*math.pi*(beta**2)*Etot)

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
dp_sep, t_phi_ns = tdefs.separatrix(phi_s, eta, T, harmonic, npts, a1, a2)
phase_sep = 2*math.pi*(t_phi_ns)/(1e9*T)


t_bucket = np.linspace(t_phi_ns[0], t_phi_ns[-1], n_slices)
t_data = np.linspace(-0.5*T_ns, 0.5*T_ns, n_slices) #data fills rf period



print "bucket phase width [rad, deg]",bucket_phase_width,180*bucket_phase_width/math.pi
print "bucket duration calc from phase [ns] ", 

bucket_dur_ns = T_ns*(bucket_phase_width/(2*math.pi))
	

#bucket_dur_ns = t_phi_ns[-1] - t_phi_ns[0]
outbucket_dur_ns = T_ns - bucket_dur_ns
print "rf period, bucket duration ",T_ns, bucket_dur_ns, outbucket_dur_ns, 
print "fraction time out of bucket ",outbucket_dur_ns/T_ns

print "number points in bucket ",nturn_points*bucket_dur_ns/T_ns
print "number points out of bucket ",nturn_points*outbucket_dur_ns/T_ns

bucket_points = int(nturn_points*bucket_dur_ns/T_ns)

print "bucket points, out-of-bucket points ",bucket_points, nturn_points - bucket_points

binslowerignore = 3 #add zeroes before profile (earlier time)
binsupperignore = 2#add zeroes after profile (later time)
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

sync_point =  bucket_points*(t_phi_sep_ns- t_phi_ns[0])/bucket_dur_ns
print "synch_frac ",sync_point

Bdot = (V_rf*np.sin(phi_s))/(bending_radius*2*math.pi*Rnom)

nbturns = 100

#crude start point
injgap =  7e-4 #time from Hminus
icrude = index_hminus + int(round(injgap/time_interval))
t_crude = tdat_chf[icrude]
winsize = int(0.5*nturn_points)
xwindow_ini = tdat_chf[icrude-winsize: icrude+winsize]
ywindow_ini = data_chf[icrude-winsize: icrude+winsize]

find_hminus = False
if not find_hminus:
	ipeak = ywindow_ini.index(min(ywindow_ini))
lower_lim = xwindow_ini[ipeak]
upper_lim = lower_lim + nbturns*tof


print "lower, upper times ",lower_lim, upper_lim

tdat_chf_ms = 1e3*np.array(tdat_chf)

show_data = False
if show_data:
	
	plt.subplot(211)
	plt.plot(tdat_chf_ms, data_chf,'k')
	#plt.axvline(x=t_hminus, color='r')
	plt.axvline(x=1e3*lower_lim, color='r')
	plt.axvline(x=1e3*upper_lim, color='r')
	#plt.xlim(lower_lim, upper_lim)
	plt.ylabel('bunch monitor (all)')
	
	plt.ylim(-0.2, 0.2)
	plt.xlim(xmin = t_hminus*1e3 , xmax = t_hminus*1e3 + 1)
	plt.subplot(212)
	plt.plot(tdat_chf_ms, data_chf,'k')
	plt.xlim(1e3*lower_lim, 1e3*(lower_lim + 1e-5))
	plt.ylabel('bunch monitor (zoom)')
	plt.xlabel('time [ms]')
	#plt.savefig('bunchmonitor_data')
	plt.show()
	#sys.exit()




print "H- signal time ",t_hminus
print "time step ",time_interval
print "nturn_points ",nturn_points
print "no. data points ",len(data_chf)

#harmonic = 1
#omegarev = beta*c/Rnom #angular rev. frequency
#T = 2*np.pi/(harmonic*omegarev) #rf period

#pmxi, pmx, pmy, pmi, pmy_av, pmy_std = find_peaks(tdat_chf, data_chf, nturn_points, find_hminus, x_low = lower_lim, x_high = upper_lim)

#t_peak_av, peak_av1 = adata.movingaverage(pmx, pmy, win_peak_av)			


print "interval ", nturn_points

sort_data = True

if sort_data:

	nbturns_fit = 400
	upper_lim_fit = lower_lim + nbturns_fit*tof
	
	_, pmx, _, pmi_peak, _, _ = find_peaks(tdat_chf, data_chf, nturn_points, find_hminus, x_low = lower_lim, x_high = upper_lim_fit, find_maxima = True)
	
	#the minimum rather than peak signalx
	t_baseline = [tdat_chf[i] for i in pmi_peak]
	
	t_delta = [0.1*(t_baseline[i] - t_baseline[i-10]) for i in range(10, len(t_baseline))] 
	#print "t_delta ",t_delta
	t_delta_mean = sum(t_delta)/len(t_delta)
	print "mean time step ",t_delta_mean
	
	tp = np.arange(nbturns_fit)
	fitpmx = np.polyfit(tp, pmx, 1)
	tof_fit = fitpmx[0]
	
	polyfitpmx = np.poly1d(fitpmx)
	
	diff = pmx - polyfitpmx(tp)
	meandiff = sum(diff)/len(diff)
	print "meandiff ",meandiff
	ampt = max(diff) - min(diff)
	timeosc = [ampt*np.sin(0.25*2*np.pi*Qs*i) - 0.5*ampt for i in np.arange(nbturns_fit)]
	
	plt.subplot(211)
	plt.plot(pmx,'ko-')
	plt.plot(polyfitpmx(tp))
	plt.subplot(212)
	plt.plot((pmx - polyfitpmx(tp)))
	plt.plot(timeosc)
	plt.show()
	#sys.exit()
	
	#T = t_delta_mean
	#plt.plot(t_delta)
	#plt.show()
	#sys.exit()
	
	
	print "len(pmi_peak) ",len(pmi_peak)
	
		
	istart =  pmi_peak[0] - 2 #+ 5 #3#index_start1 - nturn_points_half + index_max_win
	iend = istart + nbturns_fit*nturn_points
	
	tstart = tdat_chf[istart] #- 0.5*t_delta_mean 
	#print "max_in_win ",max_in_win, data_chf[istart]
	
	t_dat_all = np.array(tdat_chf[istart:iend])
	sig_dat_all = np.array(data_chf[istart:iend])

	#use ToF to define time.
	#t_dat_adj = tstart + np.linspace(0, nbturns*tof, nbturns*n_slices)
	
	syncosc = [np.sin(2*np.pi*Qs*i) for i in np.arange(nbturns_fit)]
	

	print "tstart ", tstart
	t_turn = np.linspace(0, tof_fit, n_slices)
	t_dat_adj = []
	
	timeosc = timeosc - timeosc[0]
	print "timeosc ",timeosc
	
	for it in range(nbturns_fit):
		t_this_turn = tstart + it*tof_fit + timeosc[it]+ t_turn
		t_dat_adj  = np.concatenate((t_dat_adj, t_this_turn))
		
		#print "t_this_turn ",t_this_turn
		#if it == 2:
		#	sys.exit()

	
	
	#interpolate onto time base Trev/n_slices
	sig_interp = np.interp(t_dat_adj, t_dat_all, sig_dat_all)
	t_dat_split = np.split(t_dat_adj, nbturns_fit)
	sig_dat_split = np.split(sig_interp, nbturns_fit)	
	
	#plt.plot(t_dat_adj, sig_interp,'ro-')
	#plt.plot(t_dat_all, sig_dat_all,'ko')
	#plt.show()
	
	


	use_original_time_data = False
	data_win_sig_sel = []
	if use_original_time_data:
	
		index_turn_start_l= [istart]	
		for it in range(nbturns):
			if it > 0:
				for j1, t1 in enumerate(tdat_chf[istart:iend]):
					if  abs(tstart + it*tof - t1) < 3e-8:
					
						index_turn_start = istart +j1
						index_turn_start_l.append(index_turn_start)
					
						print "turn ",it, " diff ",(tstart + it*tof) - t1
						break
		 			
		 			
		
		for it, ist in enumerate(index_turn_start_l):
			
			t_win = np.array(tdat_chf[ist:ist+nturn_points])
			data_win =np.array(data_chf[ist:ist+nturn_points])
			#print "ist, ist+nturn_points, len(data_win) ",ist, ist+nturn_points, len(data_win)
			#print "t_win limits ",t_win[0], t_win[-1],t_win[-1]- t_win[0]
			data_win_sig =  max(data_win) - data_win
			data_win_sig_sel.append(data_win_sig)
		
			plt.plot(range(1,nturn_points+1), data_win_sig + it*0.01,'.-')
	else:
		it = 0
		for data_win in sig_dat_split:
			data_win_sig =  max(data_win) - data_win
			data_win_sig_sel.append(data_win_sig)
			plt.plot(range(1,nturn_points+1), data_win_sig + it*0.01,'.-')
			it = it + 1
			
	plt.axvline(x=binslowerignore+0.5)
	plt.axvline(x=nturn_points - binsupperignore+0.5)
	plt.show()

	print data_win_sig_sel[1]
	#plt.plot(t_bucket, data_win_sig_sel[0], 'ko-')
	#plt.show()
	#sys.exit()
	
	prof_data_a = np.array(data_win_sig_sel)
	nbturns_sel = 2#nbturns#len(prof_data_a)

	


	f1 = open('kurridata.txt', 'w')
	
	print "write ", nbturns_sel, " turns to file"
	for it in range(nbturns_sel):
		yprof = prof_data_a[it]
		
		
		yprof_norm = yprof/max(yprof)
		#plt.plot(yprof,'ko-')
		
		
		for y in yprof:
			print >>f1,y
	f1.close()
	#plt.show()
	#sys.exit()
		

	#plt.plot(tdat_chf[index_turn_list[0]:index_turn_list[-1]], data_chf[index_turn_list[0]:index_turn_list[-1]],'k-')
	#for i in index_turn_list:
	#	plt.axvline(x=tdat_chf[i])

	#for i in pmi_peak:
	#	plt.axvline(x=tdat_chf[i])
	#plt.show()

	
			
	#pmi = index_turn_list

#prof_data_l = profile_data(tdat_chf, np.array(data_chf), pmi, nturn_points)
#prof_data_a = np.array(prof_data_l)

#number of points between peaks
#i_pksep = [pmi[i] - pmi[i-1] for i in range(1, len(pmi))]

#plt.plot(i_pksep,'ko')
#plt.show()
#sys.exit()

#plt.plot(pmx, pmy, 'ko-')
#plt.show()	




from subprocess import call
run_tomo_code = True
if run_tomo_code:
	input_settings = tdefs.get_input("input_default.dat")

	
	synch_index_ideal = 4	
	binwidth_time = time_interval
	

	fname = 'input_v2.dat'

	descrip = 'KURRI'
	datafile = 'kurridata.txt'
	input_settings["numframes"] = nbturns_sel
	input_settings["numbins"] = n_slices
	input_settings["synchbin"] = synch_index_ideal#int(0.5*n_slices)
	input_settings["framebinwidth"] = binwidth_time
	input_settings["binslowerignore"] = binslowerignore
	input_settings["binsupperignore"] = binsupperignore
	input_settings['Bdot'] = Bdot
	input_settings['VRF1ref'] = V_rf
	input_settings['Rnom'] = Rnom
	input_settings['rhonom'] = bending_radius
	input_settings['gammatnom'] = gamma_tr
	input_settings['pickup'] = 1.0
	input_settings['selffields'] = 0
	tdefs.write_inputfile(fname,descrip, datafile, input_settings)


	temp_files = ['plotinfo.data','profiles.data','d001.data','jmax.data','image001.data']
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

	#dp_sep, t_phi_ns = tdefs.separatrix(Qs, phi_s, eta0, T, 1)

	max_dp_sep = max(dp_sep)

	#process main tomography output file
	fname_img = 'image001.data'
	image_data_a = read_data_file(fname_img) #1D data	
	
	len_img_data = len(image_data_a)
	ndim = int(len_img_data**0.5)

	dE_a_eV = dEbin*(np.array(range(ndim)) - y0)
	dE_a_mev = 1e-6*dE_a_eV
	

	print "phi_s_out, phi_s_0 ",phi_s_out, phi_s
	print "phi_s_out, phi_s_0 ",phi_s_out*180/math.pi, phi_s*180/math.pi

	E_tot_eV_a = Etot + dE_a_eV 
	p_a = [1e-6*(Et**2 - PROTON_MASS**2)**0.5 for Et in E_tot_eV_a]

	#translate to time, dp
	
	dp_a = (p_a - p0_mev)/p0_mev
	t_synch = 1e9*phi_s_out*T/(2*math.pi)
	t_a_ns = 1e9*dtbin*(np.array(range(ndim)) - x0) + t_synch

	print "ndim ",ndim, len(t_a_ns)
	
	print "x0 ",x0
	print "t_a_ns ",t_a_ns
	
	

	#t_a_ns = np.linspace(t_phi_ns[0], t_phi_ns[-1], n_slices)

	#split into array
	image_data_spl = np.split(image_data_a, ndim)
	#transpose to get phase on horizontal axis, energy on vertical
	image_data_tr = np.transpose(image_data_spl)

	#tomo_profile = [sum(d) for d in image_data_tr]
	tomo_time_proj = [sum(d) for d in image_data_spl]
	tomo_mom_proj = [sum(d) for d in image_data_tr]


	#plt.plot(tomo_mom_proj)
	#plt.show()


	#print "t_cen lims ",t_cen[0], t_cen[-1], t_phi_ns[0], t_phi_ns[-1]



	show_dp = True
	if show_dp:


		#if run_pyheadtail:

			#counts, xbins, ybins,_ = plt.hist2d(bunch_t_ini_ns, bunch_dp_ini, bins=60, norm=LogNorm())


		plt.subplot(211)
		
		plt.title('reconstructed')
		ctom = plt.contour(t_a_ns, dp_a, image_data_tr, 20)
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
		plt.xlim(t_phi_ns[0],t_phi_ns[-1])

		x1,x2,y1,y2 = plt.axis()
		plt.xlabel('time [ns]')
		plt.ylabel(r"$\Delta$p [MeV]")

		

		
		#plt.contour(counts, 20, extent=[xbins.min(), xbins.max(), ybins.min(), ybins.max()])
		plt.ylabel(r"$\Delta$p [MeV]")

		plt.subplot(212)

		prof_data_norm = prof_data_a[0]/max(prof_data_a[0])
		tomo_t_proj_norm = tomo_time_proj/max(tomo_time_proj)

		
		plt.plot(t_data, prof_data_norm,'ko-',linewidth=2,label='data')
		plt.plot(t_a_ns, tomo_t_proj_norm,'ro-',linewidth=2,label='tomography')
		#plt.xlim(t_bucket[0], t_bucket[-1])

		plt.axvline(x=t_bucket[0], color='gray')
		plt.axvline(x=t_bucket[-1], color='gray')
		#plt.ylim(-1.1*max_dp_sep, 1.1*max_dp_sep) 

		
		#plt.xlim(t_a_ns[0],t_a_ns[-1])
		#plt.xlim(t_phi_ns[0], t_phi_ns[-1])
		#plt.ylim(-1.1*max_dp_sep, 1.1*max_dp_sep) 
		plt.xlabel('time [ns]')
		plt.ylabel('normalised intensity')
		plt.legend()
		plt.savefig('phasespace')
		plt.tight_layout()
		plt.show()
	else:
		plt.plot(t_phi_ns, dE_sep_mev,'r') #bucket
		plt.plot(t_phi_ns, -dE_sep_mev,'r') #bucket
		plt.contour(t_a_ns, dE_a_mev, image_data_tr, 20)
		plt.xlabel('[ns]')
		plt.ylabel(r"$\Delta$E [MeV]")
		plt.savefig('phasespace')
		plt.show()










