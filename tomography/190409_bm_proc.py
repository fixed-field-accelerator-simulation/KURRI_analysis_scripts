from __future__ import division
import pylab as plt
import numpy as np
import sys
import os
import matplotlib.cm as cm
import math 
import tomo_defs_190408 as tdefs
from scipy.optimize import leastsq

dirname = "../data/2019/3/28/"
#dirname = "../data/2019/4/20190401/"

read_file_offset = False
read_rf_waveform = True
multi_channel = True

show_mountain_range = False
intensity_calc = True
show_imshow = True

start_time = 0.15e-3#0.12e-3
nturns_sort = 2000



#handle arguments passed to script
if "mountain" in sys.argv:
	show_mountain_range = True

#read the data
interactive_data_sel = True
for arg in sys.argv:	
	if 'csv' in arg:
		path = arg
		interactive_data_sel = False
		fspl = arg.split('/')
		fname_sig = fspl[-1]


if interactive_data_sel:
	fname_sig, time_data, signal_data, tdat_rf, sig_rf, t_phis_file, phis_file, offset_file_ms = tdefs.read_data_files(dirname, read_rf_waveform = read_rf_waveform, read_file_offset = read_file_offset, multi_channel = multi_channel)
else:
	if multi_channel:
		time_data, waveform_data = tdefs.read_scope(path, multi_channel = multi_channel)

		waveform_data = np.array(waveform_data)
		waveform_data_tr = np.transpose(waveform_data)
	
		tdat_rf = time_data
		sig_rf = waveform_data_tr[0]
		signal_data = waveform_data_tr[1]


#show_rf = True
rf_fft = False
if rf_fft:
	t1 = 0.29e-3
	t2 = 0.32e-3
	#t1 = 2.3e-3
	#t2 = 2.32e-3
	i1 = None
	i2 = None
	for i in range(len(tdat_rf)):
		if tdat_rf[i] > t1 and i1 == None:
			i1 = i
		elif tdat_rf[i] > t2 and i1 != None:
			i2 = i
			break
	print "i1, i2 ",i1,i2
	
	if i1 == None or i2 == None:
		print "selected times not found in rf waveform"
		sys.exit()

	tdat_rf_sel_ms = 1e3*np.array(tdat_rf[i1:i2])
	sig_rf_sel = sig_rf[i1:i2]

	fft_rf = np.fft.fft(sig_rf_sel)
	fft_rf_abs = np.abs(fft_rf)
	nfft = len(fft_rf)
	fft_rf_half = fft_rf_abs[:int(0.5*nfft)]
	imax =  np.argmax(fft_rf_half)
	print "dominant freq, imax ",imax/nfft, imax
	



	fft_sel = np.zeros(nfft, dtype = complex)
	fft_sel[0] = fft_rf[0]
	fft_sel[imax] = fft_rf[imax]
	fft_sel[-imax] = fft_rf[-imax]
 	
	fft_sel[3*imax+1] = fft_rf[3*imax+1]
	fft_sel[-3*imax-1] = fft_rf[-3*imax-1]

	plt.plot(fft_rf,'k')
	plt.plot(fft_rf_abs,'k--')
	plt.plot(fft_sel)
	plt.show()


	sig_invfft = np.fft.ifft(fft_sel)

	plt.plot(sig_rf_sel, 'ko')
	plt.plot(sig_invfft,'r')
	plt.show()
	sys.exit()
	sys.exit()

time_interval = time_data[2] - time_data[1]
time_interval_ns = 1e9*time_interval

#idelay is delay in time bins between the real RF phase seen by the beam and that seen on the scope
if round(time_interval_ns) == 2:
	idelay = -47
else:
	idelay = int(round(-47*(2/time_interval_ns)))

print "time interval per bin (ns) ",time_interval_ns
print "delay (number of bins, time[ns]) ",idelay, idelay*time_interval_ns

tof_mean, tof_l, time_indices, ttime_a = tdefs.calc_tof_rf_waveform(tdat_rf, sig_rf, nturns_sort+1, start_time)
istart = time_indices[0]	


freq_compare = False
if freq_compare:

	time_offset = 0.1143*1e-3
	nskip = 500
	ttime_a_sel = ttime_a[:-1][::nskip]
	tof_sel = np.array(tof_l)[::nskip] #ToF based on zero-crossing in waveform
	f_data = 1/tof_sel #RF frequency from waveform
	
	svkfilepath = "../sep18_exp/svk20.dat"
	ke_tab_sel = [tdefs.ke_lookup(t, svkfilepath, toffset = time_offset) for t in ttime_a_sel]
	f_tab_sel = [tdefs.f_rf_lookup(t, svkfilepath, toffset = time_offset) for t in ttime_a_sel]
	
	ttime_a_sel_ms = 1e3*ttime_a_sel
	
	plt.plot(ttime_a_sel_ms, f_tab_sel,'r-',label='tabulated')
	plt.plot(ttime_a_sel_ms, f_data,'ko',label='RF waveform')
	plt.axvline(x=time_offset*1e3, color='gray')
	plt.xlabel('time from scope trigger [ms]')
	plt.legend(loc = 'upper left')
	plt.savefig('freq_compare_rf_svk')
	plt.show()
	sys.exit()
	
	
ttime_ms = 1e3*ttime_a[:-1]
nbins = int(round(tof_mean/time_interval)) #number of time bins

interpolate_data = True 
prof_data_a = tdefs.sort_data(istart, time_data, signal_data, ttime_a, time_indices, nturns_sort, nbins, interpolate_data, nshift=idelay) 

#translate bins into phase
phase_ax = np.linspace(-180,180,nbins)

if show_mountain_range:
	tdefs.mountain_range(prof_data_a, nturns_sort, nbins, incl_xaxis = True, xaxis= phase_ax, xlabel='phase [deg]')#, extra_data = [mean_phase_data, phase_profmax])			

write_intensity_file = True
if intensity_calc:
	intensity, mean_phase_data, phase_fmax, std_prof_max = tdefs.intensity_calc_fn(prof_data_a, phase_ax, ttime_a, time_interval_ns, dirname, fname_sig, write_intensity_file = write_intensity_file)
	


min_mean = min(mean_phase_data)
index_min_mean = mean_phase_data.index(min_mean)
max_mean = max(mean_phase_data[index_min_mean:])

t_mean_1 = ttime_ms[index_min_mean]

index_max_mean = mean_phase_data.index(max_mean)
t_mean_2 = ttime_ms[index_max_mean]

nw = 31
imid = int(0.5*(nw-1))
phase_fmax_movav = tdefs.running_mean(phase_fmax,nw)
t_movav = ttime_ms[imid:-imid]

#calculate extreme points based on phase_fmax_movav
index_min_max = np.argmin(phase_fmax_movav)
index_max_max = np.argmax(phase_fmax_movav)
t_max_1 = t_movav[index_min_max]
t_max_2 = t_movav[index_max_max]

print "index_min_mean ",index_min_mean, index_max_mean
print "len phase_fmax, phase_fmax_movav ",len(phase_fmax),len(phase_fmax_movav)
print "all times are in ms and are w.r.t scope trigger"
print "time of minimum in mean, required phase jump ",t_mean_1, 20 - phase_fmax_movav[index_min_mean-imid] 
print "time of maximum in mean, required phase jump",t_mean_2, 20 - phase_fmax_movav[index_max_mean-imid] 
print "time of minimum in phase_fmax_movav, required phase jump", t_max_1,20 - phase_fmax_movav[index_min_max-imid]
print "time of maximum in phase_fmax_movav, required phase jump", t_max_2,20 - phase_fmax_movav[index_max_max-imid]

if show_imshow:
	tlim1 = 1e3*ttime_a[0]
	tlim2 = 1e3*ttime_a[-1]
	extent = [-180, 180, tlim1, tlim2]
	tdefs.imshow_plot(prof_data_a,extent,ylabel='time [ms]',xlabel='arbitrary phase [deg.]', extra_data = [mean_phase_data, phase_fmax, [t_movav, phase_fmax_movav]])
	
plot_vs_time = False
if plot_vs_time:
	plt.plot(ttime_ms, mean_phase_data, 'k',label='<phase>')
	plt.plot(ttime_ms, phase_fmax, 'm',label='phase@dist_max')
	plt.plot(t_movav, phase_fmax_movav, 'b--', linewidth=2)
	plt.axvline(x=t_max_1, color='gray',linestyle='--')
	plt.axvline(x=t_max_2, color='gray',linestyle='--')
	plt.xlabel('time [ms]')
	plt.ylabel('phase [deg]')
	plt.legend()
	plt.savefig('phase_mean_max')
	plt.show()

else:
	turn_ax = np.array(range(1, nturns_sort+1))
	turn_ax_movav = np.array(turn_ax[imid:-imid])

	mean_phase_data = np.array(mean_phase_data)
	phase_fmax = np.array(phase_fmax)

	PROTON_MASS = 938.27231e6 #938.27203e6
	mass = PROTON_MASS

	time_offset = 0.1143*1e-3
	nskip = 499
	ttime_a_sel = ttime_a[0]
	svkfilepath = "../sep18_exp/svk20.dat"
	ke_turn1 = tdefs.ke_lookup(ttime_a[0], svkfilepath, toffset = time_offset)
	phi_s_nom = 20*math.pi/180
	V_rf_nom = 4e3
	harmonic = 1
	kindex = 7.6
	alpha_c = 1/(kindex+1)
	ke_turn_a = ke_turn1 + turn_ax*V_rf_nom*np.sin(phi_s_nom)
	E_turn_a = ke_turn_a + mass
	p_turn_a = np.sqrt(E_turn_a**2 - mass**2)
	betarel_a = p_turn_a/E_turn_a
	betagam_a = p_turn_a/mass
	gamma_a = betagam_a/betarel_a
	eta_a = (alpha_c - gamma_a**(-2))
	coef_a = (abs(eta_a)/((betarel_a**2)*(2*math.pi*E_turn_a)))**0.5
	coef_rel = coef_a/coef_a[0]
	#print "ke ini,final,diff ",ke_turn_a[0], ke_turn_a[-1], ke_turn_a[-1] - ke_turn_a[0]
	#print "E final/ini ",E_turn_a[-1]/E_turn_a[0]
	#plt.plot(coef_a*np.sqrt(2100*np.cos(40*math.pi/180)))
	#plt.show()


	

	nt_fit = 500
	coef_rel_fit = coef_rel[-nt_fit:]
	data_ini = np.array(phase_fmax_movav[:nt_fit])
	data_end = np.array(phase_fmax_movav[-nt_fit:])
	#data = np.array(mean_phase_data[:nt_fit])
	guess_mean = np.mean(data_ini)
	guess_amp = 0.5*(max(data_ini) - min(data_ini))

	fit = lambda x, turn_a: x[0]*np.sin(2*math.pi*x[3]*turn_a+x[1]) + x[2]


	optimize_func = lambda x, turn_a, data_in: x[0]*np.sin(2*math.pi*x[3]*turn_a+x[1]) + x[2] - data_in
	#optimize_func = lambda x, turn_a: x[0]*np.sin(2*math.pi*x[3]*coef_rel_fit*turn_a+x[1]) + x[2] - data
	guess_tune = 0.003
	guess_phase = 0
	x = [guess_amp,guess_phase, guess_mean, guess_tune]

	res = leastsq(optimize_func, x, args=(turn_ax_movav[:nt_fit], data_ini,))
	est_amp_ini, est_phase_ini, est_mean_ini, est_tune_ini = res[0]



	print "initial 500 turns"
	print "amp, phase, tune, mean found by least squares fit ",est_amp_ini, est_phase_ini, est_tune_ini, est_mean_ini
	guess_amp = 0.5*(max(data_end) - min(data_end))
	guess_mean = np.mean(data_end)
	guess_tune = est_tune_ini
	x = [guess_amp,guess_phase, guess_mean, guess_tune]
	res = leastsq(optimize_func, x, args=(turn_ax_movav[:nt_fit], data_end,))
	est_amp_end, est_phase_end, est_mean_end, est_tune_end = res[0]
	print "final 500 turns"
	print "amp, phase, tune, mean found by least squares fit ",est_amp_end, est_phase_end, est_tune_end, est_mean_end

	mean_mean_phase_data = np.mean(mean_phase_data)
	print "mean mean_phase_data, phase_fmax_movav ",mean_mean_phase_data, np.mean(phase_fmax_movav)
	
	
	phi_s_deg = mean_mean_phase_data
	#phi_s_deg = 20
	phi_s = phi_s_deg*math.pi/180
	V_rf = V_rf_nom*np.sin(phi_s_nom)/np.sin(phi_s)

	print "phi_s, V_rf ",phi_s_deg,V_rf

	phase_m = 25*math.pi/180
	Qs_nom_a = tdefs.synch_tune(V_rf_nom, phi_s_nom, harmonic, alpha_c, ke_turn_a, mass)
	Qs_a = tdefs.synch_tune(V_rf, phi_s, harmonic, alpha_c, ke_turn_a, mass)
	
	Qs_amp_nom_a = Qs_nom_a*(1 - (phase_m**2/16)*(1+(5/3)*((np.tan(phi_s_nom))**2)))

	Qs_amp_a = Qs_a*(1 - (phase_m**2/16)*(1+(5/3)*((np.tan(phi_s))**2)))


	
	print "synch tune zero_amp, at amp ",Qs_a[0], Qs_amp_a[0]
	


	#est_tune = 0.0029
	
	phi_offset = 0
	plt.plot(turn_ax, mean_phase_data - phi_offset, 'k',label='<phase>')
	plt.plot(turn_ax, phase_fmax - phi_offset, 'm',label='phase@dist_max')
	plt.plot(turn_ax_movav, phase_fmax_movav - phi_offset, 'b--', linewidth=2)
	plt.plot(turn_ax_movav, fit([est_amp_ini, est_phase_ini, est_mean_ini, est_tune_ini], turn_ax_movav)- phi_offset,'c')
	plt.plot(turn_ax_movav, fit([est_amp_end, est_phase_end, est_mean_end, est_tune_end], turn_ax_movav)- phi_offset,'c--')
	
	plt.axhline(y=phi_s_deg)

	plt.xlabel('turn')
	plt.ylabel('phase [deg]')
	#plt.legend()
	plt.savefig('phase_fit')
	plt.show()


	plt.plot(Qs_nom_a, 'k-')
	plt.plot(Qs_amp_nom_a, 'k--')
	plt.plot(Qs_a, 'r-')
	plt.plot(Qs_amp_a, 'r--')
	plt.plot([1],[est_tune_ini],'bo')
	plt.plot([nturns_sort-1],[est_tune_end],'bo')
	plt.xlabel('turn')
	plt.ylabel('synchrotron tune')
	plt.show()
