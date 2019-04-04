from __future__ import division
import pylab as plt
import numpy as np
import sys
import os
import matplotlib.cm as cm
import math 
import tomo_defs_190404 as tdefs


dirname = "../data/2019/3/27/"

read_file_offset = False
read_rf_waveform = True
multi_channel = True

show_mountain_range = False
intensity_calc = True
show_imshow = True

start_time = 0.12e-3 #0.12e-3
nturns_sort = 500

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


#plt.plot([t*1e3 for t in tdat_rf], sig_rf)
#plt.axvline(x=2.21434, color='k')
#plt.xlim(2.21, 2.22)
#plt.show()

time_interval = time_data[2] - time_data[1]
time_interval_ns = 1e9*time_interval

#idelay is delay in time bins between the real RF phase seen by the beam and that seen on the scope
if round(time_interval_ns) == 2:
	idelay = -47
else:
	idelay = int(round(-47*(2/time_interval_ns)))

print "time interval per bin (ns) ",time_interval_ns
print "delay (number of bins, time[ns]) ",idelay, idelay*time_interval_ns


tof_mean, time_indices, ttime_a = tdefs.calc_tof_rf_waveform(tdat_rf, sig_rf, nturns_sort+1, start_time)
istart = time_indices[0]	

nbins = int(round(tof_mean/time_interval)) #number of time bins

interpolate_data = True 
prof_data_a = tdefs.sort_data(istart, time_data, signal_data, ttime_a, time_indices, nturns_sort, nbins, interpolate_data, nshift=idelay) 

phase_ax = np.linspace(-180,180,nbins)

if show_mountain_range:
	tdefs.mountain_range(prof_data_a, nturns_sort, nbins, incl_xaxis = True, xaxis= phase_ax, xlabel='phase [deg]')#, extra_data = [mean_phase_data, phase_profmax])			

write_intensity_file = True
if intensity_calc:
	intensity, mean_phase_data, phase_fmax, std_prof_max = tdefs.intensity_calc_fn(prof_data_a, phase_ax, ttime_a, time_interval_ns, dirname, fname_sig, write_intensity_file = write_intensity_file)
	

find_phase_extremes = False
if find_phase_extremes:	
	ttime_ms = 1e3*ttime_a[:-1]

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

	print "time of minimum in mean, required phase jump (raw, smooth) ",t_mean_1, 20 - phase_fmax[index_min_mean], 20 - phase_fmax_movav[index_min_mean-imid] 
	print "time of maximum in mean, required phase jump (raw, smooth)",t_mean_2, 20 - phase_fmax[index_max_mean], 20 - phase_fmax_movav[index_max_mean-imid] 
	print "time of minimum in phase_fmax_movav, required phase jump (raw, smooth)", t_max_1
	print "time of maximum in phase_fmax_movav, required phase jump (raw, smooth)", t_max_2

if show_imshow:
	tlim1 = 1e3*ttime_a[0]
	tlim2 = 1e3*ttime_a[-1]
	extent = [-180, 180, tlim1, tlim2]
	
	if find_phase_extremes:
		extra_data = [mean_phase_data, phase_fmax, [t_movav, phase_fmax_movav]]
	else:
		extra_data = [mean_phase_data, phase_fmax]

	tdefs.imshow_plot(prof_data_a,extent,ylabel='time [ms]',xlabel='arbitrary phase [deg.]', extra_data = extra_data)
	
if find_phase_extremes:	
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
