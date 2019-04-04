from __future__ import division
import pylab as plt
import numpy as np
import sys
import os
import matplotlib.cm as cm
import math 
import tomo_defs_190403 as tdefs

dirname = "../data/2019/4/20190401/"

read_file_offset = False
read_rf_waveform = True
multi_channel = True

show_mountain_range = False
intensity_calc = True
show_imshow = True

start_time = 0.12e-3
nturns_sort = 500

#set idelay to -9 if 10ns per step and to -47 if 2n

#read the data
interactive_data_sel = True
if interactive_data_sel:
	fname_sig, time_data, signal_data, tdat_rf, sig_rf, t_phis_file, phis_file, offset_file_ms = tdefs.read_data_files(dirname, read_rf_waveform = read_rf_waveform, read_file_offset = read_file_offset, multi_channel = multi_channel)

time_interval = time_data[2] - time_data[1]
time_interval_ns = 1e9*time_interval

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
	intensity, mean_phase_data, phase_profmax, std_prof_max = tdefs.intensity_calc_fn(prof_data_a, phase_ax, ttime_a, time_interval_ns, dirname, fname_sig, t_phis_file = t_phis_file, write_intensity_file = write_intensity_file)
	
	
if show_imshow:
	tlim1 = 1e3*ttime_a[0]
	tlim2 = 1e3*ttime_a[-1]
	extent = [-180, 180, tlim1, tlim2]
	tdefs.imshow_plot(prof_data_a,extent,ylabel='time [ms]',xlabel='arbitrary phase [deg.]', extra_data = [mean_phase_data, phase_profmax])


ttime_ms = 1e3*ttime_a[:-1]
#plt.subplot(211)
#plt.plot(ttime_ms, intensity)
#plt.subplot(212)
min_mean = min(mean_phase_data)
index_min_mean = mean_phase_data.index(min_mean)
max_mean = max(mean_phase_data[index_min_mean:])

time_min_mean = ttime_ms[index_min_mean]

index_max_mean = mean_phase_data.index(max_mean)
time_max_mean = ttime_ms[index_max_mean]

#print "phase, mean and max, at this time ",mean_phase_data[index_min_mean], phase_profmax[index_min_mean] 

nw = 21
imid = int(0.5*(nw-1))
phase_profmax_movav = tdefs.running_mean(phase_profmax,nw)
t_movav = ttime_ms[imid:-imid]
print "time of minimum in mean, required phase jump (raw, smooth) ",time_min_mean, 20 - phase_profmax[index_min_mean], 20 - phase_profmax_movav[index_min_mean-imid] 
print "time of maximum in mean, required phase jump (raw, smooth)",time_max_mean, 20 - phase_profmax[index_max_mean], 20 - phase_profmax_movav[index_max_mean-imid] 

 
plt.plot(ttime_ms, mean_phase_data, 'k',label='<phase>')
plt.plot(ttime_ms, phase_profmax, 'r',label='phase@dist_max')
plt.plot(t_movav, phase_profmax_movav, 'b--', linewidth=2)
plt.axvline(x=time_min_mean, color='gray',linestyle='--')
plt.axvline(x=time_max_mean, color='gray',linestyle='--')
plt.xlabel('time [ms]')
plt.ylabel('phase [deg]')
plt.legend()
plt.savefig('phase_mean_max')
plt.show()
