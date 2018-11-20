from __future__ import division
import pylab as plt
import numpy as np
import sys
import os
sys.path.insert(0, '../KURRI_analysis_scripts/')
import analysedata_180919 as adata

time_interval_ns = 4

dir1 = "../data/2018/9/20/proc/"
dir2 = "../data/2018/9/21/machida_san_foil_escape/proc/"
dir3 = "../data/2018/9/25/machida-san_2018_09_25/proc/"


#dir_sel = [dir1, dir3]
#case_sel = [[1,3], [1,2,"3_50",4,5]] #dir3
#case_sel = [[1,3], [1,2,3,4,5]]

dir_sel = [dir2]
case_sel = [[1,2,3]]

offset_09_20 = -0.038801293
offset_09_21 = 0.10296254081

#dir_sel = [dir2]
#case_sel = [[1,2,3,4,5]]
offset_l = [offset_09_20, offset_09_21]
apply_offset = True


#time offset method
fix_offset = False
use_file_offset = True
usesigmax_file_offset = True
usesigmax_offset = False

def file_read(dirname, fname, ncols = 2):
	"""Read phis file"""

	ff = open(dirname+fname,'r')
	col1 = []
	col2 = []
	if ncols > 2:
		col3 = []

	for lraw in ff:
		line = lraw.split()
		col1.append(float(line[0]))
		col2.append(float(line[1]))
		if ncols > 2:
			col3.append(float(line[2]))
	ff.close()

	if ncols == 2:
		return col1, col2
	elif ncols == 3:
		return col1, col2, col3

versus_time = True

label_l = ["20/9 linear","20/9 flattop @0.95ms ","21/9 case 1","21/9 case 2","21/9 case 3","21/9 case 4", "21/9 case 5"]

label_l = ["21/9 case 1","21/9 case 2","21/9 case 3","21/9 case 4", "21/9 case 5"]
col_l = ['k','r','b','c','m','g','y','k']
compare_intensity = False
compare_fwhm = True
i = 0
if compare_intensity or compare_fwhm:
	for dir1, casel, offset_set in zip(dir_sel, case_sel, offset_l):
		for case in casel:
						
			fname_int = "case"+str(case)+"_intensity.txt"
			time_l, int_l, fwhm_l = file_read(dir1, fname_int, ncols = 3)
			time_a = np.array(time_l)
			time_a_ms = 1e3*time_a

			fwhm_a = np.array(fwhm_l)
			
			case_spl = str(case).split("_")
			case_num = case_spl[0]
			fname_phis = "case"+str(case_num)+"_phis.txt"
			t_phis_file, phis_file = file_read(dir1, fname_phis)
			
			if use_file_offset:
				try:
					fo = open(dir1+"case"+str(case)+"_toffset.txt",'r')
					if usesigmax_file_offset:
						fo.readline()
					offset = 1e3*float(fo.readline())
				except:
					offset = offset_set
			elif fix_offset:
				offset = offset_set
			elif usesigmax_offset:
				indmax = int_l.index(max(int_l))
				offset_indmax  = time_l[indmax]			
				offset = offset_indmax - 0.05

			
			#time_a = np.array(time_l)
			
			if apply_offset:
				time_a_ms = time_a_ms - offset

			if versus_time:
				if compare_intensity:
					plt.plot(time_a_ms, int_l, col_l[i]+'.', label=label_l[i])
				elif compare_fwhm:
					plt.plot(time_a_ms, time_interval_ns*fwhm_a, col_l[i]+'.', label=label_l[i])

				if len(t_phis_file) > 1:
					for t in t_phis_file[1:]:
						plt.axvline(x = t, color = col_l[i])
			else:
				plt.plot(int_l, col_l[i]+'.', label=label_l[i])	
				
			i = i+1

	if versus_time:
		plt.ylim(ymin=0)
		plt.xlim(-0.05,xmax=1.0)
		if apply_offset:
			plt.xlabel('time with offset [ms]')
		else:
			plt.xlabel('time from scope trigger [ms]')
	else:
		plt.xlim(0,xmax=750)
		plt.xlabel('turn ')
	if compare_intensity:
		plt.ylabel('integrated turn-by-turn signal [nVs]')
	else:
		plt.ylabel('FWHM [ns]')
	plt.legend(prop={'size': 10})
	plt.savefig('compare_intensities')
	plt.show()

	sys.exit()

show_freq_at_maxintensity = False
if show_freq_at_maxintensity:
	#frequency at intensity maximum calculated from RF waveform
	f_ideal = 1.582769
	f_ideal_50us = 1.5892538946474106 
	#20/9 data
	freq_at_max_2009 = [1.59054812935, 1.59079946576]
	#21/9 data
	freq_at_max_2109 = [1.58621542229, 1.58658480398, 1.58607004716, 1.58354295928, 1.58452670016]
	#25/9 data 
	freq_at_max_2509 = [1.58351755625, 1.59024490293, 1.58806419026,  1.58925732567,  1.5855347887, 1.58998865587]

	plt.plot(freq_at_max_2009, 'mo',label='20/09')
	plt.plot(freq_at_max_2109, 'ko', label='21/09')
	plt.plot(freq_at_max_2509, 'ro', label='25/09')
	plt.axhline(y=f_ideal, color='gray',linestyle='--')
	plt.axhline(y=f_ideal_50us, color='gray',linestyle='--')	
	plt.ylabel('RF frequency at intensity max [MHz]')
	plt.legend()
	plt.savefig('freq_at_peak')
	plt.show()
	sys.exit()

compare_rf = False
i = 0
if compare_rf:
	for dir1, casel in zip(dir_sel, case_sel):
		for case in casel:
			fname_rf = "case"+str(case)+"_rf.csv"
			tdat_rf, sig_rf = adata.read_scope(dir1, fname_rf)

			tdat_rf_ms = 1e3*np.array(tdat_rf)
	
			plt.plot(tdat_rf_ms, sig_rf, col_l[i]+'-', label='case'+str(case))
			
			i = i+1

	#plt.ylim(ymin=0)
	
	plt.xlim(tdat_rf[0] + 0.4, tdat_rf[0]+0.6)
	plt.xlabel('time [s]')
	plt.ylabel('rf')
	plt.legend()
	plt.show()

