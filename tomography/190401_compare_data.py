from __future__ import division
import pylab as plt
import numpy as np
import sys
import os
sys.path.insert(0, '../KURRI_analysis_scripts/')
import analysedata_180919 as adata

time_interval_ns = 8 #4

#March/April 2019
prefix = 'tek00'
dir_27 = "../data/2019/3/27/"
dir_28 = "../data/2019/3/28/"
dir_29 = "../data/2019/3/29/"

#September 2018
#prefix = 'case'
#dir1 = "../data/2018/9/20/proc/"
#dir1_m = "../data/2018/9/20/machida-san_pattern/proc/"
#dir2 = "../data/2018/9/21/machida_san_foil_escape/proc/"
#dir3 = "../data/2018/9/25/machida-san_2018_09_25/proc/"
#dir4 = "../data/2018/9/27/JB_2018_09_27/proc/"

#dir_sel = [dir1, dir3]
#case_sel = [[1,3], [1,2,"3_50",4,5]] #dir3
#case_sel = [[1,3], [1,2,3,4,5]]

dir_sel = [dir_27, dir_28, dir_29]
case_sel  = [[9,12],[0, 2],[6]]

#dir_sel = [dir_29]
#case_sel = [[0,2,4,6]]

#dir_sel = [dir_27]
#case_sel = [[9,12,14,16]]


t_jump = 0.00221557818824
#case_sel = [[1,2,3,4,5]]
#case_sel = [[1,2,"3_50","3_44.5",4,5]]
#case_sel = [["2_43", "2_48","2_53", "2_58"]]

offset_09_20 = -0.038801293
offset_09_21 = 0.10296254081

#dir_sel = [dir2]
#case_sel = [[1,2,3,4,5]]
offset_l = [offset_09_20, offset_09_21]

offset_l = [0]*3
apply_offset = False


#time offset method
fix_offset = False
use_file_offset = True
usesigmax_file_offset = True
usesigmax_offset = False

read_phis_file = False

def file_read(dirname, fname, ncols = 2):
	"""Read phis file"""

	print "read data"
	ff = open(dirname+fname,'r')
	col1 = []
	col2 = []
	col3 = []
	col4 = []
	col5 = []
	col6 = []

	for lraw in ff:
		line = lraw.split()
		col1.append(float(line[0]))
		col2.append(float(line[1]))
		if ncols > 2:
			col3.append(float(line[2]))
		if ncols == 6:
			col4.append(float(line[3]))
			col5.append(float(line[4]))
			col6.append(float(line[5]))
	ff.close()

	if ncols == 2:
		return col1, col2
	elif ncols == 3:
		return col1, col2, col3
	elif ncols == 6:
		return col1, col2, col3, col4, col5, col6

	
versus_time = True

#label_l = ["20/9 linear","20/9 flattop @0.95ms ","21/9 case 1","21/9 case 2","21/9 case 3","21/9 case 4", "21/9 case 5"]

#label_l = ["21/9 case 1","21/9 case 2","21/9 case 3","21/9 case 4", "21/9 case 5"]
#label_l = ["25/9 case 1","25/9 case 2","25/9 case 3 (fp=50)","25/9 case 3 (fp=44.5)","25/9 case 4", "25/9 case 5"]

#label_l = ["27/9 fp = 43","27/9 fp = 48","27/9 fp = 53","27/9 fp = 58"]

label_l = ['27th (0deg)','27th (20deg)','28th (0deg)', '28th (-20deg)','29th (0deg)']


#label_l = ['-15deg','-10deg','-5deg','0deg']
#label_l = ['0deg','20deg','40deg','60deg']
col_l = ['k','r','b','c','m','g','y','k']
compare_intensity = False
compare_fwhm = False
compare_phases = True
i = 0
if compare_intensity or compare_fwhm or compare_phases:
	for dir1, casel, offset_set in zip(dir_sel, case_sel, offset_l):
		for case in casel:
			print "dir1 ",dir1, casel
			if prefix != 'case':
				if case < 10:
					fname_int = prefix+"0"+str(case)+"ALL_intensity.txt"
				else:
					fname_int = prefix+str(case)+"ALL_intensity.txt"

				print "hhe"
				time_l, int_l, fwhm_l, phmean_l, phstd_l, phmax_l = file_read(dir1, fname_int, ncols = 6)
			else:
				fname_int = prefix+str(case)+"_intensity.txt"
				time_l, int_l, fwhm_l = file_read(dir1, fname_int, ncols = 3)
			
			print "len phmean_l ",len(phmean_l)
			#time_l, int_l = file_read(dir1, fname_int)
			if i == 0:
				phmean_save = phmean_l

			time_a = np.array(time_l)
			time_a_ms = 1e3*time_a
			
			if read_phis_file:
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
					fwhm_a = np.array(fwhm_l)
					plt.plot(time_a_ms, time_interval_ns*fwhm_a, col_l[i]+'.', label=label_l[i])
				elif compare_phases:
					plt.plot(time_a_ms, phmean_l, col_l[i]+'.', label=label_l[i])
					#plt.plot(time_a_ms, phstd_l, col_l[i]+'.', label=label_l[i])

				if read_phis_file:
					if len(t_phis_file) > 1:
						it = 0
						for t in t_phis_file[1:]:
							if it == len(t_phis_file) - 2:
								style = '--'
							else:
								style = '-'

							plt.axvline(x = t, color = col_l[i], linestyle=style)
							it = it + 1
			else:
				plt.plot(int_l, col_l[i]+'.', label=label_l[i])	
				
			i = i+1

	#plt.axvline(x=9e-5, color='r')
	#plt.axvline(x=0.0626, color='r')
	if versus_time:
		
		#plt.xlim(-0.05,xmax=1.5)
		if apply_offset:
			plt.xlabel('time with offset [ms]')
		else:
			plt.xlabel('time from scope trigger [ms]')
	else:
		plt.xlim(0,xmax=750)
		plt.xlabel('turn ')


	if compare_intensity:
		plt.ylim(ymin=0)
		plt.ylabel('integrated turn-by-turn signal [nVs]')
	elif compare_fwhm:
		plt.ylim(ymin=0)	
		plt.ylabel('FWHM [ns]')
	elif compare_phases:
		#plt.xlim(-0.2+t_jump*1e3, 0.2+t_jump*1e3)
		plt.axvline(x=t_jump*1e3)
		plt.ylabel('phase [deg]')

	plt.legend(prop={'size': 8})

 	if compare_intensity:
		plt.savefig('compare_intensities')
	elif compare_fwhm:
		plt.savefig('compare_fwhm')
	elif compare_phases:
		plt.savefig('compare_phases')
	plt.show()

	ph_mean_diff = [p1-p2 for p1,p2 in zip(phmean_l,phmean_save)]
	print "max, min ph_mean_diff ",max(ph_mean_diff[1000:]), min(ph_mean_diff[1000:])
	plt.plot(time_a_ms, ph_mean_diff, 'k-')
	plt.axvline(x=t_jump*1e3,color='gray',linewidth=1)
	plt.savefig('phasediff.png')
	plt.show()


#march/april specific
ph_diff_crude = [-20, 22.06659285, 41.929237987, 62.6337568622, 74.6232637078, 87.9399546656, 90.1530632616, 77.177506003, 90.4613982657]
ph_req = [-20.77, 20, 40, 60, 80, 100, 120, 140, 160]

plt.plot([0,160],[0,160],'--', color='gray')
plt.plot(ph_req, ph_diff_crude,'ko-')
plt.xlabel('phase jump applied [deg]')
plt.ylabel('shift in mean phase [deg]')
plt.savefig('mean_phase_jump_vs_pred')
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

