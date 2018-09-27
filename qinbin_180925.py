from __future__ import division
import pylab as plt
import numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import quadrature
import sys
import os
import matplotlib.cm as cm
import analysedata_150619 as adata

#Calculate time to loss at various probe positions from bunch monitor data.
#30/9/2015 This version makes use of an updated algorithm signal_loss_time_new to find the loss time.
#7/8/2018 Option to process data from multiple probes.

#probeoffset=4300.0 
time_cut = 130e-3 #cut data after time_cut


dirname = "../../data/2015/6/24/700A/"
#dirname = "../../data/2018/9/21/QinBin/"

#The offset (in mm) is added to each probe position
if "2014" in dirname or "2015" in dirname:
	probe_offset_map =  {'F01': 4300,'F02':4300, 'F03':4300,'F05':4300,"F07":4300,'F12': 4300}
else:
	probe_offset_map =  {'F1': 4705.5, 'F5': 4708}
	

#channel_id = ["F01","F02","F03","F05","F12"] #identifier for probes to be analysed 
channel_id = ["F12"] #identifier for probes to be analysed 
#channel_id = ["F1","F5"] #identifier for probes to be analysed 

#indices_select used to select particular data files associated with each probe
indices_select = [-1] #empty list to process all files
#indices_select = range(1,5) #select files with index 1 to 4 for each probe
#indices_select = [0,1,2] #select first three files for each probe

#select channel files
ch_files_sel = adata.select_filenames(dirname, indices_select, channel_id)
print "ch_files_sel ",ch_files_sel
print "file count ",[len(f) for f in ch_files_sel]

#set parameters for high pass filter
#cut_off_freq = 0.5e6 #0.5 MHz
sigtonoise = 4 #set threshold above noise to determine final loss point
#sigtonoise reduced from 5 to 3 to cope with noise in F12 seen in June 2015 data

method = 'noise' #choose "fixed" or "noise"

#31/3/14 - map each index to probe position		
#F5_pos = [890, 870, 850, 830, 800, 750, 700, 650, 600, 550, 500, 450, 425, 400, 380, 360, 340]
#F1_pos = [870, 850, 800, 750, 700, 650, 600, 550, 500, 450, 425, 400, 380, 360, 340]
#F7_pos = [870, 850, 800, 750, 700, 650, 600, 550, 500, 450, 425, 400, 380, 360, 340, 320, 300]

scan_threshold = False
plot_data = True


write_average_channel = False #for cases where more than one channel, just write the mean at each probe position to file.
pos_from_filename = True

probe_name_l = []
probe_radius_l = []
if pos_from_filename:
	for files in ch_files_sel:
		fsplit = [f.replace(".","_").split('_') for f in files] #June 2015 data
				
		probe =[f[-3] for f in fsplit]
		fpos = np.array([f[-2] for f in fsplit])


		try:
			probe_offset_a = [probe_offset_map[p] for p in probe]
		except:
			pass
		
		probe_radius = 1e-3*(fpos.astype(int) + probe_offset_a)
		
		probe_name_l.append(probe)
		probe_radius_l.append(probe_radius)
else:
	print "please provide probe radius values manually"
	print "in the format"
	print "probe_radius_l = [[probe1_pos1, probe1_pos2, probe1_pos3],[probe2_pos1,...]...]]"
	
	probe_name_l = list(channel_id)
	sys.exit()
		
print "probe_name_l ",probe_name_l
print "probe_radius_l ",probe_radius_l


clist = ['r','b','g','c']
 
#calculation loop here
loss_time_all = []
loss_time_err_all = []
for chf in ch_files_sel:
	loss_time_ch  =[]
	loss_time_err_ch = []
	if chf != []:
		for chf1 in chf:
			index = chf.index(chf1)
			#read scope data
			tdat_chf, data_chf = adata.read_scope(dirname, chf1)
			
			icut = -1
			if time_cut != None:
				it = 0
				for t in tdat_chf:
					if t > time_cut:
						icut = it
						break
					it = it + 1
					
			
			print "process data file ",chf1
			print "cut at ",tdat_chf[icut]

			 
			#plt.plot(tdat_chf, data_chf)
			#plt.show()
			

			loss_time = adata.signal_loss_time_new(tdat_chf[:icut], data_chf[:icut], sigtonoise = sigtonoise, show_data = False, wsize_half=50)
			loss_time_ch.append(loss_time[0])
			loss_time_error = loss_time[2]
			print "loss_time, error ",loss_time[0], loss_time_error
			
			if scan_threshold:
				loss_time_scan = []
				sigtonoise_list = np.array(range(6,11))
				for sn in sigtonoise_list:
					
					loss_time1 = adata.signal_loss_time(tdat_chf, data_chf, method, sigtonoise = sn)
					loss_time_scan.append(1e6*loss_time1[0])
					#print "sig/noise, loss time ",sn, loss_time
					
				out = np.polyfit(sigtonoise_list, loss_time_scan, 1, full=True)  	
				poly = np.poly1d(out[0])
				loss_time_error = poly(1)-poly(sig2noise)
				#print "loss_time_scan ",loss_time_scan
				print "loss_time_error ",loss_time_error

			loss_time_err_ch.append(loss_time_error)

			if plot_data:
				plt.plot(tdat_chf[:icut], np.array(data_chf[:icut])+index, 'k')#, label=str(F7_pos[index]), )
				plt.axvline(x=loss_time[0])
				plt.show()
				
				
			
	
	loss_time_all.append(loss_time_ch)
	loss_time_err_all.append(loss_time_err_ch)

print "loss time all ",loss_time_all, "len ",len(loss_time_all)

if scan_threshold:
	print "loss time error ",loss_time_err_all


def write_data(fout,loss_time_tofile, err_time_tofile, probe_radius):
	ff = open(fout,"w")
	print >>ff, "probe position (mm) loss time (us) time window(us)"
	for pos, losst, errt in zip(probe_radius, loss_time_tofile, err_time_tofile): 
		print pos, 1e6*losst, 1e6*errt
		#print >>ff, pos, 1e6*losst, 1e6*errt
		ff.write("%3s %5.1f %3.1f \n" % (pos , 1e6*losst, 1e6*errt))
	ff.close()
		

if not write_average_channel:
	for loss_ch, err_ch, id, probe_radius in zip(loss_time_all, loss_time_err_all, probe_name_l, probe_radius_l):
		print "loss_ch,id  ",loss_ch, id
		if loss_ch != []:
			
			fout = 'QinBin_'+id[0]+'.txt'
		
			loss_time_tofile = loss_ch
			err_time_tofile = err_ch
			
			write_data(fout, loss_time_tofile, err_time_tofile, probe_radius)
			
		plt.errorbar(loss_ch, probe_radius, xerr=err_ch, marker='o')
	
else:
	loss_time_tofile = [sum(c)/len(c) for c in np.transpose(loss_time_all)]
	err_time_tofile = [sum([c1**2 for c1 in c])**0.5 for c in np.transpose(loss_time_err_all)]
	

plt.xlabel('time of loss')
plt.ylabel('radius')
plt.savefig('qinbin_result')
plt.show()

