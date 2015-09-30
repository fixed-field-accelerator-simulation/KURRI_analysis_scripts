from __future__ import division
import pylab as plt
import numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import quadrature
import sys
import os
import matplotlib.cm as cm
import analysedata as adata

#Calculate time to loss at various probe positions from bunch monitor data.
#30/9/2015 This version makes use of an updated algorithm signal_loss_time_new to find the loss time.

#dirname ="../data/SortByDate/2014/3/24/dmag/"	 
#dirname ="../../data/SortByDate/2014/3/25/D0_1012_D1_0_F7/" # Corrector at 0A
#dirname ="../../data/SortByDate/2014/3/25/D0_1012_D1_140_F1/" # Corrector at 140A
#dirname = "../../data/SortByDate/2014/3/27/Corr_550A/" #Corrector at 140A, 550A, 700A
#dirname = "../../data/SortByDate/2014/3/31/QinBin/F5/"
dirname = "../../data/2015_06_24/900A/"


indices_select = [0] #select particular indices from each channel, empty list means process all files
#channel_id = ["ch2","ch3"] #string to identify each channel (3/25 files)
#channel_id = ["CH2","CH3"] #string to identify each channel (3/27 files) 
channel_id = ["F12"] #identifier for channels to be analysed (3/31 files)

#select channel files
ch_files_sel = adata.select_filenames(dirname, indices_select, channel_id)



#set parameters for high pass filter
#cut_off_freq = 0.5e6 #0.5 MHz
sigtonoise = 5 #set threshold above noise to determine final loss point
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

clist = ['k','r','b','g','c']
 
#calculate loss time
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
			#print "chf ",chf
			#plt.plot(tdat_chf, data_chf)
			#plt.show()
			
			#loss_time = adata.signal_loss_time(tdat_chf, data_chf, method,  cut_off_freq = 0.5e6, sigtonoise = 6, show_data = True)
			
			loss_time = adata.signal_loss_time_new(tdat_chf, data_chf, sigtonoise = sigtonoise, show_data = False, wsize_half=50)
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

				
				show_scan = False
				if show_scan:
					#plt.subplot(211)
					plt.plot(sigtonoise_list, loss_time_scan, 'ko')
					plt.plot(range(1,11), poly(range(1,11)), 'r-')
					#plt.plot([index_trans], loss_time_scan[len(sigtonoise_list)-i], 'ro')
					#plt.subplot(212)
					plt.xlabel('signal to noise ratio')
					plt.ylabel(r'loss time ($\mu$s)')
					#plt.plot(sigtonoise_list[2:], diffchange, 'ko')
					plt.show()

				
				
			loss_time_err_ch.append(loss_time_error)
				#for t in loss_time_scan:
				#	plt.axvline(x=t)
			
			if plot_data:
				plt.plot(tdat_chf, np.array(data_chf)+index, color = clist[index])#, label=str(F7_pos[index]), )
				plt.axvline(x=loss_time[0], color = clist[index])
				plt.show()

				#plt.xlim(loss_time_scan[0], loss_time_scan[0)
			#plt.xlim(loss_time[0], loss_time[0]+100)
		
		if plot_data:
			plt.xlim(0.011,0.018)
			plt.xlabel('time (s)')
			plt.ylabel('bunch monitor signal (with offset)')
			plt.title('F7 probe signals')
			plt.legend()
			plt.savefig("F7_qinbin")
			plt.show()
			#sys.exit()

			
	
	loss_time_all.append(loss_time_ch)
	loss_time_err_all.append(loss_time_err_ch)

print "loss time all ",loss_time_all, "len ",len(loss_time_all)



if scan_threshold:
	print "loss time error ",loss_time_err_all

sys.exit()

#-------------------------Rest specific to file format used at particular time-----------------


if pos_from_filename:
	#fpos = [f.split('_')[-2] for f in ch_files_sel[0]] #assume same positions are true for other channels
	#fpos = [f.split('_')[-3] for f in ch_files_sel[0]]
	
	fsplit = [f.replace(".","_").split('_') for f in ch_files_sel[0]] #June 2015 data
	probe =[f[0] for f in fsplit]
	fpos = [f[1] for f in fsplit]
	
	#fpos = [re.split('._',f) for f in ch_files_sel[0]] #June 2015 data
	print "fpos ",fpos
	
	print "probe ",probe
	fout = 'QinBin_'+channel_id[0]+'.txt'

else:
	print "Warning: Using hardcoded probe positions!"
	F5_pos = [890, 870, 850, 830, 800, 750, 700, 650, 600, 550, 500, 450, 425, 400, 380, 360, 340]
	F1_pos = [870, 850, 800, 750, 700, 650, 600, 550, 500, 450, 425, 400, 380, 360, 340]
	F7_pos = [870, 850, 800, 750, 700, 650, 600, 550, 500, 450, 425, 400, 380, 360, 340, 320, 300]
  
	if "F1" in dirname:
		fout = 'QinBin_F1.txt'
		fpos = list(F1_pos)
	elif "F5" in dirname:
		fout = 'QinBin_F5.txt'
		fpos = list(F5_pos)
	elif "F7" in dirname:
		fout = 'QinBin_F7.txt'
		fpos = list(F7_pos)
		
	
if not write_average_channel:
	for loss_ch, err_ch in zip(loss_time_all, loss_time_err_all):
		if loss_ch != []:
			loss_time_tofile = loss_ch
			err_time_tofile = err_ch
			break
else:
	loss_time_tofile = [sum(c)/len(c) for c in np.transpose(loss_time_all)]
	err_time_tofile = [sum([c1**2 for c1 in c])**0.5 for c in np.transpose(loss_time_err_all)]
	print "loss_time_tofile ",loss_time_tofile, err_time_tofile
	

plt.errorbar([int(f) for f in fpos],loss_time_tofile, yerr=err_time_tofile, color='k',marker='o')
plt.show()


#outlier dealt with individually, 550A, 
#Corr_550A_Foil_77.0mm_ST6_-6.0A_F7_300mm_S7up_CH2.csv loss time 0.000226, (and for CH3 0.000226)

	
ff = open(fout,"w")
print >>ff, "probe position (mm) loss time (us) time error(us)"
for pos, losst, errt in zip(fpos, loss_time_tofile, err_time_tofile): 
	print pos, 1e6*losst, 1e6*errt
	#print >>ff, pos, 1e6*losst, 1e6*errt
	ff.write("%3s %5.1f %3.1f \n" % (pos , 1e6*losst, 1e6*errt))
ff.close()
