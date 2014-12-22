from __future__ import division
import numpy as np
import sys
import analysedata as adata

#Find last signal signal recorded during Qin Bin scan. 
#Set 'dirname' to directory containing files to be analysed.
#Set 'indices_select' to choose files by index (indices_select=[] selects all files)
#'channel_id' can be used to select particular channels 
#
#'cut_off_freq' determines the cut off point in the high pass filter
#Set 'method' to 'fixed' or 'noise' to determine the signal threshold in two ways. In the 'fixed' case
#         a fixed voltage threshold is used (threshold_volts in signal_loss_time) while in the 'noise' case
#         the signal is assumed to be some multiple of the rms of the noise calculated at the end of the data 

#dirname ="../data/SortByDate/2014/3/24/dmag/"	 
#dirname ="../data/SortByDate/2014/3/25/D0_1012_D1_140_F7/" 
#dirname = "../data/SortByDate/2014/3/27/Corr_700A/"
dirname = "../data/SortByDate/2014/3/31/QinBin/F7/"

#
indices_select = [] #select particular indices from each channel, empty list means process all files
#channel_id = ["CH1","CH2","CH3"] #string to identify each channel
channel_id = ["F1","F5","F7"] #identifier for channels to be analysed

#select channel files
ch_files_sel = adata.select_filenames(dirname, indices_select, channel_id)

print "ch_files_sel ",ch_files_sel

#set parameters for high pass filter
cut_off_freq = 0.5e6 #0.5 MHz
nstd = 10 #set threshold above noise to determine final loss point

method = 'noise' #choose "fixed" or "noise"

#calculate loss time
loss_time_all = []
for chf in ch_files_sel:	
	loss_time_ch  =[]
	if chf != []:
		for chf1 in chf:
			#read scope data
			tdat_chf, data_chf = adata.read_scope(dirname, chf1)
			loss_time = adata.signal_loss_time(tdat_chf, data_chf, method)
			loss_time_ch.append(loss_time)
	loss_time_all.append(loss_time_ch)

print "loss time all ",loss_time_all, "len ",len(loss_time_all)
	

#-------------------------Rest specific to file format used at particular time-----------------
#31/3/14 - map each index to probe position			
F5_pos = [890, 870, 850, 830, 800, 750, 700, 650, 600, 550, 500, 450, 425, 400, 380, 360, 340]
F1_pos = [870, 850, 800, 750, 700, 650, 600, 550, 500, 450, 425, 400, 380, 360, 340]
F7_pos = [870, 850, 800, 750, 700, 650, 600, 550, 500, 450, 425, 400, 380, 360, 340, 320, 300]


loss_time_flatten = loss_time_all[0]

if "F1" in dirname:
	fout = 'QinBin_F1.txt'
	fpos = list(F1_pos)
elif "F5" in dirname:
	fout = 'QinBin_F5.txt'
	fpos = list(F5_pos)
elif "F7" in dirname:
	fout = 'QinBin_F7.txt'
	fpos = list(F7_pos)
	
ff = open(fout,"w")
print >>ff, "probe position (mm) loss time (us)"
for pos, losst in zip(fpos, loss_time_flatten): 
	print >>ff, pos, 1e6*losst
ff.close()
