from __future__ import division
import numpy as np
import sys
import analysedata as adata
import pylab as plt
import time
import math
from scipy import optimize
from scipy.special import erf
from scipy.integrate import simps

#8/4/2016 copy beamsize_150622_forpaper.py. Try to deal with Abel transform.

#Extra settings from Suzie
plt.rcParams['font.family'] = 'serif'
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex=True)

#plt.rcParams.update({'axes.labelsize':25})
#plt.rcParams.update({'xtick.labelsize':14})
#plt.rcParams.update({'ytick.labelsize':14})

#fig=plt.figure(num=1, figsize=(1.5*6, 1.5*4), dpi=150, facecolor='w', edgecolor="k")
fig=plt.figure(num=1, figsize=(6, 4), dpi=150, facecolor='w', edgecolor="k")
#ax = fig.add_subplot(111)

#plt.rcParams.update({'axes.labelsize':12})
#plt.rcParams.update({'xtick.labelsize':14})
#plt.rcParams.update({'ytick.labelsize':14})

#Find last signal signal recorded during Qin Bin scan. 
#Set 'dirname' to directory containing files to be analysed.
#Set 'indices_select' to choose files by index (indices_select=[] selects all files)
#'channel_id' can be used to select particular channels 
#
#'cut_off_freq' determines the cut off point in the high pass filter
#Set 'method' to 'fixed' or 'noise' to determine the signal threshold in two ways. In the 'fixed' case
#         a fixed voltage threshold is used (threshold_volts in signal_loss_time) while in the 'noise' case
#         the signal is assumed to be some multiple of the rms of the noise calculated at the end of the data 


# Define function for calculating a power law
powerlaw = lambda x, amp, index: amp * (x**index)
chidist = lambda x, sigma: x/(sigma**2) * np.exp(-x**2/(2*(sigma**2)))

chidistsc = lambda x, sigma, sc: sc*(x)/(sigma**2) * np.exp(-x**2/(2*(sigma**2)))
chidistmod = lambda x, sigma, sc, mod: sc*x/(sigma**2) * np.exp(-x**2/(2*(sigma**2))) + mod*(x**2)

chiexpand = lambda x, s, a: a*(x/s**2) - a*(x**3/(2*s**4))
chiexpandm = lambda x, s, a, b: a*(x/s**2) - a*(x**3/(2*s**4)) + b*(x**2)

fitfn = lambda x, a, b: a*x + b*(x**3)
fitfn3 = lambda x, a, b, c: a*x + b*(x**2) + c*(x**3)


#fgauss = lambda r, sigma:  np.exp(-r**2/(2*(sigma**2)))
#r is outer limit of integration, a is point at which Abel is being evaluated
#fgauss_abel_int = lambda r, sigma, a: np.sqrt(2*math.pi)*sigma*np.exp(-a**2/(2*(sigma**2)))*erf(np.sqrt(r**2 - a**2)/(np.sqrt(2)*sigma))



	

	

def te_to_mom(mass, te):
	"Convert total energy to momentum"
	return math.sqrt(  te**2 - mass**2)

def ke_to_te(mass, ke):
	"Convert kinetic energy to total energy"
	return mass + ke
	
def ke_to_mom(mass, ke):
	"Convert kinetic energy to momentum"
	return te_to_mom(mass, ke_to_te(mass, ke))
	
fn = lambda t, sigma_t, mean_t: (t-mean_t)*np.exp(-(t-mean_t)**2/(2*sigma_t**2))

def find_peaks(xdat, ydat, interval, find_hminus, x_low = None, x_high = None):
    
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
        index_high = len(xdat) - 1
    
    npts = index_high - index_low + 1
    

    for i in range(int(npts/interval) + 10):
        
        my = 0.
        my_min = 100
        
        winsize = int(0.5*interval)
        

            
        if i == 0 :
            istart = index_low
            xwindow = xdat[index_low:index_low+2*winsize]
            ywindow = ydat[index_low:index_low+2*winsize]
        elif i == 1:
            istart = index_low + 2*winsize
            xwindow = xdat[index_low+2*winsize:index_low+4*winsize]
            ywindow = ydat[index_low+2*winsize:index_low+4*winsize]       
        else:
            icen = mi + interval

            #test if in range
            if icen + winsize > len(xdat):
                print "reach end of data ",icen+winsize, len(xdat)
                break
            
            xwindow = xdat[icen-winsize: icen+winsize]
            ywindow = ydat[icen-winsize: icen+winsize]
            

            istart = icen - winsize
        

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
        
        
                
        if i < 0:            
            plt.plot(xwindow,ywindow,'k.-')
            plt.plot([mx],[my],'ro')
            plt.show()
            
            


            
        pmxi.append(i)
        pmx.append(mx)
        pmy.append(my)
        pmi.append(mi)
        pmy_av.append(my_av)
        pmy_std.append(my_std)
        
        #if i > 0:
        #    xgap = pmx[-1] - pmx[-2]
        #    print "xgap ",xgap, xgap/0.008
        #    
        #    lastpeaki = mi

            
        
    return pmxi, pmx, pmy, pmi, pmy_av, pmy_std

def powerlaw_fit(xdata, ydata):
	# Power-law fitting is best done by first converting
	# to a linear equation and then fitting to a straight line.
	#  y = a * x^b
	#  log(y) = log(a) + b*log(x)
	from scipy import log10
	from scipy import optimize
	

	
	
	yerr = np.array([1]*len(ydata)) #no error data
	
	logx = log10(xdata)
	logy = log10(ydata)
	logyerr = yerr / ydata
	
	# define our (line) fitting function
	fitfunc = lambda p, x: p[0] + p[1] * x
	errfunc = lambda p, x, y, err: (y - fitfunc(p, x)) / err
	
	
	pinit = [1.0, -1.0]
	out = optimize.leastsq(errfunc, pinit, args=(logx, logy, logyerr), full_output=1)
	
	pfinal = out[0]
	covar = out[1]
	
	
	#print pfinal
	#print "covariance", covar
	
	#y = amp * x^exponent
	
	exponent = pfinal[1] #index in original
	amp = 10.0**pfinal[0]
	#print "index ",index
	#print "amp ", amp
	
	exponentErr = np.sqrt( covar[0][0] )
	ampErr = np.sqrt( covar[1][1] ) * amp
	
	#print "best fit power, error ",exponent, exponentErr
	#print "best fit factor ",pfinal[0],amp
	#print "prediction ", powerlaw(xdata, amp, exponent)
	
	chisq = np.sum(((ydata - powerlaw(xdata, amp, exponent))**2),axis=0)
	#print "chisq fit ",chisq

	
	print "amp, exponent ",amp, exponent
	
	plot_fit = True
	if plot_fit:
		##########
		# Plotting data
		##########
	
		plt.clf()
		plt.plot(xdata, powerlaw(xdata, amp, exponent),'k')     # Fit
		plt.plot(xdata, 1e-5*chidist(xdata, 0.002),'r') #theoretical exponent  
		plt.plot(xdata, ydata,'k.') 
		#plt.text(min(xdata), max(ydata), 'Ampli = %5.2f +/- %5.2f' % (amp, ampErr))
		plt.text(min(xdata), max(ydata), 'Exponent = %5.2f +/- %5.2f' % (exponent, exponentErr))
		plt.xlabel('x')
		plt.ylabel('y')
		plt.title('Best Fit Power Law')
		plt.show()

		
		
	return exponent, amp, chisq
	
	
	
#dirname ="../data/SortByDate/2014/3/24/dmag/"	 
#dirname ="../data/SortByDate/2014/3/25/D0_1012_D1_140_F7/" 
#dirname = "../data/SortByDate/2014/3/27/Corr_700A/"
dirname = "../data/2014/3/31/QinBin/F7/"
#dirname = "/mnt/data1/KURRIFFAG/2014-03-31/QinBin/F7/"


#-------------------------Probe position is specific to each experiment-----------------
#31/3/14 - map each index to probe position			
F5_pos = [890, 870, 850, 830, 800, 750, 700, 650, 600, 550, 500, 450, 425, 400, 380, 360, 340]
F1_pos = [870, 850, 800, 750, 700, 650, 600, 550, 500, 450, 425, 400, 380, 360, 340]
F7_pos = [870, 850, 800, 750, 700, 650, 600, 550, 500, 450, 425, 400, 380, 360, 340, 320, 300]

if "F1" in dirname:
	fout = 'beamsize_F1.txt'
	fpos = list(F1_pos)
elif "F5" in dirname:
	fout = 'beamsize_F5.txt'
	fpos = list(F5_pos)
elif "F7" in dirname:
	fout = 'beamsize_F7.txt'
	fpos = list(F7_pos)

	
#indices_select = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13] #select particular indices from each channel, empty list means process all files
#indices_select = [-1] 
indices_select = [] 
#channel_id = ["CH1","CH2","CH3"] #string to identify each channel
channel_id = ["F7"] #identifier for channels to be analysed


for cid in channel_id:
	files = adata.get_files(dirname, cid)
	if files != []:
		break
print "files ",files


#select channel files
#ch_files_sel = adata.select_filenames(dirname, indices_select, channel_id)
#print "ch_files_sel ",ch_files_sel
#sys.exit()

if indices_select == []:
	ch_files_sel = files
elif len(indices_select) == 1 and indices_select[0] < 0:
	ch_files_sel = files[:indices_select[0]]
else:
	ch_files_sel = [files[i] for i in indices_select]
print "ch_files_sel ",ch_files_sel




#set parameters for high pass filter
cut_off_freq = 0.5e6 #0.5 MHz
nstd = 10 #set threshold above noise to determine final loss point

method = 'noise' #choose "fixed" or "noise"

#read file to map time to momentum
rfdir="../march14_exp/rf_pattern/"
rffile="SimulatedVariableK.dat"
rfdat = np.loadtxt(rfdir+rffile, skiprows=2, usecols=(0,1,2), unpack=True)
rfdat_t=np.array(rfdat[0])*1E6
rfdat_E=np.array(rfdat[1])*1E-6
rfdat_F=np.array(rfdat[2]) #frequency in hertz
mass=938.272
rfdat_p=map(lambda x: ke_to_mom(mass, x), rfdat_E)
toffset= -93.3
#make a list of momentum values corresponding to a given set of time points
fittedpfunc=adata.findPfromt(rfdat_t,rfdat_p, order=2) #using RF data 

#plt.plot(rfdat_t, rfdat_p)
#plt.plot(rfdat_t, fittedpfunc(rfdat_t),'r-')
#plt.show()


    
def movingaverage(x, data, win):
    """calculate moving window average of 1D list"""
    ic = win
    l = 2*win+1
    data_sm = [sum(data[ic - win: ic + win+1])/l for ic in range(win, len(data) - win)]
    x_sm = [x[ic] for ic in range(win, len(data)-win)]
    
    return x_sm, data_sm
    
#win_av = 1000
win_peak_av = 23

#inactive if second_smoothing is False
second_smoothing = False
win_grad_av = 5#20


#calculate loss time
loss_time_all = []
zerograd_1_time_all = []
zerograd_2_time_all = []

duration_fwhm_all = []
duration_gradmax_zerograd_all = []
duration_gradmax_est = []

grad_max_all = []
duration_full = []
duration_full_offset = []
loss_time_full = []
loss_time_full_offset = []

distance_fwhm_mm = []
distance_full_mm = []
distance_full_offset_mm = []

distance_rms_abel_mm = []
duration_rms_abel = []
Faxis_rms_abel = []




r_check = []
distrib_check = []

p_full = [] #momentum where signal is lost
p_full_offset = []

peak_av_initial = [] #record signal at start of signal loss

show_gradient = True #Generate profile for paper

show_all_data = False

show_all_profiles = False

plot_rms_window = False
if plot_rms_window:
	win_grad_av_list = [1,2,3,4,5,6,7,8,9,10,15,20,30]
	rms_abel_l = [0.00248, 0.002497, 0.002526, 0.002526, 0.002534, 0.002543, 0.002615,0.002624,  0.002632, 0.00264, 0.002647,0.002659,0.002671]
	plt.plot(win_grad_av_list, rms_abel_l, 'ko-')
	plt.show()


fpos_out = []
#vary lower lim for each case
#lower_lim_file = [300]*13 + [100]*2

print "len file for each channel ",[len(c) for c in ch_files_sel]

#lower_lim_file = [300]*(lenf-1) + [200]
lower_lim_file = [300]*len(ch_files_sel)
print "lower_lim_file ",lower_lim_file

#assumed reference momentum, dispersion, scaling index used in translation from time to distance
pref = 147
dispref = 0.594
kref = 7.6
k2 = 1/(1+kref)

#colors = plt.cm.Reds(np.linspace(1.0, 0.0, len(ch_files_sel)))
colors = plt.cm.hot(np.linspace(0,0.7,len(ch_files_sel)))
#colors = plt.cm.coolwarm(np.linspace(0,1,len(ch_files_sel)))

#offset_file = [4e-3, 3e-3, 2e-3, 0]
#offset_file = np.linspace(4e-3,0,len(ch_files_sel))
offset_file = np.linspace(0,0,len(ch_files_sel))
print "offset_file ",offset_file


test_window_size = False
if test_window_size:
	win_peak_av_list = range(1,40,2)
	nloop = len(win_peak_av_list)
else:
	nloop = len(ch_files_sel)
#nloop = 6

#win_peak_av_list = range(1,100,5)
#win_grad_av_list = range(1,25,1)
#win_grad_av_list = [1,5,10,15, 20, 25, 30, 35, 40]
#nloop = len(win_peak_av_list)
#print "win_peak_av_list ",win_peak_av_list

#ninter_abel_list = [100, 250, 500, 750, 1000, 1500]

for iloop in range(nloop):
	
	if test_window_size:
		chf1 = ch_files_sel[0]
		win_peak_av = win_peak_av_list[iloop]
	else:
		chf1 = ch_files_sel[iloop]
	
	#
	#win_peak_av = 20
	#win_grad_av = win_grad_av_list[iloop]
	#loss_time_ch  =[]
	
	#ninter_abel = ninter_abel_list[iloop] 

	print "file to process ",chf1
	
	#read scope data
	time1 = time.time()
	tdat_chf, data_chf = adata.read_scope(dirname, chf1, tfactor=1e6)
	tstep = tdat_chf[1]-tdat_chf[0]
	print "time step ",tstep

	#plt.plot(pfromt, data_chf,'b')
	#plt.show()
	#sys.exit()
	
	time2 = time.time()
	
	loss_time, index_loss_time = adata.signal_loss_time(tdat_chf, data_chf, method)
	loss_time_all.append(loss_time)
	
	time3 = time.time()
	print "signal loss time ",loss_time
	
	lower_lim = loss_time - lower_lim_file[ch_files_sel.index(chf1)]
	upper_lim = loss_time + 1200#600
	

	time_interval = tdat_chf[2] - tdat_chf[1]
	interval = int(round(640e-3/time_interval))

	find_hminus = False
	pmxi, pmx, pmy, pmi, pmy_av, pmy_std = find_peaks(tdat_chf, data_chf, interval, find_hminus, x_low = lower_lim, x_high = upper_lim)
	t_peak_av, peak_av1 = movingaverage(pmx, pmy, win_peak_av)
	
	
	#peak_av = [-p/peak_av1[0] for p in peak_av1] #normalise
	peak_av_initial.append(peak_av1[0])
	peak_av = list(peak_av1)
	

	#translate time axis t_peak_av into distance
	r_list1_av = []
	rtot = 0
	for j, t in enumerate(t_peak_av):
		if j > 0:
			deltat = t - t_peak_av[j-1]
			p2 = fittedpfunc(t)
			p1 = fittedpfunc(t_peak_av[j-1])
			deltap = p2-p1
			disp_p1 = dispref*(p1/pref)**k2
			deltar = disp_p1*deltap/p1
			rtot = rtot + deltar
			r_list1_av.append(rtot)

			
			
	#t_peak_av2, peak_av2 = movingaverage(pmx, pmy, 10)
	#p_peak_av = map(fittedpfunc, np.array(t_peak_av)+toffset)
	
	#print "momenta at losses ",map(fittedpfunc,np.array([0.1,1.3,2.0,4.0,5.0,8.0,16.0])*1e3+ toffset)
	#plt.plot(tdat_chf, data_chf)
	#plt.plot(t_peak_av, peak_av,'ro-')
	#plt.show()
	
	#grad_av = np.gradient(peak_av)
	m = 2

	deriv_vs_r = True #calcualte dI/dr if True, otherwise dI/dt where I is the bunch monitor signal
	if deriv_vs_r:
		peak_avm = peak_av[1:]
		grad_av = [(peak_avm[i+m] - peak_avm[i-m])/(r_list1_av[i+m] - r_list1_av[i-m]) for i in range(m, len(peak_avm)-m)]
	
		t_grad_av = t_peak_av[m+1:len(peak_av)-m]
		
		
		#derivative vs time - used in figure
		grad_av_vst = [(peak_av[i+m] - peak_av[i-m])/(t_peak_av[i+m] - t_peak_av[i-m]) for i in range(m, len(peak_av)-m)]
		t_grad_av_vst = t_peak_av[m:len(peak_av)-m]
	else:
		grad_av = [(peak_av[i+m] - peak_av[i-m])/(t_peak_av[i+m] - t_peak_av[i-m]) for i in range(m, len(peak_av)-m)]
		t_grad_av = t_peak_av[m:len(peak_av)-m]
	
	if second_smoothing:
	#do without this averaging?
		t_grad_av_av, grad_av_av = movingaverage(t_grad_av, grad_av, win_grad_av)
		_, peak_av_av = movingaverage(t_peak_av, peak_av, win_grad_av)
	else:
		t_grad_av_av = t_grad_av
		grad_av_av = grad_av
		peak_av_av = peak_av
	
	#print "number of peaks ",len(pmx)
	

	grad_max = 0
	grad_max_time = 0
	t_av_grad_win = []
	data_grad_win = []
	for t, g in zip(t_grad_av,grad_av):
		t_av_grad_win.append(t)
		data_grad_win.append(g)
		if g > grad_max:
			grad_max = g
			grad_max_time = t
			index_grad_max = t_grad_av.index(t)
			
	grad_max_all.append(grad_max)
	

	duration_gradmax_est.append(1/grad_max)
	
				

	#find time when gradient is half minimum before max gradient time
	halfgrad_1_time = 0 
	zerograd_1_time = None
	mi = 10
	for j in range(index_grad_max,0,-1):
		if grad_av[j] < 0.5*grad_max and halfgrad_1_time == 0:
			halfgrad_1_time = t_grad_av[j]
			index_half_1 = j
		if j < index_grad_max - mi:
			if grad_av[j] > grad_av[j+mi] and grad_av[j] < 0.1*grad_max:
				zerograd_1_time = t_grad_av[j]
				index_zsig1 = j
				break
	

						
	#find time when gradient is half minimum after max gradient time
	#and find when signal reaches zero
	halfgrad_2_time = 0 
	zerograd_2_time = None
	mi = 10
	for j in range(index_grad_max, len(grad_av)):
		if grad_av[j] < 0.5*grad_max and halfgrad_2_time == 0:
			halfgrad_2_time = t_grad_av[j]
			index_half_2 = j
		if j > index_grad_max + mi:
			if grad_av[j] > grad_av[j-mi] and grad_av[j] < 0.1*grad_max:
				zerograd_2_time = t_grad_av[j]
				index_zsig2 = j
				p_zg_2 = fittedpfunc(zerograd_2_time)
				break
				
	zerograd_1_time_all.append(zerograd_1_time)
	zerograd_2_time_all.append(zerograd_2_time)
	
	if zerograd_1_time != None:
		peak_zerograd1 = peak_av[index_zsig1]
		p_zg_1 = fittedpfunc(zerograd_1_time)
		print "time steps ",index_zsig2 - index_zsig1			
	print "zerograd_1_time ",zerograd_1_time
	print "zerograd_2_time ",zerograd_2_time_all
	
			
			
	fnd1 = False
	for j,t in enumerate(t_grad_av_av):
		if t > zerograd_1_time and not fnd1:
			print "t, zerograd_1_time ",t, zerograd_1_time
			index_zsig1_av_av = j
			fnd1 = True
		if t > zerograd_2_time:
			print "t, zerograd_2_time ",t, zerograd_2_time
			index_zsig2_av_av = j
			break	
	

	if test_window_size:
		continue
	
	


	
	
	for it, t in enumerate(t_peak_av):
		if t > zerograd_2_time:
			print "index, time zerograd_2_time in t_peak_av ",it, t
			print "std of peak data after this time ",np.std(peak_av[it:])
			
			break
			
		
	print "p at either end of signal ",p_zg_1, p_zg_2
	print "effective delta p/p ",(p_zg_1-p_zg_2)/(0.5*(p_zg_1+p_zg_2))
	
	

	eta_1 = dispref*(p_zg_1/pref)**k2
	eta_2 = dispref*(p_zg_2/pref)**k2
	#print "dispersion ratio at either end of signal",(eta_2/eta_1)
	#print "inverse momentum ratio at either end ",p_zg_1/p_zg_2
	print "grid spacing ratio ",(eta_2/eta_1)*(p_zg_1/p_zg_2)


					
	#translate time axis into distance
	r_list1 = []
	r_list_psc1 = []
	rtot = 0
	rtot_psc = 0
	for j, t in enumerate(t_grad_av_av):
		if j > 0:
			deltat = t - t_grad_av_av[j-1]
			p2 = fittedpfunc(t)
			p1 = fittedpfunc(t_grad_av_av[j-1])
			deltap = p2-p1
			disp_p1 = dispref*(p1/pref)**k2
			deltar = disp_p1*deltap/p1
			rtot = rtot + deltar
			r_list1.append(rtot)
			
			#remove momentum dependence
			deltar_psc = disp_p1*deltap/pref
			rtot_psc = rtot_psc + deltar_psc
			r_list_psc1.append(rtot_psc)
	
	

	
			
	#r at first zerograd time
	r_zg1_t = r_list1[index_zsig1-1]
	r_list = [r-r_zg1_t for r in r_list1]
	r_list_psc = [r-r_zg1_t for r in r_list_psc1]
	
	

	
	r_zg1 = 0
	r_zg2 = r_list[index_zsig2-1]
	r_half1 = r_list[index_half_1]
	r_half2 = r_list[index_half_2]
	r_gradmax = r_list[index_grad_max-1]
	
	print "len t_peak_av ",len(t_peak_av)
	print "len t_grad_av ",index_zsig2 - index_zsig1
	print "len t_grad_av_av in window ",index_zsig2_av_av - index_zsig1_av_av
	
	#plt.plot(tdat_chf, data_chf)
	#plt.plot(t_peak_av, peak_av,'ro-')
	#plt.axvline(x=zerograd_1_time)
	#plt.axvline(x=zerograd_2_time)
	#plt.xlim(zerograd_1_time, zerograd_2_time)
	#plt.show()
	
	#plt.plot(t_grad_av_av, grad_av_av, 'ko-')
	#plt.plot(t_grad_av, grad_av, 'bo')
	#plt.axvline(x=zerograd_1_time)
	#plt.axvline(x=zerograd_2_time)
	#plt.show()
	###sys.exit()

	
	distance_fwhm_mm.append(1e3*(r_half2 - r_half1))
	#distance_full_mm.append(1e3*(r_zg2 - r_zg1))
	
	#distance_full_mm.append(1e3*r_list_psc[index_zsig2-1])

	print "delta r calculation"
	print "r_gradmax ",r_gradmax
	print "r difference ",r_zg2 - r_zg1
	#print "t ratio ",(zerograd_2_time - grad_max_time)/(grad_max_time - zerograd_1_time)
	#print "r ratio ",(r_zg2-r_gradmax)/(r_gradmax - r_zg1)
	
	offset = offset_file[ch_files_sel.index(chf1)]
	
	print "offset ",offset

	if zerograd_1_time == None:
		print "no zerograd_1_time for this case "
		continue
	
	tspt = zerograd_2_time - zerograd_1_time
	offset_to_time = tspt*offset/(r_zg2 - r_zg1)
	print "fall-off time ",tspt
	print "offset translated into time ",offset_to_time
	
	p_zg_offset = fittedpfunc(zerograd_2_time + offset_to_time)
	p_full_offset.append(p_zg_offset)
	
	distance_full_offset_mm.append(1e3*(r_zg2 + offset - r_zg1))

	
	#plt.plot(pmx, pmy_av, 'r.-')
	#plt.plot(pmx, pmy_std,'k.-')
	#plt.axvline(x=zerograd_2_time)
	#plt.axvline(x=zerograd_2_time+offset_to_time)
	#plt.show()
	
	
	
	#powerlaw fit
	profile_calc = True
	fit_data = True
	if profile_calc:
		#xdata = r_list_psc[index_zsig1_av_av-1: index_zsig2-1] #profile scaled
		xdata1 = r_list[index_zsig1_av_av-1: index_zsig2-1]
		#xdata = t_grad_av_av[index_zsig1_av_av:index_zsig2]
		ydata = grad_av_av[index_zsig1_av_av:index_zsig2]
		#xdata = r_list[index_grad_max-1: index_zsig2-1]
		#ydata = grad_av_av[index_grad_max:index_zsig2]
		print "number of data points ",len(xdata1)
		
		#print "xdata ",xdata1
		
		#plt.plot(xdata1, ydata,'ko-')
		#plt.show()
		
		jzero = None
		
		find_zero_gradient = True
		if find_zero_gradient:
			for j, y in enumerate(ydata):
				if y <= 0.0 and j > index_grad_max - index_zsig1_av_av:
					jzero = j-1
					break
		else:
			jzero = index_zsig2
	
		if jzero == None:
			print "didn't find zero gradient, use final point for jzero"
			jzero = len(ydata)-1
			
				
		if jzero != None:
			xdata1 = np.array(xdata1[:jzero])
			ydata = np.array(ydata[:jzero])
		else:
			xdata1 = np.array(xdata1) 
			ydata = np.array(ydata)
		
		#xdata1 = np.array([x-xdata1[0]+1e-10 for x in xdata1]) 

		
		xdata1 = abs(xdata1 - xdata1[-1]) #x w.r.t to final point
		xdata1 = xdata1 + offset
		xdata = xdata1[::-1]#reverse x data
		ydata = ydata[::-1]#reverse y data
		pred_raw = ydata/xdata #raw profile prediction
		
		#plt.plot(xdata, ydata,'ko-')
		#plt.plot(xdata1, ydata[::-1],'ro-')
		#plt.show()
		
		#exclude first few points in profile if overlarge
		i1_raw = 0
		for index, p in enumerate(pred_raw):
			if abs(p - pred_raw[-1]) < 1.0:
				i1_raw = index
				break
				
		
		abel_transform = True
		if abel_transform:
			#Abel transform f(r) -> F(x) employing F(x1) = sum(f(r)*dr/(r**2 - x1**2)
			#where we may assume r=x along the x-axis
			
			
			test_rms_calc = False
			if test_rms_calc:
				ytest = np.random.normal(0,2,1000)
				count, bins, ignored = plt.hist(ytest, 30, normed=True)
			
				#rms of histogram counts
				rms = (sum([(c)*(b**2) for c,b in zip(count,bins)])/sum(count))**0.5
			
				print "rms ",rms
			
				plt.show()
				sys.exit()
	

			test_abel = False
			if test_abel:
				fgauss = lambda r, sigma: 1/(sigma*((2*math.pi)**0.5)) * np.exp(-r**2/(2*(sigma**2)))
				fgauss_abel_int = lambda r, sigma, a: sigma*np.exp(-a**2/(2*(sigma**2)))*erf(np.sqrt(r**2 - a**2)/(np.sqrt(2)*sigma))

				ns = 5000
				rsig = 1
				rlim = 3*rsig
				rt = np.linspace(0,rlim,ns)
				dr = rt[1]-rt[0]
				fr = fgauss(rt, rsig)
				print "limit ",erf(rlim)

				#Analytic formula
				f_abel_a =  fgauss_abel_int(rlim,rsig,rt)
				
				
				
				
				f_abel_num = []
				f_abel_trapz = []
				f_abel_trapz_wsing = []
				f_abel_simps_num = []
				for i_r,r0 in enumerate(rt[:-1]):
				
					
					fr_numeric  = fgauss(rt, rsig)*rt
					fr_integ_numeric = 2*fr_numeric[i_r+1:]/((rt[i_r+1:]**2 - r0**2)**0.5)
					
					#deal with singularity somehow
					fr_integ_numeric_wsing = np.insert(fr_integ_numeric,0, fr_integ_numeric[0])

					
					fr_sum = sum(dr*fr_integ_numeric)
					
					#fr_sum_trapz = np.trapz(fr_integ_numeric[1:]/dr,rt[1:])
					
					#plt.plot(fr_integ_numeric,'ko')
					#plt.show()
					
					
					fr_sum_trapz = np.trapz(fr_integ_numeric, dx = dr)
					fr_sum_trapz_wsing = np.trapz(fr_integ_numeric_wsing, dx = dr)
					#print "fr_sum, fr_sum_trapz ",fr_sum, fr_sum_trapz
					
					f_abel_num.append(fr_sum)
					f_abel_trapz.append(fr_sum_trapz)
					f_abel_trapz_wsing.append(fr_sum_trapz_wsing)
				
					#print "first fr ",fr_numeric[:10]
					#print "first integ ",fr_integ_numeric[:10]
				
					#print "sum ",fr_sum
					#print "trap ",
					
					#simps_int = simps(fr_integ_numeric[1:])
					#print simps_int
					#f_abel_simps_num.append(simps_int)
					
					
					#print "simpsons ",simps(fr_integ_numeric[1:])
					#print "discrepancy ",erf(rlim) - fr_sum 
					#print fgauss_abel_int(rlim, rsig, dr) - fgauss_abel_int(rlim, rsig, 0)

				#print "f_abel_num ",f_abel_num
				#plt.plot(rt,fr)
				
				plt.plot(rt, f_abel_a,'k',label='analytic')
				#plt.plot(rt[:-1], f_abel_num,'b')
				plt.plot(rt[:-1], f_abel_trapz,'r',label='numerical')
				plt.plot(rt[:-1], f_abel_trapz_wsing,'b',label='numerical + singularity term')
				#plt.plot(rt[:10], f_abel_simps_num,'g')
				
				
				plt.xlabel('x (mm) ')
				plt.ylabel(r'F(x) from Gaussian ($\sigma=1$ mm)')
				plt.legend(loc='upper right')
				plt.title('N = '+str(ns))
				plt.savefig('gauss_test')
				plt.show()
				sys.exit()
				
				
	
	
			beta_smooth = 4.6/3.65

			
			fx_data = ydata/(2*math.pi)
			
			if test_abel:
				#test!
				xdata = rt
				fx_data = fr_numeric
			
			interpolate = False
			if interpolate:
				nint = ninter_abel
				xdata_int = np.linspace(min(xdata), max(xdata), nint)
				fx_data_int = np.interp(xdata_int, xdata, fx_data)
			else:
				xdata_int = xdata
				fx_data_int = fx_data
			

			F_list = []
			F_list_trapz = []
			F_list_trapz_wsing = []
			F_n_list = []
			x_trans = []
			for j in range(len(xdata_int)-1):
				tot = 0
				tot_n = 0			
				xobs = xdata_int[j]
				#xobs_n = xdata_n[j]	
				index = j+1

				integ = []
				integ1 = []
				xinteg1 = []
				
				
				for x, fx in zip(xdata_int[j+1:len(xdata_int)-1], fx_data_int[j+1:len(xdata_int)-1]):
					
					dx = xdata_int[index+1] - x
					
					denom = (x**2 - xobs**2)**0.5
					#print "first int part ",dx, fx, 2*fx*dx/denom
					
					
					integ.append(2*fx*dx/denom)
					
					
					integ1.append(2*fx/denom)
					xinteg1.append(x)
					
					index = index + 1
				
				tot = sum(integ)
				
				
				trapz = np.trapz(integ1, x=xinteg1)
				F_list_trapz.append(trapz)
				
				if len(integ1) > 0:
					#deal with singularity
					xinteg2  = xinteg1.insert(0, 0)
					integ2  = [integ1[0]] + integ1
					trapz_wsing = np.trapz(integ2, x=xinteg2)
					F_list_trapz_wsing.append(trapz_wsing)				
				
				
				F_list.append(tot)
				#F_n_list.append(2*tot_n)
				x_trans.append(xobs)
			
	

			F_norm = [F/max(F_list) for F in F_list] 
			rms_abel_norm = (sum([(c)*(b**2) for c,b in zip(F_norm,x_trans)])/sum(F_norm))**0.5
			
			rms_abel = (sum([(c)*(b**2) for c,b in zip(F_list,x_trans)])/sum(F_list))**0.5
			rms_abel_trapz = (sum([(c)*(b**2) for c,b in zip(F_list_trapz,x_trans)])/sum(F_list_trapz))**0.5
			
			
			rms_abel_trapz_wsing = (sum([(c)*(b**2) for c,b in zip(F_list_trapz_wsing,x_trans)])/sum(F_list_trapz_wsing))**0.5
			
			distance_rms_abel_mm.append(1e3*rms_abel_trapz)
			
			
			
			for i, x in enumerate(xdata):
				if i > 0:
					if x > rms_abel and xdata[i-1] < rms_abel:
						index_xdata_rms = i

			index_t_grad_av_av_rms = index_zsig1_av_av-1 + index_xdata_rms
			duration_rms_abel1 = t_grad_av_av[index_zsig2-1]-t_grad_av_av[index_zsig1_av_av + len(xdata1)-index_xdata_rms] 
			

			print "duration_rms_abel ",duration_rms_abel1
			#plt.plot(xdata, ydata)
			duration_rms_abel.append(duration_rms_abel1)
			
			#calculate rms of radial profile from pred_raw
			pred_raw_m = pred_raw[i1_raw:]
			xdata_m = xdata[i1_raw:]
			rms_pred_raw = (sum([(c)*(b**2) for c,b in zip(pred_raw_m,xdata_m)])/sum(pred_raw_m))**0.5
			
			print "F(x) rms ",rms_abel
			print "F(x) rms trapz ",rms_abel_trapz
			print "F(x) rms trapz with sing ",rms_abel_trapz_wsing
			#print "F(x) norm rms ",rms_abel_norm
			print "rms_pred_raw ",rms_pred_raw
			
			Faxis_rms_abel.append(F_list_trapz[0])

			
		
			plot_abel_profile = False
			if plot_abel_profile:
				#plt.subplot(211)
				#plt.plot(x_trans, F_list,'ko')
				plt.plot(x_trans, F_list_trapz,'o')
				#plt.plot(rt, f_abel_a,'b-')
				#plt.plot(rt, f_abel_num,'r')
				
				plt.axvline(x=rms_abel_trapz)
				plt.ylabel(r'$\hat{F}(x)$')
				#plt.plot(x_trans, F_n_list,'r-')
				#plt.subplot(212)
				#plt.plot(xdata[i1_raw:], pred_raw[i1_raw:],'k-')
				#plt.axvline(x=rms_pred_raw)
				#plt.ylabel('f(x)')
				plt.show()
				sys.exit()
			
			

			

		
		if indices_select == []:
			print "profile calculation section expects non-empty indices_select"
			print "chf1 ",chf1, ch_files_sel.index(chf1), fpos[ch_files_sel.index(chf1)]
			#sys.exit()
			
			
			
		#write to file
		if "F1" in dirname:
			#ff = open("profile_"+str(fpos[indices_select[ch_files_sel.index(chf1)]])+"_F1.txt","w")
			ff = open("profile_"+str(fpos[ch_files_sel.index(chf1)])+"_F1.txt","w")
		if "F5" in dirname:
			ff = open("profile_"+str(fpos[ch_files_sel.index(chf1)])+"_F5.txt","w")
		if "F7" in dirname:
			#ff = open("profile_"+str(fpos[indices_select[ch_files_sel.index(chf1)]])+"_F7.txt","w")
			ff = open("profile_"+str(fpos[ch_files_sel.index(chf1)])+"_F7.txt","w")
			
			
		print >>ff, "x  f "
		print >>ff, "mm   "
		for pos, f in zip(xdata, pred_raw):
			print >> ff, '%.4f %.6f ' %(1e3*pos, f)
		ff.close()		

		
		#record distribution at some fixed scaled radius 
		r_sel = 0.0005
		isel = 0
		for r, d in zip(xdata, pred_raw):
			if r > r_sel + offset:
				r_check.append(r)
				distrib_check.append(d/max(pred_raw[isel:]))
				break
			isel = isel + 1

	
		test_degree = False
		if test_degree:
			#exponent, amp, chisq = powerlaw_fit(xdata, ydata) 
			chisq_list = []
			for degree in range(2,7):
				polyres = np.polyfit(xdata, ydata, degree, full=True)
				polyfit = np.poly1d(polyres[0])
				chisq = np.sum(((ydata - polyfit(xdata))**2),axis=0)
				chisq_list.append(chisq)
				print "polynomial fit degree, chisq ",degree,chisq
		else:
			polyres = np.polyfit(xdata, ydata, 3, full=True)
			polyfit = np.poly1d(polyres[0])			

		show_data = True
		if show_data:
			#plt.subplot(211)
			#plt.plot(xdata, ydata, 'ko')
			#plt.plot(xdata, chidistsc(xdata, xdata[-1], 1e-4),'r')
			#plt.subplot(212)
			#
			#plt.gca().set_color_cycle(colors)
			plt.plot(xdata, pred_raw, color = colors[ch_files_sel.index(chf1)])
			
			
			#plt.ylim(0,5)

			#print "polyres ",polyres 
			#print "polyfit ",polyfit
	
	
	if fit_data and profile_calc:

		
		#popt2, pcov2 = optimize.curve_fit(fitfn,xdata, ydata)
		
		popt3, pcov3 = optimize.curve_fit(fitfn3,xdata, ydata)
		fit3_chisq = np.sum(((ydata - fitfn3(xdata, popt3[0], popt3[1], popt3[2]))**2),axis=0)
		
		try:
			chi_popt, chi_pcov = optimize.curve_fit(chidistsc, xdata, ydata, p0 = [  xdata[-1], 1.94622566e-05])
			chi_chisq = np.sum(((ydata - chidistsc(xdata, chi_popt[0], chi_popt[1]))**2),axis=0)
			print "chi distribution param ",chi_popt, "residual ",chi_chisq
		except:
			pass
		
		
		
		#chimod_popt, chimod_pcov = optimize.curve_fit(chidistmod,xdata, ydata)
		#chiexp_popt, chiexp_pcov = optimize.curve_fit(chiexpand, xdata, ydata, p0 = [ xdata[-1], 1.94622566e-05])
		#chiexp_chisq = np.sum(((ydata - chiexpand(xdata, chiexp_popt[0], chiexp_popt[1]))**2),axis=0)
		##print "chi expansion param ",chiexp_popt, "residual ",chiexp_chisq
		
		
		plot_expm = False
		#try:
			#chiexpm_popt, chiexpm_pcov = optimize.curve_fit(chiexpandm, xdata, ydata, p0=[chi_popt[0], chi_popt[1], 0])
			#chiexpm_chisq = np.sum(((ydata - chiexpandm(xdata, chiexpm_popt[0], chiexpm_popt[1], chiexpm_popt[2]))**2),axis=0)
			#print "modiftied chi expansion param ",chiexpm_popt, "residual ",chiexpm_chisq
			#plot_expm = True
		#except:
			#pass
		#print "curve fit result ", popt
		
		
		
		#print "polynomial, third order param ",popt3, "residual ",fit3_chisq
		#print "chi_popt (modulated) ",chimod_popt
		
		#format time and momentum for display
		tform = "%.1f" %zerograd_2_time
		pform = "%.1f" %p_zg_2
		
		sample_check = False
		if sample_check:
			sample_data = chidistsc(xdata, chi_popt[0], chi_popt[1])
			#fit sample data with chi distribution
			samp_popt, samp_pcov2 = optimize.curve_fit(chidistsc, xdata, sample_data, p0 = [  chi_popt[0], 1.94622566e-05])
			samp_chisq = np.sum(((sample_data - chidistsc(xdata, samp_popt[0], samp_popt[1]))**2),axis=0)
			#fit sample data with Taylor expansion
			samp_exp_popt, pcov = optimize.curve_fit(chiexpand, xdata, sample_data, p0=[chiexp_popt[0], chiexp_popt[1]])
			samp_exp_chisq = np.sum(((sample_data - chiexpand(xdata, samp_exp_popt[0], samp_exp_popt[1]))**2),axis=0)		
			#fit sample data with modified Taylor expansion
			samp_expm_popt, pcov = optimize.curve_fit(chiexpandm, xdata, sample_data, p0=[chiexp_popt[0], chiexp_popt[1], 0])
			samp_expm_chisq = np.sum(((sample_data - chiexpandm(xdata, samp_expm_popt[0], samp_expm_popt[1], samp_expm_popt[2]))**2),axis=0)
			print "chi distribution param (sample data) ",samp_popt, "residual ",samp_chisq
			print "Taylor expansion param (sample data) ",samp_exp_popt, "residual ",samp_exp_chisq
			print "modiftied Taylor expansion param (sample data) ",samp_expm_popt, "residual ",samp_expm_chisq
			
			plot_sample_data = False
			if plot_sample_data:
				plt.subplot(211)
				plt.plot(xdata, sample_data, 'k.')
				plt.plot(xdata, chidistsc(xdata, samp_popt[0], samp_popt[1]), 'r-')
				plt.plot(xdata, chiexpand(xdata, samp_exp_popt[0], samp_exp_popt[1]), 'g-')
				plt.plot(xdata, chiexpandm(xdata, samp_expm_popt[0], samp_expm_popt[1], samp_expm_popt[2]), 'b-')
				plt.subplot(212)
				plt.plot(xdata, chidistsc(xdata, samp_popt[0], samp_popt[1])/xdata, 'r-')
				plt.plot(xdata, chiexpand(xdata, samp_exp_popt[0], samp_exp_popt[1])/xdata, 'g-')
				plt.plot(xdata, chiexpandm(xdata, samp_expm_popt[0], samp_expm_popt[1], samp_expm_popt[2])/xdata, 'b-')
				plt.show()
			
		plt.subplot(211)
		plt.plot(xdata, ydata, 'k.', label='data')
		#plt.plot(xdata, polyfit(xdata),'r--')
		#plt.plot(xdata, fitfn(xdata, popt2[0], popt2[1]), 'b-')
		#plt.plot(xdata, fitfn3(xdata, popt3[0], popt3[1], popt3[2]), 'g-',linewidth=2)
		plt.plot(xdata, chidistsc(xdata, chi_popt[0], chi_popt[1]),'r',linewidth=2,label=r'fitted $\chi$-distribution')
		#plt.plot(xdata, chiexpand(xdata, chiexp_popt[0], chiexp_popt[1]),'b')
		if plot_expm:
			plt.plot(xdata, chiexpandm(xdata, chiexpm_popt[0], chiexpm_popt[1], chiexpm_popt[2]),'b--',linewidth=2)
		plt.xlim(0,xdata[-1])
		plt.xlabel('r (m)')
		plt.ylabel('signal derivative')
		plt.legend()
		plt.title(r'time ($\mu$s) = '+str(tform)+'    p (MeV/c) = '+str(pform))
		plt.subplot(212)
		plt.plot(xdata[1:], pred_raw[1:], 'k.', label='data/r')
		#plt.plot(xdata[1:], polyfit(xdata[1:])/xdata[1:],'r--')
		#plt.plot(xdata, fitfn(xdata, popt2[0], popt2[1])/xdata, 'b-')
		#plt.plot(xdata, fitfn3(xdata, popt3[0], popt3[1], popt3[2])/xdata,'r-')
		plt.plot(xdata, chidistsc(xdata, chi_popt[0], chi_popt[1])/xdata,'r',linewidth=2,label=r'Gaussian distribution')
		
		if plot_expm:
			plt.plot(xdata, chiexpandm(xdata, chiexpm_popt[0], chiexpm_popt[1], chiexpm_popt[2])/xdata,'b--',linewidth=2)
		
		plt.xlim(0,xdata[-1])
		plt.ylabel('calculated f(r)')
		plt.xlabel('r (m)')
		plt.legend()
		plt.savefig('fit_dist')
		plt.show()
		sys.exit()
		
		
	#plt.plot(t_grad_av_av, grad_av_av, 'ko-')
	#plt.plot(t_grad_av, grad_av, 'bo')
	#plt.axvline(x=zerograd_1_time)
	#plt.axvline(x=zerograd_2_time)
	#plt.xlim(zerograd_1_time, zerograd_2_time)
	#plt.axvline(x=t_grad_av_av[index_zsig1 + len(xdata1)-index_xdata_rms] )
	#plt.show()
	
	compare_axes = False
	if compare_axes:
		plt.subplot(211)
		#plt.plot(t_grad_av, grad_av,'r.-')
		plt.plot(t_grad_av_av, grad_av_av,'k.-')
	
	
		#plt.plot(t_peak_av2, np.gradient(peak_av2),'k.--')
		plt.axvline(x = grad_max_time, color='k',linestyle='--',label='max gradient')
		plt.axvline(x = halfgrad_1_time, color='c')
		plt.axvline(x = halfgrad_2_time, color='c',label='FWHM')
		if zerograd_1_time != None:
			plt.axvline(x = zerograd_1_time, color='m')
		plt.axvline(x = zerograd_2_time, color='m',label='100%')
	
		#plt.plot(tsamp, ex, 'b')
		#plt.axhline(y=20e-6)
		plt.xlim(min(t_grad_av_av), max(t_grad_av_av))
		plt.xlabel(r't ($\mu$s)') 
		plt.ylabel(r'rate of signal change per $\mu$s')
		plt.subplot(212)		
		plt.plot(r_list, grad_av_av[1:],'k.-')
		plt.axvline(x=r_gradmax,color='k',linestyle='--',label='max gradient')
		plt.axvline(x = r_half1, color='c')
		plt.axvline(x = r_half2, color='c',label='FWHM')
		if zerograd_1_time != None:
			plt.axvline(x = r_zg1, color='m')
		plt.axvline(x = r_zg2, color='m',label='100%')
		plt.ylabel(r'rate of signal change per $\mu$s')
		plt.xlabel('x (m)')
		plt.xlim(min(r_list), max(r_list))
		plt.savefig('compare_axes')
		plt.show()
	

	
	
	#print "peak_zerograd1 ",peak_zerograd1
	#find rms emittance where peak falls to 39%
	
	
	duration_fwhm = halfgrad_2_time - halfgrad_1_time
	duration_fwhm_all.append(duration_fwhm)
	
	#measure of beam size may be gradient min to zero sig time (forward in time)
	duration_gradmax_zerograd = zerograd_2_time - grad_max_time
	duration_gradmax_zerograd_all.append(duration_gradmax_zerograd)
		
	
	
	
	#print "signal at grad max from window averaging ",peak_av_av[index_grad_max]
	print "beam duration estimate from grad max ",peak_av_av[index_grad_max]/grad_max
	print "FWHM duration ",duration_fwhm
	print "time from minimum gradient to zero signal (forward) ",-duration_gradmax_zerograd
	
	if zerograd_1_time != None and zerograd_2_time != None:
		loss_time_full.append(zerograd_2_time)
		loss_time_full_offset.append(zerograd_2_time + offset_to_time)
		duration_full.append(zerograd_2_time - zerograd_1_time)
		duration_full_offset.append(zerograd_2_time + offset_to_time - zerograd_1_time)
		print "time between zero signals ",zerograd_2_time - zerograd_1_time
		
		print "files index ",files.index(chf1)
		
		fpos_out.append(fpos[files.index(chf1)])

		p_full.append(p_zg_2)	
		distance_full_mm.append(1e3*(r_zg2 - r_zg1))
	
	
	tsamp = np.linspace(lower_lim,zerograd_2_time)
	ex = -(grad_max/60)*fn(tsamp, duration_gradmax_zerograd, zerograd_2_time)
	ex_scaled = [e/(zerograd_2_time - ts) for ts, e in zip(tsamp, ex)]
	
	#plt.plot(tsamp, ex)
	#plt.plot(tsamp, ex_scaled,'k-')
	#plt.axvline(x=zerograd_time,color='g')
	#plt.axvline(x=zerograd_time + duration_gradmax_zerograd,color='m')
	#plt.show()
	#sys.exit()
	
	show_signal = False
	if show_signal:
		plt.plot(tdat_chf, data_chf,'b',label='raw data')
		plt.plot(pmx, pmy,'k.',label='amplitude')
		plt.plot(t_peak_av, peak_av,'r-',label='<amplitude>')
		plt.xlim(min(pmx), max(pmx))
		plt.ylim(1.1*min(peak_av), 0.1)
		plt.xlabel(r't ($\mu$s)') 
		plt.ylabel(r'signal')
		plt.legend()
		plt.savefig("dropoff")
		plt.show()
		
	
	
	flip_x = False
	if show_gradient:
		
		sizef = 16
		dur1 = zerograd_2_time - zerograd_1_time
		tplot1 = 1e3*15.7
		tplot2 = 1e3*16.2
		ratio = dur1/(tplot2 - tplot1)
		print "ratio of duration to full time range ",ratio
		xrange1 = x_trans[-1]/ratio
		xextra = 0.5*(xrange1-x_trans[-1])
		print "x range to get same ratio ",xrange1, "extra on either side of range ",xextra
		
		
		plt.subplot(311)
		plt.plot([1e-3*t for t in tdat_chf], data_chf,'k')
		#plt.plot(pmx, pmy,'k.')
		
		plt.plot([1e-3*t for t in t_peak_av], peak_av,'r-', linewidth=2)
		
		#plt.xlim(t_peak_av[0], t_peak_av[-1])
		
		plt.xlim(15.700, 16.200)
		plt.ylim(-0.15,0.05)
		
		
		#plt.axvline(x = grad_max_time, color='k',linestyle='--')
		#plt.axvline(x = halfgrad_1_time, color='c')
		#plt.axvline(x = halfgrad_2_time, color='c')
		#if zerograd_1_time != None:
		#	plt.axvline(x = zerograd_1_time, color='m')
		#plt.axvline(x = zerograd_2_time, color='m',label='100%')
		#plt.axvline(x = zerograd_2_time + offset_to_time, color='m',linestyle='--',label='100%')
		
		plt.xlabel(r't [ms]',size=sizef) 
		plt.ylabel('Bunch monitor signal [V]',size=sizef)
		#plt.ylabel(r'normalised signal')
		
		
		show_duration = True
		if show_duration:
			plt.subplot(312)
			#plt.plot(pmx, np.gradient(pmy),'k.-')
			#plt.plot(pmx[0:-1:50], np.gradient(pmy[0:-1:50]),'k.-')
			#plt.plot(pmx[0:-1:20], np.gradient(pmy[0:-1:20]),'k.--')
		
			#plt.plot([t*1e-3 for t in t_grad_av], [1e6*g for g in grad_av],'b.-')
			plt.plot([t*1e-3 for t in t_grad_av_vst], [1e6*g for g in grad_av_vst],'b.-')
			
			if second_smoothing:
				plt.plot([t*1e-3 for t in t_grad_av_av], [1e6*g for g in grad_av_av],'b.-')
		
		
			#plt.plot(t_peak_av2, np.gradient(peak_av2),'k.--')
			#plt.axvline(x = grad_max_time, color='k',linestyle='--',label='max gradient')
			
			
			#FWHM lines
			#plt.axvline(x = halfgrad_1_time, color='k', linestyle = '--')
			#plt.axvline(x = halfgrad_2_time, color='k',linestyle = '--')
			
			if zerograd_1_time != None:
				plt.axvline(x = 1e-3*zerograd_1_time, color='k',linestyle = '--',linewidth=2)
			plt.axvline(x = 1e-3*zerograd_2_time, color='k',linestyle = '--',linewidth=2)
			
			
			#corresponding duration
			#plt.axvline(x=1e-3*t_grad_av_av[index_zsig1 + len(xdata1)-index_xdata_rms] ,color='k',linewidth=2)
			
			#plt.plot(tsamp, ex, 'b')
			#plt.axhline(y=20e-6)
			
			
			#plt.xlim(min(pmx), max(pmx))
			
			plt.xlim(15.700, 16.200)
			
			plt.xlabel(r't [ms]',size=sizef) 
			plt.ylabel(r'Signal derivative [V/s]',size=sizef)
			
			#Abel transform
			plt.subplot(313)
			
			
			plt.plot(x_trans, F_norm,'b.')
			
			
			#plt.axvline(x=x_trans[-1]-rms_abel)
			
			#plt.xlim(-xextra, x_trans[-1]+xextra)
			
			if flip_x:
				plt.xlim(x_trans[-1]+xextra, -xextra)
				plt.axvline(x = 0, color='k',linestyle = '--', linewidth=2)
				plt.axvline(x = x_trans[-1], color='k',linestyle = '--',linewidth=2)
			else:
				plt.xlim(0, x_trans[-1])
			
			plt.xlabel(r'x [mm]',size=sizef) 
			plt.ylabel(r'$\hat{F}(x)$',size=sizef)
			
			plt.ylim(0, 1)
			

			plt.axvline(x=rms_abel,color='k',linewidth=2)
			
		else:
			plt.subplot(212)
			plt.plot(r_list, grad_av_av[1:],'k.-')
			plt.axvline(x=r_gradmax,color='k',linestyle='--',label='max gradient')
			plt.axvline(x = r_half1, color='c')
			plt.axvline(x = r_half2, color='c',label='FWHM')
			if zerograd_1_time != None:
				plt.axvline(x = r_zg1, color='m')
			plt.axvline(x = r_zg2, color='m',label='100%')
			plt.ylabel(r'rate of signal change per $\mu$s',size=12)
			
			plt.xlim(min(r_list), max(r_list))
			plt.axhline(y=0,color='c')
			
		
		
		plt.legend()
		plt.savefig('gradient')
		plt.show()
		sys.exit()


		

	#scale gradient by time from zerograd
	grad_scaled = []
	t_grad_scaled = []
	for t, g in zip(t_grad_av_av, grad_av_av):
			if zerograd_2_time - t > 0:
				grad_scaled.append(g/(zerograd_2_time - t))
				t_grad_scaled.append(t)
	
	
	print "lower and upper time range ",lower_lim, upper_lim


	
	calculate_profile = False
	if calculate_profile:
		plt.subplot(211)
		plt.plot(t_grad_av, grad_av,'r.-')
		plt.plot(t_grad_av_av, grad_av_av, 'k.-')
		plt.axvline(x=zerograd_2_time,color='g')
		plt.xlim(t_grad_av[0], t_grad_av[-1])
		plt.subplot(212)
		plt.plot(t_grad_scaled, grad_scaled, 'r.-')
		plt.axvline(x=zerograd_2_time,color='g')
		plt.xlim(t_grad_av[0], t_grad_av[-1])
		plt.show()
		
	if show_all_data:
		
		#plt.plot([t - zerograd_1_time for t in t_grad_av_av], [g+3e-6*loss_time for g in grad_av_av],'k-')
		plt.plot([t - zerograd_2_time for t in t_grad_av_av], [2e3*g+p_zg_2 for g in grad_av_av],'k-')
	
#loss_time_all.append(loss_time_ch)
print "rms Abel (trap) ",distance_rms_abel_mm


#plt.show()

#plt.plot(win_grad_av_list, Faxis_rms_abel, 'ko')
#plt.ylabel('F(0)')
#plt.show()

check_interpolation = False
if check_interpolation:

	plt.xlabel('Interpolation N')
	plt.ylabel('rms of projected profile (mm) ')
	plt.plot(ninter_abel_list , distance_rms_abel_mm,'ko-')
	plt.savefig('ninterpolation')
	plt.show()
	sys.exit()

if show_all_profiles:
	#plt.plot(r_check, distrib_check, 'ko')
	#plt.ylim(0,1.0)
	plt.ylim(0,4.3)
	#plt.xlabel('r (m)')
	plt.xlabel('r scaled')
	plt.ylabel('radial profile f(r)')
	plt.savefig('profiles')
	plt.show()


if test_window_size:
	win_list = win_peak_av_list

	#plt.subplot(211)
	plt.plot(win_list, zerograd_1_time_all,'ko-')
	plt.plot(win_list, zerograd_2_time_all,'ro-')
	#plt.subplot(212)
	#plt.plot(win_list, distance_rms_abel_mm,'ko-')
	plt.savefig('win_size_check.png')
	plt.show()


	sys.exit()



if show_all_data:
	plt.ylabel('p (MeV/c)')
	plt.savefig('mountainrange')
	plt.show()

print "loss times ",loss_time_all, "len ",len(loss_time_all)

print "peak_av_initial ", peak_av_initial

print "distance_full_mm ",distance_full_mm
print " distance_full_offset_mm ", distance_full_offset_mm
print "distance_rms_abel_mm ",distance_rms_abel_mm

print "duration_rms_abel ",duration_rms_abel

print "p_full ",p_full

if len(duration_fwhm_all) > 1:
	
	plt.subplot(211)
	#plt.plot(p_full, duration_fwhm_all,'co-',label='FWHM')
	#plt.plot(loss_time_all, duration_gradmax_zerograd_all, 'ro-',label='grad max to zero')
	#plt.plot(p_full, duration_gradmax_est,'mo-',label='max gradient estimate')
	plt.plot(p_full, duration_full,'ko-',label='full')
	plt.plot(p_full, duration_rms_abel,'co-',label='rms')
	#plt.plot(loss_time_full, duration_full_offset,'ro-',label='full, offset')
	plt.ylabel(r'duration ($\mu$s)')
	plt.xlabel(r'time at loss ($\mu$s)')
	plt.legend(loc='upper center', prop={'size':8})
	plt.ylabel(r"fall-off duration [$\mu$s]")
	plt.subplot(212)
	plt.plot(p_full, distance_full_mm,'ko-',label='full')
	#plt.plot(p_full, distance_full_offset_mm,'ro-',label='full, offset')
	plt.plot(p_full, distance_rms_abel_mm,'co-',label='rms')
	#plt.plot(p_full, distance_fwhm_mm,'co-',label='FWHM')
	plt.xlabel(r'p (MeV/c)')
	plt.ylabel("beam size")
	plt.legend(loc='upper center', prop={'size':8})
	plt.ylim(0, max(distance_full_mm))
	plt.ylabel(r"beam radius [mm]")
	plt.show()

	#plt.plot(loss_time_full, grad_max_all,'bo')
	#plt.xlabel(r'time at loss ($\mu$s)')
	#plt.ylabel('gradient max')
	
	#plt.show()
	
	#plt.plot(loss_time_full,'ko')
	#plt.plot(zerograd_2_time_all, 'go')
	#plt.savefig('duration_all')
	#plt.show()




print "len fpos_out, loss_time_full ",len(fpos_out), len(loss_time_full)
print "fpos_out ",fpos_out


#plt.plot(p_full, fpos_out,'k.-')
#plt.plot(p_full_offset, fpos_out,'r.-')
#plt.show()


print "write results to file"

results = zip(fpos, loss_time_full, p_full,  distance_full_mm, distance_rms_abel_mm, duration_full, duration_rms_abel, duration_gradmax_est) 
ff = open(fout,"w")
print >>ff, "probe time    p    full, rms,  full, rms, 1/gradmax "
print >>ff, "mm     us    MeV/c   mm   mm   us     us     us "
for pos, losst, p, dfull, dfwhm, drfull, drfwhm, drg1 in results: 
	print >> ff, '%d %.2f %.2f %.2f %.2f %.2f %.2f %.2f ' %(pos, losst, p, dfull, dfwhm, drfull, drfwhm, drg1)
ff.close()

