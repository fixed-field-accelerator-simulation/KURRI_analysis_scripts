from __future__ import division
import math 
import sys
import numpy as np
import pylab as plt
import os
from scipy.optimize import minimize
from scipy.optimize import fmin

PROTON_MASS = 938.27203e6
SPEED_OF_LIGHT = 299792458

BtoEnom = lambda B, rho, Erest, q: ((q*B*rho*SPEED_OF_LIGHT)**2 + (Erest**2))**0.5
vrf = lambda phi, vrf1, vrf2, hratio, phi12: vrf1*np.sin(phi) + vrf2*np.sin(hratio*(phi-phi12))

#Hamiltonian
H_fn = lambda phi, dp, phi_s, a1, a2: a1*(dp**2) + a2*(np.cos(phi) - np.cos(phi_s) + (phi - phi_s)*np.sin(phi_s))


def running_mean(x,N):
	"""Running mean or moving window average of array x in window size N"""
	cumulsum = np.cumsum(np.insert(x,0,0))
	return (cumulsum[N:] - cumulsum[:-N])/float(N)

def running_mean_wrap(x,N):
	"""If x is periodic, calculate running mean taking into account wrapping at each end """
	
	nx = len(x)
	xtrip = np.concatenate((x,x,x))
	csum_all = running_mean(xtrip,N)
	#cspl = np.split(csum_all, 3)
	#csum = cspl[1]

	n1 = nx - int(0.5*N)
	xax1 = 0.5*N + np.arange(n1)
	xax2 = xax1[-1] + np.arange(nx)	
	csum_period = csum_all[n1:n1+nx]
	
	show_cumsum = False
	if show_cumsum:
		plt.subplot(211)
		plt.plot(xtrip)
		plt.plot(csum_all, 'r-')
		plt.plot(xax1, csum_all[:n1],'k-') 
		plt.plot(xax2, csum_period,'k-',linewidth=2) 
		plt.subplot(212)
		plt.plot(x)
		plt.plot(csum_period, 'r')
		plt.show()
	
	return csum_period

def get_floats_from_file(input_file="input_v2.dat"):

	f1 = open(input_file,'r')

	float_l = []
	for line in f1:
		if line[0] != '!':
			try:
				fl = float(line)
				float_l.append(fl)
			except:
				pass

	f1.close()


def get_input(input_file="input_v2.dat"):

	f1 = open(input_file,'r')

	float_l = []
	for line in f1:
		if line[0] != '!':
			try:
				fl = float(line)
				float_l.append(fl)
			except:
				pass

	f1.close()


	input_struct = np.zeros(1, dtype=[('numframes','i2'),('numframesignore','i2'),('numbins','i2'),('framebinwidth','f8'),('binslowerignore','i2'),('binsupperignore','i2'),('binslowerempty','i2'),('binsupperempty','i2'),('synchbin','i2'),('profilecount','i2'),('profilelength','i2'),('dtbin','f8'),('dturns','i2'),('rebin','i2'),('y0','f2'),('dEmax','f8'),('filmstart','i2'),('filmstop','i2'),('filmstep','i2'),('iterations','i2'),('npartroot','i2'),('extend_phasespace','i2'),('beamref','i2'),('machineref','i2'),('VRF1ref','f8'),('VRF1dot','f8'),('VRF2ref','f8'),('VRF2dot','f8'),('h','i2'),('hratio','f8'),('phi12','f8'),('Bref','f8'),('Bdot','f8'),('Rnom','f8'),('rhonom','f8'),('gammatnom','f8'),('Erest','f8'),('q','i2'),('selffields','i2'),('coupling','f8'),('Zovern','f8'),('pickup','f8')])
	
	input_struct['numframes'] = float_l[0]
	input_struct['numframesignore'] = float_l[1]
	input_struct['numbins'] = float_l[2]
	input_struct['framebinwidth'] = float_l[3]
	input_struct['binslowerignore'] = float_l[5]
	input_struct['binsupperignore'] = float_l[6]
	input_struct['binslowerempty'] = float_l[7]
	input_struct['binsupperempty'] = float_l[8]
	input_struct['synchbin'] = float_l[10]

	input_struct['profilecount'] = float_l[0] - float_l[1]
	input_struct['profilelength'] = math.ceil((float_l[2] - (float_l[5] + float_l[6]))/float_l[9])
	input_struct['dtbin'] = float_l[3] * float_l[6]
	input_struct['dturns'] = float_l[4]
	input_struct['rebin'] = float_l[9]
	input_struct['y0'] = input_struct['profilelength']/2
	input_struct['dEmax'] = float_l[11]
	input_struct['filmstart'] = float_l[12]
	input_struct['filmstop'] = float_l[13]
	input_struct['filmstep'] = float_l[14]
	input_struct['iterations'] = float_l[15]
	input_struct['npartroot'] = float_l[16]
	input_struct['extend_phasespace'] = float_l[17]
	input_struct['beamref'] = float_l[18]
	input_struct['machineref'] = float_l[19]
	input_struct['VRF1ref'] = float_l[20]
	input_struct['VRF1dot'] = float_l[21]
	input_struct['VRF2ref'] = float_l[22]
	input_struct['VRF2dot'] = float_l[23]
	input_struct['h'] = float_l[24]
	input_struct['hratio'] = float_l[25]
	input_struct['phi12'] = float_l[26]
	input_struct['Bref'] = float_l[27]
	input_struct['Bdot'] = float_l[28]
	input_struct['Rnom'] = float_l[29]
	input_struct['rhonom'] = float_l[30]
	input_struct['gammatnom'] = float_l[31]
	input_struct['Erest'] = float_l[32]
	input_struct['q'] = float_l[33]

	if int(float_l[34]) == 1:
		input_struct['selffields'] = True
	else:
		input_struct['selffields'] = False

	input_struct['coupling'] = float_l[35]
	input_struct['Zovern'] = float_l[36]
	input_struct['pickup'] = float_l[37]

	return input_struct


def write_inputfile(input_file, descrip, datafile, inputp):

	

	f1 = open(input_file,'w')
	print >>f1, descrip
	for i in range(10):
		print >>f1, ''
	print >>f1, "! Input data file ="
	print >>f1, datafile
	print >>f1, "! Directory in which to write all output ="
	print >>f1, "./"
	print >>f1, "! Number of frames in input data ="
	print >>f1, inputp["numframes"][0]
	print >>f1, "! Number of frames to ignore ="
	print >>f1, inputp["numframesignore"][0]
	print >>f1, "! Number of bins in each frame ="
	print >>f1, inputp["numbins"][0]
	print >>f1, "! Width (in s) of each frame bin ="
	print >>f1, inputp["framebinwidth"][0]
	print >>f1, "! Number of machine turns between frames ="
	print >>f1, inputp["dturns"][0]
	print >>f1, "! Number of frame bins before the lower profile bound to ignore ="
	print >>f1, inputp["binslowerignore"][0]
	print >>f1, "! Number of frame bins after the upper profile bound to ignore ="
	print >>f1, inputp["binsupperignore"][0]
	print >>f1, "! Number of frame bins after the lower profile bound to treat as empty"	
	print >>f1, "! at the reconstructed time = "
	print >>f1, inputp["binslowerempty"][0]
	print >>f1, "! Number of frame bins after the lower profile bound to treat as empty"	
	print >>f1, "! at the reconstructed time = "
	print >>f1, inputp["binsupperempty"][0]
	print >>f1, "! Number of frame bins to rebin into one profile bin = "
	print >>f1, inputp["rebin"][0]
	print >>f1, "! Time (in frame bins) from the lower profile bound to the synchronous phase"
	print >>f1, "! (if <0, a fit is performed) in the bunch reference frame = "
	print >>f1, inputp["synchbin"][0]
	print >>f1, "! Max energy (in eV) of reconstructed phase space (if >0) = "
	print >>f1, inputp["dEmax"][0]	
	print >>f1, "! Number of the first profile at which to reconstruct ="
	print >>f1, inputp["filmstart"][0]
	print >>f1, "! Number of the last profile at which to reconstruct ="
	print >>f1, inputp["filmstop"][0]
	print >>f1, "! Step between reconstructions ="
	print >>f1, inputp["filmstep"][0]
	print >>f1, "! Number of iterations for each reconstruction ="
	print >>f1, inputp["iterations"][0]
	print >>f1, "! Square root of the number of test particles to track per cell ="
	print >>f1, inputp["npartroot"][0]
	print >>f1, "! Flag to extend the region in phase space of map elements (if =1) ="
	print >>f1, inputp["extend_phasespace"][0]
	print >>f1, "! Reference frame for bunch parameters (synchronous phase, baseline, integral) ="
	print >>f1, inputp["beamref"][0]
	print >>f1, "! Reference frame for machine parameters (RF voltages, B-field) ="
	print >>f1, inputp["machineref"][0]
	print >>f1, "!"
	print >>f1, "! Machine and Particle Parameters:"
	print >>f1, "! Peak RF voltage (in V) of principal RF system ="
	print >>f1, inputp['VRF1ref'][0]
	print >>f1, "! and its time derivative (in V/s) ="
	print >>f1, inputp['VRF1dot'][0]
	print >>f1, "! Peak RF voltage (in V) of higher-harmonic RF system ="
	print >>f1, inputp['VRF2ref'][0]
	print >>f1, "! and its time derivative (in V/s) ="
	print >>f1, inputp['VRF2dot'][0]
	print >>f1, "! Harmonic number of principal RF system ="
	print >>f1, inputp['h'][0]
	print >>f1, "! Ratio of harmonics between RF systems ="
	print >>f1, inputp['hratio'][0]
	print >>f1, "! Phase difference (in radians of the principal harmonic) between RF systems ="
	print >>f1, inputp['phi12'][0]
	print >>f1, "! Dipole magnetic field (in T) ="
	print >>f1, inputp["Bref"][0]
	print >>f1, "! and its time derivative (in T/s) ="
	print >>f1, inputp["Bdot"][0]
	print >>f1, "! Machine radius (in m) ="
	print >>f1, inputp["Rnom"][0]
	print >>f1, "! Bending radius (in m) ="
	print >>f1, inputp["rhonom"][0]
	print >>f1, "! Gamma transition ="
	print >>f1, inputp["gammatnom"][0]
	print >>f1, "! Rest mass (in eV/c**2) of accelerated particle ="
	print >>f1, inputp["Erest"][0]
	print >>f1, "! Charge state of accelerated particle ="
	print >>f1, inputp["q"][0]
	print >>f1, "!"
	print >>f1, "! Space Charge Parameters:"
	print >>f1, "! Flag to include self-fields in the tracking (if =1) ="
	print >>f1, "0"
	print >>f1, "! Geometrical coupling coefficient = "
	print >>f1, inputp['coupling'][0]
	print >>f1, "! Reactive impedance (in Ohms per mode number) over a machine turn ="
	print >>f1, inputp['Zovern'][0]
	print >>f1, "! Effective pick-up sensitivity (in digitizer units per instantaneous Amp) ="
	print >>f1, inputp['pickup'][0]



write_input_param = False

if write_input_param:
	inputp = get_input()
	print inputp.dtype.names

	print "dtbin ",inputp['dtbin']
	print "dturns ",inputp['dturns']

	print "profilelength ",inputp['profilelength']
	print "rebin ",inputp['rebin']
	print inputp['y0']
	print inputp['dEmax']
	print "filmstart, filmstop, filmstep ",inputp['filmstart'],inputp['filmstop'],inputp['filmstep']
	print "iterations ",inputp["iterations"]
	print "beam/machine ref ",inputp['beamref'], inputp['machineref']
	print "VRF1ref, VRF1dot ",inputp['VRF1ref'],inputp['VRF1dot']
	print "VRF2ref, VRF2dot ",inputp['VRF2ref'],inputp['VRF2dot']
	print "h, hratio, phi12 ",inputp["h"], inputp["hratio"],inputp['phi12']
	print "Bref, Rnom, rhonom ",inputp['Bref'], inputp['Bdot'], inputp['Rnom'], inputp['rhonom']
	print "gamma transition, rest mass, q ",inputp['gammatnom'], inputp['Erest'],inputp['q']
	print "selffields, coupling, Zovern", inputp['selffields'], inputp["coupling"], inputp["Zovern"]



def get_plotinfo(input_file="plotinfo.data"):

	f1 = open(input_file,'r')

	float_l = []
	for line in f1:
		if line[0] == ' ':
			#print line
			lspl = line.split("=")
			try:
				
				fl = float(lspl[1])
				float_l.append(fl)
			except:
				pass

	f1.close()


	input_struct = np.zeros(1, dtype=[('profilecount','i2'),('profilelength','i2'),('dtbin','f8'),('dEbin','f8'),('eperimage','f8'),('xat0','f8'),('yat0','f8'),('phi_s','f8')])


	input_struct['profilecount'] = float_l[0]
	input_struct['profilelength'] = float_l[1]
	input_struct['dtbin'] = float_l[2]
	input_struct['dEbin'] = float_l[3]
	input_struct['eperimage'] = float_l[4]
	input_struct['xat0'] = float_l[5]
	input_struct['yat0'] = float_l[6]
	input_struct['phi_s'] = float_l[10]

	return input_struct



def set_param(inputp):
	"""Calculate some parameters"""

	Erest = inputp['Erest'][0]
	gammat = inputp['gammatnom'][0]
	Rnom = inputp['Rnom'][0]
	E0 =  BtoEnom(inputp['Bref'][0], inputp['rhonom'][0], Erest, inputp['q'][0]) #total energy
	gamma0 = E0/Erest
	eta0 = 1/(gamma0**2) - 1/(gammat**2) 
	beta0 = (1 - (1/gamma0**2))**0.5
	omegarev0 = beta0*SPEED_OF_LIGHT/Rnom
	Vpeak = inputp['VRF1ref'][0]
	
	return E0, gamma0, eta0, beta0, omegarev0, Vpeak

def find_phi0(inputp):

	phi_a = np.linspace(0, 2*math.pi, 1000)
	vrf_a = vrf(phi_a, inputp['VRF1ref'][0],inputp['VRF2ref'][0],inputp['hratio'][0],inputp['phi12'][0])

	Rnom = inputp['Rnom'][0]
	Bdot = inputp['Bdot'][0]
	rhonom = inputp['rhonom'][0]
	target = 2*math.pi*Rnom*rhonom*Bdot

	if target == 0:
		phi0 = 0
	else:
		i = 0
		for v in vrf_a:
			if i > 0:
				if (v - target)*(vrf_a[i-1] - target) < 0:
					phi0 = 0.5*(phi_a[i] + phi_a[i-1])
					break
			i = i + 1

		#print phi0
		#plt.plot(phi_a, vrf_a, 'ko')
		#plt.axhline(y=target)
		#plt.show()
		#sys.exit()
	return phi0


def phi_u_fn(phi_u, phi_s, H):

	fn = abs(np.cos(phi_u) + phi_u*np.sin(phi_s) + np.cos(phi_s) - (math.pi - phi_s)*np.sin(phi_s))  
	
	return fn

def turning_point(phi_s, H=0, phi_guess=0):
	"""phi_u and pi-phi_s are the separatrix turning points of H=0. Find phi_u for a given phi_s"""

	#phi_a = np.linspace(-2*math.pi, 2*math.pi, 1000)
	
	#plt.plot(phi_a, phi_u_fn(phi_a, phi_s, H))	
	#plt.plot(phi_a, phi_u_fn2(phi_a, phi_s, H))
	#phi_guess= phi_s
	res = fmin(phi_u_fn,phi_guess, args = (phi_s, H))
	#print "min found at ", res
	#plt.axvline(x = res['x'][0])
	#plt.show()
	#sys.exit()

	phi_u = res[0]

	#minimize(phi_u_fn,0, args = (phi_s))
	#phi_u = res['x'][0]

	return phi_u


def turning_point2(phi_s, a1, a2, N):

	phase_test = np.linspace(-2*math.pi, 2*math.pi, N)
	H_test = H_fn(phase_test, 0, phi_s, a1, a2)

	i = 0
	sw = True
	phi_u = None
	for H, ph in zip(H_test	, phase_test):
		if i > 0:
			if ph > math.pi - phi_s and sw:
				H_sep = H
				sw = False
			if not sw:
				if H > H_sep:
					phi_u = ph
					break

		i = i + 1

	return phi_u




def separatrix(phi_s, eta0, T, h, npts, a1, a2):
	"""calculate separatrix in dp, phi space. 
		phi_s: synchchronous phase
		eta0: phase slip
		T: RF period
		h: harmonic number
		bamp: bucket amp
		"""

	if phi_s == 0.0:
		phi_u = -math.pi
	else:
		phi_u =  turning_point(phi_s, phi_guess = phi_s)
		#phi_u =  turning_point2(phi_s, a1, a2, 1000)	

	#phi_a: array of phase values at which to evaluate separatrix
	if phi_u > 0:
		if phi_s > 0.5*math.pi:
			phi_a = np.linspace((math.pi - phi_s), phi_u, npts)
		else:
			phi_a = np.linspace(-(math.pi - phi_s), phi_u, npts)
	else:
		phi_a = np.linspace(phi_u, (math.pi - phi_s), npts)

	bamp = (abs(a2/a1))**0.5

	signc = -np.sign(2*np.cos(phi_s) - (math.pi - phi_s-phi_s)*np.sin(phi_s))

	dp_sep = []
	phi_real = []
	for phi in phi_a:

		dpsq = -signc*(bamp**2)*(np.cos(phi) + np.cos(phi_s) - (math.pi - phi-phi_s)*np.sin(phi_s))

		if dpsq > 0:
			dp = dpsq**0.5
			dp_sep.append(dp)
			phi_real.append(phi)

	dp_sep = np.array(dp_sep)
	phi_real_a = np.array(phi_real)
	t_phi_ns = 1e9*T*(phi_real_a)/(2*math.pi) 


	return dp_sep, t_phi_ns





def contour(phi_s, eta0, T, h, npts, a1, a2, Hamil):
	"""calculate any contour in dp, phi space on given Hamiltonian. 
		nu_s: synch tune
		phi_s: synchchronous phase
		eta0: phase slip
		T: RF period
		h: harmonic number
		a1, a2: coef in Hamiltonian for dp and phase terms
		"""


	if phi_s == 0.0:
		phi_u = -math.pi
	else:
		phi_l =  turning_point(phi_s, H=Hamil, phi_guess = -0.8)
		phi_u =  turning_point(phi_s, H=Hamil, phi_guess = phi_s)
		

	#phi_u = phi_u + 0.5*math.pi
	print "phi_l, phi_u, H ",phi_l, phi_u, Hamil
	print "math.pi - phi_s ",-(math.pi - phi_s)

	#phi_a: array of phase values at which to evaluate separatrix
	 
	if phi_u > 0:
		phi_a = np.linspace(phi_l, phi_u, npts)
	else:
		phi_a = np.linspace(phi_u, phi_l, npts)

	print "Hamil ",Hamil


	ratio = a2/a1
	signc = -np.sign(2*np.cos(phi_s) - (math.pi - phi_s-phi_s)*np.sin(phi_s))
	print "a2/a2 ",a2/a1

	dp_cont = []
	phi_real = []
	for phi in phi_a:
		#original form for opposite sign eta
		#dpsq = 0.5*(bh**2)*(np.cos(phi) + np.cos(phi_s) - (math.pi - phi - phi_s)*np.sin(phi_s))

		dpsq = -((Hamil/a1) - ratio*signc*(np.cos(phi) - np.cos(phi_s) + (phi-phi_s)*np.sin(phi_s))) 

		if dpsq > 0:
			dp = dpsq**0.5
			dp_cont.append(dp)
			phi_real.append(phi)

	print "number of points on contour ",len(phi_real)
	dp_cont = np.array(dp_cont)

	#plt.plot(phi_a, dp_sep)
	#plt.ylim(ymin=0)
	#plt.show()
	#sys.exit()

	phi_real_a = np.array(phi_real)
	t_phi_ns = 1e9*T*(phi_real_a)/(2*math.pi) 


	return dp_cont, t_phi_ns



def calc_tof_bm_data(lower_lim, points_per_turn, nturns_fit):
	"""Calculate TOF by fitting the time of successive peaks"""
	
	upper_lim_fit = lower_lim + nturns_fit*tof_guess
	
	find_hminus = False
	_, pmx, _, i_peak_l, _, _ = adata.find_peaks(tdat_chf, data_chf, points_per_turn, find_hminus, x_low = lower_lim, x_high = upper_lim_fit, find_maxima = True)
	
	tp = np.arange(nturns_fit)
	npoly = 1
	fitpmx = np.polyfit(tp, pmx, npoly)
	tof_fit_mean = fitpmx[npoly - 1]
	tof_fit_ns = 1e9*tof_fit_mean
	polyfitpmx = np.poly1d(fitpmx)
	tof_fit_a = polyfitpmx(tp)

	tof_diff = pmx - polyfitpmx(tp)
	#fit_tof_diff = np.polyfit(tp, tof_diff, 2)
	#polyfitdiff = np.poly1d(fit_tof_diff)
	#fitdiff_a = polyfitdiff(tp) - polyfitdiff[0]
	
	#ampt = max(tof_diff) - min(tof_diff)
	#timeosc = [ampt*np.sin(0.25*2*np.pi*Qs*i) - 0.5*ampt for i in np.arange(nturns_fit)]
	#timeosc = timeosc - timeosc[0]
	
	if show_time_fit:
		pmx_a = np.array(pmx)
		plt.subplot(211)
		plt.plot(1e6*pmx_a,'k-',linewidth=2,label='data')
		plt.plot(1e6*polyfitpmx(tp),'r',label='fit')
		plt.ylabel('turn-by-turn times [us]')
		plt.title('slope = '+str(tof_fit_ns)+ ' ns')
		plt.legend(loc = 'upper left')
		plt.subplot(212)
		plt.plot(1e6*(pmx_a - polyfitpmx(tp)))
		#plt.plot(1e6*polyfitdiff(tp) ,'r')
		plt.xlabel('turn')
		plt.ylabel('difference [us]')
		plt.savefig('time_fit')
		plt.show()
	
	#print (tof_fit+1e9*fitdiff_a[0])/time_interval , (tof_fit+1e9*fitdiff_a[-1])/time_interval
			
	return tof_fit_mean, i_peak_l, tof_fit_a
		
	
#find zero crossing times in rf data
def calc_tof_rf_waveform(tdat_rf, sig_rf, nturns_fit, istart):
		
	time_interval = tdat_rf[1] - tdat_rf[0]

	zero_cross_ind = []
	t_zero_cross = []
	plot_zero_calc = False
	for ind in range(istart, len(tdat_rf)):
		if sig_rf[ind-1] > 0.0 and sig_rf[ind] < 0.0:
			slope = (sig_rf[ind] - sig_rf[ind-1])
			offset = sig_rf[ind] - slope
			zerot_fac = -offset/slope
			
			if plot_zero_calc:
				xp = np.linspace(0,1,10)
				yp = slope*xp + offset
				plt.plot([0,1],[sig_rf[ind-1],sig_rf[ind]],'ko')
				plt.plot(xp, yp, 'r-')
				plt.axvline(x=zerot_fac)	
				plt.show()

			t_zero_cross.append(tdat_rf[ind-1] + zerot_fac*time_interval)
			zero_cross_ind.append(ind)

		if len(zero_cross_ind) == nturns_fit:
			break
	
	t_zero_cross = np.array(t_zero_cross)
	
	zero_cross_sep = [zero_cross_ind[i] - zero_cross_ind[i-1] for i in range(1, len(zero_cross_ind))]
	
	sig_zero_cross_sep = [sig_rf[zero_cross_ind[i]] - sig_rf[zero_cross_ind[i]-1] for i in range(1, len(zero_cross_ind))]
	t_zero_cross_sep = [t_zero_cross[i] - t_zero_cross[i-1] for i in range(1, len(t_zero_cross))]
	tof_mean = sum(t_zero_cross_sep)/len(t_zero_cross_sep)


	return tof_mean, zero_cross_ind, t_zero_cross 

