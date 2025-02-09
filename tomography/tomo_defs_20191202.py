from __future__ import division
import math 
import sys
import numpy as np
import pylab as plt
import os
from scipy.optimize import minimize, fmin
from scipy.signal import butter, lfilter, freqz, filtfilt
from scipy.integrate import quad
from scipy import ndimage
from subprocess import call
import matplotlib.cm as cm

PROTON_MASS = 938.27203e6
SPEED_OF_LIGHT = 299792458

BtoEnom = lambda B, rho, Erest, q: ((q*B*rho*SPEED_OF_LIGHT)**2 + (Erest**2))**0.5
vrf = lambda phi, vrf1, vrf2, hratio, phi12: vrf1*np.sin(phi) + vrf2*np.sin(hratio*(phi-phi12))

#Hamiltonian
H_fn = lambda phi, dp, phi_s, a1, a2: a1*(dp**2) + a2*(np.cos(phi) - np.cos(phi_s) + (phi - phi_s)*np.sin(phi_s))
H_fn_gen_pot = lambda n, phi, phi_s, coef, ang: coef*(np.cos(n*phi + ang) - np.cos(n*phi_s + ang) + (phi - phi_s)*np.sin(n*phi_s + ang))

H_fn_gen = lambda phi, dp, phi_s, a1, a2, n_l, coef_l, ang_l: a1*(dp**2) + sum([H_fn_gen_pot(n, phi, phi_s, a2*coef, ang) for n, coef, ang in zip(n_l, coef_l, ang_l)])


#synchrotron tune
Qs_fn = lambda phi0, V0, h, beta, E0, eta: np.sqrt(np.abs(eta*np.cos(phi0)) * h* V0 / (2 * math.pi * beta**2 * E0))
#dEbin 
dEbin_fn = lambda phi0, V0, h, beta, E0, eta, dtbin, omega0: beta*np.sqrt(E0*V0*np.cos(phi0)/(2*math.pi*h*abs(eta)))*dtbin*omega0

#Functions used to find synchrotron period vs amp
H0_minus_pot = lambda phi, phi_s, coef, H0: H0 -  coef*(np.cos(phi) - np.cos(phi_s) + (phi - phi_s)*np.sin(phi_s))

H0_minus_pot_gen = lambda phi, phi_s, a2, n_l, coef_l, ang_l, H0: H0 - sum([H_fn_gen_pot(n, phi, phi_s, a2*coef, ang) for n, coef, ang in zip(n_l, coef_l, ang_l)])

T_int_fn = lambda phi, phi_s, coef, H0, int_coef: 1/(-int_coef*H0_minus_pot(phi, phi_s, coef,H0))**0.5
J_int_fn = lambda phi, phi_s, coef, H0, int_coef: (-(4/int_coef)*H0_minus_pot(phi, phi_s, coef,H0))**0.5

T_int_fn_gen = lambda phi, phi_s, a2, n_l, coef_l, ang_l, H0, int_coef: 1/(-int_coef*H0_minus_pot_gen(phi, phi_s, a2, n_l, coef_l, ang_l, H0))**0.5


def read_scope(path, multi_channel = False, tfactor=1):
	"""Read scope data. The script deals with two specific scope formats - 
	that used in November 2013, and in March 2014. The two are distinguished by looking
	for the string DP04104 that appears in the later case"""

	print 'reading ', path

	f_data = open(path, 'r')
	
	tdat = []
	ydat = []

	n = 0
	for line in f_data:
    
		lspl = line.split(',')

		if n == 0:
			if lspl[1][:7] == 'DPO4104':
				nskip = 21
				xi = 0
				yi = 1
			elif len(lspl) == 2:
				nskip = 0
				xi = 0
				yi = 1
			else:
				nskip = 6 + 1
				xi = 3
				yi = 4

		if n >= nskip:
			if line !='\r\n':

				tdat.append(tfactor*float(lspl[xi]))

				if multi_channel:
					ydat.append([float(l) for l in lspl[yi:]])
				else:
					ydat.append(float(lspl[yi]))
			else:
				break

		n += 1

	f_data.close()
   
 	return tdat, ydat

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
	res = fmin(phi_u_fn, phi_guess, args = (phi_s, H))
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


def phi_cross(phi_g, phi_s, a2, H, phi_upper, n_l, coef_l, ang_l, phase_mid=None):

	"""Solution constrained to be less than phi_upper """
	#fn = abs(H_fn_gen_pot(1, phi_g, phi_s, a2, 0) - H)

	fn = abs(H_fn_gen(phi_g, 0, phi_s, 0, a2, n_l, coef_l, ang_l) - H)

	if phase_mid != None:
		phi_lower = phase_mid
	else:
		phi_lower = phi_s

	if phi_g < phi_lower:
		fn = 1e6
	elif phi_g > phi_upper:
		fn = 1e6

	return fn

def phase_map_left_to_right(phase_left, phi_s, a1, a2, phase_sep_max, n_l, coef_l, ang_l, phase_mid = None):
	"""given phase axis crossing points less than phi_s, find crossing points on the same Hamiltonian contours on the right"""

	
	#H_phase_old = H_fn(phase_left, 0, phi_s, a1, a2)

	H_phase = H_fn_gen(phase_left, 0, phi_s, a1, a2, n_l, coef_l, ang_l)

	phase_out = []
	phase_left_sol = []
	for H0, phase_left in zip(H_phase, phase_left):
		res = fmin(phi_cross, phi_s, args = (phi_s, a2, H0, phase_sep_max, n_l, coef_l, ang_l, phase_mid))
		H_tol = abs(H0 - H_fn_gen(res[0], 0, phi_s, a1, a2, n_l, coef_l, ang_l))

		if H_tol < 1e-3:
			phase_left_sol.append(phase_left)
			phase_out.append(res[0])

	phase_out = np.array(phase_out)
	phase_left_sol = np.array(phase_left_sol)

	return phase_out, phase_left_sol

def phase_extrema(phi_s):
	"""For a given phi_s find max and min phases along separatrix"""
	#calculate phase extrema along separtrix
	phi_l1 = math.pi - phi_s

	if phi_s == 0.0:
		phi_l2 = -math.pi
	else:
		#guess turning point 
		#may break down if phi_s is close to 90 degrees
		if phi_s < 50*math.pi/180:
			phi_g = 0
		elif phi_s >= 50*math.pi/180 and phi_s < 0.5*math.pi:
			phi_g = 80*math.pi/180
		else:
			phi_g = math.pi

		phi_l2 =  turning_point(phi_s, phi_guess = phi_g)
		
	return phi_l1, phi_l2


def synchrotron_period_vs_amplitude(phase_left, phi_s, a1, a2, phase_sep_max, harmonic, omegarev, eta, n_l, coef_l, ang_l):
	"""given phase axis crossing points less than phi_s, find crossing points on the same Hamiltonian contours on the right
		Integrate fixed Hamiltonian minus the potential term between these phase extremes
		A similar integral is evaluated to find the action J"""

	phase_right, phase_left = phase_map_left_to_right(phase_left, phi_s, a1, a2, phase_sep_max, n_l, coef_l, ang_l)

	#H_phase = H_fn(phase_left, 0, phi_s, a1, a2)
	H_phase = H_fn_gen(phase_left, 0, phi_s, a1, a2, n_l, coef_l, ang_l)

	int_coef = abs(2*harmonic*omegarev*eta)

	Ts_amp_l = []
	J_l = []
	phase_sol = []
	for ph_left, ph_right, H1 in zip(phase_left, phase_right, H_phase):

			#print "ph_left, ph_right, H1 ",ph_left, ph_right, H1
			#ensure extreme points start just inside crossing point avoid singularity in integral
			tol = 1e-9
			diff_lim = -1e-11

			if len(n_l) == 1:
				while  H0_minus_pot(ph_left, phi_s, a2, H1) >= 0:
					ph_left = ph_left + tol
					
				while  H0_minus_pot(ph_left, phi_s, a2, H1) < diff_lim:
					ph_left = ph_left - tol
			
				while  H0_minus_pot(ph_right, phi_s, a2, H1) >= 0:
					ph_right = ph_right - tol

				while  H0_minus_pot(ph_right, phi_s, a2, H1) < diff_lim:
					ph_right = ph_right + tol

			else:
				#print H0_minus_pot_gen(ph_left, phi_s, a2, n_l, coef_l, ang_l, H1)
				
				il = 0 
				while H0_minus_pot_gen(ph_left, phi_s, a2, n_l, coef_l, ang_l, H1) >= 0 and il < 10:
					ph_left = ph_left + tol
					il = il + 1
				#print H0_minus_pot_gen(ph_left, phi_s, a2, n_l, coef_l, ang_l, H1)
				
				#while H0_minus_pot_gen(ph_left, phi_s, a2, n_l, coef_l, ang_l, H1) < diff_lim:
				#	ph_left = ph_left + tol

				#print H0_minus_pot_gen(ph_left, phi_s, a2, n_l, coef_l, ang_l, H1)
				#sys.exit()
				il = 0 
				while H0_minus_pot_gen(ph_right, phi_s, a2, n_l, coef_l, ang_l, H1) >= 0 and il < 10:
					ph_left = ph_left - tol
					il = il + 1
				#while H0_minus_pot_gen(ph_right, phi_s, a2, n_l, coef_l, ang_l, H1) < diff_lim:
				#	ph_left = ph_left + tol
		
			
			#synchrotron period
			if len(n_l) == 1:
				integral = quad(T_int_fn, ph_left, ph_right, args=(phi_s, a2, H1, int_coef))[0]
			else:
				integral = quad(T_int_fn_gen, ph_left, ph_right, args=(phi_s, a2, n_l, coef_l, ang_l, H1, int_coef))[0]
			period_s = 2*integral
			
			print H0_minus_pot_gen
			
			ph_a = np.linspace(ph_left, ph_right, 100)
			
			
			#plt.plot(ph_a, H0_minus_pot_gen(ph_a, phi_s, a2, n_l, coef_l, ang_l, H1))
			#plt.plot(ph_a, T_int_fn_gen(ph_a, phi_s, a2, n_l, coef_l, ang_l, H1, int_coef))
			#plt.show()
			
			
			#print "period_s ",period_s

			#action
			J = (1/(2*math.pi))*quad(J_int_fn, ph_left, ph_right, args=(phi_s, a2, H1, int_coef))[0]

			Ts_calc = omegarev*period_s/(2*math.pi)
			print Ts_calc
			#sys.exit()

			try:
				b = Ts_calc + 1
				Ts_amp_l.append(Ts_calc)
				phase_sol.append(ph_right)
			except:
				pass
				
			J_l.append(J)

	Ts_amp_a = np.array(Ts_amp_l)
	J_a = np.array(J_l)
	phase_sol = np.array(phase_sol)

	return Ts_amp_a, J_a, phase_sol

def separatrix(phi_s, eta0, h, npts, a1, a2, T=None):
	"""calculate separatrix in dp, phi space. 
		phi_s: synchchronous phase
		eta0: phase slip
		T: RF period (optional)
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
	phi_sep = np.array(phi_real)

	if T != None:
		t_phi_ns = 1e9*T*(phi_sep)/(2*math.pi) 
	else:
		t_phi_ns = None

	return dp_sep, phi_sep, t_phi_ns




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
def calc_tof_rf_waveform(tdat_rf, sig_rf, nturns, start_time):
	"""Calculate time of zero-crossing of RF signal. Calulate mean TOF based on this
       tdat_rf: time data from RF waveform
		sig_rf:  RF signal
        nturns: number of turns at which to calculate zero-crossing
        start_time: desired start time. The script will look for the first zero-crossing later than start_time"""


	i = 0
	i_st = None
	for t in tdat_rf:
		if t > start_time:
			i_st = i
			break
		i = i + 1

	if i_st == None:
		print "time ",start_time, "not found"
		sys.exit()

	time_interval = tdat_rf[1] - tdat_rf[0]

	zero_cross_ind = []
	t_zero_cross = []
	plot_zero_calc = False
	for ind in range(i_st, len(tdat_rf)):

		if sig_rf[ind-1] >= 0.0 and sig_rf[ind] < 0.0:

			#skip spurious double crossing of zero
			if len(zero_cross_ind) > 0:
				if ind - zero_cross_ind[-1] <= 2:
					continue
				
			if sig_rf[ind-1] != 0.0:
				slope = (sig_rf[ind] - sig_rf[ind-1])
				offset = sig_rf[ind] - slope
				zerot_fac = -offset/slope
			else:
				zerot_fac = 0.0
			
			if plot_zero_calc:
				xp = np.linspace(0,1,10)
				yp = slope*xp + offset
				plt.plot([0,1],[sig_rf[ind-1],sig_rf[ind]],'ko')
				plt.plot(xp, yp, 'r-')
				plt.axvline(x=zerot_fac)	
				plt.show()

			#t_zero_cross.append(tdat_rf[ind-1])
			t_zero_cross.append(tdat_rf[ind-1] + zerot_fac*time_interval)
			zero_cross_ind.append(ind)


		if len(zero_cross_ind) == nturns:
			break
	
	t_zero_cross = np.array(t_zero_cross)
	
	zero_cross_sep = [zero_cross_ind[i] - zero_cross_ind[i-1] for i in range(1, len(zero_cross_ind))]
	
	sig_zero_cross_sep = [sig_rf[zero_cross_ind[i]] - sig_rf[zero_cross_ind[i]-1] for i in range(1, len(zero_cross_ind))]
	tof_l = [t_zero_cross[i] - t_zero_cross[i-1] for i in range(1, len(t_zero_cross))]
	
	if tof_l == []:
		print "zero crossing times in RF waveform not found! Please check the Rf waveform."
		sys.exit()
	else:
		tof_mean = sum(tof_l)/len(tof_l)

	#plt.plot(tdat_rf, sig_rf)
	#plt.plot(t_zero_cross, np.zeros(len(t_zero_cross)),'ro')
	#plt.xlim(1e-3,0.00101)
	#plt.show()
	#sys.exit()

	return tof_mean, tof_l, zero_cross_ind, t_zero_cross


def show_data(time_indices, time_data, signal_data, tdat_rf, sig_rf, nturns, shift=0):

	tdat_rf_ms = 1e3*np.array(tdat_rf) 
	sig_rf_sc = 0.01*np.array(sig_rf)
	time_data_ms  = 1e3*np.array(time_data)
	time_interval = time_data[1] - time_data[0]

	tpad = 2e-4

	nt = 10

	if nt > len(time_indices):
		nt = len(time_indices) - 1
	tfirst = tdat_rf_ms[time_indices[0]] -tpad
	tlast =  tdat_rf_ms[time_indices[nt]] + tpad

	plt.subplot(211)
	plt.plot(tdat_rf_ms, sig_rf_sc, 'r')
	plt.plot(time_data_ms, signal_data,'k')
	
	#plt.axvline(x=t_hminus, color='r')
	
	
	#plt.xlim(0.5,1.2)
	#plt.xlim(0, time_data_ms[-1])

	#plt.ylim(-0.2, 0.2)
	plt.ylim(-0.05, 0.05)
	plt.ylabel('bunch monitor (all)')
	plt.axvline(x=tdat_rf_ms[time_indices[0]])
	

	plt.subplot(212)
	plt.plot(time_data_ms + 1e3*time_interval*shift, signal_data,'k')
	plt.plot(tdat_rf_ms , sig_rf_sc, 'r')
	plt.xlim(tfirst, tlast)
	plt.ylim(-0.1, 0.1)

	plt.axvline(x=tdat_rf_ms[time_indices[0]])
	plt.ylabel('bunch monitor (zoom)')
	plt.xlabel('time [ms]')
	plt.savefig('bunchmonitor_data')
	plt.show()

	sys.exit()


def list_files(path):
	"""list all files in data directory"""
	def sorted_ls(path):
		mtime = lambda f: os.stat(os.path.join(path, f)).st_mtime
		return list(sorted(os.listdir(path), key=mtime))

	def first_chars(x):
		return(x[:40])

	sort_by_name = True
	if sort_by_name:
		files = sorted(os.listdir(path), key = first_chars)
	else:
		files = sorted_ls(path)

	return files


def list_csv_files(path):
	"""list csv files in data directory"""

	files = list_files(path)

	filter_string = '' #look for csv files with this string in file name
	csv_files = []
	icsv = 0 
	for i1, f1 in enumerate(files):
		if 'csv' in f1 and filter_string in f1:
				
			csv_files.append(f1)
			icsv = icsv + 1

	return csv_files

def read_bmonitor_data(path, multi_channel = False):
	"""select and read bunch monitor data"""

	csv_files = list_csv_files(path)
	print "data directory ",path
	print "the following csv files were found"
	#files =  os.listdir(dirname)
	for icsv, csvf in enumerate(csv_files):
		print "index ",icsv," file ", csvf

	rawin1 =  raw_input("select index of data file  ")
	try:
		i_file_sel = int(rawin1)
	except:
		print "input is not integer"
		sys.exit()

	fname_sig = csv_files[i_file_sel]
	time_data, signal_data = read_scope(path+fname_sig, multi_channel = multi_channel)
	
	return time_data, signal_data, fname_sig

def read_rf_data(path, fname_bm_sig=None):
	"""select and read rf waveform data
		if supplied with fname_bm_sig, it will try to find the rf file that pairs with the selected bunch monitor file swapping "rf" for "wbeam" """

	csv_files = list_csv_files(path)

	if fname_bm_sig != None:
		
		rf_file_guess = fname_bm_sig.replace("wbeam","rf")

		if rf_file_guess in csv_files:
			i_rf_file_sel = csv_files.index(rf_file_guess)
			print "found following matching rf file ",rf_file_guess	
	else:
		rawin2 = raw_input("select RF data file (enter to select "+ str(i_bmfile_sel + 1) + ")  ")
			
		try:
			i_rf_file_sel = int(rawin2)
		except:
			if rawin2 == '':
				i_rf_file_sel = i_bmfile_sel + 1
			else:
				i_rf_file_sel = None
				print "no valid RF file selected"
					
	if i_rf_file_sel != None:
		fname_rf = csv_files[i_rf_file_sel]
		tdat_rf, sig_rf = read_scope(path, fname_rf)

	return tdat_rf, sig_rf


def read_phis_file_fn(path, fname_bm_sig=None):
	"""Read text file specifying phi_s settings"""

	files = list_files(path)

	def read_file(path, fname):
		"""Read phis file"""
		ff = open(path+fname,'r')
		t_phis_file = []
		phis_file = []
		for lraw in ff:
			
			line = lraw.split()
			t_phis_file.append(float(line[0]))
			phis_file.append(float(line[1]))
		ff.close()

		return t_phis_file, phis_file

	if fname_bm_sig[:4] == 'case':
			fname_phis = 'case'+fname_bm_sig[4]+'_phis.txt'
			print "look for phis file named ",fname_phis
			if fname_phis in files:
				t_phis_file, phis_file = read_file(path,fname_phis)

				print "t_phis_file ",t_phis_file
				print "phis_file ",phis_file
		
			else:
				print "phis_file not found"
				sys.exit()
	else:
		t_phis_file = [0]
		phis_file = [0]

	return t_phis_file, phis_file

def read_offset_time_fn(path, fname_bm_sig=None):

	case_spl = fname_bm_sig.split("_")
	case_num = case_spl[0]
		
	if case_spl[1] != 'wbeam.csv':
		offset_file_name = path+str(case_num)+"_"+case_spl[1]+"_toffset.txt"
	else:
		offset_file_name = path+str(case_num)+"_toffset.txt"
		
	print "look for file ",offset_file_name
	try:
		fo = open(offset_file_name,'r')
		#fo.readline()
		offset_file_ms = 1e3*float(fo.readline())
	except:
		print "offset file not found"
		pass
		
	return offset_file_ms

def read_data_files(path, col_bm_rf, read_rf_waveform = True, read_phis_file = True, read_file_offset=True, multi_channel = False):
	"""Read bunch monitor data and, if chosen, the rf waveform, phis file and time offset files
	   If multi_channel = True, then more than one channel of data is read from the scope
	   col_bm_rf: column index of bunch monitor and rf data"""

	if multi_channel:
		time_data, waveform_data, fname_sig = read_bmonitor_data(path, multi_channel = True)

	
		waveform_data = np.array(waveform_data)
		waveform_data_tr = np.transpose(waveform_data)
		
		print "waveform_data_tr ",waveform_data_tr.shape
		

		tdat_rf = time_data
		sig_rf = waveform_data_tr[col_bm_rf[1]]
		
		#sig_rf = waveform_data_tr[2]
		#signal_data = waveform_data_tr[1] #March, April 2019
		
		#signal_data = waveform_data_tr[2] #June 20, 2019, first dataset
		#signal_data = waveform_data_tr[1] #June 20, 2019, other datasets
		
		signal_data = waveform_data_tr[col_bm_rf[0]]
		#print "waveform_data shape ",waveform_data.shape
		#print "time step ",time_data[1] - time_data[0]

		
		#plt.subplot(211)
		#plt.plot(time_data, waveform_data_tr[0],'k-')
		#plt.subplot(212)
		#plt.plot(time_data, waveform_data_tr[1],'r-')
		#plt.show()

		#plt.subplot(211)
		#plt.plot(time_data, waveform_data_tr[0],'k-')
		#plt.xlim(2e-4,2.1e-4)
		#plt.subplot(212)
		#plt.plot(time_data, waveform_data_tr[1],'r-')
		#plt.xlim(2e-4,2.1e-4)
		#plt.show()

		#sys.exit()
	else:
		time_data, signal_data, fname_sig = read_bmonitor_data(path)

		#read rf waveform
		if read_rf_waveform:
			tdat_rf, sig_rf  = read_rf_data(path, fname_bm_sig=fname_sig)

	#read phis file which contains the phis settings for data
	t_phis_file = None
	if read_phis_file:
		t_phis_file, phis_file = read_phis_file_fn(path, fname_bm_sig=fname_sig)

	#read time offset between rf and beam data 
	if read_file_offset:
		offset_file_ms = read_offset_time_fn(path, fname_bm_sig=fname_sig)
	else:
		offset_file_ms = None

	return fname_sig, time_data, signal_data, tdat_rf, sig_rf, t_phis_file, phis_file, offset_file_ms



#svk file from Tom. rf burst starts -0.5ms, i.e. before injection.
#format is index, time[s], kinetic energy [eV], rf freq [Hz], mean radius [m]
def read_svk_file(svkfilepath):
	#svkfile = "../sep18_exp/svk20.dat"
	ff = open(svkfilepath, 'r')
	
	time_limit = 17e-3
	
	time_svk = []
	ke_svk = []
	rf_freq_svk = []
	for indf ,line in enumerate(ff):
		
		lspl = line.split()
		
		time_svk.append(float(lspl[1]))
		ke_svk.append(float(lspl[2]))
		rf_freq_svk.append(float(lspl[3]))
		
		if time_svk[-1] > time_limit:
			break
	
	return time_svk, ke_svk, rf_freq_svk
	
	
	
def ke_lookup(time, svkfilepath, toffset = 0):
	""" Find kinetic energy at specified time in svk file """
	i = 0 
	ke = None
	
	time_svk, ke_svk, rf_freq_svk = read_svk_file(svkfilepath)
	
	if time_svk[0] > time - toffset:
		print "time less than minimum in table : ",time-toffset
		print "use initial energy in table"
		return ke_svk[0]
		
	for t in time_svk:
		if t > time - toffset:
			ke = ke_svk[i]	
			break	
		i = i + 1
	
	return ke
	
def f_rf_lookup(time, svkfilepath, svk_data=None, toffset = 0):
	""" Find rf frequency at specified time in svk file
		svk_data list of lists in format [time_svk, ke_svk, f_rf_svk]  """
	i = 0 
	
	if svk_data == None:
		time_svk, ke_svk, rf_freq_svk = read_svk_file(svkfilepath)
	else:
		time_svk = svk_data[0]
		rf_freq_svk = svk_data[2]

	rf_freq = None
	for t in time_svk:
		if t > time - toffset:
			rf_freq = rf_freq_svk[i]	
			break	
		i = i + 1
	
	return rf_freq


def time_offset_calc(path, ttime_a, t_phis_file, svkfilepath, svk_data, fname_bm_sig, write_offset_file, plot_time_offset, fit_freq = False):
	"""Find the time when the RF frequency crosses the ideal frequency at 11 MeV (given by Uesugi-san's table)
       Select the crossing point where the RF frequency is increasing, not decreasing.
 
       ttime_a: time at the start of each turn
	   t_phis_file: phi_s programme
	   svk_data: svk table"""


	tof_a = [ttime_a[i] - ttime_a[i-1] for i in range(1, len(ttime_a))]
	freq_a = np.array([1/t for t in tof_a])
	freq_MHz_a = 1e-6*freq_a

	#read svk table (from Uesugi-san)
	time_svk, ke_svk, rf_freq_svk = svk_data
	fzero_ideal = f_rf_lookup(0, svkfilepath, svk_data=svk_data) #frequency at zero time in table
	f_zero_ideal_MHz = 1e-6*fzero_ideal

	if fit_freq:
		#polynomial fit to frequency data
		npfit = 3
		nturns_fit = min(2000, len(ttime_a))
	
		ffit_coef = np.polyfit(ttime_a[1:nturns_fit], freq_MHz_a[0:nturns_fit-1], npfit)
		ffit = np.poly1d(ffit_coef)
		
		#find offset between time where frequency reaches fzero_ideal and t=0
		offset =  (f_zero_ideal_MHz - ffit_coef[-1])/ffit_coef[npfit-1]
		
	else:
		i1 = 1
		
		for f in freq_MHz_a[1:]:
			if (f - f_zero_ideal_MHz) > 0 and (freq_MHz_a[i1-1] - f_zero_ideal_MHz) < 0 :
				offset = ttime_a[i1]
				break
			i1 = i1 + 1
		
	offset_ms = 1e3*offset	
	
	if write_offset_file:

		case_spl = fname_bm_sig.split("_")

		case_num = case_spl[0]
		print "case_spl ",case_spl
		print "case_num ", case_num

		if case_spl[1] != 'wbeam.csv':
			offset_file_name = path+str(case_num)+"_"+case_spl[1]+"_toffset.txt"
		else:
			offset_file_name = path+str(case_num)+"_toffset.txt"

		print "offset_file_name ",offset_file_name
		
		fo = file(offset_file_name,'w')
		#fo = file(dirname+offset_file,'w')
		print >>fo, offset


	t_rf_table = np.linspace(0, ttime_a[-1], 100)
	if plot_time_offset:
		plt.plot(t_rf_table, [1e-6*f_rf_lookup(t, svkfilepath, svk_data=svk_data) for t in t_rf_table],'k-',label='lookup table')
			
		plt.plot(ttime_a[1:], freq_MHz_a, 'b.',label='RF data (1/ToF)')
		
		if fit_freq:
			plt.plot(ttime_a[:nturns_fit], ffit(ttime_a[:nturns_fit]),'r-', label='fit')

		plt.axhline(y=f_zero_ideal_MHz, color='gray',linestyle='--',linewidth=2,label='ideal inj freq')
		#plt.axvline(x=ttime_a[i1], color='gray', linestyle='--')
		#plt.axvline(x=ttime_a[i1+2500], color='gray',linestyle='--')

		plt.xlabel('time [s]')
		plt.ylabel('frequency [MHz]')
		plt.axvline(x= offset, color='r',label='calc. offset')
			
		plt.legend(loc = 'upper right')
		#plt.ylim(ymin = 1.575)
		#plt.ylim(ymax = 5)
		plt.title('offset (ms) '+str(offset_ms)[:5])
		plt.savefig('frequency_offset')
		plt.show()

	
	return offset_ms


def find_symmetry_point_shape(data):
	"""Look for point where data on either side is most similar"""

	ld = len(data)
	ldh = int(0.5*ld)

	if ld%2 == 1:
		oddcorr = 1
	else:
		oddcorr = 0

	
	chisq_l = []
	for d in range(ldh):
		#indh1 = range(d, ldh+d)
		#indh2 = range(ldh+d, ld+d)
		indh1 = range(ldh+d, ld+d - oddcorr)
		indh2 = range(d, ldh+d)
		half1 = data.take(indh1, mode='wrap')
		half2 = data.take(indh2, mode='wrap')

		half2rev = half2[::-1]		
		halfdiff = half1 - half2rev

		chisq = np.sum(halfdiff**2)
		chisq_l.append(chisq)
					
	chisq_min = min(chisq_l)
	isym = chisq_l.index(chisq_min)
	
	return isym, chisq_min


def find_symmetry_point_integral(data, plot_data = False):
	"""Look for point where data on either side is most similar"""

	ld = len(data)
	ldh = int(0.5*ld)

	if ld%2 == 1:
		oddcorr = 1
	else:
		oddcorr = 0

	chisq_l = []
	for d in range(ldh):
		#indh1 = range(d, ldh+d)
		#indh2 = range(ldh+d, ld+d)
		indh1 = range(ldh+d, ld+d - oddcorr)
		indh2 = range(d, ldh+d)
		half1 = data.take(indh1, mode='wrap')
		half2 = data.take(indh2, mode='wrap')

		half1_sum = np.sum(half1)
		half2_sum = np.sum(half2)
		halfdiff = abs(half1_sum - half2_sum)

		chisq = halfdiff

		#if plot_data:
		#	plt.clf()
		#	print "d, halfdiff ",d, halfdiff	
		#	plt.plot(half1)
		#	plt.plot(half2)
		#	plt.show()


		chisq_l.append(chisq)
					
	chisq_min = min(chisq_l)
	isym = chisq_l.index(chisq_min)

	#check if symmetry point is at foot of profile
	if isym < ldh:
		if data[isym] < data[isym + ldh]:
			isym = isym + ldh 
			#print "adjusted isym ",isym


	
	if plot_data:
		print "isym ",isym

		plt.clf()
		plt.plot(data)
		plt.show()
		

		plt.plot(chisq_l)
		plt.show()
		sys.exit()
	
	
	return isym, chisq_min


def calc_sync_point(prof_data_a, nwin, sync_bin, show_int_prof = False, show_sym_plots = False):
	"""Look for time point in data that may correspond to the synchronous phase
	   This is done by finding the point in the integrated profie where the difference between the integral on one side and the other is minimised
		i.e. Find phi that minismises diff = |Int[f,{phi_lower,phi}] - Int[f,{phi, phi_upper}]| where phi_lower,phi_upper are on the separatrix.

		prof_data_a: input turn-by-turn bunch monitor data.
		nwin: average over this number of turns when calculating chi-squared to decide if solution converges.
		sync_bin: where we expect the synchronous time bin.
		
		show_int_prof: show selected integrated bunch monitor data (integrated over turns).
		show_sym_plots: show result of symmetry calculation.

		"""

	nturns = prof_data_a.shape[0]
	nslices = prof_data_a.shape[1]

	data_integ_0 = np.sum(prof_data_a, axis=0)
	imax_col = np.argmax(data_integ_0)
	imin_col = np.argmin(data_integ_0)
			
	ld = len(data_integ_0)

	chisq_sym_turns = []
	isym_shape_list = []
	isym_integ_list = []
			
	turnhalfwin = 50 #moving window half size
	turn_list = np.arange(turnhalfwin, nturns)

	
	intensity = []
	k = 0
	col_l = ['k','r','b','m']
	for i1 in turn_list:

		data_int = np.sum(prof_data_a[i1-turnhalfwin: i1+turnhalfwin], axis=0)
		
		#old method
		#isym_shape, _ = find_symmetry_point_shape(data_int*(nwin/i1))
		#isym_shape_list.append(isym_shape)

		#isym_shape, _ = find_symmetry_point_shape(data_int)

		isym_integ, chisq_min = find_symmetry_point_integral(data_int, plot_data = False)			
		chisq_sym_turns.append(chisq_min)
		isym_integ_list.append(isym_integ)

		indsh = np.arange(ld) + isym_integ - sync_bin
		data_integ_shift = data_int.take(indsh, mode='wrap')
		intensity.append(np.trapz(data_integ_shift))

		if show_int_prof:
			if i1 in [turn_list[0], turn_list[int(0.3333*len(turn_list))], turn_list[int(0.6666*len(turn_list))], turn_list[len(turn_list)-1]]:
				plt.plot(data_integ_shift, col_l[k])
				plt.axvline(x=isym_integ, color=col_l[k])
				#plt.axvline(x=isym_shape, color=col_l[k],linestyle='--')

				k = k + 1

	
	if show_int_prof:	
		plt.show()

			
	chisq_sym_turns = np.array(chisq_sym_turns)
	range_chisq_win_l = []
	turn_chisq_check = nwin + np.arange(len(turn_list) - nwin)
	for d in range(len(turn_list) - nwin):
		chisq_win = chisq_sym_turns.take(np.arange(d, nwin+d))
		range_chisq_win = max(chisq_win) - min(chisq_win)
		range_chisq_win_l.append(range_chisq_win)


	threshold = 10
	ithres = None
	for i, r in enumerate(range_chisq_win_l):
		if r < threshold:
			ithres = i
			break

	if ithres == None:
		print "threshold not reached, used final value of symmetry point"
		ithres = len(range_chisq_win_l) - 1
	else:
		print "converges at turn ",turn_chisq_check[ithres]

	isym_integ_av = sum(isym_integ_list)/len(isym_integ_list)
	isym_final =  int(isym_integ_av)
	idelay =  isym_final - sync_bin #- int(0.5*nslices))	

	print "symmetry point av",isym_integ_av
	print "sync_bin ",sync_bin
	print "idelay ",idelay
			
	
	if show_sym_plots:

		indsh = np.arange(ld) + idelay
		data_integ_shift = data_integ_0.take(indsh, mode='wrap')

		plt.plot(data_integ_0,'ko', label='original')
		plt.plot(data_integ_shift,'b', label='shifted')
		plt.title("delay "+str(idelay))
		plt.axvline(x = isym_final)
		plt.axvline(x = sync_bin,color='gray',linestyle='--',label='sync phase')
		plt.legend()
		plt.savefig('shift_integ')
		plt.show()


		plt.subplot(211)
		plt.plot(turn_list, isym_integ_list, 'k')
		#plt.plot(turn_list, isym_shape_list, 'r',label='shp')
		plt.axhline(y=isym_integ_av)
		plt.xlabel('turn number')	
		plt.ylabel('symmetry point')
		#plt.plot(data_int)
		#plt.plot(data_integ_shift)
		plt.legend()
		plt.subplot(212)
		plt.plot(turn_list, intensity, 'ko-')
		plt.xlabel('turn')
		plt.ylabel('intensity')
		plt.savefig('symcalc')
		plt.show()

	return idelay


def butter_lowpass(cutoff, fs, order=5):
	nyq = 0.5 * fs
	normal_cutoff = cutoff / nyq
	b, a = butter(order, normal_cutoff, btype='low', analog=False)
	return b, a


def butter_lowpass_filter(data, cutoff, fs, order=5):
	b, a = butter_lowpass(cutoff, fs, order=order)
	y = lfilter(b, a, data)
	y = filtfilt(b, a, data)
	return y

def filter_data(data, time_step):
	"""Apply low pass Butterworth filter"""

	# Filter requirements.
	order = 20
	fs =  1/(time_step)#30.0       # sample rate, Hz
	cutoff = 30e6 #3.667  # desired cutoff frequency of the filter, Hz

	print "apply filter"

	# Get the filter coefficients so we can check its frequency response.
	b, a = butter_lowpass(cutoff, fs, order)

	# Plot the frequency response.
	w, h = freqz(b, a)#, worN=8000)

	freq = 0.5*fs*w/np.pi
	freq_MHz = freq*1e-6

	# Demonstrate the use of the filter.
	# First make some data to be filtered.
	T = 1.0             # seconds
	n = len(data) #int(T * fs)     # total number of samples
	t_ms = 1e3*np.linspace(0, time_step*n, n, endpoint=False)
	# "Noisy" data.  We want to recover the 1.2 Hz signal from this.
	#data = np.sin(1.2*2*np.pi*t) + 1.5*np.cos(9*2*np.pi*t) \
	#	    + 0.5*np.sin(12.0*2*np.pi*t)

	# Filter the data, and plot both the original and filtered signals.
	y = butter_lowpass_filter(data, cutoff, fs, order)

	plot_filter = False
	if plot_filter:
		plt.subplot(2, 1, 1)
		plt.plot(freq_MHz, np.abs(h), 'b')
		plt.plot(cutoff*1e-6, 0.5*np.sqrt(2), 'ko')
		plt.axvline(cutoff*1e-6, color='k')
		#plt.xlim(0, 0.5*fs)
		plt.xlim(0, 0.5*freq_MHz[-1])
		plt.title("Lowpass Filter Frequency Response")
		plt.xlabel('Frequency [MHz]')
		plt.grid()	

		plt.subplot(2, 1, 2)
		plt.plot(1e3*t_ms, data, 'b-', label='data')
		plt.plot(1e3*t_ms, y, 'r-', linewidth=2, label='filtered data')
		plt.xlabel('time [us]')
		plt.xlim(2000.4, 2001.1)
		plt.grid()
		plt.legend()

		plt.subplots_adjust(hspace=0.35)
		plt.savefig('filter_lowpass')
		plt.show()


		fft_prof_sig = np.abs(np.fft.fft(data))
		fft_prof_fil = np.abs(np.fft.fft(y))

		f_res = 1/(time_step*n)#resolution freq
		#xaxis = np.linspace(0,1,len(fft_prof_sig))
		xaxis = 1e-6*f_res*np.arange(len(fft_prof_sig))
	
		#plt.plot(fft_prof_flat,'k')
		plt.plot(xaxis,fft_prof_sig,'k', label='raw')
		plt.plot(xaxis,fft_prof_fil,'r-',label='filtered')
		plt.xlim(0,0.5*xaxis[-1])
		plt.ylim(0,5000)
		plt.xlabel('frequency [MHz]')
		plt.ylabel('|FFT|')
		plt.legend(loc = 'upper right')
		plt.savefig('filter_fft')
		plt.show()
		sys.exit()


	return y


def calc_waveform_harmonics(tdat_rf, sig_rf, t_zerocross):
	"""Calculate harmonic content of RF waveform at selected time"""
 
	index_time = 0
	print "t_zerocross ",t_zerocross[index_time]

	np_rf = 1
	t1 = t_zerocross[index_time]
	t2 = t_zerocross[index_time+np_rf]
	
	i1 = None
	i2 = None
	for i in range(len(tdat_rf)):
		if tdat_rf[i] > t1 and i1 == None:
			i1 = i-1
		elif tdat_rf[i] > t2 and i1 != None:
			i2 = i
			break

	tdat_rf_sel_ms = 1e3*np.array(tdat_rf[i1:i2])
	phase_ax1 = np.linspace(-180,180,len(tdat_rf_sel_ms))

	sig_rf_sel = sig_rf[i1:i2]
	t_rf_sel = tdat_rf[i1:i2]

	fft_rf = np.fft.fft(sig_rf_sel)
	fft_rf_abs = np.abs(fft_rf)
	nfft = len(fft_rf)
	fft_rf_half = fft_rf_abs[:int(0.5*nfft)]
	imax =  np.argmax(fft_rf_half)
	print "dominant freq, imax ",imax/nfft, imax

	#plt.plot(t_rf_sel, sig_rf_sel)
	#plt.show()

	hsing = [np_rf]
	hmult = [np_rf,2*np_rf,3*np_rf]
	fft_sing = np.zeros(nfft, dtype = complex)
	fft_mult = np.zeros(nfft, dtype = complex)

	fft_sing[0] = fft_rf[0]
	fft_mult[0] = fft_rf[0]
	for isel in hsing: 
		fft_sing[isel] = fft_rf[isel]
		fft_sing[-isel] = fft_rf[-isel]

	for isel in hmult: 
		fft_mult[isel] = fft_rf[isel]
		fft_mult[-isel] = fft_rf[-isel]

	habs_mult = [fft_rf_abs[i] for i in hmult]
	hnorm_mult = [fft_rf_abs[i]/fft_rf_abs[hsing[0]] for i in hmult]
	ang_mult = [np.angle(fft_rf[i]) for i in hmult]
 	#ang_h2 = np.angle(fft_rf[2])
	#ang_h3 = np.angle(fft_rf[3])

	print "norm. strenght of each harmonic ",hnorm_mult
	print "angle at each harmonic ",ang_mult


	#print "ang h3 - h1 ",ang_h3 - ang_h1

	#plt.plot(fft_rf,'k')
	plt.subplot(211)
	plt.bar(hmult, habs_mult)
	plt.ylabel('|FFT|')
	#plt.xlim(0,5)
	plt.subplot(212)
	plt.plot(hmult,ang_mult, 'ko')
	plt.ylim(-math.pi, math.pi)
	plt.xlabel('harmonic')
	plt.ylabel('phase [rad]')
	plt.savefig('fft_rf')
	plt.show()

	sig_invfft = np.fft.ifft(fft_sing)
	sig_invfft2 = np.fft.ifft(fft_mult)

	show_waveform = True
	if show_waveform:
		
		plt.plot(phase_ax1, sig_rf_sel, 'k-', linewidth=2, label = 'RF waveform')
		plt.plot(phase_ax1, sig_invfft,'r', label='h=1', linewidth=2)
		plt.plot(phase_ax1, sig_invfft2,'r--', label='h=1,2,3', linewidth=2)

		plt.legend()
		plt.savefig('waveform')
		plt.show()

	return

def sort_data(istart, time_data, signal_data, turn_times_a, time_indices, nturns, nslices, interpolate_data = True, nshift = 0, npfac=1):
	"""Divide bunch monitor data into individual turns. There is the option to interpolate onto regular time points. 
	istart: the index of the first time point
	turn_times_a: times at the start of each turn
	time_indices: indices associated with turn_times_a
	nturns: limit on the number of turns to sort
	nslices: number of time points per turn
	interpolate_data: interpolate onto a regular grid delimited by turn_times_a if True
	nshift: shift time by this number of points to deal with offset of rf signal (if shift_time is True)
	
	
"""
	time_interval = time_data[2] - time_data[1]
	#iend = istart + nturns*points_per_turn
	 
	#use ToF to define time.
	#t_dat_adj = tstart + np.linspace(0, nturns*tof, nturns*nslices)

	if npfac > 1 and not interpolate_data:
		print "must interpolate data if npfac != 1"
		sys.exit()

	#interpolate onto time base Trev/nslices
	t_dat_adj = []
	for it in range(1, nturns+1):
		time1 = turn_times_a[it-1]
		time2 = turn_times_a[it]
		t_oneturn_data = np.linspace(time1, time2 ,nslices) #fitdiff_a[it]
		t_dat_adj  = np.concatenate((t_dat_adj, t_oneturn_data))

	t_dat_adj = t_dat_adj + time_interval*nshift
	t_adj_last = t_dat_adj[-1]

	i1 = 0
	for i1, t1 in enumerate(time_data):
		if t1 >= t_adj_last:
			index_t_last = i1
			break	

	
	t_dat_orig = np.array(time_data[istart:index_t_last])
	sig_dat_orig = np.array(signal_data[istart:index_t_last])


	interpolate_data = True
	if interpolate_data:
		sig_interp = np.interp(t_dat_adj, t_dat_orig, sig_dat_orig)
		sig_dat_split = np.split(sig_interp, nturns)
	else:	
		sig_dat_split = np.split(np.array(signal_data[istart:istart + nturns*nslices]), nturns)


	zeropos = False
	data_win_sig_sel = []
	for data_win in sig_dat_split[:nturns]:


		if zeropos:
			data_win_sig = []		
			for num in data_win:
				if num <=0:
					data_win_sig.append(-num)
				else:
					data_win_sig.append(0)
		else:
			data_win_sig =  max(data_win) - np.array(data_win)
		
		data_win_sig_sel.append(data_win_sig)

	#tof = [turn_times_a[i] - turn_times_a[i-1] for i in range(1, len(turn_times_a))]
	#plt.plot(tof)
	#plt.show()
	#sys.exit()
	
	
	prof_data_a = np.array(data_win_sig_sel)	
	
	return prof_data_a



def mountain_range(data, nt, nslices, bucketlim_low=None, bucketlim_high=None, sync_bin = None, incl_xaxis = False, xaxis = None, xlabel='time bin',fmount='mountainrange', extra_data = None):
	"""Plot 2D data which has shape (nslices, nt) in mountain range style
		extra_data - 1D data to superimpose on figure """

	#nd = data.shape
	nd = len(data)
	ndx = nslices

	#establish range
	max_data_t = np.array([max(data[it]) for it in range(nt)])
	min_data_t = np.array([min(data[it]) for it in range(nt)])
	range_data_t = max_data_t - min_data_t
	max_range = np.max(range_data_t)
	
	ampfac = 2000
	
	#plt.plot(range_data_t)
	#plt.show()

	lw_a = np.linspace(0.8, 0.0, nt)
	lw_a = 0.4 - np.array(range(nt))/6000
	

	fig = plt.figure(figsize=(8, 8), facecolor='w')
	ax = plt.subplot(111, frameon=False)
	
	if not incl_xaxis:
		xaxis = range(1,ndx+1)

	for it in range(nt):
		ax.plot(xaxis, ampfac*data[it] + it,'k',lw = lw_a[it])

	ax.set_ylabel('turn number')

	if bucketlim_low != None:
		plt.axvline(x=bucketlim_low+0.5,color='gray')

	if bucketlim_high != None:
		plt.axvline(x=bucketlim_high+0.5,color='gray')

	if sync_bin != None:
		plt.axvline(x=sync_bin, color='r', linestyle='--')

	if extra_data != None:
		col_l = ['r','m','b']
		i = 0
		for ex_data in extra_data:
			plt.plot(ex_data, np.array(range(nt)), col_l[i])
			i = i + 1

	#show turn number of right axis
	ax.set_xlabel(xlabel)		
	ax.set_ylabel('turn number')
	plt.savefig(fmount)
	plt.show()
	
	return


def imshow_plot(data,extent, hline_l = [],  hline_l_2 = [], ylabel=None,xlabel=None, extra_data = None, interactive=False):
	"""data is in the form of a 2d array. 
	Extent defines the limits [vmin, vmax, hmin, hmax]
	hline_l list of vertical points where horizontal lines should be drawn """
	import matplotlib.cm as cm

	nt = data.shape[0]
	tax = np.linspace(extent[2],extent[3],nt)

	fig, ax1 = plt.subplots()
	ax1 = plt.gca()

	ax1.imshow(data, origin='lower', aspect='auto', extent = extent, cmap = cm.jet)

	ax1.set_xlim(-180, 180)
	ax1.set_ylim(extent[2],extent[3])
	if ylabel != None:
		ax1.set_ylabel(ylabel)

	if xlabel != None:
		ax1.set_xlabel(xlabel)

	if extra_data != None:
		col_l = ['k','m','b']
		i = 0
		for ex_data in extra_data:
			try:
				ax1.plot(ex_data, tax ,col_l[i])
			except:
				ax1.plot(ex_data[1], ex_data[0] ,col_l[i])
			i = i + 1

	for y in hline_l_2:
		ax1.axhline(y=tax[y],color='w',linestyle='-')
	
	for y in hline_l:
		ax1.axhline(y=tax[y],color='w',linestyle='--')

	#show turn number on right axis
	ax2 = ax1.twinx()
	ax2.set_ylim(1, nt)
	ax2.set_ylabel('turn number')

	
	plt.savefig('imshow_prof')

	#if interactive mode ion() is activated, figure will stay open until ioff() is encountered
	if interactive:
		plt.ion()

	plt.show()


def fwhm_calc(prof):
	"""calculate FWHM of individual profile"""

	diff = max(prof) - min(prof)
	nbins = len(prof)
	HM = diff/2
	pos_extremum = prof.argmax()
	try:
		nearest_above = (np.abs(prof[pos_extremum:-1] - HM)).argmin()
		nearest_below = (np.abs(prof[0:pos_extremum] - HM)).argmin()
			
		fwhm = pos_extremum+nearest_above - nearest_below
		fwhm_deg = 360*fwhm/nbins
	except:
		fwhm = 0

	return fwhm_deg



def intensity_calc_fn(data, phase_ax, ttime_a, tof_l, time_interval_ns, dirname, fname_sig, sym_phase, filter_data, write_intensity_file=False):
	#calculate intensity and write to file
	#intensity is calculated by integrating the signal turn-by-turn.
	
	plot_intensity = False
	show_fwhm = True
	
	if write_intensity_file:
		#construct offset file name
		if "case" in fname_sig:
			fspl =  fname_sig.split("_")

		else:
			fspl =  fname_sig.split(".")


		foffset = ''
		for str1 in fspl[:-1]:
			foffset = foffset + str1+"_"
			
		if filter_data:
			intensity_file = foffset + "intensity_filt.txt"
		else:
			intensity_file = foffset + "intensity.txt"
		print "intensity ",intensity_file

		fi = file(intensity_file,'w')
		#fi = open('intensity.txt','w')

	#calculate phase at which data is maximum for each turn
	max_prof_index = np.argmax(data, axis=1)
	phase_profmax = [phase_ax[m] for m in max_prof_index]
	#mean_phase_data = [sum([p1*f for p1,f in zip(phase_ax, prof)])/sum(prof) for prof in data]
	mean_phase_data = sym_phase
	std_phase_data = [np.sqrt(sum([f*(p1-meanp)**2 for p1,f in zip(phase_ax, prof)])) for prof,meanp in zip(data,mean_phase_data)]

	tturn_ms = 1e3*ttime_a[:-1]
	intensity_raw = []
	fwhm_raw = []
	i1 = 0
	
	nbins = data.shape[1]

	for prof_data, tof,  t1, phmean, phstd, phmax in zip(data, tof_l, ttime_a[:-1], mean_phase_data, std_phase_data, phase_profmax):

		t_int_ns = 1e9*tof/nbins		
		#int1 = time_interval_ns*np.trapz(prof_data)
		int1 = t_int_ns*np.trapz(prof_data)
		intensity_raw.append(int1)
	
		#calculate FWHM
		diff = max(prof_data) - min(prof_data)
		HM = diff/2
		pos_extremum = prof_data.argmax()
		try:
			nearest_above = (np.abs(prof_data[pos_extremum:-1] - HM)).argmin()
			nearest_below = (np.abs(prof_data[0:pos_extremum] - HM)).argmin()
			
			FWHM = pos_extremum+nearest_above - nearest_below
			fwhm_raw.append(FWHM_deg)
		except:
			FWHM = 0
			fwhm_raw.append(0)

		FWHM_deg = 360*FWHM/nbins		
	
		show_intensity_fwhm_example = False
		if show_intensity_fwhm_example:
			if t1 > 0.5e-3 and FWHM*time_interval_ns > 200: 

				#plt.plot(np.abs(prof_data[0:pos_extremum] - HM))
				#plt.show()
				print "t1 ",t1
				print "FWHM points, time [ns] ",FWHM, FWHM*time_interval_ns

				time_ax = time_interval_ns*np.arange(len(prof_data))
				plt.plot(time_ax, prof_data,'k')
				plt.axvline(x=nearest_below*time_interval_ns, color='r')
				plt.axvline(x=(pos_extremum+nearest_above)*time_interval_ns, color='r')
				plt.ylabel('signal [V]')
				plt.xlabel('time [ns]')
				plt.title('integral [nVs] '+str(int1)[:5]+ ' FWHM [ns] '+str(time_interval_ns*FWHM))

				plt.savefig('integral_fwhm.png')								
				plt.show()
				#sys.exit()


		i1 = i1 + 1
		if write_intensity_file:
			print >> fi, repr(t1), repr(int1), repr(FWHM_deg), repr(phmean),repr(phstd),repr(phmax)

	if write_intensity_file:
		fi.close()	
	
	max_intensity = max(intensity_raw)
	imax = intensity_raw.index(max_intensity)

	t_imax = ttime_a[imax]
	t_imax_ms = 1e3*t_imax	
	fwhm_raw_a = np.array(fwhm_raw)


	if plot_intensity:
		fig = plt.figure()

		print "time interval ",time_interval_ns
		print "minimum FWHM ",min(fwhm_raw)*time_interval_ns

		if show_fwhm:
			plt.subplot(211)
		else:
			plt.xlabel('time [ms]')

	
		plt.plot(1e3*ttime_a[:-1],intensity_raw,'k.')
		plt.ylabel('intensity [nVs]')
	
		plt.ylim(ymin=0)
		#plt.xlim(t_imax_ms,1e3*ttime_a[-1])
	
		#plt.axvline(x=t_imax_ms, color='m')

		print "mean of mean phase data ",sum(mean_phase_data)/len(mean_phase_data)
	
		if show_fwhm:
			plt.subplot(212)
			#plt.plot(1e3*ttime_a[:-1],time_interval_ns*fwhm_raw_a,'k.')
			plt.plot(1e3*ttime_a[:-1], mean_phase_data, label='mean')
			plt.plot(1e3*ttime_a[:-1], std_phase_data, label='std')
			#plt.ylabel('FWHM [nVs]')
			#plt.xlim(t_imax_ms,1e3*ttime_a[-1])
			#plt.ylim(ymin=0, ymax = 350)

			plt.ylabel('phase [deg]')
			#plt.axvline(x=t_imax_ms, color='m')

			if t_phis_file != None:
				for t in t_phis_file:
					if t != 0.0:
						plt.axvline(x=t+t_imax_ms, color='r')

			plt.legend()
		#plt.axvline(x=0.18)
		
		plt.savefig('intensity')
		plt.show()

	

	return intensity_raw, phase_profmax, std_phase_data, fwhm_raw


def synch_tune(V_rf, phi_s, harmonic, alpha_c, ke, mass):

	Etot = ke + mass
	p0 = np.sqrt(Etot**2 - mass**2)
	betarel = p0/Etot
	betagam = p0/mass
	gamma = betagam/betarel
	eta = (alpha_c - gamma**(-2))
	
	Qs = Qs_fn(phi_s, V_rf, harmonic, betarel, Etot, eta)

	return Qs


def longitudinal_param(V_rf, phi_s, harmonic, alpha_c, rho, ke, tof, nslices, nturns, mass):
	"""harmonic: harmonic of RF
	   tof: mean ToF
	   nslices: data points per RF period"""

	T_rf = tof/harmonic#2*np.pi/(harmonic*omegarev) #rf period
	T_rf_ns = 1e9*T_rf
	omega_rev = 2*math.pi/T_rf
	omega_rf = omega_rev*harmonic
	dphasebin = 2*math.pi/nslices
	dtbin_phase = dphasebin/omega_rf

	Etot = ke + mass
	p0 = np.sqrt(Etot**2 - mass**2)
	betarel = p0/Etot
	betagam = p0/mass
	gamma = betagam/betarel
	eta = (alpha_c - gamma**(-2))
	
	Qs = Qs_fn(phi_s, V_rf, harmonic, betarel, Etot, eta)
	Ts = 1/Qs
	yfac = np.sqrt(1 - 0.5*(math.pi - 2*phi_s)*np.tan(phi_s))
	bh = ((2*Qs)/(harmonic*abs(eta))*yfac) #bucket height
		 
	print "kinetic energies at selected times ",ke
	print "eta ",eta
	print "rev freq ",omega_rev/(2*math.pi)
	print "rel. beta at selected times ",betarel
	print "eta ",eta
	print "Qs ",Qs
	print "Ts ",Ts
	print "p [MeV/c] ",1e-6*p0
	print "bucket height ",bh


	#set Rnom to ensure consistency between omega_rf in code and omega_rf specified here
	Rnom =  betarel*SPEED_OF_LIGHT/omega_rf
	Bnom = np.sqrt((Etot**2 - mass**2)/(SPEED_OF_LIGHT*rho)**2)

	#bucket phase extrema calculation
	if phi_s == 0:
		phi_l1 = math.pi
		phi_l2 = -math.pi
	else:
		phi_l1, phi_l2 = phase_extrema(phi_s)

	if phi_l1 < phi_l2:
		phi_lower = phi_l1
		phi_upper = phi_l2
	else:
		phi_lower = phi_l2
		phi_upper = phi_l1

	sync_frac =  (phi_s-phi_lower)/(phi_upper-phi_lower)

	bucket_phase_width = abs(phi_l1 - phi_l2)
	bucket_dur_ns = T_rf_ns*(bucket_phase_width/(2*math.pi))
	outbucket_dur_ns = T_rf_ns - bucket_dur_ns

	synch_index = int(round(nslices*(phi_s - phi_l2)/(2*math.pi)))

	nslices_active = int(nslices*bucket_dur_ns/T_rf_ns) #number of slices in bucket
	nslices_outb = nslices - nslices_active
	print "bucket points, out-of-bucket points, bucket points/all points ",nslices_active, nslices_outb, nslices_active/nslices

	if phi_s == 0:
		binslz = 0
		binsuz = 0
	else:
		binslz = abs(int(nslices*(phi_lower+math.pi)/(2*math.pi)))# nslices_outb# #nslices_outb
		binsuz =   nslices - (binslz + nslices_active)

	ph_ax = np.linspace(-180,180,nslices)

	npts = 1000
	a1 = 0.5*harmonic*omega_rf*eta
	a2 = omega_rf*V_rf/(2*math.pi*(betarel**2)*Etot)
	dp_sep, phase_sep, t_sep_ns = separatrix(phi_s, eta, harmonic, npts, a1, a2, T=T_rf)

	t_phi_l2 = T_rf_ns*(phi_l2/(2*math.pi)) #time at phi_l2
	t_sep_ns = np.insert(t_sep_ns, 0, t_phi_l2)
	dp_sep = np.insert(dp_sep,0,0)
	phase_sep = np.insert(phase_sep, 0, phi_l2) 

	phase_sep_deg = (180/math.pi)*phase_sep
	sync_bin = binslz + int(nslices_active*sync_frac)
	
	#Hamiltonian on separatrix
	Hsep = H_fn(math.pi - phi_s, 0, phi_s, a1, a2)
	#time adjustment
	tadj = T_rf_ns*(phi_l2/(2*math.pi)) - outbucket_dur_ns 

	return Ts, dtbin_phase, omega_rf, Rnom, Bnom, phase_sep_deg, dp_sep, bh, sync_bin, synch_index, binslz, binsuz, tadj, Hsep, phi_l1, phi_l2

def longitudinal_turn(V_rf, phi_s, harmonic, alpha_c, ke, Rnom, nturns, dtbin_phase, mass):

	Egain_turn = np.arange(nturns)*V_rf*np.sin(phi_s)
	Ekin_turn = ke + Egain_turn
	Etot_turn = Ekin_turn + mass
	p0_turn = np.sqrt(Etot_turn**2 - mass**2)
	beta_turn = p0_turn/Etot_turn
	betagam_turn = p0_turn/mass
	gamma_turn = betagam_turn/beta_turn
	eta_turn = (alpha_c - gamma_turn**(-2))
	omega_rf_turn = beta_turn*SPEED_OF_LIGHT/Rnom #note this may not be consistent with real variation as Rnom is fixed here

	#calculate dEbin as it is done in CERN code
	dEbin_c_turn = dEbin_fn(phi_s, V_rf, harmonic, beta_turn, Etot_turn, eta_turn, dtbin_phase, omega_rf_turn)

	return Etot_turn, dEbin_c_turn, p0_turn, eta_turn, beta_turn


def file_read(filen, sel_col=0):
	""" filen: file path.
		sel_col: indext of column in data file to read. Default 0"""
	
	f1 = open(filen,'r')
	data_l = []
	for data in f1:
		dataspl = data.split()
		data_l.append(float(dataspl[sel_col]))
	data_a = np.array(data_l)

	return data_a



def filter_isolated_cells(array, struct):
	""" Return array with completely isolated single cells removed
	:param array: Array with completely isolated single cells
	:param struct: Structure array for generating unique regions
	:return: Array with minimum region size > 1
	"""

	nonzero_indices = np.nonzero(array)
	
	nd = len(array)
	array_ones = np.zeros((nd, nd))
	array_ones[nonzero_indices] = 1

	filtered_array = np.copy(array)
	
	id_regions, num_ids = ndimage.label(array_ones, structure = struct)
	
	id_sizes = np.array(ndimage.sum(array_ones, id_regions, range(num_ids + 1)))
	
	#area_mask = (id_sizes == 1)
	area_mask = (id_sizes <= 3)
	filtered_array[area_mask[id_regions]] = 0
	
	
	return filtered_array


def run_tomo_code(tomo_outdir, input_param):
 	"""Prepare input file for tomography code. Run tomography code"""

	nturns, nslices, synch_index, dtbin_phase, binslz, binsuz, recon_start, recon_stop, recon_step, Bnom, Bdot, V_rf, V_rf2, harmonic, hratio, phi12, Rnom, rho, gamma_tr = input_param

	input_settings = get_input("input_default.dat")

	fname = 'input_v2.dat'

	descrip = 'KURRI'
	datafile = 'rawdata.txt'
	input_settings["numframes"] = nturns
	input_settings["numbins"] = nslices
	input_settings["synchbin"] = synch_index #int(0.5*nslices)
	input_settings["framebinwidth"] = dtbin_phase	
	input_settings["binslowerignore"] = binslz
	input_settings["binsupperignore"] = binsuz
	input_settings["filmstart"] = recon_start
	input_settings["filmstop"] = recon_stop
	input_settings["filmstep"] = recon_step
	input_settings['Bref'] = Bnom
	input_settings['Bdot'] = Bdot
	input_settings['VRF1ref'] = V_rf
	input_settings['VRF2ref'] = V_rf2
	input_settings['h'] = harmonic
	input_settings['hratio'] = hratio
	input_settings['phi12'] = phi12
	input_settings['Rnom'] = Rnom
	input_settings['rhonom'] = rho
	input_settings['gammatnom'] = gamma_tr
	input_settings['pickup'] = 1.0
	input_settings['selffields'] = 0
	write_inputfile(fname,descrip, datafile, input_settings)

	temp_files = ['plotinfo.data','profiles.data','d001.data','d100.data','jmax.data','image001.data','image100.data']


	#run tomography code
	call(["tomo"])

	files = os.listdir(os.curdir)
	for fn in files:
		if fn[-4:] == 'data':
			os.rename(fn,tomo_outdir+"/"+fn)

	return

