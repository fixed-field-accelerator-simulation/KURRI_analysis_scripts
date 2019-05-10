from __future__ import division
import sys
import numpy as np
import math
from scipy.constants import m_p, c, e
from PyHEADTAIL.particles.particles import Particles
import PyHEADTAIL.particles.generators as generators
from PyHEADTAIL.trackers.transverse_tracking import TransverseMap
from PyHEADTAIL.trackers.simple_long_tracking import RFSystems, LinearMap
import PyHEADTAIL.cobra_functions.stats as st
import matplotlib.pyplot as plt
from PyHEADTAIL.monitors.monitors import BunchMonitor
from PyHEADTAIL.particles.slicing import UniformBinSlicer, UniformChargeSlicer
from matplotlib.colors import LogNorm
from subprocess import call
import os
#import tomo_defs_181101 as tdefs
import tomo_defs_190510 as tdefs
import time
import matplotlib.animation as animation
from scipy.optimize import leastsq
from scipy.optimize import fmin
from scipy.integrate import quad

#main settings
phi_s_deg_nom = 20
phi_s_nom = phi_s_deg_nom*math.pi/180
n_slices = 128 #number of time points her turn
binslowerignore = 0#5 #add zeroes before profile (earlier time)
binsupperignore = 0#1 #add zeroes after profile (later time)
nturns= 1000 #2000#300#100
dturns = 1
V_rf = 4e3 # in [V]

recon_start = 1
recon_stop = 100
recon_step = 25

animate = False

phi_foil = lambda V, phi_s, deltaE: np.arcsin((V*np.sin(phi_s) + deltaE)/V)

input_settings = tdefs.get_input("input_default.dat")
plt.switch_backend('TkAgg')	


gam2pMev =  1e-6*m_p * (c**2)/e

# general simulation parameters
#n_particles = 100000 #200000
n_particles = 1
n_segments = 1

# machine parameters
inj_alpha_x = 0#-1.2
inj_alpha_y = 0#15
inj_beta_x = 0.78#5.9 # in [m]
inj_beta_y = 3.75#5.7 # in [m]
Qx = 0.3045*12
Qy = 0.115*12
kindex = 7.6
alpha_c = 1/(kindex+1)
gamma_tr = (1/alpha_c)**0.5
alpha_c_arr = np.array([alpha_c]) #PyHeadtail expects alpha_c in array form

harmonic = 1
#pipe_radius = 5e-2

#Bdot=0 # in T/s
bending_radius=0.489396551697 # in m
Rnom = 4.41
circumference = Rnom*2*np.pi
phi_offset = 0 #np.arcsin(bending_radius*circumference*Bdot/V_rf) # measured from aligned focussing phase (0 or pi)

# beam parameters
Ekin = 11e6 #1.4e9 # in [eV]
intensity = 0#8e12
epsn_x = 0#1e-3*5e-6 # in [m*rad]
epsn_y = 0#1e-3*5e-6 # in [m*rad]
epsn_z = 0.05 # 4pi*sig_z*sig_dp (*p0/e) in [eVs]

bunch_duration = 250e-9 #this may later be modified if a matched distribution is created

tomo_outdir = "tomo_output"
files = os.listdir(os.curdir)
if tomo_outdir not in files:
	print "create output directory "+tomo_outdir
	os.mkdir(tomo_outdir)

def mom_to_ke(mom, mass):
	"""mom in eV/c
	gives results in eV
	"""
	ke = np.sqrt(mom**2 + mass**2) - mass

	return ke

def ke_to_mom(ke, mass):
	"""ke in eV
	gives results in eV/c
	"""
	mom = np.sqrt(((ke+mass)**2 - mass**2))

	return mom

m_p_ev = m_p*(c**2)/e #938.2720814e6

gamma = 1 + e * Ekin / (m_p * c**2)
beta = np.sqrt(1 - gamma**-2)
print('rel beta: ' + str(beta))
eta = alpha_c - gamma**-2
print('eta: ' + str(eta))
if eta < 0:
    phi_offset = np.pi - phi_offset

Etot = gamma * m_p * c**2 / e

p0 = np.sqrt(gamma**2 - 1) * m_p * c
p0_mev = 1e-6*c*p0/e #convert to MeV/c
ke_inj = mom_to_ke(1e6*p0_mev, m_p_ev)
ke_threshold = 11.42e6 #foil crossing threshold
#set synchronous phase
deltaE_foil = 760
phi_s_foil_deg = (180/math.pi)*phi_foil(V_rf, (math.pi/180)*phi_s_deg_nom, deltaE_foil)

include_foil = False #model as energy loss and shift in synchronous phase
if ke_inj < ke_threshold and include_foil:
	phi_s_deg = phi_s_foil_deg	
else:
	phi_s_deg = phi_s_deg_nom	

phi_s = phi_s_deg*math.pi/180 #0.8*2*math.pi/3 #specify synch phase. With single RF, work out Bdot	

Qs = np.sqrt(np.abs(eta*np.cos(phi_s)) * harmonic* V_rf / (2 * np.pi * beta**2 * Etot))
Qs_amp_fn = lambda phase_m, phi_s, Qs: Qs*(1 - (phase_m**2/16)*(1+(5/3)*((np.tan(phi_s))**2)))


print "Qs, Ts ",Qs, 1/Qs

beta_z = np.abs(eta) * circumference / (2 * np.pi * Qs)
print('beta_z: ' + str(beta_z))


sigma_z_0= bunch_duration*beta*c

turn_period = circumference / (beta * c)
Bdot = (V_rf*np.sin(phi_s))/(bending_radius*2*math.pi*Rnom)
p_increment_0= e*bending_radius*Bdot*turn_period

Bdot_nom = (V_rf*np.sin(phi_s_nom))/(bending_radius*2*math.pi*Rnom)
p_increment_nom = e*bending_radius*Bdot_nom*turn_period


# calculations
print "energy gain at synch phase ",V_rf*np.sin(phi_s)

omegarev = beta*c/Rnom #angular rev. frequency

sin_phi_s = bending_radius*2*math.pi*Rnom*Bdot/V_rf
#phi_s = np.arcsin(sin_phi_s)
print "Bdot, sin_phi_s, phi_s, arcsin(sin_phi_s)",Bdot, sin_phi_s, phi_s, np.arcsin(sin_phi_s)
print "Bdot_nom ",Bdot_nom

print "eta, phi_offset ",eta,phi_offset
print "tof ",2*math.pi*Rnom/(beta*c)

# BETATRON
# Loop on number of segments and create the TransverseSegmentMap
# for each segment.
s = np.arange(0, n_segments + 1) * circumference / n_segments
alpha_x = inj_alpha_x * np.ones(n_segments)
beta_x  = inj_beta_x * np.ones(n_segments)
D_x     = np.zeros(n_segments)
alpha_y = inj_alpha_y * np.ones(n_segments)
beta_y  = inj_beta_y * np.ones(n_segments)
D_y     = np.zeros(n_segments)

T = 2*np.pi/(harmonic*omegarev) #rf period
print "T [us] ",T*1e6

yfac = np.sqrt(1 - 0.5*(math.pi - 2*phi_s)*np.tan(phi_s))
bh = ((2*Qs)/(harmonic*abs(eta))*yfac)
print "bucket height ",bh

def frequency_vs_phase(data, Qs, sync_bin):
	"""FFT of each time bin
	data: 2D data to to be analyzed
	Qs: synchrotron tune
	sync_bin: synchronous time bin"""
	

	nt =  data.shape[0]
	xr = np.arange(nt)/nt
	data_flat = data.flatten()

	fmin = 0.002
	fmax = 0.01
	fch = 0.004
	indmin = None
	indmax = None
	indch= None
	for x1 in xr:
		if x1 > fmin and indmin == None:
			indmin = np.where(xr == x1)[0][0]
		if x1 > fmax and indmax == None:
			indmax = np.where(xr == x1)[0][0]
		if x1 > fch and indch == None:
			indch = np.where(xr == x1)[0][0]
	print "indmin, indmax ",indmin, indmax

		
	phase_deg = np.array([i*360/n_slices for i in range(n_slices)])
	phase_offset = phase_deg[sync_bin - 1] - phi_s_deg
	phase_deg = phase_deg - phase_offset
	phase_m_deg = phase_deg - phi_s_deg
	phase_m = phase_m_deg*math.pi/180

	i_plot = [25,45]

	maxfft_l = []
	fftabs_l = []
	for i in range(n_slices):
		fft_prof_tbin = np.fft.fft(data[:,i])
		fftabs_prof_tbin = np.abs(fft_prof_tbin)
		fftabs_l.append(fftabs_prof_tbin)
		#fftabs_l.append(fft_prof_tbin)

		maxfft1 = max(fftabs_prof_tbin[indmin:indmax])
		maxfft_l.append(maxfft1)

		#if i%int(0.25*n_slices) == 0:
		if i in i_plot:
			plt.subplot(211)
			plt.plot(data[:,i],'o-')
			plt.title("bin "+str(phase_m_deg[i]))
			plt.subplot(212)
			plt.plot(xr, fftabs_prof_tbin)
			#plt.xlim(fmin,0.01)
			plt.xlim(0,0.01)
			plt.ylim(0, max(maxfft_l))
			plt.show()

	fftabs_a = np.transpose(np.array(fftabs_l))
	clev = np.linspace(0, max(maxfft_l), 20)
	print "fftabs_a shape, n_slices, nt ",fftabs_a.shape, n_slices, nt

	
	H_phase = H_fn(phase_deg*math.pi/180, 0, phi_s, a1, a2)
	#plt.plot(phase_deg, H_phase)
	#plt.show()
	
	#J = Qs*phase_m**2/(2*harmonic*abs(eta))
	Qs_amp = Qs_amp_fn(phase_m, phi_s, Qs)
	
	
	print "phase at sync bin ",phase_deg[sync_bin - 1]
	plt.plot([0],[Qs],'bo')
	for i in i_plot:
		plt.axvline(x=phase_m_deg[i])
		plt.axvline(x=phase_m_deg[i])

	plt.plot(phase_m_deg, Qs_amp)
	plt.axhline(y=xr[indch])
	plt.contour(phase_m_deg, xr,fftabs_a, levels = clev)
	plt.ylim(fmin, fmax)
	plt.show()

	plt.plot(phase_deg, fftabs_a[indch],'k')
	plt.plot(phase_deg, fftabs_a[2*indch],'r')
	plt.show()	
	#plt.plot(phase_m, Qs_amp)
	#plt.show()
	sys.exit()


	specgram = False
	if specgram:
		from scipy import signal
		widths = np.arange(1,1000)
		#cwtmatr = signal.cwt(prof_data_flat, signal.ricker, widths)
		#plt.plot(prof_data_flat)
		print prof_data_flat[:10]

		dt = time_interval_ns
		Fs = n_slices# the sampling frequency
		print "Fs ",Fs
		nfft = 2**9
		#plt.imshow(cwtmatr, extent=[-1, 1, widths[-1], 1], origin='lower', aspect='auto')
		plt.specgram(prof_data_flat, NFFT=nfft, Fs=Fs, noverlap = int(0.5*nfft))
		plt.show()
		
		
foil_rough_calc = False
if foil_rough_calc:
	print "phi_s original, effective at phi_s = 20 deg ",phi_s_deg, phi_s_foil_deg

	phi_s_deg_a = np.linspace(0,50,100)
	phi_s_rad_a = (math.pi/180)*phi_s_deg_a
	phi_s_rad_foil_a = phi_foil(V_rf, phi_s_rad_a, deltaE_foil)
	phi_s_deg_foil_a = (180/math.pi)*phi_s_rad_foil_a

	plt.plot(phi_s_deg_a, phi_s_deg_foil_a, 'k', linewidth=2)
	plt.xlabel(r'nominal $\phi_s$ [deg]')
	plt.ylabel(r'effective $\phi_s$ [deg]')
	plt.ylim(ymin=0)
	plt.show()


	ref_track = True
		
	nturns = int(350/(1e6*T))
	turns_a = np.arange(1, nturns+1)
	time_us_a = 1e6*T*turns_a
	
	cav_delE_foil = V_rf*np.sin(phi_s_foil_deg*math.pi/180)
	cav_delE_nofoil = V_rf*np.sin(phi_s_deg_nom*math.pi/180)		
	ke = ke_inj
	ke_ref_l = []
	
	for it in range(nturns):

		if ke < ke_threshold:
			ke = ke + cav_delE_foil
			ke = ke - deltaE_foil
		else:
			ke = ke + cav_delE_nofoil	
		ke_ref_l.append(1e-6*ke)

	ke_ref_a = np.array(ke_ref_l)
	mom_ref_a = ke_to_mom(ke_ref_a*1e6, m_p_ev)
	mom_high = mom_ref_a*(1+bh)
	mom_low = mom_ref_a*(1-bh)
	
	ke_high = 1e-6*mom_to_ke(mom_high, m_p_ev)
	ke_low = 1e-6*mom_to_ke(mom_low, m_p_ev)

	mom_bh_up = mom_high - mom_ref_a
	mom_bh_down = -mom_low + mom_ref_a
	ke_bh_up = ke_high - ke_ref_a
	ke_bh_down = -ke_low + ke_ref_a

	ke_lim1_index = None
	ke_lim2_index = None
	for i1, ke1 in enumerate(ke_high):
		if ke1 > ke_threshold*1e-6:
			ke_lim1_index = i1
			break

	for i1, ke1 in enumerate(ke_low):
		if ke1 > ke_threshold*1e-6:
			ke_lim2_index = i1
			break
	print "time limit indices ",ke_lim1_index, ke_lim2_index

	if ke_lim2_index != None:
		crossing_time = time_us_a[ke_lim2_index] - time_us_a[ke_lim1_index]
		print "foil crossing time, turns ",crossing_time, int(crossing_time/(1e6*T))
	
	plt.plot(time_us_a, ke_ref_l, 'k', linewidth=2)
	plt.plot(time_us_a, ke_high, 'm', linewidth=2)
	plt.plot(time_us_a, ke_low, 'm', linewidth=2)
	plt.axhline(y = ke_threshold*1e-6, color='r')

	if ke_lim1_index != None:
		plt.axvline(x=time_us_a[ke_lim1_index], color='gray')
	if ke_lim2_index != None:
		plt.axvline(x=time_us_a[ke_lim2_index], color='gray')
	plt.ylabel('kinetic energy')
	plt.xlabel('time [us]')
	plt.show()

	#plt.plot(ke_bh_up)
	#plt.plot(ke_bh_down)
	#plt.show()

	sys.exit()



#Hamiltonian
H_fn = lambda phi, dp, phi_s, a1, a2: a1*(dp**2) + a2*(np.cos(phi) - np.cos(phi_s) + (phi - phi_s)*np.sin(phi_s))# + np.random.uniform(0,200)
H_fn_gen_pot = lambda n, phi, phi_s, coef, ang_off: coef*(np.cos(n*phi + ang_off) - np.cos(n*phi_s + ang_off) + (phi - phi_s)*np.sin(n*phi_s + ang_off))
H_fn_gen_pot2 = lambda phi, phi_s, coef, ang_off,n: coef*(np.cos(n*phi + ang_off) - np.cos(n*phi_s + ang_off) + (phi - phi_s)*np.sin(n*phi_s + ang_off))

H_fn_gen = lambda phi, dp, phi_s, a1, a2, n_l, coef_l, ang_l: a1*(dp**2) + sum([H_fn_gen_pot(n, phi, phi_s, a2*coef, ang) for n, coef, ang in zip(n_l, coef_l, ang_l)])

potential = lambda phi, phi_s, a2: a2*(np.cos(phi) - np.cos(phi_s) + (phi - phi_s)*np.sin(phi_s))
Y_factor = lambda phi_s : abs(np.cos(phi_s) - 0.5*(math.pi - 2*phi_s)*np.sin(phi_s))**0.5



omega_rf = omegarev*harmonic
a1 = 0.5*harmonic*omega_rf*eta
a2 = omega_rf*V_rf/(2*math.pi*(beta**2)*Etot)

bucket_amp = (V_rf/(math.pi*(beta**2)*Etot*harmonic*abs(eta)))**0.5
bucket_height = (2**0.5)*bucket_amp*Y_factor(phi_s)
print "beta, eta ",beta, eta
print "sin(phi_s), Y_fac ",np.sin(phi_s), Y_factor(phi_s)
print "bucket_amp, bucket_height ",bucket_amp, bucket_height
print "a2/a2 ratio ",a2/a1


#calculate points on separatrix
npts = 1000
dp_sep, phase_sep, t_sep_ns = tdefs.separatrix(phi_s, eta, harmonic, npts, a1, a2)
#phase_sep = t_sep_ns*(2*math.pi)/(1e9*T)
phase_sep_deg = (180/math.pi)*phase_sep



#sys.exit()

#plt.plot(t_sep_ns, dp_sep)
#plt.show()
#sys.exit()


ph_a = np.linspace(phase_sep[0], phase_sep[-1], n_slices)
#ph_a = np.linspace(phase_sep[0]-0.5, phase_sep[-1]+0.5, n_slices)

ph_a_deg = (180/math.pi)*ph_a
for i1, ph1 in enumerate(ph_a):
	if ph1 > phi_s:
		synch_index_ideal = i1
		print "synch_index ",synch_index_ideal
		break

binwidth_phase = abs(phase_sep[0] - phase_sep[-1])/(n_slices-1)
binwidth_time = binwidth_phase*T/(2*math.pi)

check_H = False
if check_H:
	phase_test = np.linspace(-2*math.pi, 2*math.pi, 1000)
	H_dp_zero = H_fn(phase_test, 0, phi_s, a1, a2)

	i = 0
	sw = True
	for H, ph in zip(H_dp_zero, phase_test):
		if i > 0:
			if ph > math.pi - phi_s and sw:
				H_sep = H
				print H, ph
				sw = False
			if not sw:
				if H > H_sep:
					ph_turn = ph
					print "ph_turn ",ph_turn
					break
 
			

		i = i + 1
	plt.plot(phase_test, H_dp_zero)
	plt.axvline(x = phase_sep[0])
	plt.axvline(x = phase_sep[-1])	
	plt.axvline(x= phi_s,color = 'r')
	plt.axhline(y = 0, color='k')

	plt.show()

if phi_s > 0.5*math.pi:
	a2 = -a2

#Hamiltonian on separatrix
Hsep = H_fn(math.pi - phi_s, 0, phi_s, a1, a2)

print "Hsep ",Hsep

#for phi_s1 in np.linspace(0, math.pi, 10):
#	print "phi_s, potential ",phi_s1, potential(math.pi-phi_s1, phi_s1, a2),potential(phi_s1, phi_s1, a2)

plot_H_contours = False
if plot_H_contours:

	
	print "H at synch point ",H_fn(phi_s, 0, phi_s, a1, a2)
	print "H at separatrix approx ",H_fn(phase_sep[0], 0, phi_s, a1, a2)

	Hnorm_a = np.linspace(0,1,10)

	for Hn in Hnorm_a:

		Hamil = Hn*Hsep#-1*a1

		dp_c, t_c = tdefs.contour(phi_s, eta, T, harmonic, npts, a1, a2, Hamil)
		phase_c = t_c*(2*math.pi)/(1e9*T)
		

		plt.plot(phase_c, dp_c, 'b')
		plt.plot(phase_c, -dp_c, 'b')

	plt.plot(phase, dp_sep, 'r')
	plt.plot(phase, -dp_sep, 'r')


	#plt.axhline(y=bucket_height,color='r')
	plt.show()

def phi_cross(phi_g, phi_s, H, phi_upper):

	fn = abs(H_fn_gen_pot(1, phi_g, phi_s, a2, 0) - H)
	
	if phi_g < phi_s:
		fn = 1e6
	elif phi_g > phi_upper:
		fn = 1e6

	return fn


generate_dist = True
if generate_dist:

	run_pyheadtail = False
	
	print "generate distribution at phi_s [deg] ",phi_s*180/math.pi
	print "nominal phi_s [deg]",phi_s_nom*180/math.pi
	synch_index = synch_index_ideal 
	
	vfac = 1.0

	calculate_H_period = True
	if calculate_H_period:
		ntest = 100
		#ph_a_left = np.linspace(phi_s -0.01 , phi_s-0.1, ntest)
		#ph_a_left = np.linspace(phi_s  -0.09727273 , phi_s-0.01, ntest)

		ph_a_test = np.linspace(phase_sep[0], phase_sep[-1],ntest)
		ph_a_test_deg = (180/math.pi)*ph_a_test
		
		#print "ph_a_left ",ph_a_left
		#sys.exit()
		n_l = [1,3]
		coef_l = [1, 0.1]
		ang_l = [0, 1.2]
		H_ph_gen = H_fn_gen(ph_a_test, 0, phi_s, a1, a2, n_l, coef_l, ang_l)
		#plt.plot(ph_a_test, H_ph_gen)
		#plt.show()

		phase_mid = None
		if len(n_l) == 1:
			ph_a_left = np.linspace(phase_sep[0] , phi_s-1e-3, ntest)
			H_phaseaxis_left = H_fn(ph_a_left, 0, phi_s, a1, a2)
			
		else:
			phase_min_H = ph_a_test[np.argmax(H_ph_gen)]
			ph_a_left = np.linspace(phase_sep[0] , phase_min_H - 1e-4, ntest)

			phase_mid = phase_min_H
			
		#print "phase_min_H ",phase_min_H, (180/math.pi)*phase_min_H

		#H_phaseaxis = H_fn(ph_a_test, 0, phi_s, a1, a2)

		
		ph_a_right, ph_a_left = tdefs.phase_map_left_to_right(ph_a_left, phi_s, a1, a2, phase_sep[-1], n_l, coef_l, ang_l, phase_mid = phase_mid)
		
		ph_a_left_deg = (180/math.pi)*ph_a_left
		
		H_phaseaxis_left = H_fn_gen(ph_a_left, 0, phi_s, a1, a2, n_l, coef_l, ang_l)
		
		print "len ph_a_left, H ",len(ph_a_left), len(H_phaseaxis_left)

		ph_a_right_deg = (180/math.pi)*ph_a_right

	
		show_ph_search = False
		if show_ph_search:
			plt.plot(ph_a_test_deg, H_ph_gen)
			plt.plot(ph_a_left_deg, H_phaseaxis_left, 'r') 
			plt.plot(ph_a_right_deg, H_phaseaxis_left, 'mo-') 
			plt.axvline(x=phi_s_deg, color='gray')
			plt.axvline(x=phase_sep_deg[0], color='gray')
			plt.axvline(x=phase_sep_deg[-1], color='gray')
			plt.show()
			sys.exit()

		
		Ts_amp_a, J_a, ph_right_sol = tdefs.synchrotron_period_vs_amplitude(ph_a_left, phi_s, a1, a2, phase_sep[-1], harmonic, omegarev, eta, n_l, coef_l, ang_l)

		ph_right_sol_deg = (180/math.pi)*ph_right_sol

		#integrate between left and right turning points for each value of H

		H0_minus_pot = lambda phi, phi_s, coef, H0: H0 -  coef*(np.cos(phi) - np.cos(phi_s) + (phi - phi_s)*np.sin(phi_s))
		#H0_minus_pot = lambda phi, phi_s, coef, ang_off, n, H0: H0 -  coef*(np.cos(n*phi + ang_off) - np.cos(n*phi_s + ang_off) + (phi - phi_s)*np.sin(n*phi_s + ang_off))
		H_int_fn = lambda phi, phi_s, coef, H0, int_coef: 1/(-int_coef*H0_minus_pot(phi, phi_s, coef,H0))**0.5
		


		
		Qs_amp_a = 1/Ts_amp_a

		
		#ph_hat_a = ph_a_right_deg - phi_s_deg #right crossing point
		ph_hat_a = ph_right_sol_deg #- phi_s_deg
		Qs_amp_pred = Qs_amp_fn(ph_hat_a*math.pi/180, phi_s, Qs)

		plt.plot(ph_hat_a, Qs_amp_a, 'ko',label='num. integration')
		plt.plot(ph_hat_a, Qs_amp_pred, 'r',linewidth=2,label='formula')
		#plt.xlim(0,10)
		#plt.ylim(200,210)
		plt.ylim(ymin=0)
		plt.xlabel(r'$\hat{\phi}$ [deg]')
		plt.ylabel('synchrotron tune')
		plt.legend()
		plt.savefig('synch_tune_amp')
		plt.show()

		
		


		sys.exit()

	check_nonlin_fn = True
	if check_nonlin_fn:
		h3_angle = 1.2
		ntest = 5000
		ph_a_test = np.linspace(-math.pi, math.pi,ntest)
		ph_a_test_deg = (180/math.pi)*ph_a_test

		
		gnorm = -np.sin(ph_a_test) 			
		g = -(np.sin(ph_a_test) + 0.1*np.sin(3*ph_a_test + h3_angle))
		i1 = 0
		for g1 in g:
			if i1 > 0:
				if g1 < 0 and g[i1-1] >= 0:
					ic = i1
					
	
				if g1 < -np.sin(phi_s) and g[i1-1] >= -np.sin(phi_s):
					isyn = i1
					break

			i1 = i1 + 1





		print "ic ",ic, isyn
		print "crossing phase ",ph_a_test[ic], (180/math.pi)*ph_a_test[ic]
		print "sync phase ",ph_a_test[isyn], (180/math.pi)*ph_a_test[isyn]
		plt.axhline(y=-np.sin(20*math.pi/180))
		plt.plot(ph_a_test_deg, -np.sin(ph_a_test), 'k')
		plt.plot(ph_a_test_deg, g, 'r')
		plt.savefig('pot_fn')
		plt.show()
		sys.exit()

		
	
	foil_boundary = False
	if not foil_boundary:
		dp_bin_ini = np.linspace(-vfac*bucket_height,vfac*bucket_height, n_slices)
		ph_mesh, dp_mesh = np.meshgrid(ph_a, dp_bin_ini, sparse = True)
		#Hgrid = H_fn(ph_mesh, dp_mesh, phi_s, a1, a2)
		Hgrid = H_fn(ph_mesh, dp_mesh, phi_s, a1, a2) #+ H_fn_gen_pot(3, ph_mesh, phi_s, 0.1*a2, h3_angle)
	else:
		dp_bin_ini_upper = np.linspace(0,bucket_height, n_slices)
		dp_bin_ini_lower = np.linspace(-bucket_height,0, n_slices)
		ph_mesh_upper, dp_mesh_upper = np.meshgrid(ph_a, dp_bin_ini_upper, sparse = True)
		ph_mesh_lower, dp_mesh_lower = np.meshgrid(ph_a, dp_bin_ini_lower, sparse = True)
		Hgrid_upper = H_fn(ph_mesh_upper, dp_mesh_upper, phi_s, a1, a2)
		Hgrid_lower = H_fn(ph_mesh_lower, dp_mesh_lower, phi_s_nom, a1, a2)
		Hgrid = Hgrid_upper + Hgrid_lower

		

	#binwidth_time = 9.19e-9
	#print "time per bin ",binwidth_time
	#flip sign of a1 in case of large phi_s?
	
	
	print "Hsep ",Hsep
	if Hsep > 0:
		Hgrid_inv = Hsep - Hgrid
	else:
		Hgrid_inv = -Hsep + Hgrid
		
	#plt.plot(Hgrid_inv[0],'r')
	#plt.plot(Hgrid[int(0.5*n_slices)],'k')
	#plt.plot(Hgrid_inv[int(0.5*n_slices)],'k')
	#plt.plot(Hgrid_inv[15],'g')
	#plt.plot(Hgrid_inv[-1],'r')
	#plt.show()
	#sys.exit()

	if foil_boundary:
		print ph_mesh_upper.shape, Hgrid_upper.shape
		
		plt.contour(ph_a_deg, dp_bin_ini_upper, Hgrid_upper, 10)
		plt.contour(ph_a_deg, dp_bin_ini_lower, Hgrid_lower, 10)
		plt.xlim(ph_a_deg[0], ph_a_deg[-1])
		plt.show()
		
		#sys.exit()

	#exclude points outsize separatrix
	Hexcl1_ind = np.where(Hgrid_inv < 0)
	#Hexcl1_ind = np.where(Hgrid_inv > 100)
	#Hexcl1_ind = np.where(-Hgrid < Hsep)
	#print "Hlarge_ind ",Hlarge_ind

	Hgrid_inv[Hexcl1_ind] = 0
	
	#exclude another region of Hamiltonian
	exclude_cond2 = False
	if exclude_cond2:
		Hexcl2_ind = np.where(Hgrid_inv < 0)
		#Hexcl2_ind = np.where(Hgrid_inv > 250)
		Hgrid_inv[Hexcl2_ind] = 0

	
	#gdist =  np.sqrt(Hgrid_inv)
	gdist = Hgrid_inv
	#gdist = Hgrid_inv**2
	#gdist = Hgrid_inv*(np.exp(-Hgrid_inv/Hsep) - 1)
	
	#plt.plot(gdist[int(0.5*n_slices)],'k')
	#plt.ylabel('g')
	#plt.show()

	hist_ini = np.sum(gdist, axis=0)
	hist_dp_ini = np.sum(gdist, axis=1)

	t_phi_s = phi_s*(1e9*T)/(2*math.pi)
	t_gen_ns = ph_a*(1e9*T)/(2*math.pi)
	
	#bunch_t_ini_ns = np.tile(t_gen_ns, (n_slices, 1))
	#print  bunch_t_ini_ns
	#sys.exit()
	


	show_hist = False
	if show_hist:
		phase_mean_wgt =  sum([p1*t1 for p1,t1 in zip(ph_a_deg, hist_ini)])/sum(hist_ini)
		print "mean phase ", phase_mean_wgt
		#plt.subplot(211)
		plt.plot(ph_a_deg, hist_ini, 'ko-')
		plt.axvline(x = phi_s_deg)
		plt.axvline(x= phase_mean_wgt, color='g')
		#plt.subplot(212)
		#plt.plot(dp_bin_ini, hist_dp_ini, 'ko-')
		plt.xlabel('phase [deg]')
		plt.ylabel('intensity')
		plt.savefig('histogram_gen')
		plt.show()
	
	#longitudinal phase space
	#plt.contourf(ph_a_deg, dp_bin_ini, Hgrid, 20,cmap = 'RdGy')
	
	
	plt.contourf(ph_a_deg, dp_bin_ini, -gdist, 20,cmap = 'RdGy')
	plt.plot(phase_sep_deg, dp_sep, 'r')
	plt.plot(phase_sep_deg, -dp_sep, 'r')
	#plt.colorbar()
	plt.axvline(x = phi_s_deg, color='gray', linestyle = '--')
	#plt.ylim(-1.1*bucket_height, 1.1*bucket_height)
	plt.xlim(-65, 180)
	plt.xlabel('phase [deg]')
	plt.ylabel('dp/p')
	
	plt.savefig('phasespace_gen')
	plt.show()
	
	
	sys.exit()
	

	fb = open('bunch_1D.txt','w')
	for i in range(nturns):
			for n1 in hist_ini:
				print >>fb, n1
	fb.close()

else:
	run_pyheadtail = True


def exec_pyheadtail(z_bin, modify_bunch = False, z_shift = 0, plot_tracking = False, sel_xaxis = 'phase'):
	"""modify_bunch: If True, make a change to the bunch generated by PyHeadtail"""

	set_bin_edges = True
	use_slicer = False
	
	rfsystems = RFSystems(circumference, [harmonic], [V_rf], [phi_offset],
		                  alpha_c_arr, gamma, p_increment=p_increment_0, charge=e, mass=m_p)
	rfbucket = rfsystems.get_bucket(gamma=gamma)

	

	# Generate the particle distribution
	use_particle_generator = False
	if use_particle_generator: 
		bunch = generators.ParticleGenerator(n_particles, intensity, e, m_p, circumference, gamma, 
	#		                    distribution_x=generators.gaussian2D(epsn_x), alpha_x=inj_alpha_x, beta_x=inj_beta_x,
	#		                    distribution_y=generators.gaussian2D(epsn_y), alpha_y=inj_alpha_y, beta_y=inj_beta_y,
				                distribution_z=generators.RF_bucket_distribution(rfbucket, sigma_z=sigma_z_0,margin=0)).generate()
	else:
		fake_beta_z = abs(eta)*circumference / (2 * np.pi * Qs)
		bunch = generators.generate_Gaussian6DTwiss(n_particles, intensity, e, m_p, circumference, gamma, inj_alpha_x, inj_alpha_y, inj_beta_x, inj_beta_y, fake_beta_z, epsn_x, epsn_y, epsn_z) 

	# Modify the particle distribution
	if modify_bunch:
		debunched = False
		hollow = False
		shift_bunch_z = True

		dp_cut = 0.003
		amp_cut = 0.4 #applied in debunched and hollow cases
		#z_shift = 0# 5

		zm = []
		dpm = []
		for z, dp in zip(bunch.z, bunch.dp):

			if debunched:
				if dp < dp_cut and dp > -dp_cut:
					dpm.append(dp)
					zm.append(z)
			elif hollow:		
				amp = ((z/max_z)**2 + (dp/max_dp)**2)**0.5
				if amp > amp_cut:
					dpm.append(dp)
					zm.append(z)
			elif shift_bunch_z:
				dpm.append(0)
				zm.append(z_shift)
				
		bunch.z = np.array(zm)
		bunch.dp = np.array(dpm)
	


	if plot_tracking:

		if animate:
			plt.ion()#fix ordering

		fig = plt.figure(1)

		plotevery = 25
		plot_vaxis_p  = False #set vertical axis to p rather than dp
		#sel_xaxis = 'phase' #'z' or 'phase' 

	p_ref_track = []
	ke_ref_track = []
	bunch_z_track = []
	prof_data_track = []
	track_count = 0
	for iturn in np.arange(0, nturns):

		if include_foil:
			if ke_ref_track[-1] > ke_threshold and not foil_escape:
				print "reference kinetic energy above foil crossing threshold"
				print "p increment old and new ",p_increment_0, p_increment_nom
				rfsystems = RFSystems(circumference, [harmonic], [V_rf], [phi_offset],
		                  alpha_c_arr, gamma, p_increment=p_increment_nom, charge=e, mass=m_p)
				foil_escape = True

		# track the particles (no tracking on turn 0)
		if iturn > 0:
			rfsystems.track(bunch)

		p_ref = c*bunch.p0/e
		p0_now = 1e-6*p_ref 
		ke_ref = mom_to_ke(p_ref, m_p_ev) #ref. ke energy			
		p_ref_track.append(p_ref)
		ke_ref_track.append(ke_ref)

		if include_foil and not foil_escape:
			
			ke_ref -= E_loss
			p0_foil = ke_to_mom(ke_ref, m_p_ev)
			bunch.p0 = p0_foil*e/c

			#print "bunch.p0 ", bunch.p0, bunch.beta
			#print "max, min bunch.dp ",max(bunch.dp), min(bunch.dp)

			update_bunch = False
			if update_bunch:				
				bunch_deltap = bunch.dp*p_ref
				bunch_p = p_ref + bunch_deltap
				bunch_ke = mom_to_ke(bunch_p, m_p_ev) 	


				bunch_ke -= E_loss
				
				bunch_p_foil = ke_to_mom(bunch_ke, m_p_ev)

				p_ref_delta = p_ref - p0_foil
				bunch_dp_foil = (bunch_p_foil-p0_foil)/p0_foil
				#bunch_fixdp_foil = bunch.dp*(p_ref/p0_foil)


		bucket = rfsystems.get_bucket(gamma=bunch.gamma)			
		bunch_z_track.append(bunch.z[0])

		if set_bin_edges:
			hist, bin_edges_z = np.histogram(bunch.z, bins = z_bins_ideal, density=False)
		else:
			if use_slicer:
				slice_set = slicer.slice(bunch) 
				z_bin = slice_set.z_centers
				hist = slice_set.n_macroparticles_per_slice
				t_bin =  1e9*(z_bin*(harmonic*T/circumference))	
				phase_bin = 2*math.pi*t_bin/(T*1e9)	
			else:				
				hist, bin_edges_z = np.histogram(bunch.z, bins = n_slices, density=False)


		#fix ordering
		if z_bin[0] > z_bin[-1]:
			z_bin = z_bin[::-1]
			hist = hist[::-1]
			
		if iturn%dturns == 0:
			prof_data_track.append(hist[::-1])
			track_count = track_count + 1

			if dturns > 1:
				print "print at turn ",i+1

		if plot_tracking and iturn%plotevery == 0:
			# monitor the particles


			# plot the RF bucket envelope
			z = np.linspace(*bucket.interval, num=100)
			dp = np.linspace(-1.2*bh, 1.2*bh, num=100)
			ZZ, DPP = np.meshgrid(z, dp)
			HH = bucket.hamiltonian(ZZ, DPP)

			phase_grid = z*(-harmonic*2*np.pi/circumference)
			phph, DPP = np.meshgrid(phase_grid, dp)

			if sel_xaxis == 'phase':
				xlabel = 'phase'
				bunch_phase = bunch.z*(-harmonic*2*np.pi/circumference)
				xaxis = bunch_phase*180/math.pi
				binax = z_bin*(-harmonic*2*np.pi/circumference)*180/math.pi
				xcont = phph*180/math.pi
			elif sel_xaxis == 'z':
				xlabel = 'z [m]'
				xaxis = bunch.z
				binax = z_bin
				xcont = ZZ

	
			if plot_tracking:

				if plot_vaxis_p:
					#bunch_p = bunch.dp + p0_now
					bunch_p = bunch.dp*p_ref
					DPP = DPP + p0_now

				# plot the particles in phase space
				#plt.subplot(211)

				if plot_vaxis_p:
					plt.plot(bunch.z, bunch_p, 'o')
				else:
					#plt.plot(bunch.z, bunch.dp, 'o')
					plt.plot(xaxis, bunch.dp, 'o')

				if iturn == 0:
					plt.contour(xcont, DPP, HH, levels=[0], colors='magenta')
					#plt.xlim(xcont[0], xcont[-1])
					

				if plot_vaxis_p:
					plt.ylabel('p [MeV/c]')
				else:	
					plt.ylim(-bh, bh)
					plt.ylabel('Delta p/p')
				#plt.subplot(212)

				#plt.plot(binax, hist)

				#if i == 0:
				#	plt.xlim(min(binax), max(binax))
				#	plt.xlabel(xlabel)
				
				if animate:
					plt.show()
					plt.pause(0.1)
				# plt.cla()

	if plot_tracking:
		
		if sel_xaxis == 'z':
			plt.xlabel('z[m]')
		elif sel_xaxis == 'phase':
			plt.xlabel('phase [deg]')
		plt.savefig('tracking')
		plt.show()
		
		if animate:
			plt.ioff()

	prof_data_a = np.array(prof_data_track)

	return prof_data_a, bunch_z_track

def mountain_range(prof_a, sel_xaxis):
	phase_bin =	-z_bin[::-1]*(harmonic*360/circumference)
	fig = plt.figure(figsize=(8, 8), facecolor='w')
	ax = plt.subplot(111, frameon=False)
	it = 0
	
	print "len prof_a ",len(prof_a)
	for prof_data in prof_a:

		lw = 0.5 - (it / 3000)
		#ax.plot(range(1,n_slices+1), prof_data + 0.1*it,'k', lw=lw)
		if sel_xaxis == 'z':
			ax.plot(z_bin, prof_data[::-1] + 0.1*it,'k', lw=lw)
			plt.xlabel('z[m]')
		else:
			ax.plot(phase_bin, prof_data + 0.1*it,'k', lw=lw)
			plt.xlabel('phase [deg]')
		it = it + 1

	plt.show()

def freq_analysis(bunch_z_track):

	print 'freq_analysis'
	#data = index_max_turn 
	z_data = np.array(bunch_z_track)
	phase_data = -z_data*harmonic*360/circumference

	z_data_min1 = min(z_data[:400])
	z_data_max1 = max(z_data[:400])
	phase_data_min1 = min(phase_data[:400])
	phase_data_max1 = max(phase_data[:400])

	turn_a = np.arange(len(z_data))
	zlinfit = np.poly1d(np.polyfit(turn_a, z_data,3))
	zfit_a = zlinfit(turn_a)
	#print "min/max z data first 400 turns ",z_data_min1, z_data_max1
	#print "min/max phase data first 400 turns ",phase_data_min1, phase_data_max1
	plt.plot(phase_data,'ko-')
	plt.axhline(y=phase_data_min1)
	plt.axhline(y=phase_data_max1)
	plt.show()
	
	plt.plot(turn_a, z_data - zfit_a,'ko')
	plt.xlim(xmax=400)
	#plt.plot(turn_a, zfit_a,'r-')
	plt.show()
	
	
	z_data_mod = z_data - zfit_a
	fft_data = np.fft.fft(z_data_mod)
	fft_max_abs = np.abs(fft_data)

	#plt.plot(fft_max_abs)
	#plt.show()
	fft_max_index = np.argmax(fft_max_abs[1:int(0.5*nturns)])+1
	fft_period = nturns/fft_max_index
	fft_tune = 1/fft_period
	fft_phase = np.angle(fft_data[fft_max_index])

	#print "max FFT, period, turns, resolution ",fft_max_index, fft_period, nturns, 1/nturns

	nxlim = 2*nturns/200

		
	#plt.plot(fft_data)	
	#plt.xlim(0,nxlim)
	#plt.show()

	
	#print "fft_tune, resolution ",fft_tune, 1/nturns
	#print "fft_period, range within error ",fft_period, 1/(fft_tune + 1/nturns), 1/(fft_tune - 1/nturns)
	#print "len bunch_z_track ",len(bunch_z_track)

	g_mean =  0.5*(max(z_data) + min(z_data))
	g_amp =  0.5*(max(z_data) - min(z_data))
	g_phase = fft_phase
	g_freq = fft_tune
	

	t = np.arange(nturns)

	optimize_func = lambda x: x[0]*np.sin(2*math.pi*x[1]*t+x[2]) + x[3] - z_data_mod
	#optimize_func = lambda x: g_amp*np.sin(2*math.pi*x[0]*t+x[1]) + g_mean - data

	est_amp, est_freq, est_phase, est_mean = leastsq(optimize_func, [g_amp, g_freq, g_phase, g_mean])[0]
	#est_freq, est_phase = leastsq(optimize_func, [g_freq, g_phase])[0]	

	
	#print "g_amp ",g_amp
	#print "guess, optimised mean ",g_mean, est_mean
	#print "fft_phase, est_phase ", fft_phase, est_phase
	#print "fit freq, period ",est_freq, 1/est_freq

	data_guess = g_amp*np.sin(2*math.pi*g_freq*t+g_phase + 0.5*math.pi) + g_mean
	data_fit = est_amp*np.sin(2*math.pi*est_freq*t+est_phase) + est_mean

	show_fit = True
	if show_fit:
		plt.plot(z_data_mod,'ko-')
		plt.plot(data_guess,'r')
		plt.plot(data_fit,'g')
		plt.show()
		0
	

	return fft_tune, est_freq



run_pyheadtail = True
run_freq_analysis = True
calc_freq_phase = False
if run_pyheadtail:
	nrun = 2
	
	#specify z
	#z_shift_a = np.linspace(-1.6, 13, nrun)
	#phase_ini_a = z_shift_a*harmonic*360/circumference
	
	#specify phase
	phase_ini_a = np.linspace(20, 0.6*180*(math.pi - phi_s)/math.pi, nrun)
	z_shift_a = -circumference*phase_ini_a/(360*harmonic)
	

	fft_tune_l = []
	fit_tune_l = []
	for irun in np.arange(nrun):

		z_shift = z_shift_a[irun]

		#calculate bins in z
		zlim1 = -max(phase_sep)*circumference/(2*math.pi*harmonic)
		zlim2 = -min(phase_sep)*circumference/(2*math.pi*harmonic)
		z_bins_ideal = np.linspace(zlim1, zlim2, n_slices+1)
		z_bin = 0.5*np.array([(z_bins_ideal[i]+z_bins_ideal[i-1]) for i in range(1, n_slices+1)])
		
		#x axis of longitudinal phase space 'z' or 'phase'
		sel_xaxis = 'phase'	
	
		print "expected phase_ini, z_ini ", phase_ini_a[irun], z_shift_a[irun]
		prof_data_a, bunch_z_track = exec_pyheadtail(z_bin, modify_bunch = True, z_shift = z_shift, plot_tracking = True, sel_xaxis = sel_xaxis)
		#sys.exit()
		
		prof_data_norm_a = prof_data_a/np.max(prof_data_a)

		#mountain_range(prof_data_norm_a, sel_xaxis)
	
		if run_freq_analysis:
			fft_tune, fit_tune = freq_analysis(bunch_z_track)
			print "fft_tune, fit_tune ",fft_tune, fit_tune

			fft_tune_l.append(fft_tune)
			fit_tune_l.append(fit_tune)

		if calc_freq_phase:
			frequency_vs_phase(prof_data_a, Qs, synch_index_ideal)

if run_freq_analysis:
	fw = open('freq_vs_amp.txt','w')
	
	print >>fw, "phi_s ",phi_s, " nturns ",nturns
	
	for z, ph, tune1, tune2 in zip(z_shift_a, phase_ini_a, fft_tune_l, fit_tune_l):
		print >>fw, z, ph, tune1, tune2	
	fw.close()

	print "resolution ",1/nturns
	plt.plot(phase_ini_a, fft_tune_l, 'ko-')
	plt.plot(phase_ini_a, fit_tune_l, 'ro-')
	plt.show()

	fft_tune_fac = [fft_tune_l[0]/f for f in fft_tune_l]
	fit_tune_fac = [fit_tune_l[0]/f for f in fit_tune_l]
	plt.plot(phase_ini_a, fft_tune_fac,'ko-')
	plt.plot(phase_ini_a, fft_tune_fac,'ro-')
	plt.show()

sys.exit()
#results at different z_shift (phis=20 deg)
z_shift_l = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
z_halfrange_l = [1.61, 0.58, 0.48, 1.46, 2.41, 3.33, 4.23, 5.09, 5.91, 6.7, 7.4, 8.14]

period_l = [213.7, 211.8, 212.1, 213.5, 216.6, 221.5, 228.3, 237.6, 249, 264.8, 289.7, 339.7]

plt.plot(z_halfrange_l, period_l,'ko')
plt.xlabel('z_range/2')
plt.ylabel('synchrotron period')
plt.savefig('period_vs_amp')
plt.show()



for prof_data in prof_data_a:
	for i in range(binslowerignore):
		print >>fb, "0"
	for n1 in prof_data:
		print >>fb, n1
	for i in range(binsupperignore):
		print >>fb, "0"
fb.close()



run_tomo_code = True

if run_tomo_code:

	fname = 'input_v2.dat'
	descrip = 'KURRI'
	datafile = 'bunch_1D.txt'
	
	if run_pyheadtail and track:
		input_settings["numframes"] = track_count
	else:
		input_settings["numframes"] = nturns
		
	input_settings["numbins"] = n_slices + binslowerignore + binsupperignore
	input_settings["synchbin"] = synch_index_ideal#int(0.5*n_slices)
	input_settings["framebinwidth"] = binwidth_time
	input_settings["binslowerignore"] = binslowerignore
	input_settings["binsupperignore"] = binsupperignore
	input_settings["filmstart"] = recon_start
	input_settings["filmstop"] = recon_stop
	input_settings["filmstep"] = recon_step
	input_settings['dturns'] = dturns
	input_settings['Bdot'] = Bdot
	input_settings['VRF1ref'] = V_rf
	input_settings['Rnom'] = Rnom
	input_settings['rhonom'] = bending_radius
	input_settings['gammatnom'] = gamma_tr
	
	input_settings['pickup'] = 1.0
	input_settings['selffields'] = 0
	tdefs.write_inputfile(fname,descrip, datafile, input_settings)


	#temp_files = ['plotinfo.data','profiles.data','d001.data','jmax.data','image001.data']
	#for file1 in temp_files:
	#	if os.path.isfile(file1):
	#		os.remove(file1)
	#		print "delete ",file1

	time_cerncode0 = time.time()
	call(["tomo"])
	time_cerncode1 = time.time()	
	print "Tomography code execution time ",time_cerncode1-time_cerncode0

	files = os.listdir(os.curdir)
	for fn in files:
		if fn[-4:] == 'data':
			os.rename(fn,tomo_outdir+"/"+fn)



def read_data_file(filen):
	
	f1 = open(filen,'r')
	data_l = []
	for data in f1:
		data_l.append(float(data))
	data_a = np.array(data_l)

	return data_a


show_tomo_result = True
if show_tomo_result:

	#read output parameters from plotinfo.data
	#plotinfo = tdefs.get_plotinfo()
	plotinfo = tdefs.get_plotinfo(input_file=tomo_outdir+"/plotinfo.data")
	x0 = plotinfo['xat0'][0]
	y0 = plotinfo['yat0'][0]
	dEbin =  plotinfo['dEbin'][0]
	dtbin =  plotinfo['dtbin'][0]
	phi_s_out = plotinfo['phi_s'][0] #synch phase

	print "dtbin ",dtbin

	fname_prof = 'profiles.data'
	profiles_norm_raw = read_data_file(tomo_outdir+'/'+fname_prof)
	print "nturns ",nturns
	print "len profiles_norm_raw ",len(profiles_norm_raw)

	profiles_norm_a = np.split(profiles_norm_raw, nturns)

	show_profile_out = False
	if show_profile_out:

		fname_prof = 'profiles.data'
		fname_raw = 'bunch_1D.txt'

		bunchdata= read_data_file(fname_raw)
		
					
			
			

		nprof = len(bunchdata)/n_slices
		print len(profiles_a)

		print "nprof ",nprof

		raw_spl = np.split(bunchdata, nprof)
		prof_spl = np.split(profiles_a, nprof)

		#plt.plot(profiles_a)
		for raw, prof in zip(raw_spl, prof_spl):
			plt.plot(raw,'k')
			plt.plot(prof,'ro')
		plt.show()

	

	inputp = tdefs.get_input()
	#E0, gamma0, eta0, beta0, omega_rev0 = tdefs.set_param(inputp)
	E0, gamma0, eta0, beta0, omega_rev0, Vpeak = tdefs.set_param(inputp)
	T = 2*np.pi/omega_rev0 #rf period

	#dp_sep, t_phi_ns = tdefs.separatrix(Qs, phi_s, eta0, T, 1)

	max_dp_sep = max(dp_sep)

	recon_id = np.arange(recon_start, recon_stop+1, recon_step)
	
	ndim = n_slices
	image_scaled = False #scale image data by measured profile integral
	image_all_a = []
	image_std_l = []
	image_max_l = []
	#chisq_l = []
	iloop = 0 
	for irc in recon_id:
			print "recon at turn ",irc

			fname_img = tomo_outdir+'/image000'[:-len(str(irc))]+ str(irc) + '.data'
			print "fname_img ",fname_img

			image_data_a = read_data_file(fname_img) #1D data
			image_data_trow = np.split(image_data_a, ndim) #each row corresponds to time bin
			image_data_prow = np.transpose(image_data_trow) #each row corresponds to dp/p bin

			if image_scaled:
				intensity_scale.append(intensity_totomo[irc-1])
				image_std_l.append(np.std(image_data_prow)*intensity_totomo[irc-1])
				image_max_l.append(np.amax(image_data_prow)*intensity_totomo[irc-1])
			else:
				image_std_l.append(np.std(image_data_prow))
				image_max_l.append(np.amax(image_data_prow))

			tomo_time_proj = [sum(d) for d in image_data_trow]
			tomo_mom_proj = [sum(d) for d in image_data_prow]

			image_all_a.append(image_data_prow)		

			if iloop == 0:
				tomo_time_proj_a = np.array([tomo_time_proj])
			else:
				tomo_time_proj_a = np.vstack((tomo_time_proj_a, tomo_time_proj))

			iloop = iloop + 1


	#process main tomography output file
	#fname_img = 'image001.data'
	#image_data_a = read_data_file(fname_img) #1D data	
	
	#len_img_data = len(image_data_a)
	#ndim = int(len_img_data**0.5)

	dE_a_eV = dEbin*(np.array(range(ndim)) - y0)
	dE_a_mev = 1e-6*dE_a_eV
	
	E0 = Ekin + m_p_ev #Total energy (including rest mass)

	print "phi_s_out, phi_s_0 ",phi_s_out, phi_s
	print "phi_s_out, phi_s_0 (deg) ",(180/math.pi)*phi_s_out, (180/math.pi)*phi_s
	E_tot_eV_a = E0 + dE_a_eV 
	p_a = [1e-6*(Etot**2 - m_p_ev**2)**0.5 for Etot in E_tot_eV_a]

	#translate to time, dp
	dp_a = (p_a - p0_mev)/p0_mev
	t_synch = 1e9*phi_s_out*T/(2*math.pi)
	t_tomo_ns = 1e9*dtbin*(np.array(range(ndim)) - x0) + t_synch

	ph_tomo = (2*math.pi)*t_tomo_ns/(1e9*T)
	ph_tomo_deg = (180/math.pi)*ph_tomo

	#split into array
	#image_data_spl = np.split(image_data_a, ndim)
	#transpose to get phase on horizontal axis, energy on vertical
	#image_data_tr = np.transpose(image_data_spl)

	#tomo_profile = [sum(d) for d in image_data_tr]
	#tomo_time_proj = [sum(d) for d in image_data_spl]
	#tomo_mom_proj = [sum(d) for d in image_data_tr]

	#read in image files


cont_limit = max(image_max_l)



if animate:


	#plt.subplot(211)
	levels = np.linspace(0, cont_limit, 30)
	#plt.title('reconstructed turn '+sttomo_time_proj_lr(irc)+ ' time '+time_str+' ms')
	fig, ax = plt.subplots()
	cont = ax.contour(t_tomo_ns, dp_a, image_all_a[0], 40)
	def animate(i):	
		ax.clear()	
		if image_scaled:
			cont = ax.contour(t_tomo_ns, dp_a, intensity_scale[i]*image_all_a[i], levels = levels)
		else:
			cont = ax.contour(t_tomo_ns, dp_a, image_all_a[i], levels = levels)

		ax.plot(t_sep_ns, dp_sep,'r') #bucket
		ax.plot(t_sep_ns, -dp_sep,'r') 
		ax.set_xlabel('time [ns]')
		ax.set_ylabel(r"$\Delta$p/p")
		#ax.set_xlim(1.01*t_sep_ns[0],1.01*t_sep_ns[-1])
		ax.set_ylim(1.01*bh, -1.01*bh)
		#cont.set_zdata(image_all_a[1])
		return cont


	Nt = len(image_all_a)
	anim = animation.FuncAnimation(fig, animate, frames=Nt, interval=500, blit=False)
	anim.save('tomo.gif', dpi=80, writer='imagemagick')
	#anim.save('tomo.mp4', fps=30, extra_args=['-vcodec', 'libx264'])
	#anim.save("tomo.mp4", writer=writer)
	plt.show()

else:
	
	plot_surface = False
	plot_contour = True
	
	if plot_surface:
		from mpl_toolkits.mplot3d import Axes3D
		t_mesh, dp_mesh = np.meshgrid(t_tomo_ns, dp_a)

		print "len intensity_totomo ",len(intensity_totomo)
		print "image max over all profiles ",max(image_max_l)

		ipl = 0
		for irc in recon_id:
		
			#print "image sum ",np.sum(image_all_a[ipl])

			
			#print "measured profile integral ",intensity_totomo[irc-1]
			
			if image_scaled:
				image = image_all_a[ipl]*intensity_totomo[irc-1]
			else:
				image = image_all_a[ipl]

			fig = plt.figure()
			ax = fig.gca(projection='3d')			
			plt.title('reconstructed turn '+str(irc)+ ' time '+time_str+' ms')
			ax.plot_surface(t_mesh, dp_mesh, image, cmap="seismic") #or coolwarm

			ax.set_xlabel('time [ns]')
			ax.set_ylabel(r"$\Delta$p/p")
			
			
			ax.set_zlim(0, cont_limit)
			ipl = ipl + 1	

			plt.show()

		


	print "recon_id ",recon_id
	ipl = 0
	for irc in recon_id:

		if image_scaled:
			image = image_all_a[ipl]*intensity_totomo[irc-1]
		else:
			image = image_all_a[ipl]
		
		show_imshow_spy = False
		if show_imshow_spy:
			plt.subplot(211)
			plt.imshow(image, vmax = cont_limit)
			plt.subplot(212)
			plt.spy(image)
			plt.show()
			sys.exit()
		

		#levels = np.linspace(0, 0.008, 20)
		levels = np.linspace(0, cont_limit, 30)
		print "plot phasespace at turn ",irc
		
		plt.subplot(211)

		#plt.subplot(211)
		#plt.title('reconstructed turn '+str(irc)+ ' time '+time_str+' ms')
		plt.title('reconstructed turn '+str(irc))
		if plot_contour:
			ctom = plt.contour(t_tomo_ns, dp_a, image, levels = levels)
			#ctom = plt.contour(t_tomo_ns, dp_a, image)
			#plt.colorbar(ctom)
		
		#print "ctom ",ctom.levels	
		#sys.exit()

		#plt.contourf(t_gen_ns, dp_bin_ini, Hgrid_inv, ctom.levels)
		#plt.contourf(t_gen_ns, dp_bin_ini, Hgrid_inv)

		plt.plot(t_sep_ns, dp_sep,'r') #bucket
		plt.plot(t_sep_ns, -dp_sep,'r') #bucket

		#if run_pyheadtail:
		#	plt.contour(-TT, DPP, HH, levels=[0], colors='magenta')
		#counts, xbins, ybins,_ = plt.hist2d(t_tomo_ns, dp_a, bins=60, norm=LogNorm())
		#plt.xlim(-300,300)
		plt.xlim(1.01*t_sep_ns[0],1.01*t_sep_ns[-1])
		plt.ylim(1.01*bh, -1.01*bh)

		plt.axvline(x=t_tomo_ns[int(x0)])

		x1,x2,y1,y2 = plt.axis()

		plt.xlabel('time [ns]')
		plt.ylabel(r"$\Delta$p/p")

		#plt.contour(counts, 20, extent=[xbins.min(), xbins.max(), ybins.min(), ybins.max()])
		#plt.ylabel(r"$\Delta$p [MeV]")

		plt.subplot(212)

		#t_turn_data_ns_shift = t_turn_data_ns #- 20# t_shift_ns

		#plt.plot(t_turn_data_ns, prof_data_norm,'ko-',linewidth=2,label='data')
		#plt.plot(t_turn_data_ns, prof_data_norm_a[irc-1],'k.-',linewidth=2,label='data')

		plt.plot(t_tomo_ns, profiles_norm_a[irc-1],'ko-',linewidth=2,label='data')

		#plt.plot(t_tomo_ns, prof_data_norm,'ko-',linewidth=2,label='data')

		#print "len t_tomo_ns ",len(t_tomo_ns)
		#print "tomo_time_proj_a shape ",tomo_time_proj_a.shape
		
		plt.plot(t_tomo_ns, tomo_time_proj_a[ipl],'r-',linewidth=2,label='tomography')
		
		plt.axvline(x=min(t_sep_ns), color='gray')
		plt.axvline(x=max(t_sep_ns), color='gray')

		if binslowerignore:
			plt.axvline(x=t_turn_data_ns[binslowerignore] + 0.1*t_data_step, color='gray', linestyle='--')
		if binsupperignore > 0:
			plt.axvline(x=t_turn_data_ns[n_slices - binsupperignore]-0.1*t_data_step, color='gray', linestyle='--')		

		plt.axvline(x=t_synch, color='gray', linestyle = '--')
		#plt.ylim(-1.1*max_dp_sep, 1.1*max_dp_sep) 

		plt.xlim(1.01*t_sep_ns[0],1.01*t_sep_ns[-1])
		#plt.xlim(t_tomo_ns[0],t_tomo_ns[-1])
		#plt.xlim(t_sep_ns[0], t_sep_ns[-1])
		#plt.ylim(-1.1*max_dp_sep, 1.1*max_dp_sep) 
		plt.xlabel('time [ns]')
		plt.ylabel('normalised intensity')
		plt.legend(loc = 'upper left')
		plt.savefig('phasespace_'+str(irc))
		plt.tight_layout()
		plt.show()

		ipl = ipl + 1





