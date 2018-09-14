from __future__ import division
# import required libraries
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
import tomo_defs as tdefs


#main settings
phi_s_deg_nom = 20
phi_s_nom = phi_s_deg_nom*math.pi/180
n_slices = 10
binslowerignore = 5 #add zeroes before profile (earlier time)
binsupperignore = 1 #add zeroes after profile (later time)
nbturns=100#300#100
dturns = 1
V_rf = 4e3 # in [V]

phi_foil = lambda V, phi_s, deltaE: np.arcsin((V*np.sin(phi_s) + deltaE)/V)

input_settings = tdefs.get_input("input_default.dat")
plt.switch_backend('TkAgg')	


gam2pMev =  1e-6*m_p * (c**2)/e

# general simulation parameters
n_particles = 200000
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


Qs = np.sqrt(np.abs(eta) * V_rf / (2 * np.pi * beta**2 * Etot))
print('Qs: ' + str(Qs))
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

bh = ((2*Qs)/(harmonic*abs(eta)))
print "bucket height ",bh


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
dp_sep, t_phi_ns = tdefs.separatrix(phi_s, eta, T, harmonic, npts, a1, a2)
phase_sep = t_phi_ns*(2*math.pi)/(1e9*T)
phase_sep_deg = (180/math.pi)*phase_sep

#plt.plot(t_phi_ns, dp_sep)
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

run_pyheadtail = False

generate_dist = False
if generate_dist:
	
	print "generate distribution at phi_s [deg] ",phi_s*180/math.pi
	print "nominal phi_s ",phi_s_nom*180/math.pi
	synch_index = synch_index_ideal 
	
	vfac = 1.0

	foil_boundary = False
	if not foil_boundary:
		dp_bin_ini = np.linspace(-vfac*bucket_height,vfac*bucket_height, n_slices)
		ph_mesh, dp_mesh = np.meshgrid(ph_a, dp_bin_ini, sparse = True)
		Hgrid = H_fn(ph_mesh, dp_mesh, phi_s, a1, a2)
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
	plt.plot(Hgrid[int(0.5*n_slices)],'k')
	#plt.plot(Hgrid_inv[int(0.5*n_slices)],'k')
	#plt.plot(Hgrid_inv[15],'g')
	#plt.plot(Hgrid_inv[-1],'r')
	plt.show()
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

	show_hist = True
	if show_hist:
		#plt.subplot(211)
		plt.plot(t_gen_ns, hist_ini, 'ko-')
		plt.axvline(x = t_phi_s)
		#plt.subplot(212)
		#plt.plot(dp_bin_ini, hist_dp_ini, 'ko-')
		plt.xlabel('time [ns]')
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
	plt.xlabel('phase [deg]')
	plt.ylabel('dp/p')
	
	plt.savefig('phasespace_gen')
	plt.show()
	
	
	#sys.exit()
	

	fb = open('bunch_1D.txt','w')
	for i in range(nbturns + 1):
			for n1 in hist_ini:
				print >>fb, n1
	fb.close()

else:
	run_pyheadtail = True


#

if run_pyheadtail:

	

	fake_beta_z = abs(eta)*circumference / (2 * np.pi * Qs)
	epsn_z = 0.1*epsn_z

	# Define RF systems
	rfsystems = RFSystems(circumference, [harmonic], [V_rf], [phi_offset],
		                  alpha_c_arr, gamma, p_increment=p_increment_0, charge=e, mass=m_p)
	rfbucket = rfsystems.get_bucket(gamma=gamma)

	# Generate the particle distribution
	use_particle_generator = True
	if use_particle_generator: 
		bunch = generators.ParticleGenerator(n_particles, intensity, e, m_p, circumference, gamma, 
	#		                    distribution_x=generators.gaussian2D(epsn_x), alpha_x=inj_alpha_x, beta_x=inj_beta_x,
	#		                    distribution_y=generators.gaussian2D(epsn_y), alpha_y=inj_alpha_y, beta_y=inj_beta_y,
				                distribution_z=generators.RF_bucket_distribution(rfbucket, sigma_z=sigma_z_0,margin=0)).generate()
	else:
		bunch = generators.generate_Gaussian6DTwiss(n_particles, intensity, e, m_p, circumference, gamma, inj_alpha_x, inj_alpha_y, inj_beta_x, inj_beta_y, fake_beta_z, epsn_x, epsn_y, epsn_z) 


	print "synch phase (deg) ",(180/math.pi)*rfsystems.phi_s(bunch.gamma, bunch.charge)


	voltage_mismatch = False
	if voltage_mismatch:
		# Voltage mismatch
		Vfactor=10
		#Vfactor=0.3

		rfsystems = RFSystems(circumference, [harmonic], [V_rf*Vfactor], [phi_offset],
				              alpha_c_arr, gamma, p_increment=p_increment_0, charge=e, mass=m_p)



	#print getattr(bunch, 'mean_x')

	#stats_to_store = [
	#            'mean_x', 'mean_xp', 'mean_y', 'mean_yp', 'mean_z', 'mean_dp',
	#            'sigma_x', 'sigma_y', 'sigma_z', 'sigma_dp', 'epsn_x', 'epsn_y',
	#            'epsn_z', 'macroparticlenumber' ]


	print "#bunch sigma_z, sigma_dp, macroparticle number ",bunch.sigma_z(), bunch.sigma_dp(),bunch.macroparticlenumber

	#bunch monitor saves above stats_to_store to file
	#bunchmonitor = BunchMonitor('bunch',10)

	#Alternatively UniformChargeSlicer divides beams so each bin has almost equal charge

	#slice_set = UniformBinSlicer(n_slices).slice(bunch) 
	 
	#print "dir(bunch) ",dir(bunch)

	max_dp = max(bunch.dp)
	min_dp = min(bunch.dp)
	max_z = max(bunch.z)
	min_z = min(bunch.z)

	modify_bunch = False
	if modify_bunch:

		zm = []
		dpm = []

		debunched = False
		dp_cut = 0.003

		amp_cut = 0.4
		for z, dp in zip(bunch.z, bunch.dp):

			if debunched:
				if dp < dp_cut and dp > -dp_cut:
					dpm.append(dp)
					zm.append(z)
			else:		
				amp = ((z/max_z)**2 + (dp/max_dp)**2)**0.5
				if amp > amp_cut:
					dpm.append(dp)
					zm.append(z)

		bunch.z = np.array(zm)
		bunch.dp = np.array(dpm)

	

	bucket = rfsystems.get_bucket(gamma=bunch.gamma)
	print "bucket.interval ",bucket.interval
	
	#Calculate histogram in z, dp
	set_bin_edges = True
	use_slicer = False
	if set_bin_edges:
		zlim1 = -max(phase_sep)*circumference/(2*math.pi*harmonic)
		zlim2 = -min(phase_sep)*circumference/(2*math.pi*harmonic)
		z_bins_ideal = np.linspace(zlim1, zlim2, n_slices+1)
		
		hist_ini, bin_edges_z = np.histogram(bunch.z, bins = z_bins_ideal, density=False)
		
		z_bin_ini = 0.5*np.array([(bin_edges_z[i]+bin_edges_z[i-1]) for i in range(1, n_slices+1)])	
	
	else:
		if use_slicer:
			slicer = UniformBinSlicer(n_slices=n_slices)
			fb = open('bunch_1D.txt','w')
			slice_set = slicer.slice(bunch)	
			z_bin_ini = slice_set.z_centers
			hist_ini = slice_set.n_macroparticles_per_slice
		else:
			hist_ini, bin_edges_z = np.histogram(bunch.z, bins = n_slices, density=False)
					

	if set_bin_edges or not use_slicer:
		z_bin_ini = 0.5*np.array([(bin_edges_z[i]+bin_edges_z[i-1]) for i in range(1, n_slices+1)])

	#fix ordering
	if z_bin_ini[0] > z_bin_ini[-1]:
		z_bin_ini = z_bin_ini[::-1]
		hist_ini = hist_ini[::-1]

	t_bin_ini =  1e9*(-z_bin_ini*(harmonic*T/circumference))	
	
	print "z_bin_ini limits ",z_bin_ini[0], z_bin_ini[-1]
	#slice using numpy histogram

	#momentum bins
	dp_bins_ideal = np.linspace(-bh, bh, n_slices+1)
	hist_dp_ini, bin_edges_dp = np.histogram(bunch.dp, bins = dp_bins_ideal, density=False)
	dp_bin_ini = 0.5*np.array([(bin_edges_dp[i]+bin_edges_dp[i-1]) for i in range(1, n_slices+1)])  
	

	phase_bin_ini = 2*math.pi*t_bin_ini/(T*1e9)

	#binwidth_time = (t_bin_ini[1] - t_bin_ini[0])*1e-9
	#binwidth_time = 9e-9


	for i1, ph1 in enumerate(phase_bin_ini):
		if ph1 < phi_s:
			synch_index = i1
			print "synch_index (PyHeadtail) ",synch_index
			break




	show_hist = False
	if show_hist:	
		
		hist_vs_phase = True
		if hist_vs_phase:
			plt.subplot(211)
			#plt.plot(z_bin_ini, hist_ini, 'k.-')
			plt.plot(phase_bin_ini, hist_ini,'ko-')
			plt.axvline(phi_s, color='r')
			plt.axvline(phase_bin_ini[synch_index])

			plt.subplot(212)
			plt.plot(dp_bin_ini, hist_dp_ini, 'r.-')
		else:
			plt.plot(z_bin_ini, hist_ini, 'k.-')
			

		plt.show()
		sys.exit()
	
	#print dir(slice_set)
	#print "len bunch.dpm, bunch.z ",len(bunch.dp), len(bunch.z)
	#print slice_set.slice_positions
	
	fb = open('bunch_1D.txt','w')
	#write initial distribution to file
	for i in range(binslowerignore):
		print >>fb, "0"

	for n1 in hist_ini[::-1]:
		print >>fb, n1

	for i in range(binsupperignore):
		print >>fb, "0"


	#plt.plot(bunch.z, bunch.dp,'ko')
	#plt.plot(zm, dpm,'ro')
	#plt.show()



	print "max_dp, min_dp ",max_dp, min_dp
	print "min & max momentum in bunch ",(p0_mev*(1+min_dp)), (p0_mev*(1+max_dp))

	print "min & max bunch.z & bunch length ",min_z, max_z, max_z - min_z
	



	#plot initial phase space ()

	#z_cut = 3*bunch.sigma_z()
	#slicer = UniformBinSlicer(n_slices=64, z_cuts = (-z_cut, z_cut))
	
	print dir(bucket)


	#bucket_z_l = bucket.z_ufp[1] - bucket.z_ufp[0]


	#record initial distribution (tomography will try to recover this)
	bunch_z_ini = bunch.z
	bunch_dp_ini = bunch.dp


	# track and plot phase space ()


	#bucket coordinates for initial
	z = np.linspace(*bucket.interval, num=100)
	dp = np.linspace(-1.2*max_dp, 1.2*max_dp, num=100)
	ZZ, DPP = np.meshgrid(z, dp)
	HH = bucket.hamiltonian(ZZ, DPP)
	TT = 1e9*ZZ/(beta*c)
	

	#phase_TT =TT*(2*math.pi)/(T*1e9)

	#bunch_phase = bunch_t_ini_ns*(2*math.pi)/(T*1e9)
	bunch_phase =  bunch.z*(-harmonic*2*np.pi/circumference)*180/np.pi
	bunch_t_ini_ns =  1e9*(bunch.z*(-harmonic*T/circumference))

	print "T ",T
	print "max(bunch_phase) ",max(bunch_phase)
	#bunch_t_ini_ns = 1e9*(bunch_phase/(360))*T
	#bunch_t_ini_ns_orig = 1e9*bunch_z_ini/(beta*c)

	print "len bunch t ",len(bunch_t_ini_ns),len(bunch.dp) 

	phase_sep_pyh_rad = ZZ*(-harmonic*2*np.pi/circumference)
	phase_sep_pyh = phase_sep_pyh_rad*180/np.pi

	#plt.contour(ZZ*(-harmonic*2*np.pi/circumference)*180/np.pi, DPP, HH, levels=[0], colors='magenta')
	#plt.show()
	

	print "phi_s", phi_s
	print "ph_a[synch_index] ",ph_a[synch_index]
	print "ph[0], ph[-1] ",ph_a[0], ph_a[-1]
	print "phase at t_bin_ini[0], [-1] ",2*math.pi*t_bin_ini[0]/(T*1e9), 2*math.pi*t_bin_ini[-1]/(T*1e9)
	print "t_bin extremes ",t_bin_ini[0], t_bin_ini[-1]
	
	print "synch time, implied phase ",t_bin_ini[synch_index], 2*math.pi*t_bin_ini[synch_index]/(T*1e9)

	show_initial_distribution = False
	if show_initial_distribution:
		# plot the RF bucket envelope
		
		
		plot_type1 = False #dp vs z
		plot_type2 = True #dp vs phase
		#if both are False, dp vs time

		
		if plot_type1:
			plt.subplot(211)
			plt.plot(bunch.z, bunch.dp, 'o')
			#plt.hist2d(bunch_t_ini_ns, bunch_dp_ini, bins=60, norm=LogNorm())
			
			#plt.contour(counts, 20,extent=[xbins.min(), xbins.max(), ybins.min(), ybins.max()])

			

			plt.contour(ZZ, DPP, HH, levels=[0], colors='magenta')
			plt.ylim(-1.2*max_dp, 1.2*max_dp)
			#plt.xlim(bucke) =t.z_ufp[0], bucket.z_ufp[1])
			plt.ylabel('dp/p')
			plt.xlabel('bunch duration (ns)')
			plt.subplot(212)
			plt.plot(z_bin_ini, hist_ini)
			plt.ylabel('projection')
			plt.show()
			plt.close()

		elif plot_type2:
			
			plt.subplot(211)
			counts, xbins, ybins,_ = plt.hist2d(bunch_phase, bunch_dp_ini, bins=60, norm=LogNorm())
			#plt.plot(bunch.z, bunch.dp, 'o')
			#plt.hist2d(bunch_t_ini_ns, bunch_dp_ini, bins=60, norm=LogNorm())
			
			#plt.contour(counts, 20,extent=[xbins.min(), xbins.max(), ybins.min(), ybins.max()])

			#plt.plot(phase_sep, dp_sep, 'r')
			#plt.plot(phase_sep, -dp_sep, 'r')
			plt.contour(phase_sep_pyh, DPP, HH, levels=[0], colors='magenta')
			
			#plt.plot(bunch_phase, bunch.dp, 'o')
			plt.xlim(phase_sep[0]*180/math.pi, phase_sep[-1]*180/math.pi)

			plt.ylim(-bh, bh)
		
			#plt.xlim(bucket.z_ufp[0], bucket.z_ufp[1])
			#plt.xlim(-300, 250)
			#plt.xlim(-math.pi, math.pi)
			
			plt.ylabel('dp/p')
			plt.xlabel('bunch phase [deg]')

			
			
			plt.axvline(x = phi_s*180/math.pi, color='r')
			
			plt.subplot(212)
			plt.plot(z_bin_ini, hist_ini)
			plt.ylabel('projection')
			plt.show()
			plt.close()
		
		else:
			counts, xbins, ybins,_ = plt.hist2d(bunch_t_ini_ns, bunch_dp_ini, bins=60, norm=LogNorm())

			plt.contour(-TT, DPP, HH, levels=[0], colors='magenta')
			plt.show()

		sys.exit()

	plotevery = 25

	plot_tracking = False
	if plot_tracking:
		plt.ion()#fix ordering
		fig = plt.figure(1)

		plot_vaxis_p  = False #set vertical axis to p rather than dp

	track_count = 1
	
	
	E_loss = 760 #foil loss in eV
	foil_escape = False

	p_ref_ini = c*bunch.p0/e
	ke_ref_ini = mom_to_ke(p_ref_ini, m_p_ev)
	p_ref_track = [p_ref_ini]
	ke_ref_track = [ke_ref_ini]

	track = True
	if track:
		for i in np.arange(0, nbturns, 1):

			if include_foil:
				if ke_ref_track[-1] > ke_threshold and not foil_escape:
					print "reference kinetic energy above foil crossing threshold"
					print "p increment old and new ",p_increment_0, p_increment_nom
					rfsystems = RFSystems(circumference, [harmonic], [V_rf], [phi_offset],
		                  alpha_c_arr, gamma, p_increment=p_increment_nom, charge=e, mass=m_p)
					foil_escape = True

			# track the particles
			rfsystems.track(bunch)

			p_ref = c*bunch.p0/e
			ke_ref = mom_to_ke(p_ref, m_p_ev) #ref. ke energy

			p0_now = 1e-6*p_ref# gam2pMev*np.sqrt(bunch.gamma**2 - 1) 
			
			p_ref_track.append(p_ref)
			ke_ref_track.append(ke_ref)

			#print "i, p_ref, delta_p ",i, p0_now, p_ref - p_ref_track[i]
			

			if include_foil and not foil_escape:
				
				#print "bunch.p0 ", bunch.p0, bunch.beta
				#print "max, min bunch.dp ",max(bunch.dp), min(bunch.dp)
				#update reference momentum
				
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

					#print "max bunch_ke ", max(bunch_ke), max(bunch_ke)-ke_ref
					#sys.exit()
					#directly apply energy loss to bunch
					bunch_ke -= E_loss
				
					bunch_p_foil = ke_to_mom(bunch_ke, m_p_ev)

					p_ref_delta = p_ref - p0_foil
					bunch_dp_foil = (bunch_p_foil-p0_foil)/p0_foil
					#bunch_fixdp_foil = bunch.dp*(p_ref/p0_foil)



			
			bucket = rfsystems.get_bucket(gamma=bunch.gamma)			


			if set_bin_edges:
				hist, bin_edges_z = np.histogram(bunch.z, bins = z_bins_ideal, density=False)
				z_bin = z_bin_ini
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
			
		
			if i%dturns == 0:
				track_count = track_count + 1
				if dturns > 1:
					print "print at turn ",i+1
				

				for i in range(binslowerignore):
					print >>fb, "0"
				for n1 in hist[::-1]:
					print >>fb, n1 #+ np.random.normal(0, 10000)
				for i in range(binsupperignore):
					print >>fb, "0"

			if i%plotevery == 0:
				# monitor the particles


				# plot the RF bucket envelope
				z = np.linspace(*bucket.interval, num=100)
				dp = np.linspace(-1.2*bh, 1.2*bh, num=100)
				ZZ, DPP = np.meshgrid(z, dp)
				HH = bucket.hamiltonian(ZZ, DPP)

				

		
				
				if plot_tracking:

					if plot_vaxis_p:
						#bunch_p = bunch.dp + p0_now
						bunch_p = bunch.dp*p_ref
						DPP = DPP + p0_now

					# plot the particles in phase space
					plt.subplot(211)

					if plot_vaxis_p:
						plt.plot(bunch.z, bunch_p, 'o')
					else:
						plt.plot(bunch.z, bunch.dp, 'o')

					if i == 0:
						#plt.contour(ZZ, DPP, HH, levels=[0], colors='magenta')
						plt.xlim(min(z_bin), max(z_bin))
						
					


					if plot_vaxis_p:
						plt.ylabel('p [MeV/c]')
					else:	
						plt.ylim(-bh, bh)
						plt.ylabel('Delta p/p')
					plt.subplot(212)

					plt.plot(z_bin, hist)
					if i == 0:
						plt.xlim(min(z_bin), max(z_bin))
						plt.xlabel('z [m]')
					
					plt.show()
					plt.pause(0.1)
					# plt.cla()

	else:
		#instead of tracking, repeat initial distribution

		for i in range(nbturns):
			for i in range(binslowerignore):
				print >>fb, "0"

			for n1 in hist_ini:
				print >>fb, n1

			for i in range(binsupperignore):
				print >>fb, "0"

	#plt.ioff()
	fb.close()


	if plot_tracking:
		plt.savefig('tracking')
		plt.ioff()
		sys.exit()


	show_ke = True
	if show_ke:
		ke_ref_track_mev = [1e-6*ke for ke in ke_ref_track]

		track_t_us = 1e6*T*np.arange(len(p_ref_track))
		plt.plot(track_t_us, ke_ref_track_mev,'k',linewidth=2)
		plt.axhline(y=11.42, color='r')
		plt.ylabel('kinetic energy [MeV]')
		plt.show()


run_tomo_code = True

if run_tomo_code:

	fname = 'input_v2.dat'
	descrip = 'KURRI'
	datafile = 'bunch_1D.txt'
	
	if run_pyheadtail and track:
		input_settings["numframes"] = track_count
	else:
		input_settings["numframes"] = nbturns + 1
		
	input_settings["numbins"] = n_slices + binslowerignore + binsupperignore
	input_settings["synchbin"] = synch_index_ideal#int(0.5*n_slices)
	input_settings["framebinwidth"] = binwidth_time
	input_settings["binslowerignore"] = binslowerignore
	input_settings["binsupperignore"] = binsupperignore
	input_settings['dturns'] = dturns
	input_settings['Bdot'] = Bdot
	input_settings['VRF1ref'] = V_rf
	input_settings['Rnom'] = Rnom
	input_settings['rhonom'] = bending_radius
	input_settings['gammatnom'] = gamma_tr
	
	input_settings['pickup'] = 1.0
	input_settings['selffields'] = 0
	tdefs.write_inputfile(fname,descrip, datafile, input_settings)


	temp_files = ['plotinfo.data','profiles.data','d001.data','jmax.data','image001.data']
	for file1 in temp_files:
		if os.path.isfile(file1):
			os.remove(file1)
			print "delete ",file1

	call(["tomo"])

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
	plotinfo = tdefs.get_plotinfo()
	x0 = plotinfo['xat0'][0]
	y0 = plotinfo['yat0'][0]
	dEbin =  plotinfo['dEbin'][0]
	dtbin =  plotinfo['dtbin'][0]
	phi_s_out = plotinfo['phi_s'][0] #synch phase

	print "dtbin ",dtbin

	show_profile_out = False
	if show_profile_out:

		fname_prof = 'profiles.data'
		fname_raw = 'bunch_1D.txt'

		bunchdata= read_data_file(fname_raw)
		profiles_a = read_data_file(fname_prof)

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
	E0, gamma0, eta0, beta0, omegarev0 = tdefs.set_param(inputp)
	T = 2*np.pi/omegarev0 #rf period

	#dp_sep, t_phi_ns = tdefs.separatrix(Qs, phi_s, eta0, T, 1)

	max_dp_sep = max(dp_sep)

	#process main tomography output file
	fname_img = 'image001.data'
	image_data_a = read_data_file(fname_img) #1D data	
	
	len_img_data = len(image_data_a)
	ndim = int(len_img_data**0.5)

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
	image_data_spl = np.split(image_data_a, ndim)
	#transpose to get phase on horizontal axis, energy on vertical
	image_data_tr = np.transpose(image_data_spl)

	#tomo_profile = [sum(d) for d in image_data_tr]
	tomo_time_proj = [sum(d) for d in image_data_spl]
	tomo_mom_proj = [sum(d) for d in image_data_tr]

	if not run_pyheadtail:
		t_bin_ini = t_gen_ns

	#print "t_cen lims ",t_cen[0], t_cen[-1], t_phi_ns[0], t_phi_ns[-1]
	
	max_dist_1D_ini = max(hist_ini)	
	hist_ini_norm = [float(n1)/max_dist_1D_ini for n1 in hist_ini] 

	show_dp = True
	if show_dp:
	
		print "n_slices ",n_slices
	
		versus_phase = True
		#plt.axvline(x=0)
		#plt.plot(t_phi_ns, dp_sep,'r') #bucket
		#plt.plot(t_phi_ns, -dp_sep,'r') #bucket

		#if run_pyheadtail:

			#counts, xbins, ybins,_ = plt.hist2d(bunch_t_ini_ns, bunch_dp_ini, bins=60, norm=LogNorm())

		show_projection = False
		
		if show_projection:
			plt.subplot(311)
		else:
			plt.subplot(211)
			
		plt.title('original')

		if run_pyheadtail:
			if versus_phase:
				counts, xbins, ybins,_ = plt.hist2d(bunch_phase, bunch_dp_ini, bins=60)
			else:
				counts, xbins, ybins,_ = plt.hist2d(bunch_t_ini_ns, bunch_dp_ini, bins=40, cmap="RdGy")
		else:
			if versus_phase:
				plt.contourf(ph_a_deg, dp_bin_ini, -gdist, 20, cmap='RdGy')
			else:
				plt.contourf(t_gen_ns, dp_bin_ini, -gdist, 20, cmap='RdGy')

		if versus_phase:
			plt.plot(phase_sep_deg, dp_sep,'r') #bucket
			plt.plot(phase_sep_deg, -dp_sep,'r') #bucket
			plt.xlim(min(phase_sep_deg), max(phase_sep_deg))
			#plt.xlabel('phase [deg]')	
		else:
			plt.plot(t_phi_ns, dp_sep,'r') #bucket
			plt.plot(t_phi_ns, -dp_sep,'r') #bucket	
			plt.xlim(t_phi_ns[0],t_phi_ns[-1])

		plt.ylim(-bh, bh)
		plt.ylabel("dp/p")

		if show_projection:
			plt.subplot(312)
		else:
			plt.subplot(212)
			
		plt.title('reconstructed')
		
		if versus_phase:
			ctom = plt.contour(ph_tomo_deg, dp_a, image_data_tr, 20)
		else:
			ctom = plt.contour(t_tomo_ns, dp_a, image_data_tr, 20)
		#print "ctom ",ctom.levels	
		#sys.exit()

		#plt.contourf(t_gen_ns, dp_bin_ini, Hgrid_inv, ctom.levels)
		#plt.contourf(t_gen_ns, dp_bin_ini, Hgrid_inv)

		if versus_phase:
			plt.plot(phase_sep_deg, dp_sep,'r') #bucket
			plt.plot(phase_sep_deg, -dp_sep,'r') #bucket
			plt.xlim(min(phase_sep_deg), max(phase_sep_deg))
			plt.xlabel('phase [deg]')	
		else:
			plt.plot(t_phi_ns, dp_sep,'r') #bucket
			plt.plot(t_phi_ns, -dp_sep,'r') #bucket	
			plt.xlim(t_phi_ns[0],t_phi_ns[-1])
			plt.xlabel('time [ns]')
			

		x1,x2,y1,y2 = plt.axis()
		
		plt.ylabel("dp/p")


		if show_projection:
			plt.subplot(313)
		
			plt.plot(t_bin_ini, hist_ini_norm, 'ko-',label='input')
			plt.plot(t_tomo_ns, tomo_time_proj/max(tomo_time_proj),'ro-',label='tomography')

			#plt.xlim(t_phi_ns[0], t_phi_ns[-1])
			#plt.ylim(-1.1*max_dp_sep, 1.1*max_dp_sep) 


			plt.xlim(t_phi_ns[0],t_phi_ns[-1])
		#plt.xlim(t_phi_ns[0], t_phi_ns[-1])
		#plt.ylim(-1.1*max_dp_sep, 1.1*max_dp_sep) 

		plt.savefig('phasespace_compare')
		#plt.tight_layout()
		plt.show()
	else:
		plt.plot(t_phi_ns, dE_sep_mev,'r') #bucket
		plt.plot(t_phi_ns, -dE_sep_mev,'r') #bucket
		plt.contour(t_tomo_ns, dE_a_mev, image_data_tr, 20)
		plt.xlabel('[ns]')
		plt.ylabel(r"$\Delta$E [MeV]")
		plt.savefig('phasespace')
		plt.show()



	show_profiles = True
	if show_profiles:
		print "max(hist_dp_ini) ",max(hist_dp_ini)
		
		print "t_bin_ini min/max ",min(t_bin_ini), max(t_bin_ini)
		print "t_tomo_ns min/max ",min(t_tomo_ns), max(t_tomo_ns)
		
		
		max_hist_dp_ini  = max(hist_dp_ini)
		hist_dp_ini_norm = [h/max_hist_dp_ini for h in hist_dp_ini]

		#t_tomo_ns, dp_a t_bin_ini[::-1]
		show_norm = True
		if show_norm:
			plt.subplot(211)
			plt.plot(t_bin_ini, hist_ini_norm, 'ko-',label='original', linewidth=2)
			plt.plot(t_tomo_ns, tomo_time_proj/max(tomo_time_proj),'r-',label='tomography', linewidth=2)
			plt.xlabel('time [ns]')
			plt.ylabel('normalised intensity')
			plt.legend()
			plt.subplot(212)
			plt.plot(dp_bin_ini, hist_dp_ini_norm,'k-', linewidth=2)
			plt.plot(dp_a, tomo_mom_proj/max(tomo_mom_proj),'r-',label='tomography',linewidth=2)
			plt.xlim(-bh, bh)
			plt.xlabel('dp/p')
			plt.ylabel('normalised intensity')
			plt.tight_layout()
			plt.savefig('histcompare')
		else:
			plt.subplot(211)
			plt.plot(t_bin_ini, hist_ini, 'ko-',label='input')
			plt.plot(t_tomo_ns, tomo_time_proj,'ro-',label='tomography')
			plt.xlabel('time [ns]')
			plt.legend()
			plt.subplot(212)
			plt.plot(dp_bin_ini, hist_dp_ini,'ko-')
			plt.plot(dp_a, tomo_mom_proj,'ro-',label='tomography')
			plt.xlabel('dp/p')
			plt.savefig('histcompare')		
		
	plt.show()

