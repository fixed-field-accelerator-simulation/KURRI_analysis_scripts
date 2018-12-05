from __future__ import division
import sys
import numpy as np
import math
import pylab as plt
from scipy.signal import correlate
from scipy.signal import fftconvolve

ph_a = np.linspace(0,2*math.pi, 100)

phase_bin_deg = 90
phase_bin = phase_bin_deg*math.pi/180

print np.sin(phase_bin)

r = 1
circum = 2*math.pi*r
t1 = r*phase_bin
t2 = r*(2*math.pi - phase_bin)


def inverse_fft(ffta, indices):
	
	nffa = len(ffta)
	fftm = np.zeros(nffa, dtype = np.complex_)
	
	#fftm[0] = ffta[0]

	for index in indices:
		fftm[index] = ffta[index]
		fftm[nffa - index] = ffta[nffa-index]

	sigm = np.fft.ifft(fftm)

	return sigm
	#plt.plot(sigm[:100])
	#plt.show()
	#sys.exit()


nt = 3000
#db = 200*0.01
db = 0.01


def generate_signal(sep_frac, bunch_type = 'delta'):

	#sep_frac = 1.0
	
	varphi_a = np.linspace(0,2*math.pi, nt)

	sep = sep_frac*math.pi
	varphi_beam = [sep, 2*math.pi-sep]

	db_m = db#*np.cos(sep) 

	#fakedist = np.abs(np.cos(varphi_a))
	fakedist = np.ones(nt)
	#plt.plot(varphi_a, fakedist)
	#plt.show()
	#sys.exit()

	sig_l = []
	sig1_l = []
	sig2_l = []

	sigma_p1 = 0.25*math.pi

	bunch_type = 'g'
	if bunch_type == 'delta':
		for varphi, dist in zip(varphi_a, fakedist):
			if varphi > varphi_beam[0] - 0.5*db_m and varphi < varphi_beam[0] + 0.5*db_m:
				sig_l.append(dist)
				sig1_l.append(dist)
				sig2_l.append(0.0)
			elif varphi > varphi_beam[1] - 0.5*db_m and varphi < varphi_beam[1] + 0.5*db_m:
				sig_l.append(dist)
				sig1_l.append(0.0)
				sig2_l.append(dist)
			else:
				sig_l.append(0.0)
				sig1_l.append(0.0)
				sig2_l.append(0.0)

		sig_a = np.array(sig_l)
		sig1_a = np.array(sig1_l)
		sig2_a = np.array(sig2_l)

	else:
		#print "varphi_beam ",varphi_beam
		delphi1 = []
		for vphi in varphi_a:
			delp = abs(vphi - varphi_beam[0])
			if delp > math.pi:
				delphi1.append(2*math.pi - delp)
			else:
				delphi1.append(delp)
		
		delphi2 = []
		for vphi in varphi_a:
			delp = abs(vphi - varphi_beam[1])
			if delp > math.pi:
				delphi2.append(2*math.pi - delp)
			else:
				delphi2.append(delp)
			#else:
			#	delphi2.append(2*math.pi - vphi + varphi_beam[1])

		delphi1 = np.array(delphi1)
		delphi2 = np.array(delphi2)

		#plt.plot(delphi1,'k')
		#plt.plot(delphi2,'r')
		#plt.show()

		sig1_a = np.exp(-delphi1/(2*sigma_p1))
		sig2_a = np.exp(-delphi2/(2*sigma_p1))
		sig_a = sig1_a + sig2_a 
		#print "delphi1 , sig1 ", delphi1, sig1
	
	show_signal = False
	if show_signal:
		plt.plot(varphi_a, sig1_a, 'k-')
		plt.plot(varphi_a, sig2_a)
		plt.plot(varphi_a, sig_a, 'g')
	
		plt.show()
		sys.exit()



	#plt.plot(varphi_a, sig_l, 'ko')
	#plt.show()

	#sig_total = np.array(100*sig_l)

	#plt.plot(sig_total)
	#plt.show()

	return sig_a, sig1_a, sig2_a, db_m


#sep_frac_a = np.linspace(0,1,100)
sep_frac_a = np.linspace(0,1,100)
#sep_frac_a = np.array([0.25])
fft_p_0 = []
fft_p_1 = []
fft_p_2 = []

nc = 10


#f_p_1 = int(0.01*nc*nt)
f_p_1 = nc
f_p_2 = 2*f_p_1

db_l = []

plot_sep = False
for sepf1 in sep_frac_a:

	sig1_phase = sepf1*math.pi
	sig2_phase = 2*math.pi-sig1_phase

	sig_p, sig1_p, sig2_p, dbm = generate_signal(sepf1)

	db_l.append(dbm)

	freqr = np.linspace(0, 1, nc*nt)
	sig_total = np.tile(sig_p, nc) 
	sig1_total = np.tile(sig1_p, nc) 
	sig2_total = np.tile(sig2_p, nc)

	fft = np.fft.fft(sig_total)
	
	fftabs = np.abs(fft)
	
	#sig1_total = np.append(sig1_total, [0,0,0])
	#sig1_total = np.delete(sig1_total,[0,1,2])
	#print "len sig1_total ",len(sig1_total)

	fftsig1 = np.fft.fft(sig1_total)
	fftphase_sig1 = np.angle(fftsig1)
	fftabs_sig1 = np.abs(fftsig1)
	fftsig2 = np.fft.fft(sig2_total)
	fftphase_sig2 = np.angle(fftsig2)
	fftabs_sig2 = np.abs(fftsig2)

	sig_ifft_f1 = inverse_fft(fft, [f_p_1])
	sig_ifft_f2 = inverse_fft(fft, [f_p_2])

	sig1_ifft_f1 = inverse_fft(fftsig1, [f_p_1])
	sig2_ifft_f1 = inverse_fft(fftsig2, [f_p_1])
	sig1_ifft_f2 = inverse_fft(fftsig1, [f_p_2])
	sig2_ifft_f2 = inverse_fft(fftsig2, [f_p_2])
	ind_sig1_ifftmax = np.argmax(sig1_ifft_f1[:nt]) 
	ind_sig2_ifftmax = np.argmax(sig2_ifft_f1[:nt]) 

	fft_p_0.append(fft[0])
	fft_p_1.append(fft[f_p_1])
	fft_p_2.append(fft[f_p_2])

	#fft_correlate = fftconvolve(sig1_total, sig2_total)
	#fft_correlate_abs = np.abs(fft_correlate)
	
	#print "sig1, sig2 phase ",sig1_phase, sig2_phase
	#print "predicted FFT phase ",-sig1_phase,  2*math.pi-sig2_phase
	#print "sig1_total shape ",sig1_total.shape

	#print "FFT sig1, sig2 at first harmonic ",np.real(fftsig1[f_p_1]), np.real(fftsig2[f_p_1])
	#print "FFT amplitude sig1, sig2 at first harmonic ",fftabs_sig1[f_p_1], fftabs_sig2[f_p_1]
	#print "FFT phase sig1, sig2 at first harmonic ",fftphase_sig1[f_p_1], fftphase_sig2[f_p_1]

	#print np.cos(fftphase_sig1[f_p_1]), 10*np.cos(fftphase_sig2[f_p_1])
	
	
	maxinv=  np.max(sig1_ifft_f1[:nt])
	if plot_sep:
		plt.subplot(211)
		plt.plot(sig_p, 'k')
		#plt.plot(np.tile(sig1_p,2),'b')
		#plt.plot(np.tile(sig2_p,2),'m')
		#plt.plot(sig1_p,'b')
		#plt.plot(sig2_p,'m')


		#plt.axvline(x=ind_sig1_ifftmax)
		#plt.axvline(x=ind_sig2_ifftmax)
		#plt.plot(sig2_total,'m')
		plt.subplot(212)
		
		#plt.plot(fftabs)
		#plt.plot(fft,'k')

		plt.plot(sig_ifft_f1[:nt],'b')
		plt.plot(sig1_ifft_f1[:nt],'b--')
		plt.plot(sig2_ifft_f1[:nt],'b:')

		plt.plot(sig_ifft_f2[:nt],'r-')

		#plt.plot(sig2_ifft[:nt],'m')
		#plt.plot(maxinv*np.cos(2*math.pi*0.01*np.arange(nt) + 0.5*math.pi),'k')	

		plt.axvline(x=f_p_1)
		
		#plt.plot(fftsig1,'b')
		#plt.plot(fftsig2,'m')
		#plt.xlim(0,f_p_2)
		#plt.plot(fftsig1 + fftsig2,'ro')
		#plt.plot(fft_correlate,'r')
		plt.show()
		sys.exit()
		
	#plt.plot([0.01],[fft_p_1],'o')
	#plt.plot([0.02],[fft_p_2],'o')
if plot_sep:
	plt.show()

#plt.plot(sep_frac_a, db_l)
#plt.show()

sig1_phase_a = sep_frac_a*math.pi

#plt.plot(sep_frac_a, fft_p_0, 'ko-',label='DC')
#plt.plot(sig1_phase_a, fft_p_1/np.max(fft_p_1), 'ro-',label='f')
plt.plot(sep_frac_a, fft_p_1, 'ro-',label='f')
#plt.plot(sep_frac_a, [np.abs(f) for f in fft_p_1], 'ro--',label='f')
#plt.plot(sep_frac_a, 0.5*max(fft_p_1)*(np.cos(sig1_phase_a)+np.cos(-sig1_phase_a))) 
plt.plot(sep_frac_a, max(fft_p_1)*(np.cos(sig1_phase_a)), 'm')
plt.plot(sep_frac_a, max(fft_p_1)*(np.cos(2*sig1_phase_a)), 'c')

plt.plot(sep_frac_a, fft_p_2, 'bo-',label='2f')
#plt.plot(sep_frac_a, [np.abs(f) for f in fft_p_2], 'bo--',label='2f')

plt.legend()
plt.savefig('fft_model')
plt.show()

plt.plot(sep_frac_a, [np.abs(f) for f in fft_p_1], 'ro--',label='f')
plt.plot(sep_frac_a, [np.abs(f) for f in fft_p_2], 'bo--',label='2f')
plt.legend()
plt.show()


sys.exit()



