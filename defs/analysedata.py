#!/usr/bin/python
import matplotlib
#matplotlib.rcParams['font.family'] = 'serif'
#matplotlib.rcParams['legend.fontsize']=12
from pylab import *
import numpy

#DAVID'S READSCOPE FUNCTION
def read_scope(dirname, filename):
	"""Read scope data. The script deals with two specific scope formats - 
	that used in November 2013, and in March 2014. The two are distinguished by looking
	for the string DP04104 that appears in the later case"""
	
	filepath = dirname + filename
	print 'reading ', filepath

	f_data = open(filepath, 'r')
	n = 0
	tdat = []
	ydat = []
	
	for line in f_data:
	
		timeandamp = line.split(',')

		#identify scope
		if n == 0:
			if timeandamp[1][:7] == 'DPO4104':
				nskip = 17
				xi = 0
				yi = 1
			else:
				nskip = 6
				xi = 3
				yi = 4

		if n > nskip:
			if line !='\r\n':
				tdat.append(float(timeandamp[xi]))
				ydat.append(float(timeandamp[yi]))
			else:
				break
			
		n += 1

	f_data.close()

	
	return tdat, ydat

#DAVIDS HIGHPASS FILTER FUNCTION added 26/3/14
def highpass(signal,dt,RC):
	""" Return RC high-pass output samples, given input samples, time constant RC and time interval dt"""
	alpha = RC/(RC + dt)
	y =	 [signal[0]]
	for i in range(1,len(signal)):
		y.append(alpha*y[i-1] + alpha*(signal[i] - signal[i-1]))
	return y

#DAVIDS FUNCTION TO GET ALL THE DATA FILES IN A DIRECTORY

 
#SUZIE'S ANALYSIS FUNCTIONS
def readcsv_old(fname, nsamples=100000):
	"""csv files from the oscilloscope that was used in Nov 13"""
	
	time=[]
	volts=[]

	dat = csv2rec(fname)

	for it in range(nsamples):
		time.append(float(dat[it][3]))
		volts.append(float(dat[it][4]))

	print len(time)	 
	return time, volts

def readcsv(fname, nsamples=100000):
	dat=csv2rec(fname, skiprows=17)
		
	for it in range(nsamples):
		time.append(float(dat[it][0]))
		volts.append(float(dat[it][1]))

	print len(time)	 
	return time, volts


def getpeaks(rawdata,revfreq, startindex=0,endindex=10000):
	"""This breaks up data into chunks defined by int(1/revfreq), between given start & end points
	finds the maximum value within that window & returns arrays of peak values & positions."""
	si=startindex
	step=int(1/revfreq)
	print step
	pmax=[]
	ppos=[]

	index=si+step
	while index<endindex+1:
		peak=abs(rawdata[si:index]).max()
		pmax.append(peak)
		print peak
		ppos.append(si+ numpy.where(abs(rawdata[si:index])==peak)[0][0])
		si=si+step
		index=index+step

	return pmax, ppos


def perform_fft(data, plotswitch=0, plotname="fft_spectrum.png"):
	"""Performs a basic FFT and returns the dominant frequency in the data and the freq. vs power spectrum data"""
	t=numpy.arange(len(data))
	Y = numpy.fft.fft(data)
	n=len(Y)
	power = abs(Y[1:(n/2)])**2
	nyquist=1./2

	freq=numpy.array(range(n/2))/(n/2.0)*nyquist

	if plotswitch==1:
		plot(freq[:len(power)], power, '-')
		yscale('log')
		#xlim(0, 0.04)
		#ylim(1, 100000)
		savefig(plotname)
		clf()

	maxindex=numpy.where(power==power.max())[0][0]
	print "Max. frequency in spectrum: ", freq[maxindex+1]
	return freq[maxindex+1], freq[:len(power)], power

### This is the Gaussian data smoothing function I got from:
###	 http://www.swharden.com/blog/2008-11-17-linear-data-smoothing-in-python/
def smoothListGaussian(list,degree=5):	

	window=degree*2-1  

	weight=numpy.array([1.0]*window)  

	weightGauss=[]	

	for i in range(window):	 

		i=i-degree+1  

		frac=i/float(window)  

		gauss=1/(numpy.exp((4*(frac))**2))	

		weightGauss.append(gauss)  

	weight=numpy.array(weightGauss)*weight	

	smoothed=[0.0]*(len(list)-window)  
	
	for i in range(len(smoothed)):	

		smoothed[i]=sum(numpy.array(list[i:i+window])*weight)/sum(weight)  

	return smoothed

def smoothListTriangle(list,strippedXs=False,degree=5):	 

	weight=[]  

	window=degree*2-1  

	smoothed=[0.0]*(len(list)-window)  

	for x in range(1,2*degree):weight.append(degree-abs(degree-x))	

	w=numpy.array(weight)  

	for i in range(len(smoothed)):	

		smoothed[i]=sum(numpy.array(list[i:i+window])*w)/float(sum(w))	

	return smoothed 

def peakfinder(data, thresh=2.5):
	"""Method: take derivative of data, 
	Find sign of gradient (increasing/decreasing)
	Find where the gradient of the sign changes -> take downward crossing only (gives maxima) 
	Find the array indices where that is 'True' and add one"""

	c = (diff(sign(diff(data))) < 0).nonzero()[0] + 1 # local max
	#threshold = thresh*mean(data)
	threshold=thresh
	print 'total peaks found: ', len(c)
	print 'threshold: ', threshold
	# Next remove maxima from the list if they fall below the threshold of 2.5*mean & print how many 'peaks'
	c_keep = numpy.delete(c, ravel(where(data[c]<threshold)))
	print 'peaks found above threshold: ',len(c_keep)
	return c_keep

def extremumfinder(data, thresh=2.5):
	"""Method: take derivative of data, 
	Find sign of gradient (increasing/decreasing)
	Find where the gradient of the sign changes -> take downward crossing only (gives maxima) 
	Find the array indices where that is 'True' and add one"""

	c = (diff(sign(diff(data)))).nonzero()[0] + 1 # local max
	print 'total extrema found: ', len(c)
	threshold = thresh*mean(data)
	print 'threshold: ', threshold
	# Next remove maxima from the list if they fall below the threshold of 2.5*mean & print how many 'peaks'
	c_keep = numpy.delete(c, ravel(where(abs(data[c])<threshold)))
	print 'peaks found: ',len(c_keep)
	return c_keep

def makeplot(ffile, pfile):	   
	"""JUST MAKES A PLOT EITHER WITH TWO COLUMNS OF DATA, OR A THIRD REPRESENTING YERR VALUES
	0 = timebase
	1 = volt"""
	
	xax=[]
	yax=[]
	errs=[]
	for line in file(ffile):
		line = line.split()
		xax.append(float(line[0]))
		yax.append(float(line[1]))
		if len(line)>2:
			errs.append(float(line[2]))

	fig=figure(num=1, figsize=(6.5, 10.5), dpi=150, facecolor='w', edgecolor='k')		
	if len(errs)==len(yax):
		errorbar(numpy.array(xax), numpy.array(yax), yerr=numpy.array(errs), fmt='bo-')
	else:
		plot(xax, yax, 'k.-')
	#ylim(0,0.5)
	#xlabel(r'x [mm]', size=12)
	#ylabel(r'z [mm]',size=12)
	savefig(pfile)





