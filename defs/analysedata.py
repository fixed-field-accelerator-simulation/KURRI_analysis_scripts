#!/usr/bin/python
import matplotlib
#matplotlib.rcParams['font.family'] = 'serif'
#matplotlib.rcParams['legend.fontsize']=12
from pylab import *
import numpy as np
import os


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
    y =  [signal[0]]
    for i in range(1,len(signal)):
        y.append(alpha*y[i-1] + alpha*(signal[i] - signal[i-1]))
    return y


def get_files(directory, file_type):
    """ return files of type file_type in directory"""
    
    if not os.path.exists(directory):
        print "directory ",directory, "not found"
        sys.exit()
    
    listf = os.listdir(directory)
    
    files = []
    for file1 in listf:
        if file1[0] != '.':
            if file_type in file1:
                files.append(file1)
    
    return files


def select_filenames(directory, indices_select, channel_id):
    """Select channel files according to indices"""
    
    ch_files_list = []
    for cid in channel_id:
        f = get_files(directory, cid)
        if f !=[]:
            if indices_select:
                fsel = [f[i] for i in indices_select] 
            else:
                fsel = list(f)
            ch_files_list.append(fsel)
        else:
            ch_files_list.append([])

    return ch_files_list


def signal_loss_time(tdat,data, method, cut_off_freq = 0.5e6, threshold_volts = 0.05, sigtonoise = 10):
    """Find time of last signal.  
    The low frequency component is first removed using a high pass filter (below cut_off_freq).
    
    method can be "fixed" in which case the latest time with voltage at threshold_volts is found
    or "noise" in which case the latest time with voltage sigtonoise times the rms of the noise calculated
    in the last 100 points (where the script assumes no signal)"""
    
    #parameters for filter
    dt = tdat[1] - tdat[0]
    RC = 1/(2*np.pi*cut_off_freq)
    
    #apply high pass filter
    filter_chf = highpass(data,dt,RC)

    if method == 'fixed':
        threshold = threshold_volts
    elif method == 'noise':
        noise_level = np.std(filter_chf[-100:])
        threshold = sigtonoise*noise_level
    
    for y in reversed(filter_chf):
        if abs(y) > threshold:
            loss_time = tdat[filter_chf.index(y)]
            break
    
    return loss_time
    
        
#SUZIE'S ANALYSIS FUNCTIONS


def readset(dirname, setname):
    '''Read in a set of data where a 'set' is from a number of different probes'''
    dat=[]
    for i in range(len(setname)):
        print "Reading data from: ", setname[i]
        ffile = dirname+setname[i]
        dat.append(np.loadtxt(ffile, skiprows=1, usecols=(0,1), unpack=True))
    return dat
        

def plotset(dirname, setname):
    '''Plot a set of data where a 'set' is from a number of different probes'''

    dataset = readset(dirname, setname)
    for i in range(len(dat)):
        dat=dataset[i]
        plot(dat[0], dat[1], '.-',label=setname[i])


def findRfromt(rvals, tvals):
    '''do a polynomial fit to find R as a function of t from set of datapoints'''
    rfit = np.polyfit(tvals,rvals,4)
    #r=t/a-b/a
    #make a polynomial that gives the r value for any t value
    rpoly = np.poly1d(rfit)
    return rpoly
    
def findtfromR(tvals, rvals):
    '''do a polynomial fit to find t as a function of R from set of datapoints'''
    tfit = np.polyfit(rvals,tvals,4)
    #r=t/a-b/a
    #make a polynomial that gives the t value for any r value
    tpoly = np.poly1d(tfit)
    return tpoly
    
def findFfromt(tvals, fvals):
    '''do a polynomial fit to find f as a function of t from set of datapoints'''
    ffit = np.polyfit(tvals,fvals,5)
    #make a polynomial that gives the t value for any r value
    fpoly = np.poly1d(ffit)
    return fpoly
    
def findPfromt(pvals, tvals):
    '''do a polynomial fit to find momentum as a function of t from set of datapoints'''
    pfit = np.polyfit(pvals,tvals,3)
    usepoly = np.poly1d(pfit)
    return usepoly
    
def findEfromt(tvals, Evals):
    '''do a polynomial fit to find momentum as a function of t from set of datapoints'''
    Efit = np.polyfit(tvals,Evals,5)
    usepoly = np.poly1d(Efit)
    return usepoly
    
def findPfromR(rvals, pvals):
    pfit = np.polyfit(pvals,rvals,3)
    usepoly = np.poly1d(pfit)
    return usepoly
    
def ke_to_mom(pmass, kenergy):
    '''Convert kinetic energy to momentum'''
    te= pmass+kenergy
    mom=math.sqrt(pow(te,2) - pow(pmass,2))
    return mom
    
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
        ppos.append(si+ np.where(abs(rawdata[si:index])==peak)[0][0])
        si=si+step
        index=index+step

    return pmax, ppos


def perform_fft(data, plotswitch=0, plotname="fft_spectrum.png"):
    """Performs a basic FFT and returns the dominant frequency in the data and the freq. vs power spectrum data"""
    t=np.arange(len(data))
    Y = np.fft.fft(data)
    n=len(Y)
    power = abs(Y[1:(n/2)])**2
    nyquist=1./2

    freq=np.array(range(n/2))/(n/2.0)*nyquist

    if plotswitch==1:
        plot(freq[:len(power)], power, '-')
        yscale('log')
        #xlim(0, 0.04)
        #ylim(1, 100000)
        savefig(plotname)
        clf()

    maxindex=np.where(power==power.max())[0][0]
    print "Max. frequency in spectrum: ", freq[maxindex+1]
    return freq[maxindex+1], freq[:len(power)], power

### This is the Gaussian data smoothing function I got from:
###  http://www.swharden.com/blog/2008-11-17-linear-data-smoothing-in-python/
def smoothListGaussian(list,degree=5):  

    window=degree*2-1  

    weight=np.array([1.0]*window)  

    weightGauss=[]  

    for i in range(window):  

        i=i-degree+1  

        frac=i/float(window)  

        gauss=1/(np.exp((4*(frac))**2)) 

        weightGauss.append(gauss)  

    weight=np.array(weightGauss)*weight 

    smoothed=[0.0]*(len(list)-window)  
    
    for i in range(len(smoothed)):  

        smoothed[i]=sum(np.array(list[i:i+window])*weight)/sum(weight)  

    return smoothed

def smoothListTriangle(list,strippedXs=False,degree=5):  

    weight=[]  

    window=degree*2-1  

    smoothed=[0.0]*(len(list)-window)  

    for x in range(1,2*degree):weight.append(degree-abs(degree-x))  

    w=np.array(weight)  

    for i in range(len(smoothed)):  

        smoothed[i]=sum(np.array(list[i:i+window])*w)/float(sum(w)) 

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
    c_keep = np.delete(c, ravel(where(data[c]<threshold)))
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
    c_keep = np.delete(c, ravel(where(abs(data[c])<threshold)))
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
        errorbar(np.array(xax), np.array(yax), yerr=np.array(errs), fmt='bo-')
    else:
        plot(xax, yax, 'k.-')
    #ylim(0,0.5)
    #xlabel(r'x [mm]', size=12)
    #ylabel(r'z [mm]',size=12)
    savefig(pfile)





