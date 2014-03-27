#!/usr/bin/python
import matplotlib
from pylab import *
import numpy as np
from analysedata import *
matplotlib.rcParams['legend.fontsize']=8



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
    rfit = np.polyfit(rvals,tvals,1)
    #t=ar+b so r=t/a-b/a
    tfit=np.array([1/rfit[[0]][0], -rfit[[1]][0]/rfit[[0]][0]])
    
    #make a polynomial that gives the r value for any t value
    rpoly = np.poly1d(tfit)
    return rpoly
    
def findPfromt(pvals, tvals):
    '''do a polynomial fit to find momentum as a function of t from set of datapoints'''
    pfit = np.polyfit(pvals,tvals,4)
    usepoly = np.poly1d(pfit)
    return usepoly
    
def calcsingleDispersion(pvals, rvals):
    #central momentum
    pcentral=mean(pvals)
    dp=max(pvals)-min(pvals)
    dpp=dp/pcentral 
    #print "delta p: ", dp
    #print "delta p/p: ", dpp
    dr=max(rvals)-min(rvals)
    #print "dr: ", dr
    return dr/dpp
    
def calcfitDispersion(pvals, rvals):
    pcentral=mean(pvals)
    pfit=np.polyfit((np.array(pvals)-pcentral)/pcentral, rvals, 1) #r=ap+b, dispersion = dr/(dp/p)
    return pfit[[0]][0]
    

def ke_to_mom(pmass, kenergy):
    "Convert kinetic energy to momentum"
    te= pmass+kenergy
    mom=math.sqrt(pow(te,2) - pow(pmass,2))
    return mom

#Data directory
dir="/Users/Suzie/Physics/KURRIFFAG/MARCH14/Davidcode/COD_data_140325/"
rfdir="/Users/Suzie/Physics/KURRIFFAG/DATA/2014/2014-03-25/rf_pattern/"

#Filenames for plots
pfile = "140327_dispersion_ring.pdf"

#fig=figure(num=1, figsize=(6.5, 10.5), dpi=150, facecolor='w', edgecolor='k')   

#Filenames for sets of data
set1 =["D0_0amps_F1.txt" ,"D0_0amps_F5.txt" ,"D0_0amps_F7.txt"] #468A
set2 =["D0_140amps_F1.txt", "D0_140amps_F5.txt", "D0_140amps_F7.txt"] #468A

#For a given time range, get the radial position at each time from each probe
data1=readset(dir, set1) #gets the data
time_values=np.arange(500,1500,10)
radial_values=[]

fig=figure(num=1, figsize=(6.5, 10.5), dpi=150, facecolor='w', edgecolor='k')   
pfilen="r_fitted.pdf"

for i in range(len(data1)):
    fittedrfunc=findRfromt(data1[i][0], data1[i][1])
    plot(data1[i][1], data1[i][0], '.', data1[i][1], fittedrfunc(data1[i][1]), '-', label=str(i))
    R_values=map(fittedrfunc, time_values)
    radial_values.append(R_values)

#Get R as a function of t for the average of the three probes (ie. CO movement)
avg_radial_values=[]
max_radial_offset=[]
min_radial_offset=[]
for r in np.transpose(np.array(radial_values)):
    avg_radial_values.append(np.mean(r))
    max_radial_offset.append(np.max(r)-np.mean(r))
    min_radial_offset.append(np.mean(r)-np.min(r))

plot(time_values, avg_radial_values,'r-', label='avg.')
xlabel(r'Time [$\mu$sec]', size=12)
ylabel(r'Radius [cm]',size=12)
legend()
savefig(pfile)

    
#Convert the time to an energy then a momentum (based on RF settings file)
rffile="SimulatedVariableK.dat"
rfdat = np.loadtxt(rfdir+rffile, skiprows=2, usecols=(0,1), unpack=True)
rfdat_t=np.array(rfdat[0])*1E6
rfdat_E=np.array(rfdat[1])*1E-6

mass=938.272
rfdat_p=map(lambda x: ke_to_mom(mass, x), rfdat_E)

#make a list of momentum values corresponding to a set of time points
fittedpfunc=findPfromt(rfdat_t,rfdat_p)
P_values=map(fittedpfunc, time_values)

pfile3="rfpattern.pdf"
fig=figure(num=3, figsize=(6.5, 10.5), dpi=150, facecolor='w', edgecolor='k')
subplot(2,1,1)   
plot(rfdat_t[:4000], rfdat_p[:4000], '-')
plot(time_values, P_values, 'ro')
xlabel(r'time [$\mu$sec]', size=12)
ylabel(r'Momentum [MeV/c]',size=12)
subplot(2,1,2)
plot(rfdat_t[:4000], rfdat_E[:4000], '-')
xlabel(r'time [$\mu$sec]', size=12)
ylabel(r'Energy[MeV]',size=12)
savefig(pfile3)

#CALCULATE THE DISPERSION (TWO METHODS)
#1. just take max/min value of P and R in range & use deltaR/(deltaP/P)
avg_disp = calcsingleDispersion(P_values, np.array(avg_radial_values)/1E3)
print "Average dispersion in metres: ", avg_disp

#2. make a linear fit in the region and use the gradient
grad_disp = calcfitDispersion(P_values, np.array(avg_radial_values)/1E3)
print "fitted dispersion in metres: ", grad_disp

pfile4="p_vs_r.pdf"
fig=figure(num=4, figsize=(6.5, 10.5), dpi=150, facecolor='w', edgecolor='k')   
plot( P_values,np.array(avg_radial_values)/1E3,'.')
xlabel(r'Momentum [MeV/c]', size=12)
ylabel(r'Avg. radius [m]',size=12)
savefig(pfile4)


