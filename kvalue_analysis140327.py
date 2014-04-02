#!/usr/bin/python
import matplotlib
import glob
from pylab import *
import numpy as np
from analysedata import *
#from dispersion_analysis140327 import *
matplotlib.rcParams['legend.fontsize']=8


def calcKvalue(t_rf, f_rf, E_rf, tvals, rvals):
    '''use rvalues and tvals to get k value from rf file'''
    
    fvals=map(findFfromt(t_rf, f_rf), tvals)
    E_vals=map(findEfromt(t_rf, E_rf*1E6)/1E6, tvals) #approximate value of gamma at 60MeV
    #print E_vals
    gamma_vals=map(lambda x: (x+938.0)/938.0, E_vals)
    #print gamma_vals
    kval_array=[]
    for i in range(len(rvals)-1):
        drr=(rvals[i+1]-rvals[i])/rvals[i]
        dff=(fvals[i+1]-fvals[i])/fvals[i]
        kval_array.append(gamma_vals[i]*gamma_vals[i]*dff/drr-(1.0-gamma_vals[i]*gamma_vals[i]))
        #print dff/drr
    return kval_array


colors = ['b', 'g', 'r']

#Filenames for plots
pfile1 = "r_vs_t.pdf"  
pfile2="k_vs_t.pdf"

#Data directory
#dir="/Users/Suzie/Physics/KURRIFFAG/MARCH14/Davidcode/COD_data_140325/"
#dir="/Users/Suzie/Physics/KURRIFFAG/DATA/2014/2014-03-27/"
dir="/Users/Suzie/Physics/KURRIFFAG/MARCH14/Davidcode/Qin_Bin_31_3_31/"

rfdir="/Users/Suzie/Physics/KURRIFFAG/DATA/2014/2014-03-25/rf_pattern/"

#Filenames for sets of data
#set1=glob.glob(dir+"QinBin_F?.txt")
#print set1
set1=["QinBin_F1.txt", "QinBin_F5.txt", "QinBin_F7.txt"]
#set1 =["D0_0amps_F1.txt" ,"D0_0amps_F5.txt" ,"D0_0amps_F7.txt"] #468A
#set1=["D0_550Aamps_F1.txt" ,"D0_550Aamps_F5.txt" ,"D0_550Aamps_F7.txt"] #550A 27/3/2014
#set1 =["D0_140amps_F1.txt", "D0_140amps_F5.txt", "D0_140amps_F7.txt"] #468A
#set1=["D0_400Aamps_F1.txt", "D0_400Aamps_F5.txt", "D0_400Aamps_F7.txt"] #400A 27/3/14
#set1=["D0_700amps_F1.txt", "D0_700amps_F5.txt", "D0_700amps_F7.txt"] #27/3/2014

use_radius_range=True
use_time_range=False

#For a given time range, get the radial position at each time from each probe
data1=readset(dir, set1) #gets the data
time_range=np.arange(500,1500,1)
#rad_range=np.arange(300,420,1)


radial_values=[]
time_values=[]

probenames=["F1","F5","F7"]
proberadius=4180.0 #mm
kvalue=7.6

fig=figure(num=1, figsize=(6.5, 10.5), dpi=150, facecolor='w', edgecolor='k')   

#Lookup table of time, frequency and estimated energy (based on RF settings file)
rffile="SimulatedVariableK.dat"
rfdat = np.loadtxt(rfdir+rffile, skiprows=2, usecols=(0,1,2), unpack=True)
rfdat_t=np.array(rfdat[0])*1E6 #time in musec
rfdat_E=np.array(rfdat[1])*1E-6 #Energy in MeV
rfdat_F=np.array(rfdat[2]) #frequency in hertz

#Can also map the momentum values from the Energy values.
mass=938.272
rfdat_p=map(lambda x: ke_to_mom(mass, x), rfdat_E)

#CALCULATE THE K-VALUE BASED ON FREQUENCY FROM EACH TIME POINT
#THIS ALSO USES THE LOOKUP TABLE TO FIND E AND USE ESTIMATED GAMMA FACTOR FOR A MORE EXACT RESULT
kresult=[]
toffset= -93.3 # musec timing offset from 31/3/14 data

fig=figure(num=3, figsize=(6, 15), dpi=150, facecolor='w', edgecolor='k')

smoothed_data=False
for probe in range(len(probenames)):
    if smoothed_data==True:
        rad_range_1=np.arange(min(data1[probe][0]),max(data1[probe][0]), 1)

    else: rad_range_1=sort(np.array(data1[probe][0]))
    
    fittedtfunc=findtfromR(data1[probe][1], data1[probe][0])
    tfit_values=map(fittedtfunc+toffset, rad_range_1)
    rad_range_m=(rad_range_1+proberadius)/1E3
    
    subplot(2,1,1)
    plot(tfit_values, rad_range_m, '-', color=colors[probe])
    plot(data1[probe][1]+toffset, (np.array(data1[probe][0])+proberadius)/1E3, 'x',color=colors[probe])

    ylabel(r'radius [m]')
    subplot(2,1,2)
    kresult.append(calcKvalue(rfdat_t, rfdat_F, rfdat_E, tfit_values, rad_range_m))

    plot(tfit_values[:-1],kresult[probe], '.-', color=colors[probe],label=probenames[probe])
    legend()
    xlabel(r'time [$\mu$sec]', size=12)
    ylabel(r'measured k',size=12)
    ylim([6,9])

savefig(pfile2)

sys.exit()

