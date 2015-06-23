#!/usr/bin/python
import matplotlib
#matplotlib.rcParams['font.family'] = 'serif'
#matplotlib.rcParams['legend.fontsize']=12
from pylab import *
import numpy
from analysedata import *

fdirect="/Users/Suzie/Physics/KURRIFFAG/DATA/2013_11_13_2/"
nsamples=5000
degree=10 #This controls how smoothed the data is

###TESTING###
#fname="double_240.csv"
#rawdata = readcsv(fdirect+fname, nsamples)

vals = [220, 225, 230, 235, 240, 245, 250, 255, 260, 300, 400]

outfilen="CO_data.dat"
outf=file(outfilen, 'w')

for v in vals:
    print "Probe position: ", v
    fname ="double_"+str(v)+".csv" 
       
    rawdata=readcsv(fdirect+fname, nsamples)

    #Create Gaussian smoothed data (uses window of 'degree' & shifts by 'degree' step windowing)
    smoothdat = smoothListGaussian(rawdata[1], degree)

    ###TESTING###
    #plot(abs(numpy.array(rawdata[1])[1000:1700]), 'k-')
    #plot(abs(numpy.array(smoothdat)[1000-degree:1700-degree]), 'r-')
    #show()
    #sys.exit()

    x=linspace(0, len(smoothdat), len(smoothdat))
    data=abs(numpy.array(smoothdat))

    #NOW FIND PEAKS IN SMOOTHED DATA
    c_keep = peakfinder(data, 2.0)

    # If there are more than 10 peaks, normalise them wrt the 0th (h- injection) peak and write to file
    if len(c_keep) > 11:
        ttratio = (data[c_keep[10]]/data[c_keep[0]])
        #tterror = (data[c_keep[10]]/data[c_keep[0]]).var()
        outf.write(str(v)+' '+str(ttratio)+'\n')
        plot(x, data/data[c_keep[0]])
        plot(x[c_keep], data[c_keep]/data[c_keep[0]], "o", label="max")
        legend()
        show()
    elif len(c_keep) < 11:
         outf.write(str(v)+' 0.0  0.0 \n')

    ###TESTING###
    

    #else: print "10 turns not found"
    
  
outf.close()


