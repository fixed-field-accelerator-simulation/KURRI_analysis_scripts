#!/usr/bin/python
import matplotlib
import glob
from pylab import *
import numpy as np
from analysedata import *

matplotlib.rcParams['legend.fontsize']=8

    
def calcsingleDispersion(pvals, rvals):
    '''Use to calculate a single value of the dispersion for a set of r vs p data'''
    #central momentum
    pcentral=mean(pvals)
    dp=max(pvals)-min(pvals)
    dpp=dp/pcentral 
    #print "delta p: ", dp
    #print "delta p/p: ", dpp
    dr=max(rvals)-min(rvals)
    #print "dr: ", dr
    return p*(dr/dp)

def idealDispersionvsp(pvals, kval):
    disp=[] 
    '''NOT IMPLEMENTED!!! returns an array which contains the ideal dispersion values vs momentum with a given k value'''
    return disp
    
def idealDispersionvsr(rvals, kval):
    '''returns an array which contains the ideal dispersion values vs radius with a given fixed k value'''
    disp=np.array(rvals)/(kval+1.0)
    return disp

def calcDispersion(pvals, rvals):
    '''Use to calculate the dispersion over a range of radii'''
    disp=[]
    dpp_array=[]
    
    for i in range(len(rvals)-1):
        dr=rvals[i+1]-rvals[i]
        dp=pvals[i+1]-pvals[i]
        #pcentral=(pvals[i+1]+pvals[i])/2.0
        dpp=dp/pvals[i]
        disp.append(pvals[i]*dr/dp) #changed from dr/dpp
        dpp_array.append(dpp)
    #dispersion = dr/(dp/p)
    return disp, dpp_array


if __name__ == "__main__":
    '''This analysis code needs updating/checking!!!'''

    colors = ['b', 'g', 'r']

    #Filenames for plots
    pfile1 = "r_vs_t.pdf"
    #pfile2="rfpattern.pdf"
    pfile2="k_vs_t.pdf"
    pfile3="p_vs_r.pdf"
    pfile4="Dispersion_vs_p.pdf"
    pfile5="Dispersion_vs_r.pdf"

    #fig=figure(num=1, figsize=(6.5, 10.5), dpi=150, facecolor='w', edgecolor='k')   

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
    rad_range=np.arange(300,420,1)


    radial_values=[]
    time_values=[]

    probenames=["F1","F5","F7"]
    proberadius=4180.0 #mm
    kvalue=7.6

    rad_range_m=(rad_range+4180.0)/1E3
    #print rad_range_m

    fig=figure(num=1, figsize=(6.5, 10.5), dpi=150, facecolor='w', edgecolor='k')   

    if use_radius_range==True:
        #Get R as a function of t for the three probes (ie. CO movement)
        for i in range(len(data1)):
            fittedrfunc=findRfromt(data1[i][0], data1[i][1])
            plot(data1[i][1], data1[i][0], '.', label=str(probenames[i]))
            Rfit_values=map(fittedrfunc, time_range)
            radial_values.append(Rfit_values)

    elif use_time_range==True:
        for i in range(len(data1)):
            fittedtfunc=findtfromR(data1[i][0], data1[i][1])
            tfit_values=map(fittedtfunc, rad_range)
            time_values.append(tfit_values)
            plot(data1[i][1], data1[i][0], '.', label=str(probenames[i]))

    #also calculate the average value of CO movement
    # avg_radial_values=[]
    # max_radial_offset=[]
    # min_radial_offset=[]
    # for r in np.transpose(np.array(radial_values)):
    #   avg_radial_values.append(np.mean(r))
    #   max_radial_offset.append(np.max(r)-np.mean(r))
    #   min_radial_offset.append(np.mean(r)-np.min(r))

    #plot R vs t for 3 probes and average
    # plot(time_range, avg_radial_values,'k:', label='avg.')
    for i in range(len(probenames)):
        plot(time_range, radial_values[i], color=colors[i], label=probenames[i]+"fit")
    xlabel(r'Time [$\mu$sec]', size=12)
    ylabel(r'Radius [mm]',size=12)
    legend()
    savefig(pfile1)

    
    #Convert the time to an energy then a momentum (based on RF settings file)
    rffile="SimulatedVariableK.dat"
    rfdat = np.loadtxt(rfdir+rffile, skiprows=2, usecols=(0,1,2), unpack=True)
    rfdat_t=np.array(rfdat[0])*1E6
    rfdat_E=np.array(rfdat[1])*1E-6
    rfdat_F=np.array(rfdat[2]) #frequency in hertz

    mass=938.272
    rfdat_p=map(lambda x: ke_to_mom(mass, x), rfdat_E)

    #make a list of momentum values corresponding to a given set of radial points
    fittedpfunc=findPfromt(rfdat_t,rfdat_p) #using RF data
    P_values=map(fittedpfunc, time_range) #get a range of momenta over the given time range

    #plot time vs momentum
    fig=figure(num=3, figsize=(6, 15), dpi=150, facecolor='w', edgecolor='k')
    # subplot(2,1,1)
    # plot(rfdat_t[:4000], rfdat_p[:4000], '-')
    # plot(time_range, P_values, 'r-')
    # xlabel(r'time [$\mu$sec]', size=12)
    # ylabel(r'Momentum [MeV/c]',size=12)
    # 
    # #plot time vs energy
    # subplot(2,1,2)
    # plot(rfdat_t[:4000], rfdat_E[:4000], '-')
    # xlabel(r'time [$\mu$sec]', size=12)
    # ylabel(r'Energy[MeV]',size=12)
    # savefig(pfile2)

    #subplot(3,1,1)
    #plot(rfdat_t, rfdat_E, 'k:', label="rf programme")
    #xlabel(r'time [$\mu$sec]', size=12)
    #ylabel(r'Energy [MeV]')

    #CALCULATE THE DISPERSION

    #first adjust radial values to be in [m] and with correct offset from r=0
    radii_metres=(np.array(radial_values)+proberadius)/1E3

    dispresult=[]
    dpp=[]
    for probe in range(len(probenames)):
        dispresult.append(calcDispersion(P_values, radii_metres[i])[0])
        dpp.append(calcDispersion(P_values, radii_metres[i])[1])


    #plot the dispersion with momentum
    # fig=figure(num=5, figsize=(6.5, 10.5), dpi=150, facecolor='w', edgecolor='k') 
    # plot( P_values[1:],avg_disp,'k:', label="avg")
    # for i in range(len(probenames)):
    #   plot(P_values[1:], dispresult[i], color=colors[i], label=probenames[i])
    # xlabel(r'Momentum [MeV/c]', size=12)
    # ylabel(r'Dispersion [m]',size=12)
    # legend()
    # savefig(pfile4)

    #plot deltap/p with momentum
    # pfile4a="dpp_vs_p.pdf"
    # fig=figure(num=7, figsize=(6.5, 10.5), dpi=150, facecolor='w', edgecolor='k') 
    # for i in range(len(probenames)):
    #     plot(P_values[1:], dpp[i], color=colors[i], label=probenames[i])
    # xlabel(r'Momentum [MeV/c]', size=12)
    # ylabel(r'dp/p',size=12)
    # legend()
    # savefig(pfile4a)

    #plot dispersion with radial probe position
    fig=figure(num=6, figsize=(6.5, 10.5), dpi=150, facecolor='w', edgecolor='k')   

    #plot( avg_radial_values[1:],avg_disp,'k:', label="avg")
    #radii_m=(avg_radial_values[1:]+4180.0)/1E3


    for i in range(len(probenames)):
        plot(radii_metres[i][1:], dispresult[i], color=colors[i], label=probenames[i])
    plot(rad_range_m, idealDispersionvsr(rad_range_m, kvalue))
    xlabel(r'Radius [m]', size=12)
    ylabel(r'Dispersion [m]',size=12)
    legend()
    savefig(pfile5)

    #plot radial position with momentum
    fig=figure(num=4, figsize=(6.5, 10.5), dpi=150, facecolor='w', edgecolor='k')   
    #plot( P_values,np.array(avg_radial_values),'k:', label="avg")
    for i in range(len(probenames)):
        plot(P_values, radii_metres[i], color=colors[i], label=probenames[i])
    xlabel(r'Momentum [MeV/c]', size=12)
    ylabel(r'Radius [m]',size=12)
    legend()
    savefig(pfile3)


