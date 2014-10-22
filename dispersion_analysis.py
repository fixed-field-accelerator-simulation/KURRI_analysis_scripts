#!/usr/bin/python
import matplotlib
import glob
from pylab import *
import numpy as np
from analysedata import *
from scipy import interpolate

#Written by S. Sheehy on 27/3/2014
#Updates:
#21/10/2014: Updated to use individual ranges for each probe calculation (S. Sheehy)
#


matplotlib.rcParams['legend.fontsize']=8

    
def calcsingleDispersion(pvals, rvals):
    '''Use to calculate a single value of the dispersion for a set of r vs p data'''
    #central momentum
    pcentral=mean(pvals)
    dp=max(pvals)-min(pvals)
    dpp=dp/pcentral 
    #print "delta p: ", dp
    print "delta p/p: ", dpp
    dr=max(rvals)-min(rvals)
    print "dr: ", dr
    return p*(dr/dp)
    
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
        dpp=dp/pvals[i]
        disp.append(pvals[i]*dr/dp) #changed from dr/dpp
        dpp_array.append(dpp)
    #dispersion = dr/(dp/p)
    return disp, dpp_array

def calckfromdisp(Dval, rval):
    #gamma_val=(Eval+938.0)/938.0
    #print gamma_val
    #kval=(1-Dval)*pow(gamma_val,2)-1
    kval=(rval/Dval)-1
    #print kval
    return kval

def calcKvalue_poly(t_rf, f_rf, E_rf, tvals, rvals):
    '''use rvalues and tvals to get k value from rf file'''
    #fvals=interpolate.splev(tvals, findyfromx_spline(f_rf, t_rf), der=0)
    #E_vals=interpolate.splev(tvals, findyfromx_spline(E_rf, t_rf),der=0)
    fvals=map(findFfromt(t_rf, f_rf), tvals)
    E_vals=map(findEfromt(t_rf, E_rf*1E6)/1E6, tvals) #approximate value of gamma at 60MeV
    gamma_vals=map(lambda x: (x+938.0)/938.0, E_vals)
    #print gamma_vals
    kval_array=[]
    for i in range(len(rvals)-1):
        drr=(rvals[i+1]-rvals[i])/rvals[i]
        dff=(fvals[i+1]-fvals[i])/fvals[i]
        kval_array.append(gamma_vals[i]*gamma_vals[i]*dff/drr-(1.0-gamma_vals[i]*gamma_vals[i]))
        #print dff/drr
    return kval_array
    
def calcKvalue(t_rf, f_rf, E_rf, tvals, rvals):
    '''use rvalues and tvals to get k value from rf file'''
    #fvals=interpolate.splev(tvals, findyfromx_spline(f_rf, t_rf), der=0)
    E_vals=interpolate.splev(tvals, findyfromx_spline(E_rf, t_rf),der=0)
    fvals=map(findFfromt(t_rf, f_rf), tvals)
    
    #E_vals=map(findEfromt(t_rf, E_rf*1E6)/1E6, tvals) #approximate value of gamma at 60MeV
    gamma_vals=map(lambda x: (x+938.0)/938.0, E_vals)
    #print gamma_vals
    kval_array=[]
    for i in range(len(rvals)-1):
        drr=(rvals[i+1]-rvals[i])/rvals[i]
        dff=(fvals[i+1]-fvals[i])/fvals[i]
        kval_array.append(gamma_vals[i]*gamma_vals[i]*dff/drr-(1.0-gamma_vals[i]*gamma_vals[i]))
        #print dff/drr
    return kval_array
    
if __name__ == "__main__":
    '''This analysis code needs updating/checking!!!'''

    colors = ['b', 'g', 'r']

    #Directory for plots
    resdir="Results/"
    
    #Filenames for plots
    pfile1 ="r_vs_t.pdf"
    pfile2="rfpattern.pdf"
    #pfile2="k_vs_t.pdf"
    pfile3="p_vs_r.pdf"
    pfile4="Dispersion_vs_E.pdf"
    pfile5="Dispersion_vs_r.pdf"
    pfile7="kvalue_vs_E.pdf"
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

    #get the radial position at each time from each probe
    data1=readset(dir, set1) #gets the data
    #time_range=np.arange(500,15000,1)
    #rad_range=np.arange(300,420,1)

    time_range=[]
    rad_range=[]
    radial_values=[]
    time_values=[]

    probenames=["F1","F5","F7"]
    proberadius=4300.0 #mm
    kvalue=7.6

    fig=figure(num=1, figsize=(8.5, 6.5), dpi=150, facecolor='w', edgecolor='k')   

    #Get radius as a function of t for the three probes (ie. CO movement)
    for i in range(len(data1)):
        fittedrfunc=findRfromt(data1[i][0], data1[i][1])
        j=data1[i][0].argsort() 
        time_range=np.arange(min(data1[i][1]), max(data1[i][1]), 10)
        
        print "time range: ", min(data1[i][1]), max(data1[i][1])
        #Rfit_values=interpolate.splev(time_range, findyfromx_spline(data1[i][0][j], data1[i][1][j]), der=0)
        Rfit_values=map(fittedrfunc, time_range)
        radial_values.append(Rfit_values)
        time_values.append(time_range)
        plot(data1[i][1], data1[i][0], 'x', color=colors[i], label=probenames[i])
    #plot R vs t for each of the 3 probes
    for i in range(len(probenames)):
        plot(time_values[i], radial_values[i], color=colors[i], label=probenames[i]+"fit")
    xlabel(r'Time [$\mu$sec]', size=12)
    ylabel(r'Radius [mm]',size=12)
    legend(loc=0)
    savefig(resdir+pfile1)
    close()


    #also calculate the average value of CO movement
    # avg_radial_values=[]
    # max_radial_offset=[]
    # min_radial_offset=[]
    # for r in np.transpose(np.array(radial_values)):
    #   avg_radial_values.append(np.mean(r))
    #   max_radial_offset.append(np.max(r)-np.mean(r))
    #   min_radial_offset.append(np.mean(r)-np.min(r))

    
    #Convert the time to an energy then a momentum (based on RF settings file)
    rffile="SimulatedVariableK.dat"
    rfdat = np.loadtxt(rfdir+rffile, skiprows=2, usecols=(0,1,2), unpack=True)
    rfdat_t=np.array(rfdat[0])*1E6
    rfdat_E=np.array(rfdat[1])*1E-6
    rfdat_F=np.array(rfdat[2]) #frequency in hertz

    mass=938.272
    rfdat_p=map(lambda x: ke_to_mom(mass, x), rfdat_E)
    toffset= -93.3
    
    #make a list of momentum values corresponding to a given set of time points
    fittedpfunc=findPfromt(rfdat_t,rfdat_p) #using RF data
    
    #get a range of momenta over the given time range for each probe
    P_values=[]
    E_values=[]
    for i in range(len(probenames)):
        #P_values.append(interpolate.splev(time_values[i]+toffset, findyfromx_spline(rfdat_p, rfdat_t), der=0))
        P_values.append(map(fittedpfunc, time_values[i]+toffset))
        E_values.append(map(lambda x: mom_to_ke(mass, x), P_values[i]))
    
    #plot time vs momentum
    fig=figure(num=3, figsize=(8, 10), dpi=150, facecolor='w', edgecolor='k')
    subplot(2,1,1)
    plot(rfdat_t, rfdat_p, 'k--')
    for i in range(len(probenames)):
        plot(time_values[i], P_values[i], '-', color=colors[i], label=probenames[i])
    xlabel(r'time [$\mu$sec]', size=12)
    ylabel(r'Momentum [MeV/c]',size=12)
    # 
    #plot time vs energy
    subplot(2,1,2)
    plot(rfdat_t, rfdat_E, 'k--')
    for i in range(len(probenames)):
        plot(time_values[i], E_values[i], '-', color=colors[i], label=probenames[i])
    xlabel(r'time [$\mu$sec]', size=12)
    ylabel(r'Energy[MeV]',size=12)
    savefig(resdir+pfile2)
    close()

    #CALCULATE THE DISPERSION

    #first adjust radial values to be in [m] and with correct offset from r=0
    radii_metres=[]
    

    dispresult=[]
    dpp=[]
    for probe in range(len(probenames)):
        radii_metres.append((np.array(radial_values[probe])+proberadius)/1E3)
        dispresult.append(calcDispersion(P_values[probe], radii_metres[probe])[0])
        dpp.append(calcDispersion(P_values[probe], radii_metres[probe])[1])

    #sys.exit()
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

    fig=figure(num=4, figsize=(8, 5), dpi=150, facecolor='w', edgecolor='k') 
    for i in range(len(probenames)):
        plot(radii_metres[i][1:], dispresult[i], color=colors[i], label=probenames[i])
    #plot(rad_range_m, idealDispersionvsr(rad_range_m, kvalue))
    xlabel(r'Radius [m]', size=12)
    ylabel(r'Dispersion [m]',size=12)
    ylim([0.3,0.8])
    legend()
    savefig(resdir+pfile5)
    close()
    
    fig=figure(num=5, figsize=(8, 5), dpi=150, facecolor='w', edgecolor='k') 
    for i in range(len(probenames)):
        plot(E_values[i][1:], dispresult[i], color=colors[i], label=probenames[i])
    #plot(rad_range_m, idealDispersionvsr(rad_range_m, kvalue))
    xlabel(r'Energy [MeV]', size=12)
    ylabel(r'Dispersion [m]',size=12)
    ylim([0.3,0.8])
    legend()
    savefig(resdir+pfile4)
    close()
    
    #plot radial position with momentum
    fig=figure(num=6, figsize=(8, 5), dpi=150, facecolor='w', edgecolor='k')   
    #plot( P_values,np.array(avg_radial_values),'k:', label="avg")
    for i in range(len(probenames)):
        plot(P_values[i], radii_metres[i], color=colors[i], label=probenames[i])
    xlabel(r'Momentum [MeV/c]', size=12)
    ylabel(r'Radius [m]',size=12)
    legend()
    savefig(resdir+pfile3)
    close()
    
    #plot effective k value calculated from dispersion 
    fig=figure(num=7, figsize=(8, 5), dpi=150, facecolor='w', edgecolor='k')
    k_values=[]
    for i in range(len(probenames)):
        k_values.append(map(lambda x, y: calckfromdisp(x,y), dispresult[i], radii_metres[i][1:]))
        #plot(time_values[i][1:], k_values[i], '--', color=colors[i], label=probenames[i]+'dmethod')
    xlabel(r'time [us]', size=12)
    ylabel(r'k value',size=12)
    ylim([5.5,9])
    #legend()
    #savefig(resdir+pfile7)

    #CALCULATE THE K-VALUE BASED ON FREQUENCY FROM EACH TIME POINT
    #THIS ALSO USES THE LOOKUP TABLE TO FIND E AND USE ESTIMATED GAMMA FACTOR FOR A MORE EXACT RESULT
    kresult=[]

    #fig=figure(num=3, figsize=(6, 15), dpi=150, facecolor='w', edgecolor='k')

    for probe in range(len(probenames)):

        fittedtfunc=findtfromR(data1[probe][1], data1[probe][0]) #from experimental data 
        j=data1[probe][0].argsort() #data sorted by index of sorted radial values 
        #print data1[probe][0][j]
        #print data1[probe][1][j]
        #tfit_values=interpolate.splev(radial_values[probe], findyfromx_spline(data1[probe][1][j], data1[probe][0][j]), der=0) #use cubic spline to fit time values for finer radial points
        #tfit_values=map(lambda x: x+toffset, tfit_values) #shift t values by timing offset
        tfit_values=map(fittedtfunc+toffset, radial_values[probe]) #find t values corresponding to finer points in r
        kresult.append(calcKvalue_poly(rfdat_t, rfdat_F, rfdat_E, tfit_values, radii_metres[probe])) #calculate the k value in this range
        #print kresult[probe]
        plot(time_values[probe][1:],kresult[probe], '-', color=colors[probe],label=probenames[probe])
        
        xlabel(r'time [us]', size=12)
    	ylabel(r'effective k value',size=12)
        legend()
        #xlabel(r'time [$\mu$sec]', size=12)
    savefig(resdir+pfile7)

sys.exit()

