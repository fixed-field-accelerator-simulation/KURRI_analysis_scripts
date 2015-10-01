#!/usr/bin/python
from __future__ import division
import matplotlib
#import glob
from pylab import *
import numpy as np
from analysedata import *
from scipy import interpolate
from scipy.optimize import brentq
from scipy.optimize import newton
from scipy.optimize import curve_fit

#The qinbin.py script should first be run to calculate loss times for each probe position.
#This script uses a self consistent algorithm devised by S. Machida to calculate momenta and effective k.
#This script evolved from an earlier one by S. Sheehy. 

#Scroll down to "SCRIPT START" for instructions.

def calc_response_ratio(cell_tune, ncells, cell_kick_probe = [2.5, 0.5, 5.5]):
	#ratio of closed orbit distortion due to a single kick will reduce to ratio of cosine parts if betas the same at all probe locations
	#assuming cell tune is ring tune divided by number of cells can write cosine part of response matrix as cos((pi*q/cells)*(2*nc - cells)
	#where nc is the number of cells from probe to kick. 
	
	#This function assumes three probes
	
	num = np.cos(math.pi*cell_tune*(2*abs(cell_kick_probe[-1])-ncells)) - np.cos(math.pi*cell_tune*(2*abs(cell_kick_probe[1])-ncells))
	denom = np.cos(math.pi*cell_tune*(2*abs(cell_kick_probe[0])-ncells)) - np.cos(math.pi*cell_tune*(2*abs(cell_kick_probe[-1])-ncells))
	
	ratio = num/denom
	
	return ratio


def read_equ_file(dirname, filename):
    """Read Uesugi-san's equ file that generates the rf waveform"""
    
    filepath = dirname + filename
    print 'reading ', filepath

    f_data = open(filepath, 'r')
    
    a_coef = []
    for line in f_data:
        if line[0] == "A":
            try:
                num = int(line[1])
                spl = line.split()
                
                a_coef.append(float(spl[1]))
            except:
                pass


    f_data.close()
    
    return a_coef
    

def polynomial_fit_outlier(xdat, rdat, order=1, sigma_meas = 1e-3, threshold = 1e-3):
    """Successively remove data points from fit until the difference between the fit and the remaining measurements is always less than threshold"""
    
    for iteration in range(5):
        #print "iteration ",iteration
        out = np.polyfit(xdat, rdat, order, full=True)            
        poly_coef = out[0]
        
        poly = np.poly1d(poly_coef)
        
        
        diff = [abs(poly(p)-r) for p,r in zip(xdat, rdat)]
        #print "diff ",diff
        
        #print "chisq_i ",[(poly(x)-r)**2/(sigma_meas**2) for x,r in zip(xdat, rdat)]
        chisq = sum([(poly(x)-r)**2/(sigma_meas**2) for x,r in zip(xdat, rdat)])
        
       
        
        if max(diff) < threshold:
            print "converged on iteration ",iteration
            break
        else:
            #xdat = [xdat[index] for index, value in enumerate(diff) if value <= threshold]
            #rdat = [rdat[index] for index, value in enumerate(diff) if value <= threshold]   
            xdat = [xdat[index] for index, value in enumerate(diff) if value != max(diff)]
            rdat = [rdat[index] for index, value in enumerate(diff) if value != max(diff)]     
    
    return poly_coef, chisq, xdat, rdat
    
fn3 = lambda x, a0, a1, a2, a3: a0 + a1*x + a2*(x**2) + a3*(x**3)
fn4 = lambda x, a0, a1, a2, a3, a4: a0 + a1*x + a2*(x**2) + a3*(x**3) + a4*(x**4)
fn5 = lambda x, a0, a1, a2, a3, a4, a5: a0 + a1*x + a2*(x**2) + a3*(x**3) + a4*(x**4) +a5*(x**5)
fn6 = lambda x, a0, a1, a2, a3, a4, a5, a6: a0 + a1*x + a2*(x**2) + a3*(x**3) + a4*(x**4) +a5*(x**5) +a6*(x**6)
fn7 = lambda x, a0, a1, a2, a3, a4, a5, a6, a7: a0 + a1*x + a2*(x**2) + a3*(x**3) + a4*(x**4) +a5*(x**5) +a6*(x**6) +a7*(x**7)
fn8 = lambda x, a0, a1, a2, a3, a4, a5, a6, a7, a8: a0 + a1*x + a2*(x**2) + a3*(x**3) + a4*(x**4) +a5*(x**5) +a6*(x**6) +a7*(x**7) +a8*(x**8)

def polynomial_remove_outlier(xdat, rdat, sigma_meas, order=1, outliers = [], fittype='polyfit'):
    """Remove specified outliers before fitting"""
   
    xdat = [value for index, value in enumerate(xdat) if index not in outliers]
    rdat = [value for index, value in enumerate(rdat) if index not in outliers] 
    
    if fittype == "polyfit":
        #out = np.polyfit(xdat, rdat, order, full=True, cov = False) 
        out = np.polyfit(xdat, rdat, order, full=True) 
        poly_coef = out[0]
        perr = []
    else:
        if order == 3:
            fn = fn3
        elif order == 4:
            fn = fn4
        elif order == 5:
            fn = fn5
        elif order == 6:
            fn = fn6
        elif order == 7:
            fn = fn7
        elif order == 8:     
            fn = fn8
        else:
            print "define function for order ",order
            sys.exit()
            
        popt, pcov = curve_fit(fn, np.array(xdat), np.array(rdat), sigma = sigma_meas)
        perr = np.sqrt(np.diag(pcov))
        
        poly_coef = popt[::-1]
        
    poly = np.poly1d(poly_coef)
    chisq = sum([(poly(x)-r)**2/(sigma_meas**2) for x,r in zip(xdat, rdat)])
    
    #diff = [abs(poly(p)-r) for p,r in zip(xdat, rdat)]
    
    return poly_coef, chisq, xdat, rdat, perr
    
def gamma_from_meas(k,r,f,drdt,dfdt):
    """Calculate relativisitic gamma from radius r, frequency f data and their derivatives assuming k is known
    Return gamma and momentum p"""
    
    gamma_rel_sq = (k+1)/((r*dfdt)/(f*drdt) + 1)
     
    if gamma_rel_sq > 0:
        gamma_rel = np.sqrt(gamma_rel_sq)
    else:
        print "gamma^2 < 0!"
        sys.exit()
              
    beta_rel = np.sqrt(1-(1/gamma_rel_sq))
    beta_gamma_rel = beta_rel*gamma_rel
    p = 1e-6*PROTON_MASS*beta_gamma_rel
    
    return gamma_rel, p
    
def p_from_intmeas(p0,k,r,drdt, t):
    """Calculate relativisitic gamma from radius r and its derivatives assuming k is known
    Integrate over time from one momentum to the desired one. Return momentum p."""
    
    integrand = (k+1)*(drdt/r)
    integral = np.trapz(integrand, t)         
    p = p0*np.exp(integral)
            
    return p

def pdiff(ktry, p0, i, kg, r, f, drdt, dfdt, t):
    
    kg[i] = ktry
    
    gamma_1, p_1 = gamma_from_meas(kg[i], r[i], f[i], drdt[i], dfdt[i])
    p_2 = p_from_intmeas(p0, kg[:i+1], r[:i+1], drdt[:i+1], t[:i+1])
    
    pdiff = p_1 - p_2
    
    return pdiff 
            

def calcKvalue(fvals, gamma_vals, rvals):
    """Calculate k given everything else """

    
    kval_array=[]
    for i in range(len(rvals)-1):
        drr=(rvals[i+1]-rvals[i])/rvals[i]
        dff=(fvals[i+1]-fvals[i])/fvals[i]
        kval_array.append(gamma_vals[i]*gamma_vals[i]*dff/drr-(1.0-gamma_vals[i]*gamma_vals[i]))
        #print dff/drr
    return kval_array
    
def run_analysis(add_errors, rfit_order, write_orbit_data = False, probe_outliers = [[],[],[]]):    
    """ Fit orbit data and work out effective k at each probe and for mean r. If add_errors is True, add simulated errors to measurements for 
    error propagation study"""
    
    
    #store data in one list of lists
    radii_all = []
    radii_thisfit_all = [] #radius data used in fit (exclude outliers)
    xdata_all = [] #all time data
    xdata_thisfit_all = [] #time data used in fit (excludes outliers)
    xdata_err = [] #error bar applied to each xdata

    #Store polynomial fits
    rvstfn_all = [] #r(t) polynomial fit
    coef_all = [] #r(t) polynomial coeficients
    perr_all = [] #standard deviation of parameters from covariance matrix
    
    check_data = False #check raw derivative of data to select outliers.
    
    #scan order of polynomial fit and check chi-squared. 
    scan_order = False
    order_list = range(3,10) #
    

    
    #For each probe fit r vs t data
    #************************************************************************************************************
    for i_probe in range(len(data1)):

        j=data1[i_probe][0].argsort() #sort by probe position
         

        #r and t data
        radii_metres = 1e-3*(data1[i_probe][0][j] + proberadius) #position data w.r.t centre
        time = data1[i_probe][1][j]  + toffset
        
        #add random errors
        if add_errors:
            #radii_metres = radii_metres + np.array([np.random.normal(0, sigma_r_meas) for i in range(len(radii_metres))])
            radii_metres = radii_metres + np.array([np.random.uniform(-readoff_r_err, readoff_r_err) for i in range(len(radii_metres))])
            time_err = np.array([np.random.uniform(0, s) for s in data1[i_probe][2][j]])
            time = time + time_err
        xdata_err.append(data1[i_probe][2][j])
        
        xdata = list(time)
            
        radii_all.append(radii_metres)
        xdata_all.append(xdata)
        
        #polynomial fit
        poly_coef, chisq, xdata_out, rdat_out, perr = polynomial_remove_outlier(xdata, radii_metres, sigma_r_meas, order=rfit_order, outliers=probe_outliers[i_probe], fittype = 'polyfit')
            
        fittedr = np.poly1d(poly_coef)
            
        #print "chisq ",chisq,  "chisq/(N-nfit) ",chisq/(len(xdata_out) - rfit_order)
        #print "xdata ",xdata
        #print "xdata_out ",xdata_out
        #print "rdat_out ",rdat_out
            
        #print "all, filtered data points ",len(xdata), len(xdata_out)
        #print "polynomial ",fittedr
        #print "perr ",perr
 
        coef_all.append(poly_coef)
        perr_all.append(perr)
        radii_thisfit_all.append(rdat_out)
        xdata_thisfit_all.append(xdata_out)
        rvstfn_all.append(fittedr)

        leg_text = probenames[i_probe]
            
            
        if check_data:
            #sort data by time
            tr = zip(xdata, radii_metres)
            tr.sort()
            truz =  zip(*tr)
            
            xdata = truz[0]
            radii_metres = truz[1]
            
            all_data = False
            if all_data:
                deriv = [abs((radii_metres[i+1] - radii_metres[i])/(xdata[i+1] - xdata[i])) for i in range(len(radii_metres)-1)]
                xderiv = [0.5*(xdata[i+1] + xdata[i]) for i in range(len(radii_metres)-1)]
                derivderiv = [(deriv[i+1] - deriv[i])/(xderiv[i+1]-xderiv[i]) for i in range(len(deriv) - 1)]
            else:
                deriv = [abs((rdat_out[i+1] - rdat_out[i])/(xdata_out[i+1] - xdata_out[i])) for i in range(len(rdat_out)-1)]
                xderiv = [0.5*(xdata_out[i+1] + xdata_out[i]) for i in range(len(rdat_out)-1)]
            print "i_probe, deriv ",i_probe, deriv
            
            drdtfn_p = np.poly1d.deriv(fittedr)
            #print "xout, dout ",xout, dout
            subplot(211)
            if all_data:
                plot(xdata, radii_metres, color=colors[i_probe], marker='o',linewidth=0, label=probenames[i_probe])
            else:
                plot(xdata_out, rdat_out, color=colors[i_probe], marker='o',linewidth=0, label=probenames[i_probe])
            plot(xdata, fittedr(xdata),color=colors[i_probe])
            subplot(212)
            plot(xderiv, deriv, color=colors[i_probe], marker='o',linewidth=0, label=probenames[i_probe])
            plot(xdata, drdtfn_p(xdata),color=colors[i_probe])
            #plt.plot(xdata[1:-1], derivderiv, color=colors[i_probe], marker='o')
            #ylim(0, 0.0001)
            #yscale('log')
            ylim(0,1e-4)           

            #show()
            #sys.exit()
            
        if scan_order and not check_data:
            chisq_order = []
            for order in order_list:
                #poly_coef, chisq, xdata_out, rdat_out = polynomial_fit_outlier(xdata_out, rdat_out, sigma_r_meas,  order=order, threshold=6e-3) 
                out = np.polyfit(xdata_out, rdat_out, order, full=True)
                poly = np.poly1d(out[0])
                chisq = sum([(poly(x)-r)**2/(sigma_r_meas**2) for x,r in zip(xdata_out, rdat_out)])
                chisq_order.append(chisq/(len(xdata_out) - (order+1)))
                
            plt.plot(order_list, chisq_order, color=colors[i_probe], linestyle='-', marker='o',label=probenames[i_probe])
            plt.axhline(y=1, color='k',linestyle = '--')
    
    if scan_order:
        xlabel('order')
        ylabel('chi-squared/(N - order+1)')
        legend()
        savefig(resdir+pchisq)
        plt.show()
        sys.exit()
        
    if check_data:
        
        plt.subplot(211)
        ylabel('r (m)')
        ylim(4.55,5.25)
        legend(loc = 'upper left')
        plt.subplot(212)
        ylabel(r'$\frac{r_{i+1} -r_i}{t_{i+1} - t_i}$') 
        
        xlabel('time (ms)')
        
        savefig(resdir+pcheckdata)
        show()
        sys.exit()
    
    calc_mean_polynomial = True
    if calc_mean_polynomial:
        #<R> calculated by averaging polynomials for each probe
        rvstfn_mean = rvstfn_all[0]
        for ifn in range(1, len(data1)): 
            rvstfn_mean = rvstfn_mean + rvstfn_all[ifn]
        rvstfn_mean  = rvstfn_mean/len(data1)   
        
    else:
        #<R> calculated by fitting data points from all probes at once
        xdata_flat = [item for sublist in xdata_thisfit_all for item in sublist]
        radii_flat = [item for sublist in radii_thisfit_all for item in sublist]
        #poly_coef, chisq, xdata_out, rdat_out = polynomial_fit_outlier(xdata_flat, radii_flat, order=rfit_order, threshold=2e-3)
        print "xdata_flat, xdata_out ",len(xdata_flat), len(xdata_out)
        poly_coef, chisq, xdata_out, rdat_out = polynomial_remove_outlier(xdata_flat, radii_flat, sigma_r_meas, order=6, outliers=[])
        rvstfn_mean = np.poly1d(poly_coef)
        #plt.plot(xdata_flat, radii_flat,'ko')
        #plt.plot(xdata_out, rdat_out,'ko')
        #plt.show()
    

    drdtfn = np.poly1d.deriv(rvstfn_mean) #d<r>/dt
    drdt2fn = np.poly1d.deriv(drdtfn) #d^2 <r>/dt^2
    
    
    #Set time range for mean radius calculation
    min_xdata_com = 0
    max_xdata_com = min([max(f) for f in xdata_all])
    #print "min, max xdata values ",min_xdata_com, max_xdata_com
    nfit = 500

    t_fit = np.linspace(min_xdata_com, max_xdata_com, nfit)
    
    #plot orbit data and fit 
    plot_orbit_data = True
    if plot_orbit_data:
        for ip in range(len(data1)): 
            #choose plot region for this figure
            t_fit1 = np.linspace(0,max(xdata_thisfit_all[ip]), nfit)
             
            #plot fit
            plot(t_fit1, rvstfn_all[ip](t_fit1), '-', color=colors[ip], label=probenames[ip])
            
            #plot data
            xerrl = np.array([0]*len(xdata_all[ip]))
            xerrr = np.array(xdata_err[ip])
            errorbar(xdata_all[ip], radii_all[ip], yerr = 1e-3, xerr = [xerrl, xerrr], color=colors[ip], marker='*', linestyle='None')
            errorbar(xdata_thisfit_all[ip], radii_thisfit_all[ip], color=colors[ip], marker='o', linestyle='None')
            
            #plot(xdata_thisfit_all[ip], [d-p for d,p in zip(radii_thisfit_all[ip], rvstfn_all[ip](xdata_thisfit_all[ip]))], color=colors[ip], marker='o', linestyle='None') 
        
        #plot mean
        plot(t_fit, rvstfn_mean(t_fit), 'k-',label=r'$\left< r(t) \right>$',linewidth=2)
        #plot(t_fit, rvstfn_mean_raw(t_fit), 'k--',label='mean (fit all data)')
        xlabel(r't ($\mu$s)')
        ylabel('r (m)')
        legend(loc="upper left", prop={'size':12})
        #title('polynomial fit order: '+str(rfit_order))
        #ylim(4.55,5.25)
        #xlim(0, 16100)
        savefig(resdir+pfile1)
        show()



    #error analysis using sympy just for case of rfit_order ==6
    propagate_errors = False
    if propagate_errors and rfit_order == 6:
        import sympy

        x = sympy.Symbol('x')
        a0, a1, a2, a3, a4, a5, a6 = sympy.symbols('a0 a1 a2 a3 a4 a5 a6')
        sig_a0, sig_a1, sig_a2, sig_a3, sig_a4, sig_a5, sig_a6 = sympy.symbols('sig_a0 sig_a1 sig_a2 sig_a3 sig_a4 sig_a5 sig_a6')
        expr1 = a0 + a1*x + a2*x**2 + a3*x**3 + a4*x**4 + a5*x**5 + a6*x**6
        expr2 = sympy.diff(expr1, x)
        da0 = (sig_a0*sympy.diff(expr1/expr2, a0))**2
        da1 = (sig_a1*sympy.diff(expr1/expr2, a1))**2
        da2 = (sig_a2*sympy.diff(expr1/expr2, a2))**2
        da3 = (sig_a3*sympy.diff(expr1/expr2, a3))**2
        da4 = (sig_a4*sympy.diff(expr1/expr2, a4))**2
        da5 = (sig_a5*sympy.diff(expr1/expr2, a5))**2
        da6 = (sig_a6*sympy.diff(expr1/expr2, a6))**2
        #r/(dr/dt)
        expr_err_final = da0 +da1 + da2 + da3 + da4 + da5 + da6
        print "expr_err_final ",expr_err_final

        kerrcoef_fn_all = []
        for ip in range(len(data1)):
            afit = coef_all[ip][::-1]
            perr = perr_all[ip]
            print "afit ",afit
            print "perr ",perr
            
            #just check the radius fit works as expected
            expr1_sub = expr1.subs({a0:afit[0], a1:afit[1], a2:afit[2], a3:afit[3], a4:afit[4], a5:afit[5], a6:afit[6]})
        
            expr_sub1 = expr_err_final.subs({a0:afit[0], a1:afit[1], a2:afit[2], a3:afit[3], a4:afit[4], a5:afit[5], a6:afit[6]})
            expr_sub2 = expr_sub1.subs({sig_a0:perr[0], sig_a1:perr[1], sig_a2:perr[2], sig_a3:perr[3], sig_a4:perr[4], sig_a5:perr[5], sig_a6:perr[6]})
            print "expr_sub2 ",expr_sub2

            #convert to normal lambda functions in x
            radius_fn = sympy.lambdify(x, expr1_sub)
            kerrcoef_fn = sympy.lambdify(x, expr_sub2) #must still multiply by gamma^2*(dfdt/f)
            
            kerrcoef_fn_all.append(kerrcoef_fn)
        #plt.plot(np.sqrt(np.abs(expr_fn(t_fit))))
        #plt.show()
        

       
    #************************************************************************************************************
    #Use Shinji's algorithm to find self-consistent momenta and effective k
    #
    calc_effective_k = True
    if calc_effective_k:

        Ffn_win = rf_deriv_fn #simply use polynomial from equ file directly
        scrf = 1e-6
        
        dFdtfn = np.poly1d.deriv(Ffn_win)
        
        #pfn_win = findPfromt(rfdat_t_win, rfdat_p_win, order=2)
        
        #kg = np.linspace(6.95, 6.35, nfit)
        kg = np.linspace(7, 7, nfit)

        #evaluate functions at points in time window
        rvstfn_mean_tfit = rvstfn_mean(t_fit)
        Ffn_win_tfit = Ffn_win(scrf*t_fit)
        drdtfn_tfit = drdtfn(t_fit)
        dFdtfn_tfit = scrf*dFdtfn(scrf*t_fit)
        
        #check error analysis
        #for ip in range(3):
        #    kerr = np.sqrt(kerrcoef_fn_all[ip](t_fit)*(dFdtfn_tfit/Ffn_win_tfit)**2)
        #    plt.plot(kerr)
        #plt.show()
        #sys.exit()
        
        #plt.subplot(211)
        #plt.plot(t_fit, rvstfn_mean_tfit, 'm--')
        #plt.subplot(212)
        #plt.plot(t_fit, drdtfn_tfit, 'r-')
        #plt.show()    
        
    
        #prf_t_fit = pfn_win(t_fit) #p from rf table
    
        #initial energy assumed to be 11 MeV
        ke0 = 11e6
        mom0 = ke_to_mom(PROTON_MASS, ke0)
        beta_rel = ke_to_relativistic_beta(ke0, PROTON_MASS)
        beta_gamma_rel =  ke_to_relativistic_beta_gamma(ke0, PROTON_MASS)
        gamma0 = beta_gamma_rel/beta_rel
        #print "ke, mom, gamma0 ",ke0, mom0, gamma0

        #calculate p at each point in time window

        p0 = 1e-6*mom0 #assume injection energy at start
        #p0 = rfdat_p_win[0]
        
        #print "p0 ",p0
        
        
        
        #scaling index and momenta for <R>, i.e. R(t) averaged over probes
        #-----------------------------------------------------------------
        p1_r_mean = []
        p2_r_mean = []
        for il in range(nfit):
            
            kmin = (rvstfn_mean_tfit[il]*dFdtfn_tfit[il])/(Ffn_win_tfit[il]*drdtfn_tfit[il])

            
            if il == 0:
                #k1 = 6.95
                #k2 = 7.0
                k1 = kmin
            else:
                #k1 = sol
                k1 = kmin
            
            
            #Brentq root-finding method requires specifiy k1,k2 and pdiff must change sign at these values
            #sol = brentq(pdiff, k1, k2, args = (p0, il, kg, rvstfn_mean_tfit, Ffn_win_tfit, drdtfn_tfit, dFdtfn_tfit, t_fit))
            
            #Newton-Raphson or Sextant method
            sol = newton(pdiff, k1, args = (p0, il, kg, rvstfn_mean_tfit, Ffn_win_tfit, drdtfn_tfit, dFdtfn_tfit, t_fit))
         
            kg[il] = sol
            
            #print "iter, k, p diff ",il, sol,  pdiff(sol, p0, il, kg, rvstfn_mean_tfit, Ffn_win_tfit, drdtfn_tfit, dFdtfn_tfit, t_fit)

            gam, p_il_1 = gamma_from_meas(kg[il], rvstfn_mean_tfit[il], Ffn_win_tfit[il], drdtfn_tfit[il], dFdtfn_tfit[il])
            p_il_2 = p_from_intmeas(p0, kg[:il+1], rvstfn_mean_tfit[:il+1], drdtfn_tfit[:il+1], t_fit[:il+1])
            
            #print "k sol, min, p diff ",sol, kmin, p_il_1 - p_il_2
            
            p1_r_mean.append(p_il_1)
            p2_r_mean.append(p_il_2)

                    
        #solution k
        k_r_mean = list(kg)
        
        #print "ksol ",k_r_mean
        #print "p1 ",p1_sol
        #print "p2 ",p2_sol
        #print "rfdat table p ",prf_t_fit
        
        plot_cod_amp = False
        if plot_cod_amp:
            #cod_amp = rvstfn_all[1](t_fit)-rvstfn_all[2](t_fit) #2014 3 probe case
            cod_amp = rvstfn_all[3](t_fit)-rvstfn_all[4](t_fit)
            max_cod = max(cod_amp)
            min_cod = min(cod_amp)
            index_max = list(cod_amp).index(max_cod)
            index_min = list(cod_amp).index(min_cod)
            
            p_max_cod = p1_r_mean[index_max]
            p_min_cod = p1_r_mean[index_min]
            print "min cod, p ",min(cod_amp), p_min_cod
            print "max_cod, p ",max_cod, p_max_cod
            print "p at end ",p1_r_mean[-1]
            print "p ratio ",p_max_cod/p_min_cod
            print "cod ratio ",max_cod/min_cod
            print "expected COD at end with constant tune ",min_cod*p_min_cod/p_max_cod
            
            cod_pred = sin(math.pi*3.65)*min_cod*p_min_cod/(sin(math.pi*3.85)*p_max_cod)
            print "expected COD at end with varying tune ",cod_pred
            
            print "residual COD with varying tune ",max_cod - cod_pred 
            
            
            #ratio = 
            probe_ratio = (rvstfn_all[0](t_fit)-rvstfn_all[2](t_fit))/(rvstfn_all[1](t_fit)-rvstfn_all[0](t_fit))
            
            print "probe ratio initial, at 12ms, at end ",probe_ratio[index_min], probe_ratio[index_max], probe_ratio[-1]
            
            rerr = 0.5e-3
            nerr = 500
            rerr_norm1 = np.random.normal(0, rerr, nerr)
            rerr_norm2 = np.random.normal(0, rerr, nerr)
            rerr_norm3 = np.random.normal(0, rerr, nerr)
            
            prat_list = []
            for r1,r2,r3 in zip(rerr_norm1, rerr_norm2, rerr_norm3):
                probe_ratio = (rvstfn_all[0](t_fit)-rvstfn_all[2](t_fit) + r1 - r2)/(rvstfn_all[1](t_fit)-rvstfn_all[0](t_fit) + r3 - r1)
                prat = probe_ratio[index_min]
                prat_list.append(prat)
                
            print "error (numerical) at index_min ",np.std(prat_list)

            
                
            ratio_err = probe_ratio*((2*rerr/(rvstfn_all[0](t_fit)-rvstfn_all[2](t_fit)))**2 + (2*rerr/(rvstfn_all[1](t_fit)-rvstfn_all[0](t_fit)))**2)**0.5
            

            
            print "ratio_err inital, at 12ms, at end ",ratio_err[index_min], ratio_err[index_max],ratio_err[-1]
            
            print "ratio at tunes 3.65, 3.85 ",calc_response_ratio(3.65/12, 12), calc_response_ratio(3.85/12, 12)
            
            #calculate error in ratio calculation from measured tune uncertainty
            tune_err = 0.005            
            tune_err_norm = np.random.normal(0, tune_err, 500)
            print "ratio err (numerical) ", np.std(calc_response_ratio((3.65+tune_err_norm)/12, 12)),np.std(calc_response_ratio((3.85+tune_err_norm)/12, 12))
            #print "ratio error from tune ",calc_response_ratio_err(3.65/12, 12, 0.005/12)
            
            

            
            plt.subplot(211)
            plot(p1_r_mean, rvstfn_all[1](t_fit)-rvstfn_all[2](t_fit), '-', color='r')
            #plt.axhline(y=max_cod)
            
            #plt.axvline(x=p_max_cod)
            plt.xlabel("p [MeV/c]")
            plt.ylabel("COD (F5-F7)")
            plt.savefig("cod_amp")
            
            plt.subplot(212)
            plot(p1_r_mean, probe_ratio, '-', color='k')
            plt.show()
            
            sys.exit()
        
        
        #plt.plot(t_fit, p1_r_mean)
        #plt.show()


        #calculate effective k and momenta for each probe
        #-----------------------------------------------------------------
        k_probe_all = []
        p1_probe_all = []
        p2_probe_all = []
        #p_probe_rftab = []
        time_data_probe  = []
        r_fit_all = []
        pdat_all = []
        for ip in range(len(data1)):
            
            rvstfn = rvstfn_all[ip]
            drdtfn_probe = np.poly1d.deriv(rvstfn)
            
            
            #time_data = list(xdata_thisfit_all[ip])
            
            #time_data = np.linspace(0, max(xdata_thisfit_all[ip]), nfit)
            time_data = np.linspace(0,max(xdata_thisfit_all[ip]), nfit)
            #time_data = t_fit #use same time range as for mean r
            time_data_probe.append(time_data)
            
            r_fit = rvstfn(time_data)
            drdt_fit = drdtfn_probe(time_data)
            F_fit = Ffn_win(scrf*time_data)
            dFdt_fit = scrf*dFdtfn(scrf*time_data)
            
            r_fit_all.append(r_fit)#save for later plotting
            
            kg = np.linspace(7, 7, len(time_data))
            

            #p_probe_rftab.append(pfn_win(time_data))

            
            #k1 = 6.9
            #print "r_fit ",r_fit
            #print "drdt fit ",drdt_fit
            
            #p0 = pfn_win(xdata_thisfit_all[ip][0]) #p from rf table
            p0 = 1e-6*mom0
            
            gamma_vals = []
            
            k_probe = []
            p1_probe = []
            p2_probe = []
            
            for il in range(len(time_data)):
                
                kmin = (r_fit[il]*dFdt_fit[il])/(F_fit[il]*drdt_fit[il]) #gamma^2 = 1 condition met if k = kmin
                
                k1 = kmin         
                sol = newton(pdiff, k1, args = (p0, il, kg, r_fit, F_fit, drdt_fit, dFdt_fit, time_data))
                #print "k sol, min ",sol, kmin
                
                kg[il] = sol
                
                gam, p_il_1 = gamma_from_meas(kg[il], r_fit[il], F_fit[il], drdt_fit[il], dFdt_fit[il])
                p_il_2 = p_from_intmeas(p0, kg[:il+1], r_fit[:il+1], drdt_fit[:il+1], time_data[:il+1])                    
                
                #print "p_il_1, p_il_2 ",p_il_1, p_il_2
                
                k_probe.append(sol)
                p1_probe.append(p_il_1)
                p2_probe.append(p_il_2)
                
                gamma_vals.append(gam)
            
            p_t_fit = np.polyfit(time_data, p1_probe, 3) 

            p_t_fitp = np.poly1d(p_t_fit)
            pdat = p_t_fitp(xdata_all[ip])
            #plt.plot(time_data, p1_probe,'k-')
            #plt.plot(time_data, p_t_fitp(time_data),'ro')
            #plt.plot(xdata_all[ip], pdat,'bo')
            
            pdat_all.append(pdat)
            
            #plt.show()
            #sys.exit()
                
            #kd = calcKvalue(F_fit, gamma_vals, r_fit)
            #print "kd ",kd, "len(kd) ",len(kd)
            #print "len time_data ",len(time_data)
            
            #if ip == 1:
                #plt.plot(time_data[:-1], kd, 'ko')
                #plt.show()

                
            k_probe_all.append(k_probe)
            p1_probe_all.append(p1_probe)
            p2_probe_all.append(p2_probe)
            
                
        #dpp_est = (rfdat_p_win[-1] - rfdat_p_win[0])/rfdat_p_win[0]
        #print "estimated freq change ",((1/(gamma0**2)) - (1/(kg+1)))*dpp_est*rfdat_F_win[0], " actual ",rfdat_F_win[-1] - rfdat_F_win[0]

        plot_p_vs_t = False
        if plot_p_vs_t:
            plt.plot(t_fit, p1_r_mean,'k-')
            for ip in range(len(data1)):
                plt.plot(time_data_probe[ip], p1_probe_all[ip],color=colors[ip],label=probenames[ip])
            plt.ylabel('p data (MeV/c)')
            plt.xlabel(r'time ($\mu$s)')
            #savefig(resdir+pfile3)
            plt.show()
        #plt.show()
        
        plot_k_vs_p = True
        plot_k_vs_t = False
        if plot_k_vs_p:
            #plt.subplot(211)
            
            for ip in range(len(data1)):
                plt.plot(p1_probe_all[ip], k_probe_all[ip],color=colors[ip],label=probenames[ip])
            plt.plot(p1_r_mean, k_r_mean, 'k-', label='mean radius', linewidth=2)
            plt.ylabel('effective k')
            plt.xlim(1e-6*mom0, max([max(p) for p in p1_probe_all]))
            #plt.subplot(212)
            #plt.plot(p1_r_mean, drdtfn_tfit, 'b.-')
            plt.xlabel('p (MeV/c)')
            ylim(6, 8.5)
            plt.ylabel('effective k')
            plt.legend(loc = "upper left")
            title('polynomial fit order: '+str(rfit_order))
            savefig(resdir+pfile5)
            plt.show()
            
        if plot_k_vs_t:
            #plt.subplot(311)
            
            for ip in range(len(data1)):
                plt.plot(time_data_probe[ip], k_probe_all[ip], color=colors[ip],label=probenames[ip])
            plt.plot(t_fit, k_r_mean, 'k-',linewidth=2)
            plt.ylabel('effective k')
            #plt.subplot(312)
            #plt.plot(t_fit, drdtfn_tfit, 'b.-')
            #plt.plot(t_fit, drdt2fn(t_fit), 'g.-')
            
            #plt.plot(t_fit, (rvstfn_mean_tfit*dFdtfn_tfit)/(Ffn_win_tfit*drdtfn_tfit), 'b.-')
            #plt.ylabel('dF/dt')  
            #plt.subplot(313)
            #plt.plot(t_fit, (dFdtfn_tfit)/(drdtfn_tfit), 'm.-')

            plt.xlabel(r'time ($\mu$s)')
            #plt.ylabel('dr/dt')
            #ylim(6, 8.5)
            plt.ylabel('effective k')
            plt.legend(loc = "upper left")
            #title('polynomial fit order: '+str(rfit_order))
            savefig(resdir+pfile4)
            plt.show()
            
        plot_r_vs_p = False
        if plot_r_vs_p:
            #plt.subplot(211)
            
            for ip in range(len(data1)):
                plt.plot(p1_probe_all[ip], r_fit_all[ip],color=colors[ip],label=probenames[ip])
                plt.plot(pdat_all[ip], radii_all[ip], color=colors[ip], marker='o',linestyle="None")
            #plt.plot(p1_r_mean, k_r_mean, 'k-', label='mean radius', linewidth=2)
            #plt.ylabel('effective k')
            #plt.xlim(1e-6*mom0, max([max(p) for p in p1_probe_all]))
            #plt.subplot(212)
            #plt.plot(p1_r_mean, drdtfn_tfit, 'b.-')
            plt.xlabel('p (MeV/c)')
            #ylim(6, 8.5)
            #plt.ylabel('effective k')
            #plt.legend(loc = "upper left")
            #title('polynomial fit order: '+str(rfit_order))
            #savefig(resdir+pfile5)
            plt.xlim(143, 356)
            plt.show()
            sys.exit()
        

    if write_orbit_data:
        print "write results to file"
        for ip in range(len(data1)):
            #write data to file
            f = open('qinbin_data_proc_'+probenames[ip]+'.txt', 'w')
            f.write("%9s %5s %5s %5s %s \n" % ("time", "pos", "p", "t_err", "included"))
            f.write("%9s %5s %5s %5s \n" % ("us", "m","MeV/c", "us"))
            #f.write("us   m   us \n")
            i1 = 0
            for x, r, p, xerr_r in zip(xdata_all[ip], radii_all[ip], pdat_all[ip], xdata_err[ip]):
                f.write("%9.3f %5.3f %5.3f %5.3f %s \n" % (x , r, p, xerr_r, str(i1 not in probe_outliers[ip])))
                i1 = i1 + 1
            f.close()

            t_fit1 = np.linspace(0,max(xdata_thisfit_all[ip]), nfit)
            
            f = open('r_fit_'+probenames[ip]+'.txt', 'w')
            f.write("%9s %5s %5s \n" % ("time","pos","p"))
            f.write("%9s %5s %5s \n" % ("us", "m", "MeV/c"))
            #f.write("us   m   us \n")
            for x, r, p in zip(t_fit1, rvstfn_all[ip](t_fit1), p1_probe_all[ip]):
                f.write("%9.3f %5.3f %5.3f \n" % (x , r, p))
            f.close()
            
            f = open('r_fit_param_'+probenames[ip]+'.txt', 'w')
            print >>f,"fit coef (highest to lowest order"
            for c in coef_all[ip]:
                print >>f, str(c), 
            f.close()
            
        f = open('r_fit_mean.txt', 'w')
        f.write("%9s %5s %5s \n" % ("time", "pos","p"))
        f.write("%9s %5s %5s \n" % ("us", "m", "MeV/c"))
        #f.write("us   m   us \n")
        for x, r, p in zip(t_fit, rvstfn_mean(t_fit), p1_r_mean):
            f.write("%9.3f %5.3f %5.3f \n" % (x , r, p))
        f.close()        
            
            
    return t_fit, rvstfn_all, coef_all, p1_probe_all, time_data_probe,  k_probe_all, p1_r_mean, k_r_mean



#*****************************************************************************************

#if __name__ == "__main__":
#SCRIPT START

#Instructions
#Give location of orbit data "datadir1" and "subdir". Data will be in datadir/subdir
#Datafiles for each probe given by the list"set1". For example set1=["QinBin_F1.txt", "QinBin_F5.txt", "QinBin_F7.txt"]
#Choose a polynomial fitting order by setting the parameter "fit_order". 
#Exclude datapoints by index for each probe by setting probe_outliers.


#set up  useful stuff
PROTON_MASS = 938.27203e6

plt.rcParams.update({'axes.labelsize':18})
plt.rcParams.update({'xtick.labelsize':14})
plt.rcParams.update({'ytick.labelsize':14})
matplotlib.rcParams['legend.fontsize']=12

from sys import platform as _platform
if _platform == "linux" or _platform == "linux2":
    # linux
    #datadir1="/home/pvq65952/accelerators/KURRI_ADS/march14_exp/qinbin_results"
    #datadir1="/home/pvq65952/accelerators/KURRI_ADS/scripts/beamsize"
    datadir1="/home/pvq65952/accelerators/KURRI_ADS/paper2015/data"
    #datadir1="/home/pvq65952/accelerators/KURRI_ADS/scripts/qinbin/Suzie_150611/"
    rfdir="/home/pvq65952/accelerators/KURRI_ADS/march14_exp/rf_pattern/"
elif _platform == "darwin" or _platform == "OS X":
    # OS X
    #datadir1="/Users/pvq65952/accelerators/KURRI_ADS/march14_exp/qinbin_results"
    datadir1="/Users/pvq65952/accelerators/KURRI_ADS/scripts/qinbin"
    #datadir1 = "/Users/pvq65952/accelerators/KURRI_ADS/paper2015/data"
    rfdir="/Users/pvq65952/accelerators/KURRI_ADS/march14_exp/rf_pattern/"    
elif _platform == "win32" or _platform == "win64":
    # Windows...
    print "Windows detected. Specify data directories"
    sys.exit()
    

#fig=figure(num=1, figsize=(6.5, 10.5), dpi=150, facecolor='w', edgecolor='k')   

#select data subdirectory    
#subdir = "310314"
#subdir = "2014-03-31"
#subdir = "D0_0amps"
#subdir = "D0_140amps"
#subdir = "D0_400amps"
#subdir = "D0_550amps"
#subdir = "D0_700amps"
    
#sub directories for paper    
subdir = "2015-06-24/900A"
#subdir = "keff_p_qinbin" #31/3/2014 long-time base Qinbin data here
#subdir = "cod_qinbin"

fit_order = 3 #order of polynomial fit to orbit data subsequently used to calculate effective k.

#Directory for plots
resdir="Results/"
figformat = '.pdf'
    
#Filenames for plots
pfile1 ="r_vs_t"+figformat
pfile2="rfpattern"+figformat
pfile3="pdiff"+figformat    
pfile4="k_vs_t"+figformat
pfile5="k_vs_p"+figformat
pchisq = "chisquared"+figformat
pcheckdata = "checkdata"+figformat
pfile_kerr = "k_vs_t_werr"+figformat

datadir = datadir1+"/"+subdir+"/"

#Filenames for sets of data
if subdir == "keff_p_qinbin":
    set1=["QinBin_F1.txt", "QinBin_F5.txt", "QinBin_F7.txt"]
    probe_outliers = [[14],[16],[15,16]] #index of probe outliers after reordering 
elif subdir == "cod_qinbin":
    #uncomment "set1", "probe_outliers" pair 
    #Corrector scan
    #set1=["QinBin_400A_F1.txt", "QinBin_400A_F5.txt", "QinBin_400A_F7.txt"]
    #set1=["QinBin_550A_F1.txt", "QinBin_550A_F5.txt", "QinBin_550A_F7.txt"]
    set1=["QinBin_700A_F1.txt", "QinBin_700A_F5.txt", "QinBin_700A_F7.txt"]
    #29/6/15 update probe outliers. Exclude cases with no turns.
    #probe_outliers = [[0],[0],[]] #400 A
    #probe_outliers = [[],[0],[]] #550 A
    probe_outliers = [[],[0,1,2],[0,1]] #700A
else:
    set1=["QinBin_F01.txt", "QinBin_F02.txt", "QinBin_F03.txt","QinBin_F05.txt","QinBin_F12.txt"]
    #probe_outliers = [[23],[],[],[],[0,24,25]] # 24/6/2015, 700A case
    probe_outliers = [[],[],[],[],[]] # 24/6/2015, 900A case
    
#D1 current - not corrector!
#set1=["QinBin_0A_F1.txt", "QinBin_0A_F5.txt", "QinBin_0A_F7.txt"]
#set1=["QinBin_140A_F1.txt", "QinBin_140A_F5.txt", "QinBin_140A_F7.txt"]
#probe_outliers = [[0,1,2],[0,1,2],[]] #0A
#probe_outliers = [[0,1],[0,1],[]] #140 A


#probe_outliers = [[],[],[]]
#31/12/14 (k effective data)
   

#set1 =["D0_0amps_F1.txt" ,"D0_0amps_F5.txt" ,"D0_0a[1e-3]*len(ydata)mps_F7.txt"] #468A
#set1=["D0_550Aamps_F1.txt" ,"D0_550Aamps_F5.txt" ,"D0_550Aamps_F7.txt"] #550A 27/3/2014
#set1 =["D0_140amps_F1.txt", "D0_140amps_F5.txt", "D0_140amps_F7.txt"] #468A
#set1=["D0_400Aamps_F1.txt", "D0_400Aamps_F5.txt", "D0_400Aamps_F7.txt"] #400A 27/3/14
#set1=["D0_700amps_F1.txt", "D0_700amps_F5.txt", "D0_700amps_F7.txt"] #27/3/2014

#get the radial position at each time from each probe
data1=readset(datadir, set1) #gets the data


#read .equ file for rf waveform info
a_coef = read_equ_file(rfdir, "20140116_svk20BS.equ")
#a_coef = read_equ_file(rfdir, "20140115_svk02BS.equ")
rf_fn = np.poly1d(a_coef[::-1]+[0])
#print "rf fn ",rf_fn
rf_deriv_fn = rf_fn.deriv()

#probenames=["F1","F5","F7"] #2014 
probenames=["F01","F02","F03","F05","F12"] #2015
colors = ['b', 'r','c','g','m']
proberadius=4300.0 #mm 
#toffset= -93.28
toffset = 0
    
#Read rf table (just for comparison purposes)
#fittedpfunc, fittedFfunc, rfdat_t, rfdat_F, rfdat_E, rfdat_p = rf_convert(rfdir) 
    

 
#Error in orbit measurements. This is used in polynomial fitting. 
#If add_errors is True it is also used to work out uncertainty of the the entire analysis chain.
readoff_r_err = 1e-3
sigma_r_meas = (1/3**0.5)*readoff_r_err

#Add simulated errors to work out uncertainty in results. These are disregarded if add_errors =False
#sigma_t_meas = 90
#nrun = 5
#seed = 940289
#np.random.seed(seed)

#r_staterr_list = [np.random.normal(0, sigma_r_meas) for i in range(nrun)]


#print "data1 ",data1[0]
#NB: rf_order = 1 for COD measurements, 4 for keff measurement over greater momentum range
t_fit0, rvstfn_all0, coef_all0, p1_probe_all0, t_probe_all0, k_probe_all0, p1_r_mean0, k_r_mean0 = run_analysis(add_errors = False, rfit_order = fit_order, probe_outliers = probe_outliers, write_orbit_data = True)



sys.exit()

#Rest of script calculates dispersion, does some error analysis.

#dispersion
disp_p10 = rvstfn_all0[0](t_fit0)/(1+np.array(k_probe_all0[0]))
disp_p20 = rvstfn_all0[1](t_fit0)/(1+np.array(k_probe_all0[1]))
disp_p30 = rvstfn_all0[2](t_fit0)/(1+np.array(k_probe_all0[2]))

rvstfn_mean0 =  (rvstfn_all0[0] + rvstfn_all0[1] + rvstfn_all0[2])/3  #mean r(t)     
disp_mean0 = rvstfn_mean0(t_fit0)/(1+np.array(k_r_mean0))

#plt.plot(t_fit0, disp_p10)
#plt.plot(t_fit0, disp_p20)
#plt.plot(t_fit0, disp_p30)
#plt.show()


check_dispersion = True
if check_dispersion:
    out = np.polyfit(p1_r_mean0, rvstfn_mean0(t_fit0), 4, full=True) 
    r_p_fn = np.poly1d(out[0]) #r(p)
    drdp_fn = np.poly1d.deriv(r_p_fn)

    #plt.plot(t_fit0, p1_r_mean0)
    #plt.plot(p1_r_mean0, rvstfn_mean0(t_fit0))
    #plt.plot(p1_r_mean0, r_p_fn(p1_r_mean0), 'r')

    plt.plot(p1_r_mean0*drdp_fn(p1_r_mean0))
    #plt.plot(p1_r_mean0, k_r_mean0)
    #plt.plot(p1_r_mean0, disp_mean0)
    plt.show()
    sys.exit()


order_scan = False
if order_scan:
    order_list = range(3,9)
    
    print "order ",order_list
    for order in order_list:
        t_fit, rvstfn_all, coef_all, p1_probe_all,t_probe_all0, k_probe_all, p1_r_mean, k_r_mean = run_analysis(add_errors = False, rfit_order = order, probe_outlier = probe_outliers)
        rfn_mean = (rvstfn_all[0] + rvstfn_all[1] + rvstfn_all[2])/3
        plt.plot(t_fit, k_r_mean, linewidth =1+ 0.5*order_list.index(order), label=str(order))
    plt.legend(loc = "upper left")
    plt.show()
    sys.exit()
    
coef_p1_samp = []
coef_p2_samp = []
coef_p3_samp = []
k_p1_samp = []
k_p2_samp = []
k_p3_samp = []
p_p1_samp = []
p_p2_samp = []
p_p3_samp = []
k_r_mean_samp = []
p_r_mean_samp = []

disp_p1_samp = []
disp_p2_samp = []
disp_p3_samp = []
disp_r_mean_samp = []
for ir in range(nrun):
    t_fit, rvstfn_all, coef_all, p1_probe_all, t_probe_all0, k_probe_all, p1_r_mean, k_r_mean = run_analysis(add_errors = True, rfit_order = 4)

    coef_p1_samp.append(list(coef_all[0]))
    coef_p2_samp.append(list(coef_all[1]))
    coef_p3_samp.append(list(coef_all[2]))
    
    k_p1_samp.append(k_probe_all[0])
    k_p2_samp.append(k_probe_all[1])
    k_p3_samp.append(k_probe_all[2])
    k_r_mean_samp.append(k_r_mean)
    
    p_p1_samp.append(p1_probe_all[0])
    p_p2_samp.append(p1_probe_all[1])
    p_p3_samp.append(p1_probe_all[2])
    p_r_mean_samp.append(p1_r_mean)
    
    
    disp_p1_samp.append(rvstfn_all[0](t_probe_all0[0])/(1+np.array(k_probe_all[0])))
    disp_p2_samp.append(rvstfn_all[1](t_probe_all0[1])/(1+np.array(k_probe_all[1])))
    disp_p3_samp.append(rvstfn_all[2](t_probe_all0[2])/(1+np.array(k_probe_all[2])))
    
    rvstfn_mean =  (rvstfn_all[0] + rvstfn_all[1] + rvstfn_all[2])/3  #mean r(t)
    disp_r_mean_samp.append(rvstfn_mean(t_fit)/(1+np.array(k_r_mean)))
    #plt.subplot(211)
    #for ip in range(3):
    #    plot(t_fit, rvstfn_all[ip](t_fit), '-', color=colors[ip], label=probenames[ip])
    #plt.subplot(212)
    #for ip in range(3):
    #    plt.plot(p1_probe_all[ip], k_probe_all[ip],'--',color=colors[ip],label=probenames[ip])
    
    #plt.plot(rvstfn_mean(t_fit)/(1+np.array(k_r_mean)))
#plt.show()

    
coef_p1_t = np.transpose(coef_p1_samp)
coef_p2_t = np.transpose(coef_p2_samp) 
coef_p3_t = np.transpose(coef_p3_samp) 
#generate polynomial from mean coeficients
p1_poly_mean = np.poly1d([sum(c)/len(c) for c in coef_p1_t])
p2_poly_mean = np.poly1d([sum(c)/len(c) for c in coef_p2_t])
p3_poly_mean = np.poly1d([sum(c)/len(c) for c in coef_p3_t])

show_mean_rfit = False
if show_mean_rfit:
    plt.plot(t_fit, p1_poly_mean(t_fit))
    plt.plot(t_fit, p2_poly_mean(t_fit))
    plt.plot(t_fit, p3_poly_mean(t_fit))
    plt.show()


#mean and stddev of effective k values calculated at each probe
mean_k_p1 = [sum(c)/len(c) for c in np.transpose(k_p1_samp)]
mean_k_p2 = [sum(c)/len(c) for c in np.transpose(k_p2_samp)]
mean_k_p3 = [sum(c)/len(c) for c in np.transpose(k_p3_samp)]
std_k_p1 = [np.std(c) for c in np.transpose(k_p1_samp)]
std_k_p2 = [np.std(c) for c in np.transpose(k_p2_samp)]
std_k_p3 = [np.std(c) for c in np.transpose(k_p3_samp)]

#momentum error bars
std_p_p1 = [np.std(c) for c in np.transpose(p_p1_samp)]
std_p_p2 = [np.std(c) for c in np.transpose(p_p2_samp)]
std_p_p3 = [np.std(c) for c in np.transpose(p_p3_samp)]


#dispersion error bars
std_disp_p1 = [np.std(c) for c in np.transpose(disp_p1_samp)]
std_disp_p2 = [np.std(c) for c in np.transpose(disp_p2_samp)]
std_disp_p3 = [np.std(c) for c in np.transpose(disp_p3_samp)]

mean_k_r_mean = [sum(c)/len(c) for c in np.transpose(k_r_mean_samp)]
std_k_r_mean = [np.std(c) for c in np.transpose(k_r_mean_samp)]

std_p_r_mean = [np.std(c) for c in np.transpose(p_r_mean_samp)]
std_disp_r_mean = [np.std(c) for c in np.transpose(disp_r_mean_samp)]
#print "mean_k_p1 ",mean_k_p1
#print "std_k_p1 ",std_k_p1

#mean fit for nominial case
rvstfn_mean0 = (rvstfn_all0[0] + rvstfn_all0[1] + rvstfn_all0[2])/3
print "rvstfn_mean0 ",rvstfn_mean0 


#plt.subplot(212)
#plt.errorbar(t_fit, rvstfn_mean0(p1_r_mean0))
#plt.show()

show_dispersion = False
if show_dispersion:
    plt.errorbar(t_fit0, disp_p10, yerr = std_disp_p1)
    plt.errorbar(t_fit0, disp_p20, yerr = std_disp_p2)
    plt.errorbar(t_fit0, disp_p30, yerr = std_disp_p3)
    plt.show()

    
show_p_vs_t = False
if show_p_vs_t:
    print "std_p_r_mean ",std_p_r_mean
    print "std_p_p1 ",std_p_p1
    print "std_p_p2 ",std_p_p2
    print "std_p_p3 ",std_p_p3
    plt.errorbar(t_fit, p1_r_mean0, yerr = std_p_r_mean)
    plt.errorbar(t_probe_all0[0], p1_probe_all0[0], yerr = std_p_p1)
    plt.errorbar(t_probe_all0[1], p1_probe_all0[1], yerr = std_p_p1)
    plt.errorbar(t_probe_all0[2], p1_probe_all0[2], yerr = std_p_p1)
    plt.show()
    sys.exit()

ne = 98
ie1= 5

print "error bars (final point), mean, p1,p2,p3 ",std_k_r_mean[-1], std_k_p1[-1], std_k_p2[-1], std_k_p3[-1]
print "final k at each probe ", mean_k_p1[-1],mean_k_p2[-1],mean_k_p3[-1]
lastk_p1 = np.transpose(k_p1_samp)[-1]
mean_var_p1 = [sum(lastk_p1[0:i])/len(lastk_p1[0:i]) for i in range(1,len(lastk_p1))]
lastk_p2 = np.transpose(k_p2_samp)[-1]
mean_var_p2 = [sum(lastk_p2[0:i])/len(lastk_p2[0:i]) for i in range(1,len(lastk_p2))]
lastk_p3 = np.transpose(k_p3_samp)[-1]
mean_var_p3 = [sum(lastk_p3[0:i])/len(lastk_p3[0:i]) for i in range(1,len(lastk_p3))]
print "last k ", lastk_p1

print "mean_var p1",mean_var_p1
print "mean_var p2",mean_var_p2
print "mean_var p3",mean_var_p3
#plot(mean_var_p1,color=colors[0],marker='o')
#plot(mean_var_p2,color=colors[1],marker='o')
#plot(mean_var_p3,color=colors[2],marker='o')
#show()



write_results = True
if write_results:
    
    #k, p from mean radius
    f = open('k_r_mean.txt', 'w')
    f.write("%9s %7s %5s %5s  %6s %7s \n" % ("time", "p", "k", "error", "disp", "error"))
    f.write("%9s %7s %5s %5s  %6s %7s \n" % ("us", "MeV/c", "   ", "   ","m","m"))
    #f.write("us   m   us \n")
    for x, p, k, kerr, disp, disp_err in zip(t_fit, p1_r_mean0, k_r_mean0, std_k_r_mean, disp_mean0, std_disp_r_mean):
        f.write("%9.3f %5.3f %5.3f %5.3f %7.5f %7.5f \n" % (x , p,  k, kerr, disp, disp_err))
    f.close()
    
    #k, p from F1 radius fit 
    f = open('k_r_F1.txt', 'w')
    f.write("%9s %7s %5s %5s  %6s %7s \n" % ("time", "p", "k", "error", "disp", "error"))
    f.write("%9s %7s %5s %5s  %6s %7s \n" % ("us", "MeV/c", "   ", "   ","m","m"))
    #f.write("us   m   us \n")
    for x, p, k, kerr, disp, disp_err  in zip(t_probe_all0[0], p1_probe_all0[0], k_probe_all0[0], std_k_p1, disp_p10, std_disp_p1):
        print "disp, disp_err ",disp, disp_err
        f.write("%9.3f %5.3f %5.3f %5.3f %7.5f %7.5f \n" % (x , p,  k, kerr, disp, disp_err))
    f.close()      

    #k, p from F5 radius fit 
    f = open('k_r_F5.txt', 'w')
    f.write("%9s %5s %6s %6s %4s %5s \n" % ("time", "p", "k", "error", "disp", "error"))
    f.write("%9s %5s %3s  %3s %5s  %5s \n" % ("us", "MeV/c", "   ", "   ","m","m"))
    #f.write("us   m   us \n")
    for x, p, k, kerr, disp, disp_err  in zip(t_probe_all0[1], p1_probe_all0[1], k_probe_all0[1], std_k_p2, disp_p20, std_disp_p2):
        f.write("%9.3f %5.3f %5.3f %5.3f %7.5f %7.5f \n" % (x , p,  k, kerr, disp, disp_err))
    f.close() 

    #k, p from F7 radius fit 
    f = open('k_r_F7.txt', 'w')
    f.write("%9s %5s %6s %6s %4s %5s \n" % ("time", "p", "k", "error", "disp", "error"))
    f.write("%9s %5s %3s  %3s %5s  %5s \n" % ("us", "MeV/c", "   ", "   ","m","m"))
    #f.write("us   m   us \n")
    for x, p, k, kerr, disp, disp_err  in zip(t_probe_all0[2], p1_probe_all0[2], k_probe_all0[2], std_k_p3, disp_p30, std_disp_p3):
        f.write("%9.3f %5.3f %5.3f %5.3f %7.5f %7.5f \n" % (x , p,  k, kerr, disp, disp_err))
    f.close() 
    

plt.errorbar(t_fit, k_r_mean0, color = 'k')
plt.errorbar(t_fit0[ie1::ne], k_r_mean0[ie1::ne], yerr = std_k_r_mean[ie1::ne], ecolor='k',linestyle='None')
#plt.errorbar(p1_r_mean0, mean_k_r_mean, yerr = std_k_r_mean, color = 'k')

plt.errorbar(t_probe_all0[0], k_probe_all0[0], color=colors[0], label=probenames[0])
plt.errorbar(t_probe_all0[0][ie1::ne], k_probe_all0[0][ie1::ne], yerr= std_k_p1[ie1::ne], ecolor=colors[0],linestyle='None')

plt.errorbar(t_probe_all0[1], k_probe_all0[1], color=colors[1], label=probenames[1])
plt.errorbar(t_probe_all0[1][ie1::ne], k_probe_all0[1][ie1::ne], yerr= std_k_p2[ie1::ne], ecolor=colors[1],linestyle='None')

plt.errorbar(t_probe_all0[2], k_probe_all0[2], color=colors[2], label=probenames[2])
plt.errorbar(t_probe_all0[2][ie1::ne], k_probe_all0[2][ie1::ne], yerr= std_k_p3[ie1::ne], ecolor=colors[2],linestyle='None')

#plt.xlim(0, max([max(p) for p in t_probe_all0])+100)


xlabel(r'time ($\mu$s)')
ylabel('effective k')
xlim(0, 16100)
ylim(6, 8.8)
legend(loc = "upper left")
savefig(resdir+pfile_kerr)
show()



std_ceof_p1 = [np.std(c) for c in coef_p1_t]
std_ceof_p2 = [np.std(c) for c in coef_p2_t]
std_ceof_p3 = [np.std(c) for c in coef_p3_t]

print "std_ceof_p1 ",std_ceof_p1
print "std_ceof_p2 ",std_ceof_p2
print "std_ceof_p3 ",std_ceof_p3
    