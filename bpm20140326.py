title = 'Read scope trace'

import sys
sys.path[:0] = ['../../..']
import os

import numpy
import scipy
import scipy.optimize

#import Tkinter
import Pmw
import time
#import shutil
import string
import pylab as plt
import math

#from Tkinter import *


#
#
#
class Bpm:
    def __init__(self):
        global dirname, filenum, file_c1, file_c2, file_c3, file_c4

        dirname = '/Users/machida/Documents/FFAGOthers/KURRI/KURRI2014data/'
        filenum = '2014-03-26/TUNE_Meas/'
        curr = '-15.0A_D1000A'
        file_c1 = dirname + filenum + 'F_773.3A_Foil_86mm_ST6_' + curr + '_S11_1st_ch3.csv'
        file_c2 = dirname + filenum + 'F_773.3A_Foil_86mm_ST6_' + curr + '_S7up_ch1.csv'
        file_c3 = dirname + filenum + 'F_773.3A_Foil_86mm_ST6_' + curr + '_S7dw_ch2.csv'
       # filenum = '2014-03-24/ST5V_D950/'
       # curr = '-2.93'
       # file_c1 = dirname + filenum + 'D950_ST5_' + curr + '_ch1.csv'
       # file_c2 = dirname + filenum + 'D950_ST5_' + curr + '_ch2.csv'
       # file_c3 = dirname + filenum + 'D950_ST5_' + curr + '_ch3.csv'
       # filenum = '2014-03-24/dmag/'
       # curr = '1200'
       # file_c1 = dirname + filenum + 'D' + curr + '_ch1.csv'
       # file_c2 = dirname + filenum + 'D' + curr + '_ch2.csv'
       # file_c3 = dirname + filenum + 'D' + curr + '_ch3.csv'

        self.ch23out()

       # vertical bpm
        self.bpm = Vbpm()

       # read data
        self.bpmsig = self.bpm.possig()
       # bpmsig = self.bpm.difsig()
       # bpmsig = self.bpm.difsig_wired()
        print 'data length =', len(self.bpmsig)

        gx = []
        gy = []
        for i in range(len(self.bpmsig)):
            gx.append(self.bpmsig[i][0])
            gy.append(self.bpmsig[i][1])

        plt.xlabel('time (us)')
        plt.ylabel('bunch signal and peak (arb.)')
        plt.xlim(2., 40.)
        plt.ylim(-0.22,0.22)
        plt.plot(gx, gy, 'r-')
        plt.show()

# Ch2 out
    def ch2out(self):
        self.ch2 = Ch2()
        c2 = self.ch2.scope_ch
        self.gx2 = []
        self.gy2 = []
        for i in range(len(c2)):
            self.gx2.append(c2[i][0])
            self.gy2.append(c2[i][1])
        plt.xlabel('time (us)')
        plt.ylabel('bunch signal and peak (arb.)')
        plt.xlim(0., 40.)
        plt.plot(self.gx2, self.gy2, 'r-')
        plt.show()

# Ch3 out
    def ch3out(self):
        self.ch3 = Ch3()
        c3 = self.ch3.scope_ch
        self.gx3 = []
        self.gy3 = []
        cfac = 1.8
        for i in range(len(c3)):
            self.gx3.append(c3[i][0])
            self.gy3.append(cfac*c3[i][1])
        plt.xlabel('time (us)')
        plt.ylabel('bunch signal and peak (arb.)')
        plt.xlim(0., 40.)
        plt.plot(self.gx3, self.gy3, 'b-')
        plt.show()

# CH2+Ch3 out
    def ch23out(self):
        self.ch2out()
        self.ch3out()
        plt.xlabel('time (us)')
        plt.ylabel('bunch signal and peak (arb.)')
        plt.xlim(0., 40.)
        plt.plot(self.gx2, self.gy2, 'r-')
        plt.plot(self.gx3, self.gy3, 'b-')
        plt.show()

#
# NAFF tune
#
class Tune:
    def __init__(self):
        self.bpmraw = Bpm()
        self.naff()

    def naff(self):
        def avcal(dpt):
            n = 0
            av = 0
            for i in range(0,len(dpt)):
                if dpt[i] < 1.E9:
                    n += 1
                    av += dpt[i]
            if n > 0:
                av = av/n
            else:
                av = 1.E9

            for i in range(0,len(dpt)):
                if dpt[i] < 1.E9:
                    dpt[i] = dpt[i] - av
            return av

        def resal(q):
            twopi = 2*math.pi
            cs = 0
            sn = 0
            for i in range(0,len(dpt)):
                if dpt[i] < 1.E9:
                   # no filter
                   #cs += dpt[i]*math.cos(twopi*q*i)
                   #sn += dpt[i]*math.sin(twopi*q*i)
                   # Hanning filter
                   # cs += dpt[i]*math.cos(twopi*q*ind[i]) * 2*(math.sin(twopi*ind[i]/(2*len(dpt))))**2
                   # sn += dpt[i]*math.sin(twopi*q*ind[i]) * 2*(math.sin(twopi*ind[i]/(2*len(dpt))))**2
                    cs += dpt[i]*math.cos(twopi*q*ind[i]) * 2*(math.sin(twopi*ind[i]/(2*self.nw)))**2
                    sn += dpt[i]*math.sin(twopi*q*ind[i]) * 2*(math.sin(twopi*ind[i]/(2*self.nw)))**2
           #esm = 1./math.sqrt(math.sqrt(cs*cs + sn*sn))
            if cs*cs+sn*sn !=0.:
                esm = 1./((cs*cs + sn*sn))
            else:
                esm = 1.
            return esm

        ind = []
        dpt = []
        self.nw = 30
        for i in range(self.nw):
            ind.append(i)
            dpt.append(self.bpmraw.bpmsig[i][1])
            av = avcal(dpt)
        self.htune = scipy.optimize.fminbound(resal, 0.001, 1.0, xtol = 1.e-05, disp=3)
        if self.htune > 0.5:
            self.htune = 1. - self.htune
        print 'tune = ', self.htune

#
#  vertical BPM
#
class Vbpm:
    def __init__(self):
        self.s12 = Ch1()
        self.s7up = Ch2()
        self.s7dw = Ch3()
        self.idelay = 0
        self.cfac = 1.8

    def difsig(self):
        self.pkup = self.s7up.peak()
        self.pkdw = self.s7dw.peak()
        len_up = len(self.pkup)
        len_dw = len(self.pkdw)
        print len_up, len_dw
        if len_up < len_dw:
            len_s7 = len_up
        else:
            len_s7 = len_dw
        self.dif = []
        for i in range(len_s7):
            self.dif.append( ( self.pkup[i][0], self.pkup[i][1]-self.cfac*self.pkdw[i][1] ) )
        return self.dif

    def sumsig(self):
        self.pkup = self.s7up.peak()
        self.pkdw = self.s7dw.peak()
        len_up = len(self.pkup)
        len_dw = len(self.pkdw)
        if len_up < len_dw:
            len_s7 = len_up
        else:
            len_s7 = len_dw
        self.sum = []
        for i in range(len_s7):
            self.sum.append( ( self.pkup[i][0], self.pkup[i][1]+self.cfac*self.pkdw[i][1] ) )
        return self.sum

    def possig(self):
        self.difsig()
        self.sumsig()
        self.pkup = self.s7up.peak()
        self.pkdw = self.s7dw.peak()
        len_up = len(self.pkup)
        len_dw = len(self.pkdw)
        if len_up < len_dw:
            len_s7 = len_up - self.idelay
        else:
            len_s7 = len_dw - self.idelay
        self.pos = []
        for i in range(len_s7):
            if self.sum[i][1] != 0.:
                self.pos.append( ( self.pkup[i][0], self.dif[i][1]/self.sum[i][1] ) )
            else:
                self.pos.append( ( self.pkup[i][0], 1000000. ) )
        f_data_w = open('pksig', 'w')
        for i in range(len(self.pos)):
            f_data_w.write('%13.6e %13.6e\n' % (float(self.pos[i][0]), float(self.pos[i][1])))
        f_data_w.close()
        return self.pos

    def difsig_wired(self):
        len_up = len(self.s7up.scope_ch)
        len_dw = len(self.s7dw.scope_ch)
        if len_up < len_dw:
            len_s7 = len_up
        else:
            len_s7 = len_dw
        self.dif_wired = []
        for i in range(len_s7):
            self.dif_wired.append( ( self.s7up.scope_ch[i][0], self.s7up.scope_ch[i][1] ) )
        return self.dif_wired

    def peak(self):
        nw = 0
        fst = 0
        pmin = 1.E9
        pmax = -1.E9
        self.peakvalue = []
        for i in range(len(self.dif)):
            if (self.dif[i][0] > 2.86 + 0.64*nw) & (self.dif[i][0] < 2.86 + 0.64*(nw+1) - 0.1):
                if fst == 0:
                    basel = self.dif[i][1]
                    fst = 1
                if pmin > self.dif[i][1]:
                    pmin = self.dif[i][1]
                    pmin0 = self.dif[i][0]
                if pmax < self.dif[i][1]:
                    pmax = self.dif[i][1]
                    pmax0 = self.dif[i][0]
            if (self.dif[i][0] > 2.86 + 0.64*nw) & (self.dif[i][0] > 2.86 + 0.64*(nw+1) - 0.1):
                baser = self.dif[i][1]
                pmin = pmin - (basel+baser)/2
                pmax = pmax - (basel+baser)/2
                if abs(pmin) > abs(pmax):
                    self.peakvalue.append( (pmin0, pmin) )
                else:
                    self.peakvalue.append( (pmax0, pmax) )
                nw += 1
                fst = 0
                pmin = 1.E9
                pmax = -1.E9
        return self.peakvalue
#
#  read trace ch1 data
#
class Ch1:
    def __init__(self):
        print 'reading ', file_c1, '.'
        f_ch1_data = open(file_c1, 'r')
        n = 0
        self.scope_ch = []
        for line in f_ch1_data:
            timeandamp = line.split(',')
            n += 1
            if ((n > 18) & (n <=100000)):
                self.scope_ch.append( ( 1.E6*float(timeandamp[0]), float(timeandamp[1]) ) )
        f_ch1_data.close()
        # check
        f_data_w = open('check_ch1', 'w')
        for i in range(1,len(self.scope_ch)):
            f_data_w.write('%13.6e %13.6e\n' % (self.scope_ch[i][0], self.scope_ch[i][1]))
        f_data_w.close()

    def average(self):
        amp_av = 0.
        for i in range(0,len(self.scope_ch)):
            amp_av += self.scope_ch[i][1]
        amp_av = amp_av/len(self.scope_ch)
        print 'average of ch1 ', amp_av
        return amp_av
#
#  read trace ch2 data
#
class Ch2:
    def __init__(self):
        print 'reading ', file_c2, '.'
        f_ch2_data = open(file_c2, 'r')
        n = 0
        self.scope_ch = []
        for line in f_ch2_data:
            timeandamp = line.split(',')
            n += 1
            if ((n > 18) & (n <=100000)):
                self.scope_ch.append( ( 1.E6*float(timeandamp[0]), float(timeandamp[1]) ) )
        f_ch2_data.close()
        # check
        f_data_w = open('check_ch2', 'w')
        for i in range(1,len(self.scope_ch)):
            f_data_w.write('%13.6e %13.6e\n' % (self.scope_ch[i][0], self.scope_ch[i][1]))
        f_data_w.close()

    def average(self):
        amp_av = 0.
        for i in range(0,len(self.scope_ch)):
            amp_av += self.scope_ch[i][1]
        amp_av = amp_av/len(self.scope_ch)
        print 'average of ch2 ', amp_av
        return amp_av

    def peak(self):
        nw = 0
        fst = 0
        pmin = 1.E9
        pmax = -1.E9
        self.peakvalue = []
        for i in range(len(self.scope_ch)):
            if (self.scope_ch[i][0] > 2.86 + 0.64*nw) & (self.scope_ch[i][0] < 2.86 + 0.64*(nw+1) - 0.1):
                if fst == 0:
                    basel = self.scope_ch[i][1]
                    fst = 1
                if pmin > self.scope_ch[i][1]:
                    pmin = self.scope_ch[i][1]
                    pmin0 = self.scope_ch[i][0]
                if pmax < self.scope_ch[i][1]:
                    pmax = self.scope_ch[i][1]
                    pmax0 = self.scope_ch[i][0]
            if (self.scope_ch[i][0] > 2.86 + 0.64*nw) & (self.scope_ch[i][0] > 2.86 + 0.64*(nw+1) - 0.1):
                baser = self.scope_ch[i][1]
                pmin = pmin - (basel+baser)/2
                pmax = pmax - (basel+baser)/2
                if abs(pmin) > abs(pmax):
                    self.peakvalue.append( (pmin0, pmin) )
                else:
                    self.peakvalue.append( (pmax0, pmax) )
                nw += 1
                fst = 0
                pmin = 1.E9
                pmax = -1.E9
        f_data_w = open('pksig_ch2', 'w')
        for i in range(len(self.peakvalue)):
            f_data_w.write('%13.6e %13.6e\n' % (float(self.peakvalue[i][0]), float(self.peakvalue[i][1])))
        f_data_w.close()
        return self.peakvalue
#
#  read trace ch3 data
#
class Ch3:
    def __init__(self):
        print 'reading ', file_c3, '.'
        f_ch3_data = open(file_c3, 'r')
        n = 0
        self.scope_ch = []
        for line in f_ch3_data:
            timeandamp = line.split(',')
            n += 1
            if ((n > 18) & (n <=100000)):
                self.scope_ch.append( ( 1.E6*float(timeandamp[0]), float(timeandamp[1]) ) )
        f_ch3_data.close()
        # check
        f_data_w = open('check_ch3', 'w')
        for i in range(1,len(self.scope_ch)):
            f_data_w.write('%13.6e %13.6e\n' % (self.scope_ch[i][0], self.scope_ch[i][1]))
        f_data_w.close()

    def average(self):
        amp_av = 0.
        for i in range(0,len(self.scope_ch)):
            amp_av += self.scope_ch[i][1]
        amp_av = amp_av/len(self.scope_ch)
        print 'average of ch3 ', amp_av
        return amp_av

    def peak(self):
        nw = 0
        fst = 0
        pmin = 1.E9
        pmax = -1.E9
        self.peakvalue = []
        for i in range(len(self.scope_ch)):
            if (self.scope_ch[i][0] > 2.86 + 0.64*nw) & (self.scope_ch[i][0] < 2.86 + 0.64*(nw+1) - 0.1):
                if fst == 0:
                    basel = self.scope_ch[i][1]
                    fst = 1
                if pmin > self.scope_ch[i][1]:
                    pmin = self.scope_ch[i][1]
                    pmin0 = self.scope_ch[i][0]
                if pmax < self.scope_ch[i][1]:
                    pmax = self.scope_ch[i][1]
                    pmax0 = self.scope_ch[i][0]
            if (self.scope_ch[i][0] > 2.86 + 0.64*nw) & (self.scope_ch[i][0] > 2.86 + 0.64*(nw+1) - 0.1):
                baser = self.scope_ch[i][1]
                pmin = pmin - (basel+baser)/2
                pmax = pmax - (basel+baser)/2
                if abs(pmin) > abs(pmax):
                    self.peakvalue.append( (pmin0, pmin) )
                else:
                    self.peakvalue.append( (pmax0, pmax) )
                nw += 1
                fst = 0
                pmin = 1.E9
                pmax = -1.E9
        f_data_w = open('pksig_ch3', 'w')
        for i in range(len(self.peakvalue)):
            f_data_w.write('%13.6e %13.6e\n' % (float(self.peakvalue[i][0]), float(self.peakvalue[i][1])))
        f_data_w.close()
        return self.peakvalue



#
#  Zero crossing points
#    BPM signla in Ch2.
#
class Scope:
    def __init__(self):
        self.ch = Ch2()

    def rawdata(self):
        return self.ch.scope_ch


#
#  Top level
#
if __name__ == '__main__':
    dialog = None

    widget = Tune()
