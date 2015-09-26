import os
import sys
#import wx
try:
    import wxversion
    wxversion.select('3.0')
    import wx
except:
    import wx
    if wx.__version__[0] != '3':
        print "wx version 3.0 required, you have ",wx.__version__
import math
import numpy
import scipy
import scipy.optimize
import matplotlib
import pylab as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas


# constants  
clight = 2.99792458E8
rm0c2 = 0.938272029E0


class MyDialog(wx.Dialog):
    def __init__(
            self, parent, ID, title, pos=wx.DefaultPosition,
            size=wx.DefaultSize, style=wx.DEFAULT_FRAME_STYLE
            ):

        wx.Dialog.__init__(self, parent, ID, title, pos, size, style)
        panel = wx.Panel(self, -1)

       # input value
        if os.path.isfile('rfparam'):
           # if parameters are defined in the previous segment, use it as the starting value.
            f_r = open('rfparam', 'r')
            for line in f_r:
                pclast = line.split()[2]
                tmlast = str(float(line.split()[0])*1.e3)
                phlast = str(float(line.split()[7])*1.e0)
                if ((float(phlast) < 1.e-6) & (float(phlast) > -1.e-6)):
                    phlast = str(0.)
                vtlast = str(float(line.split()[8])*1.e6)
                if ((float(vtlast) < 1.e-6) & (float(vtlast) > -1.e-6)):
                    vtlast = str(0.)
            f_r.close()
            self.pc1in = pclast
        else:
           # for the first segment, enter default values.
            tmlast = '0.'
            phlast = '20.'
            vtlast = '4.'
            self.pc1in = '0.1441'
        self.ke1in = str((math.sqrt(rm0c2**2+float(self.pc1in)**2)-rm0c2)*1.e3)
        t1ini = tmlast
        t2ini = '30.'
        v1ini = vtlast
        v2ini = vtlast
        p1ini = phlast
        p2ini = phlast
        self.labl_pc = wx.StaticText(self, -1, "momentum [GeV/c]")
        self.labl_ke = wx.StaticText(self, -1, "k energy [MeV]")
        self.labl_tm = wx.StaticText(self, -1, "time [ms]")
        self.labl_ph = wx.StaticText(self, -1, "phis [deg]")
        self.labl_vt = wx.StaticText(self, -1, "volt [kV]")
        self.text_ctrl_pcsg1 = wx.TextCtrl(self, -1, self.pc1in, style=wx.TE_RIGHT)
        self.text_ctrl_kesg1 = wx.TextCtrl(self, -1, self.ke1in, style=wx.TE_RIGHT)
        self.text_ctrl_time1 = wx.TextCtrl(self, -1, t1ini, style=wx.TE_RIGHT)
        self.text_ctrl_phis1 = wx.TextCtrl(self, -1, p1ini, style=wx.TE_RIGHT)
        self.text_ctrl_volt1 = wx.TextCtrl(self, -1, v1ini, style=wx.TE_RIGHT)
        self.text_ctrl_time2 = wx.TextCtrl(self, -1, t2ini, style=wx.TE_RIGHT)
        self.text_ctrl_phis2 = wx.TextCtrl(self, -1, p2ini, style=wx.TE_RIGHT)
        self.text_ctrl_volt2 = wx.TextCtrl(self, -1, v2ini, style=wx.TE_RIGHT)

        button1 = wx.Button(self, wx.ID_OK)

       # layout of the small window for parameters input
        layout0 = wx.BoxSizer(wx.VERTICAL)
        layout0.Add(self.labl_pc,flag=wx.SHAPED|wx.ALIGN_CENTER)
        layout0.Add(self.text_ctrl_pcsg1,flag=wx.SHAPED|wx.ALIGN_CENTER)
        layout1 = wx.BoxSizer(wx.VERTICAL)
        layout1.Add(self.labl_ke,flag=wx.SHAPED|wx.ALIGN_CENTER)
        layout1.Add(self.text_ctrl_kesg1,flag=wx.SHAPED|wx.ALIGN_CENTER)
        layout2 = wx.BoxSizer(wx.VERTICAL)
        layout2.Add(self.labl_tm,flag=wx.SHAPED|wx.ALIGN_CENTER)
        layout2.Add(self.text_ctrl_time1,flag=wx.SHAPED|wx.ALIGN_CENTER)
        layout2.Add(self.text_ctrl_time2,flag=wx.SHAPED|wx.ALIGN_CENTER)
        layout3 = wx.BoxSizer(wx.VERTICAL)
        layout3.Add(self.labl_ph,flag=wx.SHAPED|wx.ALIGN_CENTER)
        layout3.Add(self.text_ctrl_phis1,flag=wx.SHAPED|wx.ALIGN_CENTER)
        layout3.Add(self.text_ctrl_phis2,flag=wx.SHAPED|wx.ALIGN_CENTER)
        layout4 = wx.BoxSizer(wx.VERTICAL)
        layout4.Add(self.labl_vt,flag=wx.SHAPED|wx.ALIGN_CENTER)
        layout4.Add(self.text_ctrl_volt1,flag=wx.SHAPED|wx.ALIGN_CENTER)
        layout4.Add(self.text_ctrl_volt2,flag=wx.SHAPED|wx.ALIGN_CENTER)
        layout5 = wx.BoxSizer(wx.HORIZONTAL)
        layout5.Add(layout0,flag=wx.SHAPED|wx.ALIGN_TOP)
        layout5.Add(layout2,flag=wx.SHAPED|wx.ALIGN_CENTER)
        layout5.Add(layout3,flag=wx.SHAPED|wx.ALIGN_CENTER)
        layout5.Add(layout4,flag=wx.SHAPED|wx.ALIGN_CENTER)
        layout5.Add(layout1,flag=wx.SHAPED|wx.ALIGN_TOP)
        layout6 = wx.BoxSizer(wx.VERTICAL)
        layout6.Add(layout5,flag=wx.SHAPED|wx.ALIGN_CENTER)
        layout6.Add(button1,flag=wx.SHAPED|wx.ALIGN_CENTER)

        self.SetSizer(layout6)

#       button2 = wx.Button(self, 1003, "Close Me")
#       button2.SetPosition((15, 115))
#       self.Bind(wx.EVT_BUTTON, self.OnCloseMe, button2)
#       self.Bind(wx.EVT_CLOSE, self.OnCloseWindow)
#
#   def OnCloseMe(self, event):
#       self.Close(True)
#
#   def OnCloseWindow(self, event):
#       self.Destroy()


#
# This class defines the main routine of the script
#
# 2015/02/17
# Shinji Machida
#
#class TestFrame(wx.Frame):
#   def __init__(self, parent, title):
#       wx.Frame.__init__(self, parent, title=title, size=(800,525))
#
class TestPanel(wx.Panel):
    def __init__(self, *args, **kwds):
        kwds["style"] = wx.DEFAULT_FRAME_STYLE
        wx.Panel.__init__(self, *args, **kwds)

       # split window by 1: Splitter, 2: sizer
        sw = 2
        if sw == 1:
            self.sp = wx.SplitterWindow(self)
           # right panel
            tt = [0., 1.e-3]
            vv = [0., 1.e-6]
            pp = [0., 1.]
            self.p1 = Matplot(self.sp, tt, pp, vv)
           # left panel
            self.p2 = wx.Panel(self.sp, style=wx.SUNKEN_BORDER)
            self.sp.SplitVertically(self.p2, self.p1, 400)
        else:
           # right panel
            tt = [0., 1.e-3]
            vv = [0., 1.e-6]
            pp = [0., 1.]
            self.p1 = Matplot(self, tt, pp, vv)
           # left panel
            self.p2 = wx.Panel(self, style=wx.SUNKEN_BORDER)

#      # statusbar
#       self.statusbar = self.CreateStatusBar()
#       self.statusbar.SetStatusText("Hola")

       # input field of segments
        self.labl = wx.StaticText(self.p2, -1, "# of segments:")
        self.segm = wx.TextCtrl(self.p2, -1, "1", style=wx.TE_CENTRE)

       # input order of fitting
        self.labl_f = wx.StaticText(self.p2, -1, "order of fitting (frequency):")
        self.ord_f = wx.TextCtrl(self.p2, -1, "5", style=wx.TE_CENTRE)

       # input order of fitting
        self.labl_a = wx.StaticText(self.p2, -1, "order of fitting (amplitude):")
        self.ord_a = wx.TextCtrl(self.p2, -1, "5", style=wx.TE_CENTRE)

       # button to show another panel
        b = wx.Button(self.p2, -1, "Apply")
        self.Bind(wx.EVT_BUTTON, self.OnButton, b)

       # button to quit
       #button = wx.Button(self.p2, 1003, "Make .equ file and exit", (0,100))
        button = wx.Button(self.p2, -1, "Make .equ file and exit")
        self.Bind(wx.EVT_BUTTON, self.OnCloseMe, button)
        self.Bind(wx.EVT_CLOSE, self.OnCloseWindow)

        self.k_global = 2
        rb = wx.RadioBox(
                self.p2, -1, "selection of k", wx.DefaultPosition, wx.DefaultSize,
                ['TOSCA k', 'const k'], 2, wx.RA_SPECIFY_COLS
                )
        self.Bind(wx.EVT_RADIOBOX, self.EvtRadioBox, rb)

       # layout in self.p2
        layout0 = wx.BoxSizer(wx.VERTICAL)
        layout0.Add(rb,flag=wx.SHAPED|wx.ALIGN_CENTER)
        layout0.Add(self.labl,flag=wx.SHAPED|wx.ALIGN_CENTER)
        layout0.Add(self.segm,flag=wx.SHAPED|wx.ALIGN_CENTER)
        layout0.Add(b,flag=wx.SHAPED|wx.ALIGN_CENTER)
        layout0.Add(button,flag=wx.SHAPED|wx.ALIGN_CENTER)
        layout1 = wx.BoxSizer(wx.VERTICAL)
        layout1.Add(self.labl_f,flag=wx.SHAPED|wx.ALIGN_CENTER)
        layout1.Add(self.ord_f,flag=wx.SHAPED|wx.ALIGN_CENTER)
        layout1.Add(self.labl_a,flag=wx.SHAPED|wx.ALIGN_CENTER)
        layout1.Add(self.ord_a,flag=wx.SHAPED|wx.ALIGN_CENTER)
        layout2 = wx.BoxSizer(wx.VERTICAL)
        layout2.Add(layout0,flag=wx.SHAPED|wx.ALIGN_CENTER)
        layout2.Add(layout1,flag=wx.SHAPED|wx.ALIGN_CENTER)

        self.SetSizer(layout2)
        layout2.Fit(self)
       #self.Layout()

       # layout of combined with self.p1 and self.p2
        sizer_1 = wx.BoxSizer(wx.HORIZONTAL)
      # sizer_1.Add(self.p2, flag=wx.SHAPED|wx.ALIGN_CENTER_VERTICAL)
        sizer_1.Add(self.p2, flag=wx.GROW|wx.ALIGN_CENTER_VERTICAL)
        sizer_1.Add(self.p1, flag=wx.GROW|wx.ALIGN_CENTER_VERTICAL)
        self.SetSizer(sizer_1)
        sizer_1.Fit(self)
        self.Layout()

    def EvtRadioBox(self, event):
        if event.GetInt() == 0:
            self.k_global = 2
        else:
            self.k_global = 1


    def OnButton(self, evt):
       # remove file for graphic and fitting later
        if os.path.isfile('rfparam'):
            os.remove('rfparam')

        self.nseg = int(self.segm.GetValue())
        self.nord_f = int(self.ord_f.GetValue())
        self.nord_a = int(self.ord_a.GetValue())
       # print 'number of segment', self.nseg

        self.pcsg1 = []
        self.time1 = []
        self.phis1 = []
        self.volt1 = []
        self.time2 = []
        self.phis2 = []
        self.volt2 = []
        for iseg in range(0,self.nseg):
            win = MyDialog(self, -1, "time, phis and voltage", size=(600, 120),
                      style = wx.DEFAULT_FRAME_STYLE)
            res = win.ShowModal()
           # Obtain values from another window
            if res == wx.ID_OK:
                self.pcsg1.append(float(win.text_ctrl_pcsg1.GetValue()))
                self.time1.append(float(win.text_ctrl_time1.GetValue())*1.e-3)
                self.phis1.append(float(win.text_ctrl_phis1.GetValue())*math.pi/180.)
                self.volt1.append(float(win.text_ctrl_volt1.GetValue())*1.e-6)
                self.time2.append(float(win.text_ctrl_time2.GetValue())*1.e-3)
                self.phis2.append(float(win.text_ctrl_phis2.GetValue())*math.pi/180.)
                self.volt2.append(float(win.text_ctrl_volt2.GetValue())*1.e-6)

                if self.time1[iseg] == self.time2[iseg]:
                    print 'Start and end time should be different.'
                    break

                tt = [self.time1[0], self.time2[iseg]]
                if min(self.volt1+self.volt2) > 0.:
                    vv = [0.9*min(self.volt1+self.volt2), 1.1*max(self.volt1+self.volt2)]
                else:
                    vv = [-0.1*max(self.volt1+self.volt2), 1.1*max(self.volt1+self.volt2)]
                if min(self.phis1+self.phis2) > 0.:
                    pp = [0.9*min(self.phis1+self.phis2), 1.1*max(self.phis1+self.phis2)]
                else:
                    pp = [-0.1*max(self.phis1+self.phis2), 1.1*max(self.phis1+self.phis2)]
               #print [self.time1[0], (self.volt1), (self.phis1)]
               #print [self.time2[iseg], (self.volt2), (self.phis2)]
                print 'at t0 ', self.time1[iseg], self.phis1[iseg]*180./math.pi, self.volt1[iseg]
                print 'at t1 ', self.time2[iseg], self.phis2[iseg]*180./math.pi, self.volt2[iseg]

                ph0 = self.phis1[iseg]
                dph = (self.phis2[iseg]-self.phis1[iseg])/(self.time2[iseg]-self.time1[iseg])
                vt0 = self.volt1[iseg]
                dvt = (self.volt2[iseg]-self.volt1[iseg])/(self.time2[iseg]-self.time1[iseg])
                p0_iseg = self.pcsg1[iseg]                
                t0 = self.time1[iseg]
                t1 = self.time2[iseg]
                kgl = self.k_global

               # make parameter table 'rfparam' for the segment
                rf = Varik(ph0, dph, vt0, dvt, p0_iseg, t0, t1, kgl)

#                p1_iseg = rf.getp1()
#                win.pc1in = str(p1_iseg)
#                print 'p1_iseg', p1_iseg, win.pc1in
##                rffreq = rf.rf_freq()

                self.p1 = Matplot(self, tt, pp, vv)

#                self.p1 = Matplot(self.sp, xx, yy)
#                self.p2 = wx.Panel(self.sp, style=wx.SUNKEN_BORDER)
#                self.sp.SplitVertically(self.p2, self.p1, 200)

                sizer_1 = wx.BoxSizer(wx.HORIZONTAL)
              # sizer_1.Add(self.p2, flag=wx.SHAPED|wx.ALIGN_CENTER_VERTICAL)
                sizer_1.Add(self.p2, flag=wx.GROW|wx.ALIGN_CENTER_VERTICAL)
                sizer_1.Add(self.p1, flag=wx.GROW|wx.ALIGN_CENTER_VERTICAL)
                self.SetSizer(sizer_1)
                sizer_1.Fit(self)
                self.Layout()

            win.Destroy()

        Awg_file(self.nord_f, self.nord_a)

#       button = wx.Button(win, 1003, "Close Me")
#       button.SetPosition((15, 115))
#       self.Bind(wx.EVT_BUTTON, win.OnCloseMe, button)
#       self.Bind(wx.EVT_CLOSE, win.OnCloseWindow)
#

    def OnCloseMe(self, event):
        self.Close(True)

    def OnCloseWindow(self, event):
        self.Destroy()
        sys.exit()

#       self.siButton = wx.Button(self.p2, -1, "Si", size=(40,20), pos=(10,80))
#       self.siButton.Bind(wx.EVT_BUTTON, self.Si)
#
#       self.noButton = wx.Button(self.p2, -1, "No", size=(40,20), pos=(10,60))
#       self.noButton.Bind(wx.EVT_BUTTON, self.No)
#
#   def Si(self, event):
#       time1 = self.text_ctrl_time1.GetValue()
#       volt1 = self.text_ctrl_volt1.GetValue()
#       time2 = self.text_ctrl_time2.GetValue()
#       volt2 = self.text_ctrl_volt2.GetValue()
#      # self.statusbar.SetStatusText("Si")
#       self.statusbar.SetStatusText(time1+time2)
#
#       xx = [float(time1), float(time2)]
#       yy = [float(volt1), float(volt2)]
#       self.p1 = p1(self, xx, yy)
#
#       sizer_1 = wx.BoxSizer(wx.HORIZONTAL)
#       sizer_1.Add(self.p2, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
#       sizer_1.Add(self.p1, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
#       self.SetSizer(sizer_1)
#       sizer_1.Fit(self)
#       self.Layout()
#
#   def No(self, event):
#       self.statusbar.SetStatusText("No")


class Matplot(wx.Panel):
    def __init__(self, parent, ttlim, pplim, vvlim):
#       wx.Panel.__init__(self, parent, -1, size=(600,500))
        wx.Panel.__init__(self, parent, -1)

        tt = []
        fr = []
        pc = []
        ke = []
        ba = []
        ph = []
        vt = []
        vc = []
        if os.path.isfile('rfparam'):
            f_r = open('rfparam', 'r')
            for line in f_r:
                tt.append(float(line.split()[0])*1.e3)
                fr.append(float(line.split()[1])*1.e-6)
                pc.append(float(line.split()[2])*1.e0)
                ke.append(float(line.split()[3])*1.e3)
                ba.append(float(line.split()[6])*1.e0)
                ph.append(float(line.split()[7])*1.e0)
                vt.append(float(line.split()[8])*1.e6)
                vc.append(float(line.split()[9])*1.e0)
            f_r.close()
        else:
            tt = [0., 1.]
            ttlim = [0., 0.001]
            fr = [1., 1.]
            pc = [0.1441, 0.1441]
            ke = [11., 11.]
            ba = [1., 1.]
            ph = [1., 1.]
            pplim = [0., math.pi/180.*50.]
            vt = [1., 1.]
            vvlim = [0., 1.e-6]
            vc = [1., 1.]

        self.figure = matplotlib.figure.Figure()

        self.axes = self.figure.add_subplot(221)
        self.axes.set_xlabel("time [ms]")
        self.axes.set_ylabel("phis [degree]")
        self.axes.set_xlim(x*1.e3 for x in ttlim)
        self.axes.set_ylim(x*180./math.pi for x in pplim)
        self.axes.plot(tt,ph)

        self.axes = self.figure.add_subplot(222)
        self.axes2 = self.axes.twinx()
        self.axes.set_xlabel("time [ms]")
        self.axes.set_ylabel("voltage [kv]", color='g')
        self.axes2.set_ylabel("amplitude correction", color='r')
        self.axes.set_xlim(x*1.e3 for x in ttlim)
        self.axes.set_ylim(x*1.e6 for x in vvlim)
        self.axes.plot(tt,vt,'g-')
        self.axes2.plot(tt,vc,'r-')

        self.axes = self.figure.add_subplot(223)
        self.axes.set_xlabel("time [ms]")
        self.axes.set_ylabel("bucket area [eV.s]")
        self.axes.set_xlim(x*1.e3 for x in ttlim)
#       self.axes.set_ylim(ymin=0)
        self.axes.set_ylim([0,0.3])
        self.axes.plot(tt,ba)

        self.axes = self.figure.add_subplot(224)
        self.axes2 = self.axes.twinx()
        self.axes.set_xlabel("time [ms]")
#       self.axes.set_ylabel("momentum [GeV/c]", color='g')
        self.axes.set_ylabel("kinetic energy [MeV]", color='g')
        self.axes2.set_ylabel("frequency [MHz]", color='r')
        self.axes.set_xlim(x*1.e3 for x in ttlim)
#       self.axes2.set_ylim(ymin=0)
        self.axes2.set_ylim([1.,5.])
#       self.axes.plot(tt,pc,'g-')
        self.axes.plot(tt,ke,'g-')
        self.axes2.plot(tt,fr,'r-')

        self.canvas = FigureCanvas(self, -1, self.figure)


#---------------------------------------------------------------------------

class TestPanel2(wx.Panel):
#   def __init__(self, parent, log):
    def __init__(self, parent):
#       self.log = log
        wx.Panel.__init__(self, parent, -1)

        b = wx.Button(self, -1, "Create and Show a Frame", (50,50))
        self.Bind(wx.EVT_BUTTON, self.OnButton, b)


    def OnButton(self, evt):
#       win = MyFrame(self, -1, "This is a wx.Frame", size=(350, 200),
        win = MyDialog(self, -1, "This is a wx.Frame", size=(600, 120),
                  style = wx.DEFAULT_FRAME_STYLE)
        win.Show(True)


#app = wx.App(redirect=False)
## frame = TestFrame(None, "Hello, world!")
#frame = TestFrame(None, -1, "Hello, world!")
#frame.Show()
#app.MainLoop()


class TabPanel(wx.Panel):
    """
    This will be the first notebook tab
    """
    #----------------------------------------------------------------------
    def __init__(self, parent):
 
        wx.Panel.__init__(self, parent=parent, id=wx.ID_ANY)
 
        sizer = wx.BoxSizer(wx.VERTICAL)
        txtOne = wx.TextCtrl(self, wx.ID_ANY, "")
        txtTwo = wx.TextCtrl(self, wx.ID_ANY, "")
 
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(txtOne, 0, wx.ALL, 5)
        sizer.Add(txtTwo, 0, wx.ALL, 5)
 
        self.SetSizer(sizer)


class NotebookDemo(wx.Notebook):
    """
    Notebook class
    """
 
    #----------------------------------------------------------------------
    def __init__(self, parent):
        wx.Notebook.__init__(self, parent, id=wx.ID_ANY, style=
                             wx.BK_DEFAULT
                             #wx.BK_TOP 
                             #wx.BK_BOTTOM
                             #wx.BK_LEFT
                             #wx.BK_RIGHT
                             )
 
        # Create the first tab and add it to the notebook
        tabOne = TestPanel(self)
#       tabOne.SetBackgroundColour("Gray")
        self.AddPage(tabOne, "rf parameters")
 
        # Create and add the second tab
        tabTwo = TestPanel2(self)
        self.AddPage(tabTwo, "Tab two")
 
        # Create and add the third tab
        TabThree = TestPanel2(self)
        self.AddPage(TabThree, "Tab three")
 
        self.Bind(wx.EVT_NOTEBOOK_PAGE_CHANGED, self.OnPageChanged)
        self.Bind(wx.EVT_NOTEBOOK_PAGE_CHANGING, self.OnPageChanging)
 
    def OnPageChanged(self, event):
        old = event.GetOldSelection()
        new = event.GetSelection()
        sel = self.GetSelection()
#       print 'OnPageChanged,  old:%d, new:%d, sel:%d\n' % (old, new, sel)
        event.Skip()
 
    def OnPageChanging(self, event):
        old = event.GetOldSelection()
        new = event.GetSelection()
        sel = self.GetSelection()
#       print 'OnPageChanging, old:%d, new:%d, sel:%d\n' % (old, new, sel)
        event.Skip()


class DemoFrame(wx.Frame):
    """
    Frame that holds all other widgets
    """
 
    #----------------------------------------------------------------------
    def __init__(self):
        """Constructor"""

        wx.Frame.__init__(self, None, wx.ID_ANY,
                          "KURRI FFAG",
                          size=(900,570)
                          )
        panel = wx.Panel(self)
 
        notebook = NotebookDemo(panel)
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(notebook, 1, wx.ALL|wx.EXPAND, 5)
        panel.SetSizer(sizer)
        self.Layout()
 
        self.Show()


#
# This is main physics class
#
# 2015/02/17
# Shinji Machida
#
class Varik:
    def __init__(self, ph0, dph, vt0, dvt, p0_iseg, t0, t1, kgl):
        global filenm
        filenm = 'out_scode_r2_all'

        self.ph0 = ph0
        self.dph = dph
        self.vt0 = vt0
        self.dvt = dvt
        self.p0_iseg = p0_iseg
        self.p1_iseg = self.p0_iseg
        self.t0 = t0
        self.t1 = t1
        self.kgl = kgl
        self.ck = Indxk()
        self.rf_freq()


    def rf_freq(self):
        voltcr = Volt_corr()
#        self.ph0 = ph0
#        self.dph = dph
#        self.vt0 = vt0
#        self.dvt = dvt
       # enter phis
#        phis_d_str = raw_input('synchro phase in deg. \n')
#        phis_d = float(phis_d_str)
#        phis_r = math.pi/180.0*phis_d
       #phis_r = math.pi*20./180.
        phis_r = self.ph0
       # enter selection of k
#        kselect_str = raw_input('1: k=7.645, 2: k from TOSCA smoothed, 3: k from TOSCA without smoothing \n')
#        kselect = int(kselect_str)
        kselect = self.kgl
        if kselect == 2:
            print 'TOSCA k is chosen.'
        elif kselect == 1:
            print 'const k=7.645 is chosen.'
        else:
            print 'correct k is not chosen.'
       # print 'selected k is ', kselect
       # enter voltage
#        volt_str = raw_input('voltage [kV] \n')
#        volt = float(volt_str)*1.E-6
       #volt = 4.e-6
        volt = self.vt0

       # enter p and get initial f
       #p0 = 0.1441
        p0 = self.p0_iseg
#        dpp_str = raw_input('frequency offset in term of dp/p (p0=0.1441 GeV/c). positive dp/p gives higher freq at t=0 \n')
#        dpp = float(dpp_str)
        dpp = 0.
        pinit = p0*(1. + dpp)
        f_r = open(filenm, 'r')
        self.pc = []
        self.tt = []
        for line in f_r:
            self.pc.append(float(line.split()[1]))
            self.tt.append(float(line.split()[9]))
        f_r.close()
        pcbegin = self.pc[0]
        pcend = self.pc[len(self.pc)-1-1]
       # first find frequency for p0=0.1441
        pct0 = p0
        if pct0 < pcbegin:
            print 'table does not have number below', pcbegin, '[GeV/c]', 'requested momentum is', pct0, '[GeV/c]'
            sys.exit()
        else:
            for i in range(1,len(self.pc)-1):
                if pct0 < self.pc[i]:
                    ttt0 = (self.tt[i]-self.tt[i-1])/(self.pc[i]-self.pc[i-1])*(pct0-self.pc[i-1]) + self.tt[i-1]
                    frt0 = 1./ttt0
                    break
       # print 'frequency at', pct0, 'is', frt0, '.'
       # second find frequency for p0=0.1441*()
        pct = pinit
        if pct < pcbegin:
            print 'table does not have number below', pcbegin, '[GeV/c]', 'requested momentum is', pct, '[GeV/c]'
            sys.exit()
        else:
            for i in range(1,len(self.pc)-1):
                if pct < self.pc[i]:
                    ttt = (self.tt[i]-self.tt[i-1])/(self.pc[i]-self.pc[i-1])*(pct-self.pc[i-1]) + self.tt[i-1]
                    frt = 1./ttt
                    break
       # print 'frequency at', pct, 'is', frt, '.'
        tinit = 1.E-6*ttt
       # print 'dp/p and df/f', (pct-pct0)/pct0, (frt-frt0)/frt0
       # tinit = 12*0.05264829740046983E-6
       # tinit = tinit*(1.E0 - 0.012636E0) # + 20 kHz
       # tinit = tinit*(1.E0 - 0.006318E0) # + 10 kHz
       # tinit = tinit*(1.E0 + 0.006318E0) # - 10 kHz
       # tinit = tinit*(1.E0 + 0.012636E0) # - 20 kHz

       # frequency curve directly from f vs ke
       # self.tom = Frqcy(tinit, volt, phis_r)

        pc0 = pinit
        ke0 = math.sqrt( rm0c2**2 + pc0**2 ) - rm0c2
        beta0 = (pc0/rm0c2)/((ke0+rm0c2)/rm0c2)
        tt0 = tinit
        rr_ref = beta0*clight*tt0/(2*math.pi)
        rr0 = rr_ref
        freq0 = beta0*clight/(2*math.pi*rr0)
        indexk = self.indx(kselect,pc0)
        bck0 = self.ba(beta0,ke0,volt,phis_r,rr0,indexk)
       #print 'initial bucket area', bck0*1.e9, '[eV.s]'

       # 1: BA, 2: phase, 3: volt
#        iop_str = raw_input('select condition, 1: BA, 2: phase, 3: volt \n')
#        iop = int(iop_str)
        iop = 14

        f_data_w1 = open('accfrq.dat', 'w')
        f_data_w2 = open('rfparam', 'a')
        f_data_w3 = open('accvlt.dat', 'w')
        dt = 1.E-6
#       f_w = open('kused.txt', 'w')
        for it in range(0,50000):
            tabs = it*dt
           #f_data_w1.write('%16.9e %16.9e\n' % (tabs, freq0))
            f_data_w1.write('%13.6e %13.6e\n' % (tabs, freq0))
           #f_data_w1.write(str(tabs)+' '+str(freq0)+' '+'\n')

            if iop == 1:
                dkedt = volt*math.sin(phis_r)*freq0
            elif iop == 12:
                if tabs < 2.5e-4:
                    volt = (float(volt_str)*1.E-6)/2 + tabs/2.5e-4 * (float(volt_str)*1.E-6)/2
                else:
                    volt = float(volt_str)*1.E-6
                dkedt = volt*math.sin(phis_r)*freq0
            elif iop == 13:
                if tabs < 2.5e-4:
                    volt = (float(volt_str)*1.E-6)/4 + tabs/2.5e-4 * (float(volt_str)*1.E-6)*3./4
                else:
                    volt = float(volt_str)*1.E-6
                dkedt = volt*math.sin(phis_r)*freq0
            elif iop == 14:
                volt = (self.vt0 + self.dvt*tabs)
                phis_r = (self.ph0 + self.dph*tabs)
                dkedt = volt*math.sin(phis_r)*freq0
            elif iop == 2:
                gam = 1./math.sqrt(1.-beta0**2)
                omega = beta0*clight/rr0
                eta = math.sqrt((1./(indexk+1.)-1./gam**2)**2)
                alpha = bck0/16*math.sqrt((2*math.pi*omega**2*eta)/(beta0**2*(ke0+rm0c2)*volt))
                phis_r = math.asin((1.-alpha)/(1.+alpha))
                dkedt = volt*math.sin(phis_r)*freq0
            elif iop == 3:
                gam = 1./math.sqrt(1.-beta0**2)
                omega = beta0*clight/rr0
                eta = math.sqrt((1./(indexk+1.)-1./gam**2)**2)
                alpha = (1.-math.sin(phis_r))/(1.+math.sin(phis_r))
                volt = (bck0/(16*alpha))**2*(2*math.pi*omega**2*eta)/(beta0**2*(ke0+rm0c2))
                dkedt = volt*math.sin(phis_r)*freq0
            else:
                print 'error'
           # bucket area calculation for this case
            bck0 = self.ba(beta0,ke0,volt,phis_r,rr0,indexk)

           # get voltage correction factor
            vc0 = voltcr.get_fact(freq0)
           #f_data_w3.write('%16.9e %16.9e\n' % (tabs, volt))
            f_data_w3.write('%13.6e %13.6e\n' % (tabs, volt))
           #f_data_w3.write(str(tabs)+' '+str(volt)+' '+'\n')

            f_data_w2.write('%13.6e %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e\n' % \
             (tabs+self.t0, freq0, pc0, ke0, beta0, rr0, bck0*1.e9, phis_r*180./math.pi, volt, vc0))

# update energy by dkedt, which is eVsin\phi
# translate it to momentum and beta
# calculate index k at new momentum
# calculate rr with new index k
            rr_b = rr0
            pc_b = pc0
            ke0 = ke0 + dkedt*dt
            pc0 = math.sqrt( (ke0+rm0c2)**2 - rm0c2**2 )
           # already at the top momentum
            if pc0 > pcend:
                break
            beta0 = (pc0/rm0c2)/((ke0+rm0c2)/rm0c2)
            indexk = self.indx(kselect,pc0)
            rr0 = rr_b*(pc0/pc_b)**(1/(indexk+1))
            freq0 = beta0*clight/(2*math.pi*rr0)
           # end of this segment
           # factor of 1.000000001 to make sure correct window size
            if tabs*1.000000001 >= (self.t1 - self.t0):
                print 'interval, pc0, ke0 = ', self.t1-self.t0, pc0, ke0
                self.p1_iseg = pc0
                break
#       f_w.close()
        f_data_w1.close()
        f_data_w2.close()
        f_data_w3.close()


    def getp1(self):
        return self.p1_iseg


    def indx(self,kselect,pc0):
        if kselect == 1:
            indexk = 7.645
        elif kselect == 2:
            instk = self.ck.intrp(pc0)
            indexk = instk[1]
        elif kselect == 3:
            instk = self.ck.intrp(pc0)
            indexk = instk[0]
        return indexk


    def ba(self,beta0,ke0,volt,phis_r,rr0,indexk):
       # calculate bucket area
        omega = beta0*clight/rr0
        gam = 1./math.sqrt(1.-beta0**2)
        eta = math.sqrt((1./(indexk+1.)-1./gam**2)**2)
        alpha = (1.-math.sin(phis_r))/(1.+math.sin(phis_r))
        bck0 = 16*math.sqrt(beta0**2*(ke0+rm0c2)*volt/(2*math.pi*omega**2*eta))*alpha
        return bck0


class Indxk:
    def __init__(self):
       # read table from TOSCA based tracking
        f_r = open(filenm, 'r')
        self.pc = []
        self.ke = []
        self.rr = []
        self.kk = []
        keplt = []
        pcplt = []
        rrplt = []
        for line in f_r:
            self.pc.append(float(line.split()[1]))
            self.ke.append(float(line.split()[0]))
            self.rr.append(float(line.split()[7]))
        f_r.close()
        for i in range(1,len(self.pc)-1):
            invkm1 = (self.rr[i+1]-self.rr[i-1])/(self.pc[i+1]-self.pc[i-1])*self.pc[i]/self.rr[i]
            self.kk.append(1./invkm1 - 1.)
            pcplt.append(self.pc[i])
            keplt.append(self.ke[i])
            rrplt.append(self.rr[i])

       # polynomial fit
        deg = 7 # up to pc**(deg) term
        self.p = numpy.polyfit(self.pc,self.rr,deg)
        rfit = []
        drdp = []
        self.kfit = []
        f_w = open('krawandfit', 'w')
        for i in range(1,len(self.pc)-1):
            rfitv = self.p[deg]
            drdpv = 0.
            for ideg in range(1,deg+1):
                rfitv = rfitv + self.p[deg-ideg]*self.pc[i]**ideg
                drdpv = drdpv + ideg*self.p[deg-ideg]*self.pc[i]**(ideg-1)
            rfit.append(rfitv)
            drdp.append(drdpv)
            self.kfit.append(self.rr[i]/self.pc[i]/drdpv-1.)
           # for output
            f_w.write(str(self.pc[i])+' '+str(self.ke[i])+' '+str(self.rr[i])+' '+ \
             str(self.kk[i-1])+' '+str(self.kfit[i-1])+'\n')
        f_w.close()
       #print 'index k is calculated between ', self.pc[1], '< pc <', self.pc[len(self.pc)-2], '[GeV/c]'

       # check by plot
       # p30 = numpy.poly1d(numpy.polyfit(pc,rr,deg))
##        plt.plot(pcp, rrp, 'r-', pcp, p30(pcp), 'b--')
#        plt.plot(pcplt, rrplt, 'r-', pcplt, rfit, 'b--')
#        plt.show()
##        plt.plot(kep, kk, 'b--')
##        plt.show()
#        plt.plot(pcplt, self.kk, 'r-', pcplt, self.kfit, 'b--')
#        plt.show()


    def intrp(self, pct):
       # return both raw and smoothed k value
        self.pct = pct

        if pct < self.pc[1]:
            print 'need table for smaller pc.'
            sys.exit()
        elif pct > self.pc[len(self.pc)-2]:
            print 'need table for larger pc.', pct, self.pc[len(self.pc)-2]
            sys.exit()
        else:
# p[0], p[1], p[2], ... , p[n-1], p[n]
#       k[0], k[1], ... , k[n-2]
            for i in range(2,len(self.pc)-1):
                if pct <= self.pc[i]:
                    self.instk = (self.kk[i-1]-self.kk[i-2])/(self.pc[i]-self.pc[i-1])*(self.pct-self.pc[i-1]) + self.kk[i-2]
                    self.instkfit = (self.kfit[i-1]-self.kfit[i-2])/(self.pc[i]-self.pc[i-1])*(self.pct-self.pc[i-1]) + self.kfit[i-2]
                    return self.instk, self.instkfit


class Volt_corr:
    def __init__(self):
       # check file exists
        if os.path.isfile('volt_coef'):
            f_r = open('volt_coef', 'r')
            self.p1 = []
            deg = int(f_r.readline().split()[0])
            for i in range(0,deg):
                self.p1.append(float(f_r.readline().split()[0]))
            self.fr0 = float(f_r.readline().split()[0])
            self.am0 = float(f_r.readline().split()[0])
            f_r.close()
        else:
            self.p1 = self.table()
            f_w = open('volt_coef', 'w')
            f_w.write(str(len(self.p1))+'\n')
            for deg in range(0,len(self.p1)):
                f_w.write(str(self.p1[deg])+'\n')
            f_w.write(str(self.fr0)+'\n')
            f_w.write(str(self.am0)+'\n')
            f_w.close()


    def table(self):

        deg = 5 # up to x**(deg) term

        clock=40000000

# from "20140115_svk02BS.equ"
        A0=   1582566.7294843209     
        A1=   97407274.425706580     
        A2=  -1380843783.6368339     
        A3=   32716167542.246567     
        A4=  -844170908364.90601     
        A5=  -6794756964383.6104     
        A6=   1172345739223412.2     
        A7=  -22229952186944820.     

        RDF=1.7
        D0=  0.83665564e0      *RDF
        D1= -0.13153234e3      /D0
        D2=  0.26421600e5      /D0
        D3= -0.24545073e7      /D0
        D4=  0.11180200e9      /D0
        D5= -0.20178026e10     /D0

        dT = 0.1e-3
        Ttotal = 33.2e-3
        Tacc   = 19.284309589869611e-3
        Tini   = 0.1e-3
        Tfin   = Ttotal-Tacc-Tini-dT

        size=(Tacc+dT)*clock
        size = int(size)

        timep1 = []
        voltp1 = []
        freqp1 = []
        amplp1 = []
        dtime = 1./clock
        for i in range(0,size):
            time = i*dtime
            temp_volt = math.sin(6.2831*(A0*(time-dT)+A1*(time-dT)**2+A2*(time-dT)**3+A3*(time-dT)**4+A4*(time-dT)**5+ \
             A5*(time-dT)**6+A6*(time-dT)**7+A7*(time-dT)**8))* \
             (1+D1*time+D2*time**2+D3*time**3+D4*time**4+D5*time**5)
            freq = (A0+2*A1*(time-dT)+3*A2*(time-dT)**2+4*A3*(time-dT)**3+5*A4*(time-dT)**4+ \
             6*A5*(time-dT)**5+7*A6*(time-dT)**6+8*A7*(time-dT)**7)*1.e-6
            ampl = 1+D1*time+D2*time**2+D3*time**3+D4*time**4+D5*time**5
           # pack in to array
            timep1.append(time)
            voltp1.append(temp_volt)
            freqp1.append(freq)
            amplp1.append(ampl)

# from "20140116_svk20BS.equ"
        A0=   1582564.5755691414     
        A1=   66632406.287613735     
        A2=  -646353884.88147271     
        A3=   10486198316.107014     
        A4=  -185039824855.80661     
        A5=  -1035018773031.0421     
        A6=   120905933776207.17     
        A7=  -1567779204687203.5     

        RDF=1.7
        D0=  0.82937840099633109      *RDF
        D1=  -86.421882251732328      /D0
        D2=   11669.422950503042      /D0
        D3=  -752381.57659578114      /D0
        D4=   24089167.621970635      /D0
        D5=  -305756230.44664323      /D0

        dT = 0.1e-3
        Ttotal = 33.2e-3
        Tacc   = 28.191532709927699e-3
        Tini   = 0.1e-3
        Tfin   = Ttotal-Tacc-Tini-dT

        size=(Tacc+dT)*clock
        size = int(size)

        timep2 = []
        voltp2 = []
        freqp2 = []
        amplp2 = []
        dtime2 = 1./clock
        for i in range(0,size):
            time = i*dtime
            temp_volt = math.sin(6.2831*(A0*(time-dT)+A1*(time-dT)**2+A2*(time-dT)**3+A3*(time-dT)**4+A4*(time-dT)**5+ \
             A5*(time-dT)**6+A6*(time-dT)**7+A7*(time-dT)**8))* \
             (1+D1*time+D2*time**2+D3*time**3+D4*time**4+D5*time**5)
            freq = (A0+2*A1*(time-dT)+3*A2*(time-dT)**2+4*A3*(time-dT)**3+5*A4*(time-dT)**4+ \
             6*A5*(time-dT)**5+7*A6*(time-dT)**6+8*A7*(time-dT)**7)*1.e-6
            ampl = 1+D1*time+D2*time**2+D3*time**3+D4*time**4+D5*time**5
           # pack in to array
            timep2.append(time)
            voltp2.append(temp_volt)
            freqp2.append(freq)
            amplp2.append(ampl)

       # combine two set of data together
        freqp = freqp1 + freqp2
        amplp = amplp1 + amplp2

       # polynomial fit
        self.p = numpy.polyfit(freqp,amplp,deg)
        afitp = []
        for i in range(0,len(freqp1)):
            afit = self.p[deg]
            for ideg in range(1,deg+1):
                afit = afit + self.p[deg-ideg]*freqp1[i]**ideg
            afitp.append(afit)
       # polynomial fit with constraint, x=0, y=0
        p1 = self.poly_fit_woconst(freqp,amplp,deg)
        self.fr0 = freqp[0]
        self.am0 = amplp[0]
        afitpc = []
        for i in range(0,len(freqp1)):
            afit = self.am0
            for ideg in range(1,deg+1):
                afit = afit + p1[ideg-1]*(freqp1[i]-self.fr0)**ideg
            afitpc.append(afit)
       # pack to self.p1
        self.p1 = []
        for ideg in range(0,deg):
            self.p1.append(p1[ideg])

#       plt.plot(freqp1, amplp1, 'r-', freqp2, amplp2, 'g-', freqp1, afitp, 'b--', freqp1, afitpc, 'b-')
#       plt.show()

        return self.p1


    def poly_fit_woconst(self, xin, yin, deg):
        '''Fits a linear fit of the form mx+b to the data'''
        if ((deg < 1) | (deg > 10)):
            print 'degree must be between 3 and 10'
            print 'stop'
            sys.exit()
       # subtract offset and make list numpy.array
        xin_0 = xin[0]
        yin_0 = yin[0]
        x = []
        y = []
        x[:] = [z - xin_0 for z in xin]
        y[:] = [z - yin_0 for z in yin]
        x = numpy.array(x)
        y = numpy.array(y)
       # fitting from x to x**deg
        if deg == 1:
            fitfunc = lambda params, x: params[0]*x
             #create fitting function
        elif deg == 2:
            fitfunc = lambda params, x: params[0]*x + params[1]*x**2
             #create fitting function
        elif deg == 3:
            fitfunc = lambda params, x: params[0]*x + params[1]*x**2 + params[2]*x**3
             #create fitting function
        elif deg == 4:
            fitfunc = lambda params, x: params[0]*x + params[1]*x**2 + params[2]*x**3 + params[3]*x**4
             #create fitting function
        elif deg == 5:
            fitfunc = lambda params, x: params[0]*x + params[1]*x**2 + params[2]*x**3 + params[3]*x**4 + \
             params[4]*x**5    #create fitting function
        elif deg == 6:
            fitfunc = lambda params, x: params[0]*x + params[1]*x**2 + params[2]*x**3 + params[3]*x**4 + \
             params[4]*x**5 + params[5]*x**6    #create fitting function
        elif deg == 7:
            fitfunc = lambda params, x: params[0]*x + params[1]*x**2 + params[2]*x**3 + params[3]*x**4 + \
             params[4]*x**5 + params[5]*x**6 + params[6]*x**7    #create fitting function
        elif deg == 8:
            fitfunc = lambda params, x: params[0]*x + params[1]*x**2 + params[2]*x**3 + params[3]*x**4 + \
             params[4]*x**5 + params[5]*x**6 + params[6]*x**7 + params[7]*x**8    #create fitting function
        elif deg == 9:
            fitfunc = lambda params, x: params[0]*x + params[1]*x**2 + params[2]*x**3 + params[3]*x**4 + \
             params[4]*x**5 + params[5]*x**6 + params[6]*x**7 + params[7]*x**8 + params[8]*x**9   #create fitting function
        elif deg == 10:
            fitfunc = lambda params, x: params[0]*x + params[1]*x**2 + params[2]*x**3 + params[3]*x**4 + \
             params[4]*x**5 + params[5]*x**6 + params[6]*x**7 + params[7]*x**8 + params[8]*x**9 + params[9]*x**10
        errfunc = lambda p, x, y: fitfunc(p, x) - y              #create error function for least squares fit

       # polynomial fit to find initial value
        pini = numpy.polyfit(x,y,deg)
        init_p = []
        for i in range(0,deg):
            init_p.append(pini[deg-1-i])
        init_p_numpy = numpy.array(init_p)  #bundle initial values in initial parameters

        #calculate best fitting parameters (i.e. m and b) using the error function
        p1, success = scipy.optimize.leastsq(errfunc, init_p_numpy.copy(), args = (x, y))
        f = fitfunc(p1, x)          #create a fit with those parameters

        return p1


    def get_fact(self, fr):
        self.fr = fr*1.e-6
        fact = self.am0
        for ideg in range(1,len(self.p1)+1):
            fact += self.p1[ideg-1]*(self.fr-self.fr0)**ideg
        return fact


class Awg_file:
    def __init__(self, nord_f, nord_a):
        self.nord_f = nord_f
        self.nord_a = nord_a
       # write header lines
        mf_b = self.make_file_begin()

        f_r = open('rfparam', 'r')
        icount = 0
        tm0 = -1.
        ntable = 1
        self.tmseg = []
        self.frseg = []
        self.vcseg = []
       # phase at the end of each segment
        self.pofft = []
       # phis jump at the end of each segment [deg]
        self.pjump = [0.]
        for line in f_r:
            tm1 = float(line.split()[0])
            fr1 = float(line.split()[1])
            vc1 = float(line.split()[9])
            ph1 = float(line.split()[7])
            if tm1 <= tm0:
               # write each part of segment
               #print 'awg', self.tmseg[0], self.tmseg[-1]
                pht = ph1 - ph0
                if ((pht < 1.e-6) & (pht > -1.e-6)):
                    pht = 0.
                self.pjump.append(pht)
                mf = self.make_file(ntable)
                ntable += 1
                self.tmseg = []
                self.frseg = []
                self.vcseg = []
            tm0 = tm1
            fr0 = fr1
            vc0 = vc1
            ph0 = ph1
            self.tmseg.append(tm0)
            self.frseg.append(fr0)
            self.vcseg.append(vc0)
       # write each part of segment
       #print 'awg', self.tmseg[0], self.tmseg[-1]
        self.pjump.append(0.)
        mf = self.make_file(ntable)

       # write end lines
        mf_e = self.make_file_end(ntable)


    def make_file(self, ntable):
        self.ntable = ntable
# write comment
        f_w = open('tmp_all.equ', 'a')
        f_w.write('\'***** Coefficients An for phase and Dn for voltage *****\n')
        f_w.write('Tacc'+str(ntable)+' = '+str(self.tmseg[-1]-self.tmseg[0])+'\n')
        if ntable == 1:
            f_w.write('dT = 0.5e-3\n')
        else:
            f_w.write('dT = 0.\n')
        f_w.write('size = (Tacc'+str(ntable)+'+dT)*clock'+'\n')
        f_w.close()

# frequency curve
       #deg = 5
        deg = self.nord_f
        self.p1 = self.poly_fit_woconst(self.tmseg, self.frseg, deg)
       # print self.p1
        tg = []
        frr = []
        frf = []
        for i in range(0,len(self.tmseg)):
            tg.append(self.tmseg[i])
            frr.append(self.frseg[i])
            frt = self.frseg[0]
            for ideg in range(1,deg+1):
                frt += self.p1[ideg-1]*(self.tmseg[i]-self.tmseg[0])**ideg
            frf.append(frt)
        frd = []
        for i in range(0,len(tg)):
            frd.append((frf[i]-frr[i])/frr[i])
#        plt.plot(tg, frd, 'g-')
#        plt.show()

       # convert from frequency to phase
        b1 = []
        b1.append(self.frseg[0])
        for ideg in range(1,deg+1):
            b1.append(self.p1[ideg-1]/(ideg+1))
       #print b1
        f_w = open('tmp_all.equ', 'a')
        for ideg in range(0,deg+1):
            f_w.write('A'+str(ideg)+' = '+str(b1[ideg])+'\n')
        f_w.close()        

        phase = 'A0*(time-dT)'
        for ideg in range(1,ideg+1):
            phase += '+A'+str(ideg)+'*(time-dT)^'+str(ideg+1)
        phase += '+phoff+pjump'
        phase = '2*'+str(math.pi)+'*('+phase+')'

       # phase at the end
        poff = 0.
        for ideg in range(1,ideg+1+1):
            poff += b1[ideg-1]*(self.tmseg[-1]-self.tmseg[0])**ideg
        self.pofft.append(poff)

# voltage amplitude correction
       #deg = 5
        deg = self.nord_a
        self.p2 = self.poly_fit_woconst(self.tmseg, self.vcseg, deg)
       # print self.p2
        tg = []
        vcr = []
        vcf = []
        for i in range(0,len(self.tmseg)):
            tg.append(self.tmseg[i])
            vcr.append(self.vcseg[i])
            vct = self.vcseg[0]
            for ideg in range(1,deg+1):
                vct += self.p2[ideg-1]*(self.tmseg[i]-self.tmseg[0])**ideg
            vcf.append(vct)
        vcd = []
        for i in range(0,len(tg)):
            vcd.append((vcf[i]-vcr[i])/vcr[i])
#        plt.plot(tg, vcd, 'g-')
#        plt.show()

       # repack coefficient to c1
        c1 = []
        c1.append(self.vcseg[0])
        for ideg in range(1,deg+1):
            c1.append(self.p2[ideg-1])
       #print c1
        f_w = open('tmp_all.equ', 'a')
        for ideg in range(0,deg+1):
            f_w.write('D'+str(ideg)+' = '+str(c1[ideg])+'\n')

        vfact = 'D0'
        for ideg in range(1,ideg+1):
            vfact += '+D'+str(ideg)+'*(time-dT)^'+str(ideg)
        vfact = '*('+vfact+')'

# write segment
       # phase offset from the previous segments
        phoffsum = 0.
        for itable in range(0,ntable-1):
            phoffsum += self.pofft[itable]
        f_w.write('phoff = '+str(phoffsum)+'\n')
        f_w.write('pjump = '+str(self.pjump[ntable-1])+'/360.\n')
        f_w.write('\"tmp_acc'+str(ntable)+'.wfm\"=sin('+phase+')'+vfact+'\n')
        f_w.write('\n')
        f_w.close()        


    def make_file_begin(self):
        if os.path.isfile('tmp_all.equ'):
            os.remove('tmp_all.equ')

        f_w = open('tmp_all.equ', 'w')

        f_w.write('\'**********          **********\n')
        f_w.write('\'    tmp_all.equ\n')
        f_w.write('\'\n')
        f_w.write('\'\n')
        f_w.write('\'    by Python script\n')
        f_w.write('\'    Shinji Machida\n')
        line = '\'**********          **********\n'
        f_w.write('\'**********          **********\n')
        f_w.write('\n')

        f_w.write('clock = 40000000 \' = 40 M\n')
        f_w.write('\n')

        f_w.write('Tini = 0.1e-3\n')
        f_w.write('Tfin = 0.1e-3\n')
        f_w.write('\n')

        f_w.write('\'***** Initial blank pattern *****\n')
        f_w.write('size = Tini*clock\n')
        f_w.write('\"tmp_ini.wfm\"=0.\n')
        f_w.write('\n')

        f_w.write('\'***** Final blank pattern *****\n')
        f_w.write('size = Tfin*clock\n')
        f_w.write('\"tmp_fin.wfm\"=0.\n')
        f_w.write('\n')

        f_w.close()


    def make_file_end(self, ntable):
        self.ntable = ntable
        f_w = open('tmp_all.equ', 'a')

        f_w.write('\'***** Concatinate temp file *****\n')
        f_w.write('\"tmp1.wfm\"=join(\"tmp_ini.wfm\",\"tmp_acc1.wfm\")\n')
        for itable in range(2,ntable+1):
            f_w.write('\"tmp'+str(itable)+'.wfm\"=join(\"tmp'+str(itable-1)+'.wfm\",\"tmp_acc'+str(itable)+'.wfm\")\n')
        f_w.write('\"tmp_all.wfm\"=join(\"tmp'+str(ntable)+'.wfm\",\"tmp_fin.wfm\")\n')
        f_w.write('\n')

        f_w.write('\'***** Delete temp file *****\n')
        f_w.write('delete(\"tmp_ini.wfm\")\n')
        f_w.write('delete(\"tmp_fin.wfm\")\n')
        for itable in range(1,ntable+1):
            f_w.write('delete(\"tmp_acc'+str(itable)+'.wfm\")\n')
            f_w.write('delete(\"tmp'+str(itable)+'.wfm\")\n')

        f_w.close()


    def poly_fit_woconst(self, xin, yin, deg):
        '''Fits a linear fit of the form mx+b to the data'''
        if ((deg < 1) | (deg > 10)):
            print 'degree must be between 3 and 10'
            print 'stop'
            sys.exit()
       # subtract offset and make list numpy.array
        xin_0 = xin[0]
        yin_0 = yin[0]
        x = []
        y = []
        x[:] = [z - xin_0 for z in xin]
        y[:] = [z - yin_0 for z in yin]
        x = numpy.array(x)
        y = numpy.array(y)
       # fitting from x to x**deg
        if deg == 1:
            fitfunc = lambda params, x: params[0]*x
             #create fitting function
        elif deg == 2:
            fitfunc = lambda params, x: params[0]*x + params[1]*x**2
             #create fitting function
        elif deg == 3:
            fitfunc = lambda params, x: params[0]*x + params[1]*x**2 + params[2]*x**3
             #create fitting function
        elif deg == 4:
            fitfunc = lambda params, x: params[0]*x + params[1]*x**2 + params[2]*x**3 + params[3]*x**4
             #create fitting function
        elif deg == 5:
            fitfunc = lambda params, x: params[0]*x + params[1]*x**2 + params[2]*x**3 + params[3]*x**4 + \
             params[4]*x**5    #create fitting function
        elif deg == 6:
            fitfunc = lambda params, x: params[0]*x + params[1]*x**2 + params[2]*x**3 + params[3]*x**4 + \
             params[4]*x**5 + params[5]*x**6    #create fitting function
        elif deg == 7:
            fitfunc = lambda params, x: params[0]*x + params[1]*x**2 + params[2]*x**3 + params[3]*x**4 + \
             params[4]*x**5 + params[5]*x**6 + params[6]*x**7    #create fitting function
        elif deg == 8:
            fitfunc = lambda params, x: params[0]*x + params[1]*x**2 + params[2]*x**3 + params[3]*x**4 + \
             params[4]*x**5 + params[5]*x**6 + params[6]*x**7 + params[7]*x**8    #create fitting function
        elif deg == 9:
            fitfunc = lambda params, x: params[0]*x + params[1]*x**2 + params[2]*x**3 + params[3]*x**4 + \
             params[4]*x**5 + params[5]*x**6 + params[6]*x**7 + params[7]*x**8 + params[8]*x**9   #create fitting function
        elif deg == 10:
            fitfunc = lambda params, x: params[0]*x + params[1]*x**2 + params[2]*x**3 + params[3]*x**4 + \
             params[4]*x**5 + params[5]*x**6 + params[6]*x**7 + params[7]*x**8 + params[8]*x**9 + params[9]*x**10
        errfunc = lambda p, x, y: fitfunc(p, x) - y              #create error function for least squares fit

       # polynomial fit to find initial value
        pini = numpy.polyfit(x,y,deg)
        init_p = []
        for i in range(0,deg):
            init_p.append(pini[deg-1-i])
        init_p_numpy = numpy.array(init_p)  #bundle initial values in initial parameters

        #calculate best fitting parameters (i.e. m and b) using the error function
        p1, success = scipy.optimize.leastsq(errfunc, init_p_numpy.copy(), args = (x, y))
        f = fitfunc(p1, x)          #create a fit with those parameters

        return p1


if __name__ == "__main__":
    app = wx.App(redirect=False)
#    frame = TestFrame(None, -1, "")
#    frame.Show()
    frame = DemoFrame()
    app.MainLoop()
