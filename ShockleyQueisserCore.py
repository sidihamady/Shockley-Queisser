#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ======================================================================================================
# Solar Cell Shockley-Queisser Limit Calculator
# Code written by:
#   Pr. Sidi Hamady
#   Université de Lorraine, France
#   sidi.hamady@univ-lorraine.fr
# See Copyright Notice in COPYRIGHT
# HowTo in README.md and README.pdf
# https://github.com/sidihamady/Shockley-Queisser
# http://www.hamady.org/photovoltaics/ShockleyQueisser.zip
# ======================================================================================================

# ShockleyQueisserCore.py
#   the class ShockleyQueisserCore implements the program core functionality
#   only the constructor and the calculate function are to be called from outside the class
#   example (to put in a test.py file, for instance):
#
#       #!/usr/bin/env python
#       #-*- coding: utf-8 -*-
#
#       from ShockleyQueisserCore import *
#
#       SCC = ShockleyQueisserCore(verbose = False)
#
#       SCC.calculate(
#           TargetBandgap           = 1.1,
#           TargetBandgapTop        = 0.0,
#           Temperature             = 300.0,
#           SolarConcentration      = 1.0,
#           OutputFilename          = './ShockleyQueisserOutput',
#           useGUI                  = True
#       )
#

# import as usual
import math
import numpy as np
import scipy as sp
import sys, os, time
import threading

TkFound = False
TkRet   = ''

# try to load the tkinter and matplotlib modules
# should be always installed in any Linux distribution
# (for Windows, just use some ready-to-use packages such as anaconda (https://www.anaconda.com/distribution/))
try:

    import matplotlib
    matplotlib.use('TkAgg')
    import matplotlib.pyplot as pl
    from matplotlib.backends.backend_pdf import PdfPages
    import matplotlib.backends.backend_tkagg
    from matplotlib.font_manager import FontProperties
    if sys.version_info[0] < 3:
        # Python 2.7.x
        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
        from matplotlib.backends.backend_tkagg import NavigationToolbar2Tk as NavigationToolbar2TkAgg
        import Tkinter as Tk
        import ttk
        import tkFileDialog
        import tkFont
        import tkMessageBox
    else:
        # Python 3.x
        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
        from matplotlib.backends.backend_tkagg import NavigationToolbar2Tk as NavigationToolbar2TkAgg
        import tkinter as Tk
        import tkinter.ttk as ttk
        import tkinter.filedialog as tkFileDialog
        import tkinter.font as tkFont
        import tkinter.messagebox as tkMessageBox
    # end if

    class NavigationToolbar(NavigationToolbar2TkAgg):
        """ custom Tk toolbar """
        def __init__(self, chart):
            NavigationToolbar2TkAgg.__init__(self, chart.canvas, chart.root)
            self.chart = chart
        # end __init__
        try:
            toolitems = [tt for tt in NavigationToolbar2TkAgg.toolitems if tt[0] in ('Home', 'Zoom')]
            toolitems.append(('AutoScale', 'Auto scale the plot', 'hand', 'onAutoScale'))
            toolitems.append(('Save', 'Save the plot', 'filesave', 'onSave'))
        except:
            pass
        # end try
        def onAutoScale(self):
            self.chart.onAutoScale()
        # end onAutoScale
        def onSave(self):
            self.chart.onSave()
        # end onSave
    # end NavigationToolbar

    TkFound = True

except ImportError as ierr:
    # if Tkinter is not found, just install or update python/numpy/scipy/matplotlib/tk modules
    TkRet = "\n! cannot load Tkinter:\n  " + ("{0}".format(ierr)) + "\n"
    pass
except Exception as excT:
    TkRet = "\n! cannot load Tkinter:\n  %s\n" % str(excT)
    pass
# end try

# suppress a nonrelevant warning from matplotlib
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

# calculations done in a secondary thread, not on UI
class CalculationThread(threading.Thread):
    def __init__(self, id, func):
        threading.Thread.__init__(self)
        self.id     = id
        self.func   = func
    # end __init__
    def run(self):
        self.func()
    # end run
# end CalculationThread

# the core class
class ShockleyQueisserCore(object):
    """ the Shockley-Queisser calculator core class """

    def __init__(self, verbose = True):
        """ the Shockley-Queisser calculator class constructor """

        self.name               = "Solar Cell Shockley-Queisser Limit Calculator"
        self.__version__        = "Version 1.0 Build 1811"

        # Basic constants
        self.pi                 = 3.14159265358979      # 
        self.q                  = 1.602176e-19          # elementary charge
        self.h                  = 6.626068e-34          # Planck constant
        self.c                  = 2.99792458e+8         # light speed in vacuum
        self.hc                 = 1.986445213e-25       # h c
        self.nmeV               = 1239.84207            # nm to eV 

        # AM1.5 solar spectrum file to put here (usually ASTM AM1.5 G-173)
        # Solar spectrum file name (SolarSpectrum_AM15G.txt included):
        # two columns (wavelength in nm  and irrandiance in W/m2/nm):
        # 280.0	4.7309E-23
        # 280.5	1.2307E-21
        # 281.0	5.6895E-21
        # 281.5	1.5662E-19
        # ...
        # 4000.0	7.1043E-03
        # the first two rows are skipped
        self.SolarSpectrumAMX   = './SolarSpectrum_AM15G.txt'
        self.WavelengthMin      = 300.0         # nm
        self.WavelengthMax      = 4000.0        # nm
        self.BandgapMin         = self.nmeV / self.WavelengthMax
        self.BandgapMax         = self.nmeV / self.WavelengthMin
        # the bandgap range in eV
        self.BandgapRange       = None
        # the used data delimiter (usually TAB) in the AM1.5 ASCII file
        self.DataDelimiter      = '\t'

        self.Bandgap            = None
        self.Efficiency         = None

        # Shockley-Queisser calculated parameters
        self.SQ_Efficiency      = None
        self.SQ_JSC             = None
        self.SQ_VOC             = None
        self.SQ_FF              = None
        self.SQ_Vm              = None
        self.SQ_Jm              = None
        self.SQ_Bandgap         = None
        self.SQ_Voltage         = None
        self.SQ_Current         = None
        self.SQ_Done            = False
        #

        # Target calculated parameters
        self.Target_Efficiency  = None
        self.Target_VOC         = None
        self.Target_JSC         = None
        self.Target_FF          = None
        self.Target_Vm          = None
        self.Target_Jm          = None
        self.Target_Voltage     = None
        self.Target_Current     = None
        #

        self.running            = False
        self.threadfinish       = None
        self.thread             = None
        self.actionbutton       = None
        self.timerduration      = 200       # in milliseconds

        self.report             = None

        self.useGUI             = True
        self.root               = None
        self.GUIstarted         = False
        self.PlotInitialized    = False

        self.SpectrumLoaded     = False

        self.tic                = 0.0
        self.Counter            = 0
        self.CounterMax         = 0

        # one can set verbose to False to disable printing output
        self.verbose            = verbose
        if not self.verbose:
            print("\nverbose set to False: printing output disabled")
        # end if

        return

    # end __init__

    def calculate(self, 
        TargetBandgap           = 1.1,
        TargetBandgapTop        = 0.0,
        Temperature             = 300.0,
        SolarConcentration      = 1.0,
        OutputFilename          = './ShockleyQueisserOutput',
        useGUI                  = True):
        """ the Shockley-Queisser calculator main function """

        # recalculate the Shockley-Queisser curve if changing temperature or solar concentration
        if self.SQ_Done and ((Temperature != self.Temperature) or (SolarConcentration != self.SolarConcentration)):
            self.SQ_Done = False
        # end if

        # useGUI: the calculator can be used in graphical (GUI) mode or command-line only mode. 
        #   in command-line mode (useGUI = False) the results are printed out and saved in text files.
        #   the command-line mode is useful to perform specific calculations such as multijunction solar cell efficiency
        if useGUI and (not TkFound):
            # if Tkinter is not found, just install or update python/numpy/scipy/matplotlib/tk modules
            print(TkRet)
        # end if
        self.useGUI = useGUI if TkFound else False

        # Temperature: in Kelvin (from 100 K to 700 K)
        self.Temperature            = Temperature if ((Temperature >= 100.0) and (TargetBandgap <= 700.0)) else 300.0
        self.kTeV                   = 0.02585202874091 * self.Temperature / 300.0     # in eV
        self.kTJ                    = self.kTeV * self.q                              # in J

        # TargetBandgap: Target bandgap in eV (from 0.2 eV to 6 eV): to compare with the Shockley-Queisser Limit
        self.Target_Bandgap         = TargetBandgap if ((TargetBandgap >= 0.2) and (TargetBandgap <= 6.0)) else 1.1

        # Target top bandgap: used to take into account the part of solar spectrum already absorbed (for example in a top cell).
        #   useful to calculate the overall efficiency in a multijunction solar cell.
        #   For example for double junction solar cell, follow the steps below:
        #       1. set TargetBandgap to 1.65 eV and the TargetBandgapTop to 0, and calculate the corresponding efficiency and current-voltage characteristic.
        #       2. set TargetBandgap to 0.95 eV and the TargetBandgapTop to 1.65, and calculate the corresponding efficiency and current-voltage characteristic.
        #       deduce from the previous data the overall double junction solar cell efficiency.
        #       examples are given in ShockleyQueisserTJ.py and ShockleyQueisserDJ.py.
        self.Target_Bandgap_Top     = TargetBandgapTop if ((TargetBandgapTop >= 0.2) and (TargetBandgapTop <= 6.0) and (TargetBandgapTop > TargetBandgap)) else 0.0

        # Solar concentration (1 sun to 1000 suns)
        self.SolarConcentration     = SolarConcentration if ((SolarConcentration >= 1.0) and (SolarConcentration <= 1000.0)) else 1.0

        # OutputFilename: Output file name without extension (used to save figure in PDF format if in GUI mode, and the text output data).
        #   set to None to disable.
        self.OutputFilename         = OutputFilename
        if self.OutputFilename and (not self.OutputFilename.endswith('.pdf')):
            self.OutputFilename     = self.OutputFilename + '.pdf'
        # end if

        if self.useGUI:
            # GUI mode: calculation done in a working thread
            self.startGUI()
        else:
            # command-line mode
            self.start()
        # end if

        return

    # end calculate

    def isRunning(self):
        if (self.thread is None):
            return self.running
        # end if
        threadalive = self.thread.isAlive() if (sys.version_info[0] < 3) else self.thread.is_alive()
        if (not threadalive):
            self.thread  = None
            self.running = False
        # end if
        return self.running
    # end isRunning

    def setRunning(self, running = True):
        self.running = running
        if self.actionbutton is not None:
            self.actionbutton["text"] = "Calculate"
            if self.running:
                self.actionbutton.configure(style='Red.TButton')
            else:
                self.actionbutton.configure(style='Black.TButton')
                self.actionbutton = None
            # end if
        # end if
    # end setRunning

    # init the Tkinter GUI
    def startGUI(self):

        if self.GUIstarted or (not self.useGUI):
            return
        # end if

        try:
            self.plotcount      = 3
            self.curvecount     = 4
            self.xLabel         = {}
            self.yLabel         = {}
            self.xLabel[0]      = '$Wavelength\ (nm)$'
            self.xLabel[1]      = '$Bandgap\ (eV)$'
            self.xLabel[2]      = '$Voltage\ (V)$'
            self.yLabel[0]      = '$Irradiance\ (W/m^2/nm)$'
            self.yLabel[1]      = '$Efficiency\ (\%)$'
            self.yLabel[2]      = '$Current\ (mA/cm^2)$'

            self.root = Tk.Tk()
            self.root.bind_class("Entry","<Control-a>", self.onEntrySelectAll)
            self.root.bind_class("Entry","<Control-z>", self.onEntryUndo)
            self.root.bind_class("Entry","<Control-y>", self.onEntryRedo)
            self.root.withdraw()
            self.root.wm_title(self.name)

            self.figure = matplotlib.figure.Figure(figsize=(10,8), dpi=100, facecolor='#F1F1F1', linewidth=1.0, frameon=True)

            self.figure.subplots_adjust(top = 0.9, bottom = 0.1, left = 0.09, right = 0.95, wspace = 0.25, hspace = 0.25)

            self.plot       = {}
            # figure top-left part:     solar spetrum
            self.plot[0]    = self.figure.add_subplot(221)
            self.plot[0].set_xlim(0.0,      4000.0)
            self.plot[0].set_ylim(0.0,      1.8)
            # figure top-right part:    efficiency vs bandgap
            self.plot[1]    = self.figure.add_subplot(222)
            self.plot[1].set_xlim(0.0,      4.5)
            self.plot[1].set_ylim(0.0,      35.0)
            # figure bottom-left part:  current-voltage characteristic
            self.plot[2]    = self.figure.add_subplot(223)
            self.plot[2].set_xlim(0.0,      1.2)
            self.plot[2].set_ylim(-50.0,    10.0)
            # figure bottom-right part: report
            self.plot[4]    = self.figure.add_subplot(224)
            self.plot[4].set_xticks([])
            self.plot[4].set_yticks([])

            self.line0a     = None
            self.line0b     = None
            self.line0c     = None
            self.line2a     = None
            self.line2b     = None

            spx  = 6
            spy  = 12
            spxm = 1
            parFrame = Tk.Frame(self.root)
            parFrame.pack(fill=Tk.X, side=Tk.TOP, padx=spx, pady=spx)

            self.LLabel = Tk.Label(parFrame, text=" ")
            self.LLabel.pack(fill=Tk.X, side=Tk.LEFT, expand=True, padx=(spxm, spxm), pady=spy)
            
            self.SolarConcentrationLabel = Tk.Label(parFrame, text="Solar Concentration: ")
            self.SolarConcentrationLabel.pack(side=Tk.LEFT, padx=(spxm, spxm), pady=spy)
            SolarConcentrationValidate = (parFrame.register(self.onFloatValidate), '%P')
            self.SolarConcentrationEdit = Tk.Entry(parFrame, width=7, validate="key", vcmd=SolarConcentrationValidate)
            self.SolarConcentrationEdit.pack(side=Tk.LEFT, padx=(spxm, spx), pady=spy)
            self.SolarConcentrationEdit.insert(0, ("%.1f" % self.SolarConcentration) if ((self.SolarConcentration is not None) and (self.SolarConcentration >= 1.0) and (self.SolarConcentration <= 1000.0)) else "")
            self.SolarConcentrationEdit.prev = None
            self.SolarConcentrationEdit.next = None

            self.TargetLabel = Tk.Label(parFrame, text="Target bandgap (eV): ")
            self.TargetLabel.pack(side=Tk.LEFT, padx=(spx, spxm), pady=spy)
            TargetValidate = (parFrame.register(self.onFloatValidate), '%P')
            self.TargetEdit = Tk.Entry(parFrame, width=7, validate="key", vcmd=TargetValidate)
            self.TargetEdit.pack(side=Tk.LEFT, padx=(spxm, spx), pady=spy)
            self.TargetEdit.insert(0, ("%.3f" % self.Target_Bandgap) if (self.Target_Bandgap is not None) else "")
            self.TargetEdit.prev = None
            self.TargetEdit.next = None
            self.TargetTopLabel = Tk.Label(parFrame, text="Top bandgap (> target): ")
            self.TargetTopLabel.pack(side=Tk.LEFT, padx=(spx, spxm), pady=spy)
            self.TargetTopEdit = Tk.Entry(parFrame, width=7, validate="key", vcmd=TargetValidate)
            self.TargetTopEdit.pack(side=Tk.LEFT, padx=(spxm, spx), pady=spy)
            self.TargetTopEdit.insert(0, ("%.3f" % self.Target_Bandgap_Top) if ((self.Target_Bandgap_Top is not None) and (self.Target_Bandgap_Top > 0.0) and (self.Target_Bandgap_Top < self.Target_Bandgap)) else "")
            self.TargetTopEdit.prev = None
            self.TargetTopEdit.next = None

            self.TemperatureLabel = Tk.Label(parFrame, text="Temperature (K): ")
            self.TemperatureLabel.pack(side=Tk.LEFT, padx=(spx, spxm), pady=spy)
            TemperatureValidate = (parFrame.register(self.onFloatValidate), '%P')
            self.TemperatureEdit = Tk.Entry(parFrame, width=7, validate="key", vcmd=TemperatureValidate)
            self.TemperatureEdit.pack(side=Tk.LEFT, padx=(spxm, spx), pady=spy)
            self.TemperatureEdit.insert(0, ("%.1f" % self.Temperature) if (self.Temperature is not None) else "")
            self.TemperatureEdit.prev = None
            self.TemperatureEdit.next = None

            self.btnstyle_red = ttk.Style()
            self.btnstyle_red.configure("Red.TButton", foreground="#DE0015")
            self.btnstyle_black = ttk.Style()
            self.btnstyle_black.configure("Black.TButton", foreground="black")

            self.btnCalculate = ttk.Button(parFrame, width=18, text="Calculate", compound=Tk.LEFT, command=self.onStart)
            self.btnCalculate.pack(side=Tk.LEFT, padx=spx, pady=spy)
            self.btnCalculate.configure(style="Black.TButton")
            self.root.bind('<Return>', self.onEnter)

            self.RLabel = Tk.Label(parFrame, text=" ")
            self.RLabel.pack(fill=Tk.X, side=Tk.LEFT, expand=True, padx=(spxm, spxm), pady=spy)

            self.canvas = FigureCanvasTkAgg(self.figure, master=self.root)
            self.canvas._tkcanvas.config(highlightthickness=0)

            try:
                self.toolbar = NavigationToolbar(self)
            except:
                self.toolbar = None
                pass
            # end try
            self.toolbar.pack(side=Tk.BOTTOM, fill=Tk.X)

            self.toolbar.update()
            self.canvas._tkcanvas.pack(side=Tk.LEFT, fill=Tk.BOTH, expand=1)

            self.canvas.draw()
            
            self.root.protocol('WM_DELETE_WINDOW', self.onClose)

            self.linecolor          = ['g', 'm', 'b', 'r']
            self.linestyle          = ['-', '-', '-', '-']
            self.linesize           = [1.5, 1.5, 1.5, 1.5]
            self.line               = {}
            self.scatter            = {}
            for idc in range(0, self.curvecount):
                self.line[idc]      = None
                self.scatter[idc]   = None
            # end for

            # center the window
            iw = self.root.winfo_screenwidth()
            ih = self.root.winfo_screenheight()
            isize = (1020, 680)
            ix = (iw - isize[0]) / 2
            iy = (ih - isize[1]) / 2
            self.root.geometry("%dx%d+%d+%d" % (isize + (ix, iy)))

            self.root.minsize(800, 600)

            self.fontsize = 10

            for idp in range(0, self.plotcount):
                try:
                    self.plot[idp].tick_params(axis='x', labelsize=self.fontsize)
                    self.plot[idp].tick_params(axis='y', labelsize=self.fontsize)
                except:
                    [tx.label.set_fontsize(self.fontsize) for tx in self.plot[idp].xaxis.get_major_ticks()]
                    [ty.label.set_fontsize(self.fontsize) for ty in self.plot[idp].yaxis.get_major_ticks()]
                    pass
                # end try

                self.plot[idp].set_xlabel(self.xLabel[idp], fontsize=self.fontsize)
                self.plot[idp].set_ylabel(self.yLabel[idp], fontsize=self.fontsize)
            # end for

            if (os.name == "nt"):
                self.root.iconbitmap(r'iconmain.ico')
            else:
                iconmain = Tk.PhotoImage(file='iconmain.gif')
                self.root.tk.call('wm', 'iconphoto', self.root._w, iconmain)
            # end if

            self.popmenu = Tk.Menu(self.root, tearoff=0)
            self.popmenu.add_command(label="Calculate", command=self.onStart)
            self.popmenu.add_separator()
            self.popmenu.add_command(label="Auto scale", command=self.onAutoScale)
            self.popmenu.add_separator()
            self.popmenu.add_command(label="Close", command=self.onClose)
            self.popmenu.add_separator()
            self.popmenu.add_command(label="About...", command=self.onAbout)
            self.root.bind("<Button-3>", self.onPopmenu)

            self.root.deiconify()
            self.setFocus()

            self.GUIstarted = True

            #self.start()

            self.root.mainloop()

        except Exception as excT:
            excType, excObj, excTb = sys.exc_info()
            excFile = os.path.split(excTb.tb_frame.f_code.co_filename)[1]
            strErr  = "\n! cannot initialize GUI:\n  %s\n  in %s (line %d)\n" % (str(excT), excFile, excTb.tb_lineno)
            print(strErr)
            if self.GUIstarted and self.root:
                self.root.quit()
                self.root.destroy()
            # end if
            os._exit(1)
            # never reached
            pass
        # end try

    # end startGUI

    def isIncSorted(self, arr):
        for ii in range(arr.size - 1):
            if arr[ii + 1] <= arr[ii] :
                return False
            # end if
        # end for
        return True
    # end isIncSorted

    # load the solar spectrum file
    def loadSpectrum(self):

        if self.SpectrumLoaded:
            return True
        # end if

        try:

            # load the spectral data from the input file and calculate the total power
            self.SolarSpectrumData  = np.loadtxt(self.SolarSpectrumAMX, delimiter=self.DataDelimiter, skiprows=2, usecols=(0,1))
            self.Wavelength         = self.SolarSpectrumData[:,0]                   # nm
            self.Irradiance         = self.SolarSpectrumData[:,1]                   # W/m2/nm
            # check the data consistency
            if ((len (self.Wavelength)              < 100)                      or 
                (len (self.Irradiance)              < 100)                      or 
                (len (self.Irradiance)              != len(self.Wavelength))    or 
                (not (self.Wavelength               >= 200.0).all())            or 
                (not (self.Wavelength               <= 10000.0).all())          or 
                (not (self.Irradiance               >= 0.0).all())              or 
                (not (self.Irradiance               <= 10.0).all())             or
                (not self.isIncSorted(self.Wavelength))):
                raise Exception('invalid wavelength/irradiance')
            # end if
            
            self.Energy             = self.nmeV / self.Wavelength                   # eV
            self.SolarPower         = sp.trapz(self.Irradiance, x=self.Wavelength)  # for AM1.5 solar spectrum, the total power is close to 1000 W/m2 or 100 mW/cm2
            if (self.SolarPower < 1.0) or (self.SolarPower > 10000.0):
                raise Exception('invalid total power density')
            # end if

            self.WavelengthMin      = self.Wavelength[0] + 10.0
            self.WavelengthMax      = self.Wavelength[len(self.Wavelength) - 1] - 10.0
            self.BandgapMin         = self.nmeV / self.WavelengthMax    # in eV
            self.BandgapMax         = self.nmeV / self.WavelengthMin    # in eV
            # the bandgap range in eV
            self.BandgapRange       = np.arange(self.BandgapMin, self.BandgapMax, (self.BandgapMax - self.BandgapMin) / 500.0)
            self.CounterMax         = len(self.BandgapRange)

            self.SpectrumLoaded     = True
            self.SQ_Done            = False

            return True

        except Exception as excT:
            excType, excObj, excTb = sys.exc_info()
            excFile = os.path.split(excTb.tb_frame.f_code.co_filename)[1]
            strErr  = "\n! cannot load the solar spectrum data:\n  %s\n  in %s (line %d)\n" % (str(excT), excFile, excTb.tb_lineno)
            print(strErr)
            if self.GUIstarted and self.root:
                self.root.quit()
                self.root.destroy()
            # end if
            os._exit(1)
            # never reached
            pass
        # end try

    # end loadSpectrum

    # Planck distribution
    def PlanckDistribution(self, Energy):
        aPlanck = np.array([])
        for aE in Energy:
            tE      = aE * self.q       # convert E from eV to J
            tP      = (2.0 * self.pi / ((self.h ** 3.0) * (self.c ** 2.0))) * (tE ** 2) / (math.exp(tE / self.kTJ) - 1.0)
            aPlanck = np.append(aPlanck, -self.q * tP)
        # end for
        return aPlanck
    # end PlanckDistribution

    # calculate the efficiency (and other photovoltaic parameters) for a given bandgap
    # Theory by W. Shockley and H. J. Queisser in Journal of Applied Physics 32 (1961)
    def calculateEfficiency(self, Bandgap, BandgapTop = 0.0):

        try:

            if (Bandgap < self.BandgapMin) or (Bandgap > self.BandgapMax):
                raise Exception("invalid bandgap: %.3f" % Bandgap)
            # end if
            aLambda         = self.nmeV / Bandgap
            aLambdaLow      = 0.0
            CutSpectrum     = (BandgapTop > Bandgap) and (BandgapTop >= self.BandgapMin) and (BandgapTop <= self.BandgapMax)
            if CutSpectrum:
                aLambdaLow  = self.nmeV / BandgapTop
            # end if
            aWavelength     = np.copy(self.SolarSpectrumData[:,0])
            aWavelength     = aWavelength[(aWavelength >= aLambdaLow) & (aWavelength <= aLambda)]
            aEnergy         = self.nmeV / aWavelength
            indexT          = np.where(np.logical_and(self.SolarSpectrumData[:,0] >= aLambdaLow, self.SolarSpectrumData[:,0] <= aLambda))
            aIrradiance     = self.SolarConcentration * np.take(self.SolarSpectrumData[:,1], indexT, axis=0)
            aFlux           = aIrradiance * aWavelength * 1e-9 / self.hc
            aJSC            = self.q * sp.trapz(aFlux, x=aWavelength)          # Short-Circuit Current in A/m2
            aPlanck         = self.PlanckDistribution(aEnergy)
            aJ0             = self.q * sp.trapz(aPlanck, x=aEnergy)            # Dark Current in A/m2
            aVOC            = self.kTeV * math.log((aJSC / aJ0) + 1.0)         # Open-Circuit Voltage in V
            aVstep          = aVOC / 500.0
            aVoltage        = np.arange(0.0, aVOC + aVstep, aVstep)
            aCurrent        = np.array([])
            aVm             = 0.0
            aJm             = 0.0
            aPm             = 0.0
            for aV in aVoltage:
                aJ          = -aJSC + (aJ0 * (math.exp(aV / self.kTeV) - 1.0))
                aCurrent    = np.append(aCurrent, aJ)
                if (aJ < 0.0) and (math.fabs(aJ * aV) > aPm):
                    aPm     = math.fabs(aJ * aV)
                    aVm     = aV
                    aJm     = aJ
                    aFF     = aPm / (aJSC * aVOC)
                # end if
            # end for
            aEff            = -100.0 * np.min(aCurrent * aVoltage) / (self.SolarPower * self.SolarConcentration)
            return (aEff, aVOC, aJSC, aFF, aVm, aJm, aVoltage, aCurrent, None)

        except Exception as excT:

            excType, excObj, excTb = sys.exc_info()
            excFile = os.path.split(excTb.tb_frame.f_code.co_filename)[1]
            strErr  = "\n! %s\n  in %s (line %d)\n" % (str(excT), excFile, excTb.tb_lineno)
            return (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, np.array([]), np.array([]), strErr)
            # never reached
            pass

        # end try

    # end calculateEfficiency

    def monitorCalculation(self):
        running = self.isRunning()
        try:
            if not running:
                self.setRunning(running = False)
                if self.threadfinish is not None:
                    self.threadfinish()
                    self.threadfinish = None
                # end if
                if self.GUIstarted:
                    self.btnCalculate["text"] = "Calculate"
                # end if
                return
            # end if
            if self.root:
                if self.GUIstarted and (self.CounterMax > 0) and (not self.SQ_Done):
                    tPC = (100 * self.Counter) / self.CounterMax
                    self.btnCalculate["text"] = "Calculate" if (tPC > 99) else ("Calculate (%.0f %%)" % tPC)
                # end if
                self.root.after(self.timerduration if ((self.timerduration >= 100) and (self.timerduration <= 1000)) else 200, self.monitorCalculation)
            # end if
        except Exception as excT:
            pass
        # end try
    # end monitorCalculation

    # start calculations of the efficiency vs bandgap curve
    def start(self):

        if self.isRunning():
            return False
        # end if

        if self.useGUI and self.GUIstarted:
            # GUI mode

            # get the input parameters
            strT = self.SolarConcentrationEdit.get().strip("\r\n\t")
            try:
                fT = float(strT)
                if (fT >= 1.0) and (fT <= 1000.0):
                    # recalculate the Shockley-Queisser curve if changing temperature or solar concentration
                    if self.SQ_Done and (fT != self.SolarConcentration):
                        self.SQ_Done = False
                    # end if
                    self.SolarConcentration = fT
                else:
                    self.SolarConcentrationEdit.delete(0, Tk.END)
                    self.SolarConcentrationEdit.insert(0, "%.1f" % self.SolarConcentration)
                # end if
            except ValueError:
                pass
            # end try

            strT = self.TargetEdit.get().strip("\r\n\t")
            try:
                fT = float(strT)
                if (fT >= 0.2) and (fT <= 6.0):
                    self.Target_Bandgap = fT
                else:
                    self.TargetEdit.delete(0, Tk.END)
                    self.TargetEdit.insert(0, "%.3f" % self.Target_Bandgap)
                # end if
            except ValueError:
                pass
            # end try

            strT = self.TargetTopEdit.get().strip("\r\n\t")
            try:
                fT = float(strT)
                if (fT >= 0.2) and (fT <= 6.0) and (fT > self.Target_Bandgap):
                    self.Target_Bandgap_Top = fT
                else:
                    self.TargetTopEdit.delete(0, Tk.END)
                    self.TargetTopEdit.insert(0, "%.3f" % self.Target_Bandgap_Top)
                # end if
            except ValueError:
                pass
            # end try

            strT = self.TemperatureEdit.get().strip("\r\n\t")
            try:
                fT = float(strT)
                if (fT >= 100.0) and (fT <= 700.0):
                    # recalculate the Shockley-Queisser curve if changing temperature or solar concentration
                    if self.SQ_Done and (fT != self.Temperature):
                        self.SQ_Done = False
                    # end if
                    self.Temperature    = fT                                                # Temperature in K
                    self.kTeV           = 0.02585202874091 * self.Temperature / 300.0       # kT in eV
                    self.kTJ            = self.kTeV * self.q                                # kT in J
                else:
                    self.TemperatureEdit.delete(0, Tk.END)
                    self.TemperatureEdit.insert(0, "%.1f" % self.Temperature)
                # end if
            except ValueError:
                pass
            # end try

            # start calculations
            self.actionbutton = self.btnCalculate
            self.setRunning(running = True)
            self.thread = CalculationThread(id=1, func=self.run)
            self.threadfinish = self.updatePlot
            self.thread.start()
            self.monitorCalculation()

        else:
            # Command-line mode (or GUI initialization step)

            if  ((self.Target_Bandgap       < 0.2)          or 
                 (self.Target_Bandgap       > 6.0)          or 
                 (self.Temperature          < 100.0)        or 
                 (self.Temperature          > 700.0)        or 
                 (self.SolarConcentration   < 1.0)          or 
                 (self.SolarConcentration   > 1000.0)):
                strErr = "\n! cannot perform calculations:\n  invalid input parameters\n"
                if self.useGUI and (self.report is not None):
                    self.report.set_color('red')
                    self.report.set_text(strErr)
                    self.canvas.draw()
                # end if
                print(strErr)
            else:
                if self.useGUI:
                    self.actionbutton = self.btnCalculate
                    self.setRunning(running = True)
                    self.thread = CalculationThread(id=1, func=self.run)
                    self.threadfinish = self.updatePlot
                    self.thread.start()
                    self.monitorCalculation()
                else:
                    SQ_Done = self.SQ_Done
                    self.setRunning(running = True)
                    self.run()
                    self.setRunning(running = False)
                    self.doSave(self.OutputFilename)
                    if self.verbose:
                        if not SQ_Done:
                            treport = ("Shockley-Queisser limit for T = %.1f K and SC = %.1f sun%s\nBandgap ; Efficiency ; JSC ; VOC ; Fill Factor" % (self.Temperature, self.SolarConcentration, "" if (self.SolarConcentration <= 1.0) else "s")) + ("\n%.3f eV ; %05.2f %% ; %06.3f mA/cm2 ; %05.3f V ; %05.3f %%" % (self.SQ_Bandgap, self.SQ_Efficiency, 0.1 * self.SQ_JSC, self.SQ_VOC, 100.0 * self.SQ_FF))
                            print("\n======================================================================\n" + treport + "\n======================================================================\n")
                        # end if
                        treport = ("Target for T = %.1f K and SC = %.1f sun%s\nBandgap ; Efficiency ; JSC ; VOC ; Fill Factor" % (self.Temperature, self.SolarConcentration, "" if (self.SolarConcentration <= 1.0) else "s")) + ("\n%.3f eV ; %05.2f %% ; %06.3f mA/cm2 ; %05.3f V ; %05.3f %%" % (self.Target_Bandgap, self.Target_Efficiency, 0.1 * self.Target_JSC, self.Target_VOC, 100.0 * self.Target_FF))
                        print("\n----------------------------------------------------------------------\n" + treport + "\n----------------------------------------------------------------------\n")
                    # end if
                # end if
            # end if
        # end if

     # end start

    # calculate the efficiency vs bandgap curve
    def run(self):

        try:
            if self.verbose:
                print("\ncalculating...")
            # end if

            self.Counter = 0

            # load the spectral data
            self.loadSpectrum()
            if not self.SpectrumLoaded:
                if self.verbose:
                    print("\ndone.")
                # end if
                return False
            # end if

            # to determine the calculation duration
            ticT = time.time()

            # Shockley-Queisser calculated parameters
            if not self.SQ_Done:

                self.Bandgap            = np.array([])
                self.Efficiency         = np.array([])

                self.SQ_Efficiency      = 0.0
                self.SQ_JSC             = 0.0
                self.SQ_VOC             = 0.0
                self.SQ_FF              = 0.0
                self.SQ_Bandgap         = 0.0
                self.SQ_Voltage         = np.array([])
                self.SQ_Current         = np.array([])

                for aGap in self.BandgapRange:
                    (aEff, aVOC, aJSC, aFF, aVm, aJm, aVoltage, aCurrent, strRet) = self.calculateEfficiency(aGap, 0.0)
                    if strRet is not None:
                        if self.useGUI and self.GUIstarted and (self.report is not None):
                            self.report.set_color('red')
                            self.report.set_text(strRet)
                            self.canvas.draw()
                        # end if
                        if self.verbose:
                            print(strRet)
                            print("\ndone.")
                        # end if
                        return False
                    # endif
                    if aEff > self.SQ_Efficiency:
                        self.SQ_Efficiency  = aEff
                        self.SQ_JSC         = aJSC
                        self.SQ_VOC         = aVOC
                        self.SQ_FF          = aFF
                        self.SQ_Vm          = aVm
                        self.SQ_Jm          = aJm
                        self.SQ_Bandgap     = aGap
                        self.SQ_Voltage     = np.copy(aVoltage)               # in V
                        self.SQ_Current     = 0.1 * np.copy(aCurrent)         # in mA/cm2
                    # end if
                    self.Bandgap            = np.append(self.Bandgap, aGap)
                    self.Efficiency         = np.append(self.Efficiency, aEff)

                    self.tic                = float(time.time() - ticT)
                    self.Counter           += 1

                # end for

                self.SQ_Done = True

            # end if

            # calculate the efficiency for the target bandgap value
            (self.Target_Efficiency, self.Target_VOC, self.Target_JSC, self.Target_FF, self.Target_Vm, self.Target_Jm, self.Target_Voltage, self.Target_Current, strRet) = self.calculateEfficiency(self.Target_Bandgap, self.Target_Bandgap_Top)
            if strRet is not None:
                if self.verbose:
                    print (strRet)
                    print("\ndone.")
                # end if
                return False
            # endif
            self.Target_Current = 0.1 * np.copy(self.Target_Current)         # convert from A/m2 to mA/cm2

            self.datax      = {}
            self.datax[0]   = self.Wavelength
            self.datax[1]   = self.Bandgap
            self.datax[2]   = self.SQ_Voltage
            self.datax[3]   = self.Target_Voltage
            self.datay      = {}
            self.datay[0]   = self.Irradiance
            self.datay[1]   = self.Efficiency
            self.datay[2]   = self.SQ_Current
            self.datay[3]   = self.Target_Current

            self.tic = float(time.time() - ticT)

            if self.verbose:
                print("\ndone. elapsed time = %.3f sec." % self.tic)
            # end if

            return True

        except Exception as excT:

            excType, excObj, excTb = sys.exc_info()
            excFile = os.path.split(excTb.tb_frame.f_code.co_filename)[1]
            strErr  = "\n! %s\n  in %s (line %d)\n" % (str(excT), excFile, excTb.tb_lineno)
            if self.useGUI and (self.report is not None):
                self.report.set_color('red')
                self.report.set_text(strErr)
                self.canvas.draw()
            # end if
            if self.verbose:
                print(strErr)
            # end if
            return False
            # never reached
            pass

        # end try

    # end run

    def setFocus(self):
        if (not self.useGUI) or (not self.root):
            return
        # end if
        self.root.attributes('-topmost', 1)
        self.root.attributes('-topmost', 0)
        self.root.after(10, lambda: self.root.focus_force())
    # end setFocus

    # plot the efficiency vs bandgap curve
    def updatePlot(self):

        if (not self.useGUI):
            return
        # end if

        try:
 
            if not self.PlotInitialized:

                idp = 0
                for idc in range(0, self.curvecount):
                    self.line[idc], = self.plot[idp].plot(np.array([]), np.array([]), self.linestyle[idc], linewidth=self.linesize[idc], zorder=4)
                    self.line[idc].set_color(self.linecolor[idc])
                    if (idp < (self.plotcount - 1)):
                        idp += 1
                    # end if
                # end for

                self.scatter[0], = self.plot[1].plot(np.array([]),      np.array([]),   'b^', zorder=4, label=" ")
                self.scatter[0].set_markerfacecolor('b')
                self.scatter[0].set_markeredgecolor('b')
                self.scatter[0].set_markersize(5)

                self.scatter[1], = self.plot[1].plot(np.array([]),      np.array([]),   'ro', zorder=4, label=" ")
                self.scatter[1].set_markerfacecolor('r')
                self.scatter[1].set_markeredgecolor('r')
                self.scatter[1].set_markersize(5)

                self.scatter[2], = self.plot[2].plot(np.array([]),      np.array([]),   'b^', zorder=4, label=" ")
                self.scatter[2].set_markerfacecolor('b')
                self.scatter[2].set_markeredgecolor('b')
                self.scatter[2].set_markersize(5)

                self.scatter[3], = self.plot[2].plot(np.array([]),      np.array([]),   'ro', zorder=4, label=" ")
                self.scatter[3].set_markerfacecolor('r')
                self.scatter[3].set_markeredgecolor('r')
                self.scatter[3].set_markersize(5)

                self.plot[0].legend(['ASTM AM1.5 G-173 (1 sun)'], loc='upper right', fontsize='x-small')
                self.plot[1].legend(numpoints=1, fontsize='x-small', loc='best')
                self.plot[2].legend(numpoints=1, fontsize='x-small', loc='best')

                self.line0a = self.plot[0].axvline(x=0,                 ymin=0, ymax=1, linewidth=1, color='b',     linestyle='-.')
                self.line0b = self.plot[0].axvline(x=0,                 ymin=0, ymax=1, linewidth=1, color='r',     linestyle='-.')
                self.line0c = self.plot[0].axvline(x=0,                 ymin=0, ymax=1, linewidth=1, color='olive', linestyle='-.')

                self.line2a = self.plot[2].axhline(y=0,                 xmin=0, xmax=1, linewidth=2, color='k')
                self.line2b = self.plot[2].axvline(x=0,                 ymin=0, ymax=1, linewidth=2, color='k')

                self.plot[4].set_xticks([])
                self.plot[4].set_yticks([])
                (tleft, tright) = self.plot[4].get_xlim()
                (tbottom, ttop) = self.plot[4].get_ylim()
                afont = FontProperties()
                tfont = afont.copy()
                tfont.set_style('normal')
                tfont.set_weight('bold')
                tfont.set_size('x-small')
                self.report = self.plot[4].text(0.5 * (tleft + tright), 0.5 * (tbottom + ttop), "", horizontalalignment='center', verticalalignment='center', clip_on=True, fontproperties=tfont, color='black', transform=self.plot[4].transAxes)

                for idp in range(0, self.plotcount):
                    self.plot[idp].get_xaxis().set_visible(True)
                    self.plot[idp].get_yaxis().set_visible(True)
                # end for

                self.PlotInitialized = True

            # end if

            idp = 0
            for idc in range(0, self.curvecount):
                self.line[idc].set_xdata(self.datax[idc])
                self.line[idc].set_ydata(self.datay[idc])
            # end for

            self.scatter[0].set_xdata(self.SQ_Bandgap)
            self.scatter[0].set_ydata(self.SQ_Efficiency)
            self.scatter[0].set_label("Max   : %05.2f %% for %.3f eV" % (self.SQ_Efficiency, self.SQ_Bandgap))

            self.scatter[1].set_xdata(self.Target_Bandgap)
            self.scatter[1].set_ydata(self.Target_Efficiency)
            self.scatter[1].set_label("Target: %05.2f %% for %.3f eV" % (self.Target_Efficiency, self.Target_Bandgap))

            self.scatter[2].set_xdata(self.SQ_Vm)
            self.scatter[2].set_ydata(0.1 * self.SQ_Jm)
            self.scatter[2].set_label("Max   : %05.2f %% for %.3f eV" % (self.SQ_Efficiency, self.SQ_Bandgap))

            self.scatter[3].set_xdata(self.Target_Vm)
            self.scatter[3].set_ydata(0.1 * self.Target_Jm)
            self.scatter[3].set_label("Target: %05.2f %% for %.3f eV" % (self.Target_Efficiency, self.Target_Bandgap))

            self.plot[0].legend(['ASTM AM1.5 G-173 (1 sun)'], loc='upper right', fontsize='x-small')
            self.plot[1].legend(numpoints=1, fontsize='x-small', loc='best')
            self.plot[2].legend(numpoints=1, fontsize='x-small', loc='best')

            self.line0a.set_xdata(self.nmeV / self.SQ_Bandgap)
            self.line0b.set_xdata(self.nmeV / self.Target_Bandgap)
            if (self.Target_Bandgap_Top > self.Target_Bandgap):
                self.line0c.set_xdata(self.nmeV / self.Target_Bandgap_Top)
            else:
                self.line0c.set_xdata((self.nmeV / self.BandgapMax) + 10.0)
            # end if

            treport = ("Bandgap (for Max)   : %.3f eV\n\n" % self.SQ_Bandgap) + ("Bandgap (for Target): %.3f eV" % self.Target_Bandgap) + "\n\nEfficiency ; JSC ; VOC ; Fill Factor" + ("\n\nMax   : %05.2f %% ; %06.3f mA/cm2 ; %05.3f V ; %05.3f %%" % (self.SQ_Efficiency, 0.1 * self.SQ_JSC, self.SQ_VOC, 100.0 * self.SQ_FF)) + ("\n\nTarget: %05.2f %% ; %06.3f mA/cm2 ; %05.3f V ; %05.3f %%" % (self.Target_Efficiency, 0.1 * self.Target_JSC, self.Target_VOC, 100.0 * self.Target_FF))
            self.report.set_color('black')
            self.report.set_text(treport)

            for idp in range(0, self.plotcount):
                self.plot[idp].relim()
                self.plot[idp].autoscale()
            # end for

            self.setFocus()
            self.canvas.draw()

        except Exception as excT:

            excType, excObj, excTb = sys.exc_info()
            excFile = os.path.split(excTb.tb_frame.f_code.co_filename)[1]
            strErr  = "\n! cannot plot data:\n  %s\n  in %s (line %d)\n" % (str(excT), excFile, excTb.tb_lineno)
            if self.useGUI and (self.report is not None):
                self.report.set_color('red')
                self.report.set_text(strErr)
                self.canvas.draw()
            # end if
            if self.verbose:
                print(strErr)
            # end if
            pass

        # end try

    # end updatePlot

    def onFloatValidate(self, sp):
        try:
            if (not sp):
                return True
            # end if
            slen = len(sp)
            if slen < 1:
                return True
            # end if
            if (slen > 12):
                sp = sp[:12]
                slen = len(sp)
            #
            if (slen <= 12):
                try:
                    spr = sp + '0'
                    float(spr)
                except ValueError:
                    if any(c not in '0123456789-+.e' for c in sp):
                        return False
                    # end if
                    pass
                # end try
                self.TargetEdit.prev = sp
                return True
            # end if
            return False
        except:
            return True
        # end try
    # end onFloatValidate

    def onEnter(self, tEvent):
        return self.start()
    # end onEnter

    def onStart(self):
        return self.start()
    # end onStart

    def doSave(self, strFilename, savePDF = False):

        if (not strFilename) or self.isRunning():
            return
        # end if

        try:

            if savePDF and self.useGUI and self.GUIstarted:
                # save figure in PDF format
                pdfT = PdfPages(strFilename)
                pdfT.savefig(self.figure)
                pdfT.close()
                # and in PNG format
                strPNG = os.path.splitext(strFilename)[0]
                strPNG = strPNG + '.png'
                pl.savefig(strPNG, dpi=600)
            # end if

            # save output data in text format
            strF = os.path.splitext(strFilename)[0]
            fileEff = strF + '_Efficiency.txt'          # efficiency vs bandgap
            fileJVM = strF + '_JV_Max.txt'              # current-voltage characteristic corresponding to the maximum efficiency
            fileJVT = strF + '_JV_Target.txt'           # current-voltage characteristic corresponding to the target bandgap
            np.savetxt(fileEff, np.c_[self.Bandgap, self.Efficiency],               fmt='%.4f\t%.6f', delimiter=self.DataDelimiter, newline='\n', header='Bandgap (eV)\tEfficiency (%)')
            np.savetxt(fileJVM, np.c_[self.SQ_Voltage, self.SQ_Current],            fmt='%.4f\t%.6f', delimiter=self.DataDelimiter, newline='\n', header=("Max Efficiency: %05.2f %% for bandgap = %.3f eV\n" % (self.SQ_Efficiency, self.SQ_Bandgap)) + 'Voltage (V)\tCurrent (mA/cm2)')
            np.savetxt(fileJVT, np.c_[self.Target_Voltage, self.Target_Current],    fmt='%.4f\t%.6f', delimiter=self.DataDelimiter, newline='\n', header=("Target Efficiency: %05.2f %% for bandgap = %.3f eV\n" % (self.Target_Efficiency, self.Target_Bandgap)) + ("Max Efficiency: %05.2f %% for bandgap = %.3f eV\n" % (self.SQ_Efficiency, self.SQ_Bandgap)) + 'Voltage (V)\tCurrent (mA/cm2)')

        except Exception as excT:

            strErr = "\n! cannot save output data:\n  %s\n" % str(excT)
            if self.useGUI and self.GUIstarted and (self.report is not None):
                self.report.set_color('red')
                self.report.set_text(strErr)
                self.canvas.draw()
            # end if
            if self.verbose:
                print(strErr)
            # end if
            return False
            # never reached
            pass

        # end try

    # end doSave

    def onSave(self):
        if (not self.useGUI) or self.isRunning():
            return
        # end if
        fileopt = {}
        fileopt['defaultextension'] = '.pdf'
        fileopt['filetypes'] = [('PDF files', '.pdf')]
        fileopt['initialfile'] = self.OutputFilename
        fileopt['parent'] = self.root
        fileopt['title'] = 'Save figure'
        pdfFilename = tkFileDialog.asksaveasfilename(**fileopt)
        if pdfFilename:
            self.OutputFilename = pdfFilename
            self.doSave(self.OutputFilename, savePDF = True)
        # end if
    # end onSave

    def onAutoScale(self):
        if (not self.useGUI) or self.isRunning():
            return
        # end if
        for idp in range(0, self.plotcount):
            self.plot[idp].relim()
            self.plot[idp].autoscale()
        # end for
        self.canvas.draw()
    # end onAutoScale

    def onEntryUndo(self, event):
        if not self.useGUI:
            return
        # end if
        try:
            if event.widget.prev is not None:
                event.widget.next = event.widget.get()
                strT = event.widget.prev
                idx = event.widget.index(Tk.INSERT)
                event.widget.delete(0, Tk.END)
                event.widget.insert(0, strT)
                event.widget.prev = strT
                event.widget.icursor(idx + 1)
        except:
            pass
        # end try
    # end onEntryUndo

    def onEntryRedo(self, event):
        if not self.useGUI:
            return
        # end if
        try:
            if event.widget.next is not None:
                idx = event.widget.index(Tk.INSERT)
                strT = event.widget.prev
                event.widget.delete(0, Tk.END)
                event.widget.insert(0, event.widget.next)
                event.widget.prev = strT
                event.widget.icursor(idx + 1)
        except:
            pass
        # end try
    # end onEntryRedo

    def onEntrySelectAll(self, event):
        if not self.useGUI:
            return
        # end if
        try:
            event.widget.select_range(0, Tk.END)
        except:
            pass
        # end try
    # end onEntrySelectAll

    def onAbout(self):
        if not self.useGUI:
            return
        # end if
        tkMessageBox.showinfo(self.name,
                             (self.name                                                         +
                              "\n"                                                              +
                              self.__version__                                                  +
                              "\nCopyright(C) 2018-2019 Pr. Sidi OULD SAAD HAMADY \n"           +
                              "Université de Lorraine, France \n"                               +
                              "sidi.hamady@univ-lorraine.fr \n"                                 +
                              "https://github.com/sidihamady/Shockley-Queisser \n"              +
                              "http://www.hamady.org/photovoltaics/ShockleyQueisser.zip \n"     +
                              "Under MIT license \nSee Copyright Notice in COPYRIGHT"),
                              parent=self.root)
    # end onAbout

    def onPopmenu(self, event):
        if (not self.useGUI):
            return
        # end if
        try:
            self.popmenu.post(event.x_root, event.y_root)
        except:
            pass
        # end try
    # end onPopmenu

    def onClose(self):
        if (not self.useGUI) or (self.root is None):
            return
        # end if
        try:
            if self.isRunning():
                tkMessageBox.showinfo(self.name, "Please wait until calculations done.", parent=self.root)
                return
            # end if
            if not tkMessageBox.askyesno(self.name, "Close " + self.name + "?", default=tkMessageBox.NO, parent=self.root):
                return
            # end if
            self.root.quit()
            self.root.destroy()
            self.root = None
        except:
            pass
        # end try
    # end onClose

# end ShockleyQueisserCore class
