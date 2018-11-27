#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ======================================================================================================
# Solar Cell Shockley-Queisser Limit Calculator
# Code written by:
#   Sidi Hamady
#   Université de Lorraine, France
#   sidi.hamady@univ-lorraine.fr
# See Copyright Notice in COPYRIGHT
# HowTo in README.md and README.pdf
# https://github.com/sidihamady/Shockley-Queisser
# http://www.hamady.org/photovoltaics/ShockleyQueisser.zip
# ======================================================================================================

# ShockleyQueisserTJ.py
#   an example on how to use the command-line mode to perform specific calculations
#   such as multijunction solar cell efficiency:
#   calculates the efficiency of a triple junction solar cell and plots current-voltage characteristics

# import the program core class
from ShockleyQueisserCore import *

# create an instance of the core class
SCC = ShockleyQueisserCore()

# set the triple junction solar cell parameters
aTargetBandgap      = np.array([1.9,    1.4,        0.9])      # bandgaps
aTargetBandgapTop   = np.array([0,      1.9,        1.4])      # top bandgap for each cell (0 for the top cell)
aTargetLabel        =          ["Top",  "Middle",   "Bottom"]
aTargetLen          = len(aTargetBandgap)
#

# initialize the output parameters
aJSC                = np.array([])
aVOC                = np.array([])
aFF                 = np.array([])
aV                  = {}
aJ                  = {}
jj                  = 0
jjx                 = 0
aVOCx               = 0.0
#

for TB,TBT in zip(aTargetBandgap, aTargetBandgapTop):

    # set useGUI to False to use the command-line mode
    SCC.calculate(
            TargetBandgap           = TB,
            TargetBandgapTop        = TBT,
            Temperature             = 300.0,
            SolarConcentration      = 1.0,
            OutputFilename          = None,
            useGUI                  = False
            )

    aJSC    = np.append(aJSC, SCC.Target_JSC)
    aVOC    = np.append(aVOC, SCC.Target_VOC)
    aFF     = np.append(aFF,  SCC.Target_FF)

    aV[jj]  = np.copy(SCC.Target_Voltage)
    aJ[jj]  = np.copy(SCC.Target_Current)
    if (jj == 0) or ((jj > 0) and (SCC.Target_VOC < aVOCx)):
        aVOCx   = SCC.Target_VOC
        jjx     = jj
     # end if
    jj += 1

# end for

# get the multijunction solar cell current-voltage characteristic
aVx = np.array([])
aJx = np.array([])
for ii in range(0, len(aV[jjx])):
    tV = 0.0
    tJ = 0.0
    for jj in range(0, aTargetLen):
        # sum voltage
        tV += aV[jj][ii]
        # take the min current (J is negative, from -JSC to 0)
        if (jj == 0) or ((jj > 0) and (aJ[jj][ii] > tJ)):
            tJ = aJ[jj][ii]
        # end if
    # end for
    aVx = np.append(aVx, tV)
    aJx = np.append(aJx, tJ)
# end for

aPm     = np.min(aJx * aVx)                                 # nominal power in mW/cm2
aPsolar = 0.1 * SCC.SolarPower * SCC.SolarConcentration     # solar power, to convert from W/m2 to mW/cm2
aEff    = -100.0 * aPm / aPsolar                            # efficiency in percentage
treport = ("Triple junction solar cell efficiency = %.3f %%\n with a nominal power of %.3f mW/cm2" % (aEff, -aPm))
print("\n----------------------------------------------------------------------\n" + treport + "\n----------------------------------------------------------------------\n")

# plot the current-voltage characteristics
fig         = pl.figure(figsize=(10, 6), dpi=100, facecolor='#F1F1F1', linewidth=1.0, frameon=True)
fig.canvas.set_window_title('Triple Junction Solar Cell')
ax          = fig.add_subplot(111)
linestyle   = ['-', '-', '-', '-']
linecolor   = ['m', 'b', 'r', 'olive']
linesize    = [1.5, 1.5, 1.5, 2.0]
for jj in range(0, aTargetLen):
    tline,  = ax.plot(aV[jj], aJ[jj], linestyle[jj], linewidth=linesize[jj], label=aTargetLabel[jj])
    tline.set_color(linecolor[jj])
# end for
tline,      = ax.plot(aVx, aJx, linestyle[aTargetLen], linewidth=linesize[aTargetLen], label="Triple Junction")
tline.set_color(linecolor[aTargetLen])
leg         = pl.legend(loc='upper center', ncol=4, frameon=1, fontsize='medium')
legframe    = leg.get_frame()
legframe.set_facecolor('white')
legframe.set_edgecolor('black')
leg.get_frame().set_alpha(0.5)
pl.axhline(y=0, xmin=0, xmax=1, linewidth=2, color='k')
pl.axvline(x=0, ymin=0, ymax=1, linewidth=2, color='k')
ax.set_xlabel('$Voltage\ (V)$',         fontsize=14)
ax.set_ylabel('$Current\ (mA/cm^2)$',   fontsize=14)
pl.title('Triple Junction Solar Cell Current-Voltage Characteristic')
pl.show()
#