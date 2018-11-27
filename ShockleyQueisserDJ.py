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

# ShockleyQueisserDJ.py
#   an example on how to use the command-line mode to perform specific calculations
#   such as multijunction solar cell efficiency:
#   calculates and plots the efficiency of a double junction solar cell with respect to the top and bottom bandgap

# import the program core class
from ShockleyQueisserCore import *

# create an instance of the core class (and disable printing output)
SCC = ShockleyQueisserCore(verbose = False)

aGapTop = np.arange(0.25, 3.25, 0.05)
aGapBot = np.copy(aGapTop)
aEff    = np.zeros((len(aGapBot), len(aGapTop)))

# calculation takes about 2 minutes on a quad-core 2.9 GHz processor for 1800 bandgap pairs
print("\ncalculating...")
ticT = time.time()

GapTopMax   = 0.0
GapBotMax   = 0.0
EffMax      = 0.0

for ii in range(0, len(aGapBot)):

    for jj in range(ii + 1, len(aGapTop)):

        GapT = aGapTop[jj]
        GapB = aGapBot[ii]

        aV      = {}
        aJ      = {}
        aVOCx   = 0.0
        cc      = 0
        ccx     = 0
        TBT     = 0.0
        for TB in (GapT, GapB):
            # set useGUI to False to use the command-line mode
            SCC.calculate(
                    TargetBandgap           = TB,
                    TargetBandgapTop        = TBT,
                    Temperature             = 300.0,
                    SolarConcentration      = 1.0,
                    OutputFilename          = None,
                    useGUI                  = False
                    )
            aV[cc]  = np.copy(SCC.Target_Voltage)
            aJ[cc]  = np.copy(SCC.Target_Current)
            if (cc == 0) or ((cc > 0) and (SCC.Target_VOC < aVOCx)):
                aVOCx   = SCC.Target_VOC
                ccx     = cc
            # end if
            cc += 1
            TBT = TB
        # end for

        try:
            # get the double junction solar cell current-voltage characteristic
            aVx = np.array([])
            aJx = np.array([])
            for nn in range(0, len(aV[ccx])):
                tV = 0.0
                tJ = 0.0
                for cc in range(0, 2):
                    # sum voltage
                    tV += aV[cc][nn] if (nn < len(aV[cc])) else aV[cc][len(aV[cc]) - 1]
                    # take the min current (J is negative, from -JSC to 0)
                    if (cc == 0) or ((cc > 0) and (aJ[cc][nn] > tJ)):
                        tJ = aJ[cc][nn] if (nn < len(aV[cc])) else aJ[cc][len(aV[cc]) - 1]
                    # end if
                # end for
                aVx = np.append(aVx, tV)
                aJx = np.append(aJx, tJ)
            # end for
            aPm             = np.min(aJx * aVx)                                 # nominal power in mW/cm2
            aPsolar         = 0.1 * SCC.SolarPower * SCC.SolarConcentration     # solar power, to convert from W/m2 to mW/cm2
            aEff[jj][ii]    = -100.0 * aPm / aPsolar                            # efficiency in percentage
            if (aEff[jj][ii] > EffMax):
                EffMax      = aEff[jj][ii]
                GapTopMax   = GapT
                GapBotMax   = GapB
            # end if
        except:
            pass
        # end try

    # end for

# end for

print("\ndone. elapsed time = %.3f sec." % float(time.time() - ticT))
treport = ("Double junction max efficiency = %.2f %% obtained for bandgaps = (%.3f ; %.3f) eV" % (EffMax, GapTopMax, GapBotMax))
print("\n----------------------------------------------------------------------\n" + treport + "\n----------------------------------------------------------------------\n")

fig = pl.figure(figsize=(10, 6), dpi=100, facecolor='#F1F1F1', linewidth=1.0, frameon=True)
fig.canvas.set_window_title('Double Junction Solar Cell')
ax  = fig.add_subplot(111)
pl.contourf(aGapBot, aGapTop, aEff, 4 * len(aGapBot))
cb  = pl.colorbar()
pl.plot(GapBotMax, GapTopMax, 'w+')
cb.ax.set_ylabel('Efficiency (%)')
ax.set_ylabel('Top Bandgap (eV)')
ax.set_xlabel('Bottom Bandgap (eV)')
pl.title('Double Junction Solar Cell Efficiency vs Bandgaps')
pl.show()