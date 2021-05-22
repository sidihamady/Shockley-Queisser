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

# ShockleyQueisserCurve.py
#   calculates and plots the efficiency of a PN junction solar cell with respect to the bandgap

# import the program core class
from ShockleyQueisserCore import *

# create an instance of the core class
SCC = ShockleyQueisserCore(verbose = False)

print("\ncalculating...")

# set useGUI to False to use the command-line mode
SCC.calculate(
        TargetBandgap           = 1.1,
        TargetBandgapTop        = 0.0,
        Temperature             = 300.0,
        SolarConcentration      = 1.0,
        OutputFilename          = None,
        useGUI                  = False
        )

# plot the efficiency-bandgap curve
fig         = pl.figure(figsize=(10, 6), dpi=100, facecolor='#F1F1F1', linewidth=1.0, frameon=True)
fig.canvas.set_window_title('Shockley-Queisser Efficiency vs Bandgap')
ax          = fig.add_subplot(111)
tline,      = ax.plot(SCC.Bandgap, SCC.Efficiency, '-', linewidth=3.0)
tline.set_color('b')
ax.set_xlabel('$Bandgap\ (eV)$',        fontsize=16)
ax.set_ylabel('$Efficiency\ (\%)$',     fontsize=16)
pl.xticks(fontsize=16)
pl.yticks(fontsize=16)
pl.title('Shockley-Queisser Efficiency vs Bandgap')
pdfT = PdfPages('ShockleyQueisserCurve.pdf')
pdfT.savefig(fig)
pdfT.close()
np.savetxt('ShockleyQueisserCurve.txt', np.c_[SCC.Bandgap, SCC.Efficiency], fmt='%.6g\t%.6g', delimiter=SCC.DataDelimiter, newline='\n', header='Shockley-Queisser Efficiency vs Bandgap\nBandgap (eV)\tEfficiency (%)')
pl.savefig('ShockleyQueisserCurve.png', dpi=600)
pl.grid(True)
pl.xlim([0.5, 3])
pl.ylim([5, 34])
pl.show()
#