#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ======================================================================================================
# Solar Cell Shockley-Queisser Limit Calculator
# Code written by:
#   Sidi Hamady
#   Université de Lorraine, France
#   sidi.hamady@univ-lorraine.fr
# See Copyright Notice in COPYRIGHT
# HowTo in README.md
# https://github.com/sidihamady/Shockley-Queisser
# http://www.hamady.org/photovoltaics/ShockleyQueisser.zip
# ======================================================================================================

# ShockleyQueisser.py
#   implements the program interface used to start he calculator
#       execute ShockleyQueisser.py from the command line prompt by typing:
#           python -u ShockleyQueisser.py
#       or by double clicking on it (depending on the operating system settings)
#       or from within your editor, if possible.
#       in the graphical interface, change the parameters you want and press 'Calculate'.

# import the program core class in ShockleyQueisserCore.py
from ShockleyQueisserCore import *

# 1. create an instance of the core class
SCC = ShockleyQueisserCore(verbose = False)

# 2. change the parameters and calculate
SCC.calculate(
        TargetBandgap           = 1.1,                          # Target bandgap in eV (from 0.2 eV to 6 eV): to compare with the Shockley-Queisser Limit
        TargetBandgapTop        = 0.0,                          # Target top bandgap: used to take into account the part of solar spectrum already absorbed (for example in a top cell).
                                                                #   useful to calculate the overall efficiency in a multijunction solar cell.
                                                                #   For example for double junction solar cell, follow the steps below:
                                                                #       1. set TargetBandgap to 1.65 eV and the TargetBandgapTop to 0, and calculate the corresponding efficiency and current-voltage characteristic.
                                                                #       2. set TargetBandgap to 0.95 eV and the TargetBandgapTop to 1.65, and calculate the corresponding efficiency and current-voltage characteristic.
                                                                #       deduce from the previous data the overall double junction solar cell efficiency.
                                                                #       examples are given in ShockleyQueisserTJ.py and ShockleyQueisserDJ.py.
        Temperature             = 300.0,                        # Temperature in Kelvin (from 100 K to 700 K)
        SolarConcentration      = 1.0,                          # Solar concentration (1 sun to 1000 suns)
        OutputFilename          = './ShockleyQueisserOutput',   # Output file name without extension (used to save figure in PDF format if using GUI, and the text output data). Set to None to disable.
        useGUI                  = True                          # The calculator can be used in graphical (GUI) mode or command-line only mode. In command-line mode (useGUI = False) the results are printed out and saved in text files.
                                                                #   the command-line mode is useful to perform specific calculations such as multijunction solar cell efficiency (see example in README.md).
        )

# main output data:
# SCC.SQ_Efficiency         SCC.SQ_VOC          SCC.SQ_JSC          SCC.SQ_FF           SCC.SQ_Voltage (array)          SCC.SQ_Current (array)
# SCC.Target_Efficiency     SCC.Target_VOC      SCC.Target_JSC      SCC.Target_FF       SCC.Target_Voltage (array)      SCC.Target_Current (array)