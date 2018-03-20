#!/usr/bin/python
"""
This script calls DataReducer to reduce NMH/data/pid_result_XXX.root to analysis format.
 
Usage:
   ./caller.py

"""

from ROOT import *

gROOT.ProcessLine(".L DataReducer.C+")
gROOT.ProcessLine("DataReducer dr")
gROOT.ProcessLine("dr.Loop()")

