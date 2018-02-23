#!/usr/local/bin/python -i
from ROOT import *

gROOT.ProcessLine(".L DataReducer.C+")
gROOT.ProcessLine("DataReducer dr")
gROOT.ProcessLine("dr.Loop()")
#gROOT.ProcessLine("dr.CreateSimpleClone()")
