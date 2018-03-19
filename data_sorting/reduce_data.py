#!/usr/bin/python
from ROOT import *

gROOT.ProcessLine(".L DataReducer.C+")
gROOT.ProcessLine("DataReducer dr")
gROOT.ProcessLine("dr.Loop()")

