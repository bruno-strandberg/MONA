#!/usr/bin/env python
"""
This script can be used to check that RestoreParity worked correctly.
 
Usage:
    data_sorting.py -p PID_FILE [-s SUMMARY_FILES_DIR]
    data_sorting.py -h                                                                     
                                                                                          
Option:
    -p PID_FILE           PID file from ECAP, e.g. NMHDIR/data/pid_result_10Apr2018.root
    -s SUMMARY_FILES_DIR  PID_FILE split to summary files by scripts in data_sorting/ [default: $NMHDIR/data/mc_end/data_atmnu]
    -h --help         Show this screen

"""

from docopt import docopt
from ROOT import *
import os

#*****************************************************************************
def summarize():
    """
    Function that summarizes the numbers of summary files.
    """
    
    flavors  = ["elec", "muon", "tau"]
    inters   = ["NC", "CC"]
    energies = ["1-5", "3-100"]

    for f in flavors:
        for inter in inters:
            for en in energies:
                if ( inter == "NC" and f != "elec" ): continue
                if (en == "1-5" and f == "tau"): continue

                outlist = os.popen( "ls $NMHDIR/data/mc_end/data_atmnu/*{0}*{1}*{2}* 2>/dev/null".format(f, inter, en) ).read().split()
                outname = "{0} {1} {2}".format(f, inter, en)
                print( "{0}\t: {1} files".format(outname, len(outlist)) )

#*****************************************************************************

args = docopt(__doc__)
gSystem.Load("$NMHDIR/common_software/libnmhsoft.so")
gROOT.ProcessLine(".L $NMHDIR/data_sorting/DataReducer.C+")

summarize()

#=========================================
# define comparison histograms
#=========================================
gandalf_EvsTh_pid = TH2D("gandalf_EvsTh_pid","gandalf_EvsTh_pid",50,0,100,40,-1,1)
shower_EvsTh_pid  = TH2D("shower_EvsTh_pid" ,"shower_EvsTh_pid" ,50,0,100,40,-1,1)
gandalf_EvsTh_sum = TH2D("gandalf_EvsTh_sum","gandalf_EvsTh_sum",50,0,100,40,-1,1)
shower_EvsTh_sum  = TH2D("shower_EvsTh_sum" ,"shower_EvsTh_sum" ,50,0,100,40,-1,1)

#=========================================
# loop over events in PID tree, fill histogram
#=========================================
pidfile   = TFile(args['-p'],"READ")
pidtree   = pidfile.Get("PID")
pp        = DataReducer(pidtree)

print("Looping PID tree:")
for i in range( 0, pp.fChain.GetEntries() ):
    if ( i % 100000 == 0 ): print("Event {}".format(i))
    pp.fChain.GetEntry(i)
    if (pp.is_neutrino == 0): continue
    gandalf_EvsTh_pid.Fill(pp.gandalf_energy_corrected, pp.gandalf_dir_z)
    shower_EvsTh_pid.Fill(pp.dusj_energy_corrected, pp.dusj_dir_z)

#=========================================
# loop over summary events, fill histogram
#=========================================
summaryfiles = os.popen( "ls {}/summary_*".format(args['-s']) ).read().split()

print("Looping over summary files:")
for fnr,sf in enumerate(summaryfiles):
    if (fnr%100==0): print("File {}".format(fnr) )
    sp = SummaryParser(sf)
    for i in range( sp.fChain.GetEntries() ):
        sp.fChain.GetEntry(i)
        gandalf_EvsTh_sum.Fill(sp.gandalf_energy_nu, sp.gandalf_dir_z)
        shower_EvsTh_sum.Fill(sp.shower_energy_nu, sp.shower_dir_z)
        
fout = TFile("fout.root","RECREATE")
gandalf_EvsTh_pid.Write()
shower_EvsTh_pid.Write()
gandalf_EvsTh_sum.Write()
shower_EvsTh_sum.Write()
fout.Close()
