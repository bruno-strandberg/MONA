#!/usr/bin/python
"""
This script can be used to call the EffMass.C macro.
 
Usage:
    EM_caller.py
    EMH_caller -h                                                                     
                                                                                          
Option:                                                                                   
    -h --help         Show this screen
"""

import sys
import os
from docopt import docopt

#*****************************************************************************
nmhdir        = os.environ['NMHDIR']
summary_dir   = nmhdir + "/data/mc_end/data_atmnu/"    #summary files dir
gseagen_dir   = nmhdir + "/data/mc_start/data_atmnu/"  #gseagen files on sps dir

#*****************************************************************************
def get_missing_files():
    """
    Returns the list of missing EffMhists files.
    """

    summary_files = os.popen( "ls {}".format(summary_dir) ).read().split()
    effmh_files   = os.popen( "ls output/" ).read().split()
    missing_files = []

    for sf in summary_files:
        name = "EffMhists" + sf[ sf.find("_"): ]

        if (name not in effmh_files):
            missing_files.append( name )

    return missing_files
            
#*****************************************************************************
def summarize():

    flavors  = ["elec", "muon", "tau"]
    inters   = ["NC", "CC"]
    energies = ["1-5", "3-100"]

    for f in flavors:
        for inter in inters:
            for en in energies:
                if ( inter == "NC" and f != "elec" ): continue
                if (en == "1-5" and f == "tau"): continue

                outlist = os.popen( "ls output/*{0}*{1}*{2}* 2>/dev/null".
                                    format(f, inter, en) ).read().split()
                outname = "{0} {1} {2}".format(f, inter, en)
                print( "{0}\t: {1} files".format(outname, len(outlist)) )
                
#*****************************************************************************

def call_effmass():

    flavors  = ["elec", "muon", "tau"]
    inters   = ["NC", "CC"]

    for f in flavors:
        for inter in inters:
            if (inter == "NC" and f != "elec"): continue

            if (inter == "NC"):
                print("DEBUG: No NC files, skipping.")
                continue

            combname = "combined_output/EffMhists_{0}_{1}.root".format(f, inter)
            outname  = "combined_output/EffMass_{0}_{1}.root".format(f, inter)

            os.system( "hadd {0} output/*{1}*{2}*".format(combname, f, inter) )
            os.system("root -b -q 'EffMass.C({0}{1}{2},{0}{3}{2})'".format('"',combname,'"', outname))

#*****************************************************************************
if __name__ == "__main__":
    missing = get_missing_files()
    summarize()
    call_effmass()

    #for m in missing:
    #    print m
