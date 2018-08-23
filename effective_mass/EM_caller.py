#!/usr/bin/python
"""
This script can be used to call the EffMass.C macro.
 
Usage:
    EM_caller [--missing] [--nohadd]
    EM_caller -h                                                                     
                                                                                          
Option:        
    --missing         Only print summary info and missing EffMhists.C outputs, if any
    --nohadd          Do not re-combine EffMhists outputs, but use the ones already in combined_output/
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

def call_effmass(nohadd):

    flavors  = ["elec", "muon", "tau"]
    inters   = ["NC", "CC"]

    for f in flavors:
        for inter in inters:
            if (inter == "NC" and f != "elec"): continue

            combname = "combined_output/EffMhists_{0}_{1}.root".format(f, inter)
            outname  = "combined_output/EffMass_{0}_{1}.root".format(f, inter)

            # if hadd of EffMhists outputs not requested, check that the combined
            # output exists in combined_output/ dir
            if (nohadd):
                if ( not os.path.isfile(combname) ):
                    raise ValueError("File {} missing".format(combname) )
            else:
                os.system( "hadd {0} output/*{1}*{2}*".format(combname, f, inter) )


            os.system("root -b -q 'EffMass.C+({0}{1}{2},{0}{3}{2})'".format('"',combname,'"', outname))

#*****************************************************************************

def copy_to_datadir():

    filenames = ["combined_output/EffMass_elec_CC.root",
                 "combined_output/EffMass_muon_CC.root",
                 "combined_output/EffMass_tau_CC.root",
                 "combined_output/EffMass_elec_NC.root"]

    copy_data = raw_input("Update files in NMH/data/eff_mass/ (typically Y) [Y/N] ?")
    
    if copy_data == "Y":

        print("Updating files in NMH/data/eff_mass/")
        for f in filenames:
            os.system( "cp {} $NMHDIR/data/eff_mass/".format(f) )

    else:
        print("Not updating files in NMH/data/eff_mass/, if required perform manually.")

#*****************************************************************************
if __name__ == "__main__":
    args = docopt(__doc__)

    missing = get_missing_files()
    summarize()

    if ( args['--missing'] ):
        for m in missing:
            print m
    else:
        call_effmass(args['--nohadd'])
        copy_to_datadir()
