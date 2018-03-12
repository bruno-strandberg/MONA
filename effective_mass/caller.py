#!/usr/bin/python
"""                                                                                                 

This script can be used to call the EffectiveMass.C macro.
 
Usage:
    caller -l RUNNR_LOW -u RUNNR_UP -f FLAVOUR... -i INTERACTION -e E_START... -v VEFF_OPT
    caller -h                                                                     
                                                                                          
Option:                                                                                   
    -l RUNNR_LOW      Lowest run number
    -u RUNNR_UP       Highest run number
    -f FLAVOUR        Neutrino flavour, 0 - e, 1 - mu, 2 - tau, may select several
    -i INTERACTION    Interaction type, 0 - nc, 1 - cc
    -e E_START        Energy start region, 1 or 3, may select several
    -v VEFF_OPT       Effective volume calculation option, 0 - interaction volume, 1 - can, 2 - Agen
    -h --help         Show this screen

"""

import sys
import os
import math
from docopt import docopt

#*****************************************************************
# some global directories
summary_dir = "/sps/km3net/users/bstrand/data/NMH/data/mc_end/data_atmnu/"
gseagen_dir = "/sps/km3net/users/bstrand/data/NMH/data/mc_start/data_atmnu/"

#*****************************************************************

def execute_effmass_calc(args):

    flavours     = { 0: "elec", 1: "muon", 2: "tau" }
    interactions = { 0: "NC", 1: "CC" }
    energies     = { 1: "1-5", 3: "3-100" }

    for flavour in args['-f']:

        for energy in args['-e']:

            f = flavours[ int(flavour) ]
            i = interactions[ int(args['-i']) ]
            e = energies[ int(energy) ]
            
            if (f == "tau"):
                print("Tau parsing not yet supported.")
                raise Exception("Exiting.")
            
            for fnr in range( int(args['-l']), int(args['-u']) + 1 ):
                
                summaryf = "{0}summary_{1}-{2}_{3}GeV_{4}.root".format(summary_dir, f, i, e, fnr)
                gseagenf = "{0}gSeaGen_{1}-{2}_{3}GeV_{4}.root".format(gseagen_dir, f, i, e, fnr)
            
                if ( not os.path.isfile(summaryf) ):
                    print("File {} missing, continuing.".format(summaryf))
                    continue
            
                if ( not gseagen_file_exists(gseagenf, f, i, e, fnr) ):
                    print ("File {} missing, continuing.".format(gseagenf))
                    continue
            
                exec_cmd = get_execution_cmd(summaryf, gseagenf, flavour, args['-i'], energy, 
                                             fnr, args['-v'])
            
                os.system( exec_cmd )
        

#*****************************************************************

def gseagen_file_exists(gseagenf, flav, interact, en, fnr):

    filefound = False

    if ( os.path.isfile(gseagenf) ):
        filefound = True
    
    #try to copy the relevant file over from irods, untar and remove the tar
    else:
        runspertar = 50
        fnr_max = int( math.ceil( float(fnr)/runspertar ) )*50
        fnr_min = fnr_max - 50 + 1
        irodsdir = '/in2p3/km3net/mc/atm_neutrino/KM3NeT_ORCA_115_23m_9m/v1.0/gSeaGen/'
        tarname = "{0}-{1}_{2}GeV_{3}-{4}.root.tar.gz".format(flav, interact, en, fnr_min, fnr_max)
        curdir = os.getcwd()

        os.chdir(gseagen_dir)
        os.system( "iget {0}{1} -P".format(irodsdir, tarname) )
        os.system("tar xvzf {}".format(tarname) )
        os.system("rm {}".format(tarname) )
        os.chdir(curdir)

        if ( os.path.isfile(gseagenf) ):
            filefound = True
        else:
            filefound = False

    return filefound

#*****************************************************************
        
def get_execution_cmd(sname, gname, flavour, interaction, estart, runnr, veff_o):
    
    syscmd  = "root -b -q 'EffectiveMass.C+("
    syscmd += '"' + sname + '",'
    syscmd += '"' + gname + '",'
    syscmd += str(flavour) + ","
    syscmd += str(interaction) + ","
    syscmd += str(estart) + ","
    syscmd += str(runnr) + ","
    syscmd += str(veff_o)
    syscmd += ")'"

    return syscmd

#*****************************************************************

if __name__ == "__main__":

    args = docopt(__doc__)
    execute_effmass_calc(args)
