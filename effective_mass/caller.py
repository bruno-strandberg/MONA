#!/usr/bin/python
"""
This script can be used to call the EffMhists.C macro.
 
Usage:
    caller -l RUNNR_LOW -u RUNNR_UP -f FLAVOUR... -i INTERACTION -e E_START... [-c ATMMU_CUT] [-v VEFF_OPT] [--rvol RVOL] [--zmin ZMIN] [--zmax ZMAX] [--combine] [--combstr=<cs>]
    caller -h                                                                     
                                                                                          
Option:                                                                                   
    -l RUNNR_LOW      Lowest run number
    -u RUNNR_UP       Highest run number
    -f FLAVOUR        Neutrino flavour, 0 - e, 1 - mu, 2 - tau, may select several
    -i INTERACTION    Interaction type, 0 - nc, 1 - cc
    -e E_START        Energy start region, 1 or 3, may select several
    --combine         Combine the output to one file in combined_output. Flavors are kept separate.
    --combstr=<cs>    Identifier string added to the combined output file [default: ]
    -h --help         Show this screen

    ==================FINE TUNING BELOW========================================

    -c ATMMU_CUT      PID cut to reject atmospheric muons [default: 0.05]
    -v VEFF_OPT       Vgen option, 0 - interaction volume, 1 - can, 2 - custom [default: 1]
    --rvol RVOL       Radius of custom volume [default: 0]
    --zmin ZMIN       Z minimum of custom volume [default: 0]
    --zmax ZMAX       Z maximum of custom volume [default: 0]

"""

import sys
import os
import math
from docopt import docopt

#*****************************************************************
# some global directories
cwd = os.getcwd()
nmhdir = os.environ['NMHDIR']
summary_dir = nmhdir + "/data/mc_end/data_atmnu/"
gseagen_dir = nmhdir + "/data/mc_start/data_atmnu/"
irodsdir    = '/in2p3/km3net/mc/atm_neutrino/KM3NeT_ORCA_115_23m_9m/v1.0/gSeaGen/'
irodsfiles  = []

#*****************************************************************

def execute_effmass_calc(args):
    """This function calls the rootscript EffMhists.C.

    args is an argument dictionary created by docopt.
    """
    
    flavours     = { 0: "elec", 1: "muon", 2: "tau" }
    interactions = { 0: "NC", 1: "CC" }
    energies     = { 1: "1-5", 3: "3-100" }

    for flavour in args['-f']:

        effmass_files = []
        
        for energy in args['-e']:

            f = flavours[ int(flavour) ]
            i = interactions[ int(args['-i']) ]
            e = energies[ int(energy) ]
                        
            for fnr in range( int(args['-l']), int(args['-u']) + 1 ):
                
                summaryf = "{0}summary_{1}-{2}_{3}GeV_{4}.root".format(summary_dir, f, i, e, fnr)
                gseagenf = "{0}gSeaGen_{1}-{2}_{3}GeV_{4}".format(gseagen_dir, f, i, e, fnr)
                effmassf = "output/EffMhists_{0}-{1}_{2}GeV_{3}.root".format(f, i, e, fnr)
                
                if ( not os.path.isfile(summaryf) ):
                    print("File {} missing, continuing.".format(summaryf))
                    continue
            
                # gSeaGen file extension determined by the function below
                gsgexists, extension = gseagen_file_exists(gseagenf, f, i, e, fnr)
                gseagenf += extension
                if ( not gsgexists ):
                    print ("File {} missing, continuing.".format(gseagenf))
                    continue

                exec_cmd = get_execution_cmd(summaryf, gseagenf, flavour, args['-i'], energy, fnr,
                                             args['-c'], args['-v'],
                                             args['--rvol'], args['--zmin'], args['--zmax'])

                os.system( exec_cmd )
                effmass_files.append(effmassf)

            #end loop over file numbers

        #end loop over energies
        
        #combine runs (if there are runs to combine)
        if ( args['--combine'] and len(effmass_files) > 0 ):

            comb_name = "combined_output/EffMhists_{0}-{1}{2}.root".format(f, i, args['--combstr'])

            syscmd = "hadd {} ".format(comb_name)
            for outf in effmass_files:
                syscmd += outf + " "

            os.system(syscmd)
            os.system("root -b -q 'EffMass.C({0}{1}{2})'".format('"',comb_name,'"'))
        
    #end loop over flavors

#*****************************************************************

def gseagen_file_exists(gseagenf, flav, interact, en, fnr):
    """
    This function checks whether the gSeaGen file exists in ../data/mc_start/...

    If the file does not exist, it tries to fetch it from IRODS. Returns true if file
    exists or file sucessfully retrieved, returns the file extension.
    """

    global irodsfiles
    extension = ""

    # if file exists locally (either in .root or .evt) return true

    if ( os.path.isfile(gseagenf+".root") ):
        extension = ".root"
        return True, extension

    if ( os.path.isfile(gseagenf+".evt") ):
        extension = ".evt"
        return True, extension
    
    # first time populate the list with all files on IRODS

    if len(irodsfiles) == 0:
        irodsfiles = os.popen( "ils {}".format(irodsdir) ).read().split()[1:]
        irodsfiles.append("end_of_files") #to avoid repeating if IRODS is not available
        
    # check if files exist in either .root or .evt format on IRODS

    runspertar = 50
    fnr_max = int( math.ceil( float(fnr)/runspertar ) )*50
    fnr_min = fnr_max - 50 + 1
    tar1 = "{0}-{1}_{2}GeV_{3}-{4}.root.tar.gz".format(flav, interact, en, fnr_min, fnr_max)
    tar2 = "{0}-{1}_{2}GeV_{3}-{4}.evt.tar.gz".format(flav, interact, en, fnr_min, fnr_max)
    tarname = ""

    if tar1 in irodsfiles:
        tarname = tar1
        extension = ".root"
    elif tar2 in irodsfiles:
        tarname = tar2
        extension = ".evt"
    else:
        print( "Cannot find file {} in IRODS".format(gseagenf) )
        return False, extension

    # copy over the file from IRODS and untar

    curdir = os.getcwd()
    os.chdir(gseagen_dir)
    os.system( "iget {0}{1} -P".format(irodsdir, tarname) )
    os.system("tar xvzf {}".format(tarname) )
    os.system("rm {}".format(tarname) )
    os.chdir(curdir)

    # check if file found

    return os.path.isfile(gseagenf+extension), extension


#*****************************************************************
        
def get_execution_cmd(sname, gname,
                      flavour, interaction, estart, runnr,
                      atmmu_cut, veff_o,
                      rvol, zmin, zmax):
    """
    This function creates the execution command for the root script.
    """

    
    syscmd  = "root -b -q 'EffMhists.C+("
    syscmd += '"' + sname + '",'
    syscmd += '"' + gname + '",'
    syscmd += str(flavour) + ","
    syscmd += str(interaction) + ","
    syscmd += str(estart) + ","
    syscmd += str(runnr) + ","
    syscmd += str(atmmu_cut) + ","
    syscmd += str(veff_o) + ","
    syscmd += str(rvol) + ","
    syscmd += str(zmin) + ","
    syscmd += str(zmax)
    syscmd += ")'"

    return syscmd

#*****************************************************************

if __name__ == "__main__":

    args = docopt(__doc__)
    execute_effmass_calc(args)
