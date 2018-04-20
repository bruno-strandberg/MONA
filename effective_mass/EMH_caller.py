#!/usr/bin/python
"""
This script can be used to call the EffMhists.C macro.
 
Usage:
    EMH_caller -l RUNNR_LOW -u RUNNR_UP -f FLAVOUR... -i INTERACTION... -e E_START... (--local | --farm) [--nmin NMIN] [-c ATMMU_CUT] [-v VEFF_OPT] [--rvol RVOL] [--zmin ZMIN] [--zmax ZMAX]
    EMH_caller -h                                                                     
                                                                                          
Option:                                                                                   
    -l RUNNR_LOW      Lowest run number
    -u RUNNR_UP       Highest run number
    -f FLAVOUR        Neutrino flavour, 0 - e, 1 - mu, 2 - tau, may select several
    -i INTERACTION    Interaction type, 0 - nc, 1 - cc, may select several
    -e E_START        Energy start region, 1 or 3, may select several
    --local           Run locally
    --farm            Run on the farm
    --nmin NMIN       Minimum nr of files to analyse per farm job [default: 190]
    -h --help         Show this screen

    ==================FINE TUNING BELOW========================================

    -c ATMMU_CUT      PID cut to reject atmospheric muons [default: 1.0]
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
nmhdir        = os.environ['NMHDIR']
summary_dir   = nmhdir + "/data/mc_end/data_atmnu/"    #summary files dir
gseagen_dir   = nmhdir + "/data/mc_start/data_atmnu/"  #gseagen files on sps dir
irodsdir      = '/in2p3/km3net/mc/atm_neutrino/KM3NeT_ORCA_115_23m_9m/v1.0/gSeaGen/' #gsg@irods
irodsfiles    = os.popen( "ils {}".format(irodsdir) ).read().split()[1:] #list of gsg@irods
to_be_fetched = []                 #list of gsgfiles set to be fetched from IRODS (for farming)
corruptfiles  = ['elec-CC_1-5GeV_251-300.root.tar.gz'] #list of corrupt files

#*****************************************************************

def execute_effmass_calc(args):
    """This function calls the rootscript EffMhists.C.

    args is an argument dictionary created by docopt.
    """
    
    flavours     = { 0: "elec", 1: "muon", 2: "tau" }
    interactions = { 0: "NC", 1: "CC" }
    energies     = { 1: "1-5", 3: "3-100" }
    cmds   = []
    nfiles = 0
    njobs  = 0

    for flavour in args['-f']:
      for energy in args['-e']:
        for interaction in args['-i']:

          f = flavours[ int(flavour) ]
          i = interactions[ int(interaction) ]
          e = energies[ int(energy) ]
                  
          for fnr in range( int(args['-l']), int(args['-u']) + 1 ):
            
            summaryf = "{0}summary_{1}-{2}_{3}GeV_{4}.root".format(summary_dir, f, i, e, fnr)
            gseagenf = "{0}gSeaGen_{1}-{2}_{3}GeV_{4}".format(gseagen_dir, f, i, e, fnr)
            effmassf = "output/EffMhists_{0}-{1}_{2}GeV_{3}.root".format(f, i, e, fnr)
            
            if ( not os.path.isfile(summaryf) ):
                if args['--local']:
                    print("File {} missing, continuing.".format(summaryf))
                continue
            
            # gSeaGen file extension determined by the function below
            gsgexists, extension, gsg_cmds = gseagen_file_exists(gseagenf, f, i, e, fnr)
            gseagenf += extension
            cmds.extend(gsg_cmds)

            if ( not gsgexists ):
                print ("File {} missing, continuing.".format(gseagenf))
                continue

            cmds.append( get_execution_cmd(summaryf, gseagenf, flavour, interaction, 
                                           energy, fnr, args['-c'], args['-v'], 
                                           args['--rvol'], args['--zmin'], args['--zmax']) )
            
            # if run locally execute all commands and clear list
            if args['--local']:
                for cmd in cmds:
                    os.system(cmd)
                cmds = list()
            
            # running on the farm
            else:
            
                # send to farm when there are more than 100 files to be analysed and
                # when file nr is 50, 100, ... This condition avoids conflicts for fetching
                # from IRODS
                if (fnr % 50 == 0) and nfiles >= int(args['--nmin']):
                    submit_farm_job(cmds, njobs)
                    cmds = list()
                    nfiles = 0
                    njobs += 1
                else:
                    nfiles += 1
                
            #end loop over file numbers

        #end loop over energies
                
    #end loop over flavors

    # if there are commands left, execute them locally
    for cmd in cmds:
        os.system(cmd)

#*****************************************************************

def gseagen_file_exists(gseagenf, flav, interact, en, fnr):
    """
    This function checks whether the gSeaGen file exists in ../data/mc_start/...
    If the file does not exist, it looks for it in IRODS. Returns true if file
    exists locally or on IRODS, false otherwise. Returns file extension. Returns the system
    command to be executed to copy the file from IRODS.
    """

    global to_be_fetched
    extension = ""
    syscmd = []

    # if file exists locally (either in .root or .evt) return true

    if ( os.path.isfile(gseagenf+".root") ):
        extension = ".root"
        return True, extension, syscmd

    if ( os.path.isfile(gseagenf+".evt") ):
        extension = ".evt"
        return True, extension, syscmd
            
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
        return False, extension, syscmd

    # check that file is not corrupt, add to to_be_fetched list

    if (tarname in corruptfiles):
        print( "File {} marked as corrupt in EMH_caller.py".format(tarname) )
        return False, extension, syscmd

    if (tarname in to_be_fetched):
        return True, extension, syscmd

    # create and return the system command
    syscmd.append( "timeout 600 iget {0}{1} {2} -P\n".format(irodsdir, tarname, gseagen_dir) )
    syscmd.append( "tar xvzf {0}{1} -C {0}\n".format(gseagen_dir, tarname) )
    syscmd.append( "rm {0}{1}\n".format(gseagen_dir, tarname) )
    to_be_fetched.append( tarname )

    return True, extension, syscmd

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

def submit_farm_job(cmds, njobs):

    """
    This function creates a script and sends it to the farm.
    """

    os.system("mkdir -p tmp/")
    scriptn = "tmp/farm_job_{}.sh".format(njobs)
    scriptf = open(scriptn, "w")
    scriptf.write('#!/bin/bash\n\n')

    cwd = os.getcwd()
    scriptf.write( "cd /sps/km3net/users/bstrand/Code/Jpp/dev_rootv6\n" )
    scriptf.write( "source setenv.sh\n" )
    scriptf.write( "cd {}\n".format(cwd) )
    scriptf.write( "source /afs/in2p3.fr/home/b/bstrand/software/cern/root_v6.12.06/bin/thisroot.sh\n" )
    scriptf.write( "source /usr/local/shared/bin/irods_env.sh\n" )
    scriptf.write( "source /sps/km3net/users/bstrand/Code/NMH/setenv.sh\n\n" )
    
    for cmd in cmds:
        scriptf.write(cmd+"\n")

    scriptf.close()

    os.system( "qsub -P P_km3net -l sps=1 -o {0} -e {0} {1}".format(cwd+"/tmp", scriptn) )

#*****************************************************************

if __name__ == "__main__":

    args = docopt(__doc__)
    execute_effmass_calc(args)
