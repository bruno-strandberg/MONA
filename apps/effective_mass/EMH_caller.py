#!/usr/bin/python
"""
This script can be used to call the EffMhists.C macro for processing many files with default settings (useful for ORCA analyses). For each file in SUMMARYDIR, it will look for a file in GSGDIR; if the file is not found an error is raised. If the gSeaGen file is found, the script sets up a call to EffMhists program.
 
Usage:
    EMH_caller --gsgdir GSGDIR --summarydir SUMMARYDIR [--odir OUTPUTDIR] [--rmin RMIN] [--rmax RMAX] [--selstr SELSTR] [--nmin NMIN] (--local | --farm) 
    EMH_caller -h                                                                     
                                                                                          
Option:
    --gsgdir GSGDIR             Local gSeaGen directory
    --summarydir SUMMARYDIR     Local summary files directory
    --odir OUTPUTDIR            Output dir [default: output/]
    --rmin RMIN                 Lowest run number [default: 1]
    --rmax RMAX                 Highest run number [default: 10000]
    --selstr SELSTR             Selection string for files in SUMMARYDIR, e.g. '*muon-CC*' (note the quotation marks!) [default: ]
    --nmin NMIN                 Minimum nr of files to analyse per farm job [default: 190]
    --local                     Run locally
    --farm                      Run on the farm
    -h --help                   Show this screen

"""

from docopt import docopt
args = docopt(__doc__)

import sys
import os
import math

#*****************************************************************
# some global directories
nmhdir        = os.environ['NMHDIR']
summary_dir   = os.path.abspath(args['--summarydir'])
gseagen_dir   = os.path.abspath(args['--gsgdir'])

#*****************************************************************

def execute_effmass_calc(args):
    """This function calls the rootscript EffMhists.C.

    args is an argument dictionary created by docopt.
    """

    summaryfiles = os.popen( "ls {}".format(args['--summarydir']) ).read().split()
    gseagenfiles = os.popen( "ls {}".format(args['--gsgdir']) ).read().split()

    flavors      = ['elec', 'muon', 'tau']
    interactions = ['NC','CC']

    for sf in summaryfiles:
        
        # extract flavor, NC/CC, energy range and file number from summary file name
        # and look for the corresponding gSeaGen file

        flav = inter = erange = fnr = ''

        flav  = [f for f in flavors if f in sf]
        inter = [i for i in interactions if i in sf]
        erange = sf[ sf[0:sf.index("GeV")].rfind('_')+1 : sf.index("GeV") ]
        fnr = sf[ sf.rfind('_') : sf.rfind('.')+1 ]

        if ( len(flav) != 1 or len(inter) != 1 or flav[0] == '' or inter[0] == '' or erange == '' or fnr == '' ):
            raise Exception("Flavor and interaction extraction failed!")
        
        gsgfile = [g for g in gseagenfiles if (flav[0] in g and inter[0] in g and erange in g and fnr in g)]

        if (len(gsgfile) != 1):
            raise Exception( "Could not find gSeaGen file for summary file {}".format(sf) )

    quit()

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
                                           energy, fnr, args['-c'], args['-n'], args['-v'], 
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
                      atmmu_cut, noise_cut, veff_o,
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
    syscmd += str(noise_cut) + ","
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
    scriptf.write( "cd {0}\n".format( os.environ['JPP_DIR'] ) )
    scriptf.write( "source setenv.sh\n" )
    scriptf.write( "cd {0}\n".format(cwd) )
    scriptf.write( "source {0}/bin/thisroot.sh\n".format( os.environ['ROOTSYS'] ) )
    scriptf.write( "source /usr/local/shared/bin/irods_env.sh\n" )
    scriptf.write( "source {0}/setenv.sh\n\n".format( os.environ['NMHDIR'] ) )
    
    for cmd in cmds:
        scriptf.write(cmd+"\n")

    scriptf.close()

    os.system( "qsub -P P_km3net -l sps=1  -l vmem=2G -l ct=5:00:00 -o {0} -e {0} {1}".format(cwd+"/tmp", scriptn) )

#*****************************************************************

if __name__ == "__main__":

    args = docopt(__doc__)
    execute_effmass_calc(args)
