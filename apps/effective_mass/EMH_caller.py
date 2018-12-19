#!/usr/bin/python
"""
This script can be used to call the EffMhists program for processing many files with default settings (useful for ORCA analyses). For each file in SUMMARYDIR, it will look for a file in GSGDIR; if the file is not found the summary file is ignored. If the gSeaGen file is found, the script sets up a call to EffMhists program.
 
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
import time

#*****************************************************************
# some global directories
nmhdir        = os.environ['NMHDIR']
summary_dir   = os.path.abspath(args['--summarydir'])
gseagen_dir   = os.path.abspath(args['--gsgdir'])

#*****************************************************************

def execute_effmass_calc(args):
    """This function calls the function EffMhists

    args is an argument dictionary created by docopt.
    """

    #===========================================================
    # create the executable command for each selected file in summary dir
    #===========================================================

    summaryfiles = os.popen( "ls {}/{}".format(args['--summarydir'], args['--selstr']) ).read().split()
    gseagenfiles = os.popen( "ls {}".format(args['--gsgdir']) ).read().split()

    flavors      = ['elec', 'muon', 'tau']
    interactions = ['NC','CC']

    jobcmds   = []

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
            print ("WARNING! Could not find gSeaGen file for summary file {}".format(sf))
            continue

        sumf = os.path.abspath(args['--summarydir']) + "/" + sf[ sf.rfind('/')+1 : ]
        gsgf = os.path.abspath(args['--gsgdir']) + "/" + gsgfile[0]

        jobcmds.append( get_execution_cmd( sumf, gsgf, os.path.abspath(args['--odir']) ) )

    #===========================================================
    # now bundle the commands into farm scripts
    #===========================================================
    
    cmdsperjob = []
    jobscripts = []

    for i,c in enumerate(jobcmds):

        cmdsperjob.append(c)

        if ( len(cmdsperjob) == int(args['--nmin']) or i == ( len(jobcmds)-1 ) ):
            jobscripts.append( create_farm_script(cmdsperjob, len(jobscripts)) )
            cmdsperjob = []

    #===========================================================
    # finally execute the jobs locally or on the farm
    #===========================================================
            
    for js in jobscripts:

        if args['--local']:
            syscmd = "bash {}".format(js)
        else:
            syscmd = "qsub -P P_km3net -l sps=1 -l ct=5:00:00 -o {0} -e {0} {1}".format(os.getcwd()+"/tmp", js)

        #print syscmd
        os.system(syscmd)

#*****************************************************************
        
def get_execution_cmd(summaryfile, gseagenfile, outputdir):
    """
    This function creates the execution command for the application.
    """
    
    syscmd  = os.getcwd() + "/./EffMhists "
    syscmd += "-d {} ".format(outputdir)
    syscmd += "-g {} ".format(gseagenfile)
    syscmd += "-s {} ".format(summaryfile)

    return syscmd

#*****************************************************************

def create_farm_script(cmds, jobnr):
    """
    This function creates a script to be sent to the farm
    """

    cwd = os.getcwd()
    os.system("mkdir -p tmp/")
    scriptn = "tmp/farm_job_{}_timestamp{}.sh".format(jobnr,time.time())

    scriptf = open(scriptn, "w")

    scriptf.write('#!/bin/bash\n\n')
    scriptf.write( "cd {0}\n".format( os.environ['JPP_DIR'] ) )
    scriptf.write( "source setenv.sh\n" )
    scriptf.write( "cd {0}\n".format(cwd) )
    scriptf.write( "source {0}/bin/thisroot.sh\n".format( os.environ['ROOTSYS'] ) )
    scriptf.write( 'export OSCPROBDIR="{}"\n'.format(os.environ['OSCPROBDIR']) )
    scriptf.write( "source {0}/setenv.sh\n\n".format( os.environ['NMHDIR'] ) )
    
    for cmd in cmds:
        scriptf.write(cmd+"\n")

    scriptf.close()

    return scriptn

#*****************************************************************

if __name__ == "__main__":

    args = docopt(__doc__)
    execute_effmass_calc(args)
