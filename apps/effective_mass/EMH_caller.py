#!/usr/bin/python
"""
This script can be used to call the EffMhists program for processing many files with default settings (useful for ORCA analyses). For each file in SUMMARYDIR, it will look for a file in GSGDIR; if the file is not found the summary file is ignored. If the gSeaGen file is found, the script sets up a call to EffMhists program.
 
Alternatively, the script can take a directory with trigger files as argument, in which case the gSeaGen directory is ignored. It will look up a trigger file per each summary file and extract the input gSeaGen file(s) corresponding to the summary file by using JPringMeta of Jpp. It will check whether it can access/find the corresponding gSeaGen files, if that is the case it will create and call commands for EffMhists execution.

Usage:
    EMH_caller --summarydir SUMMARYDIR [--gsgdir GSGDIR] [--trigdir TRIGDIR] [--odir OUTPUTDIR] [--rmin RMIN] [--rmax RMAX] [--selstr SELSTR] [--nmin NMIN] (--local | --lyonfarm | --nikheffarm) 
    EMH_caller -h                                                                     
                                                                                          
Option:
    --summarydir SUMMARYDIR     Local summary files directory
    --gsgdir GSGDIR             Local gSeaGen directory [default: ]
    --trigdir TRIGDIR           Local trigger files directory [default: ]
    --odir OUTPUTDIR            Output dir [default: output/]
    --rmin RMIN                 Lowest run number [default: 1]
    --rmax RMAX                 Highest run number [default: 10000]
    --selstr SELSTR             Selection string for files in SUMMARYDIR, e.g. '*muon-CC*' (note the quotation marks!) [default: ]
    --nmin NMIN                 Minimum nr of files to analyse per farm job [default: 190]
    --local                     Run locally
    --lyonfarm                  Run on the lyon farm
    --nikeffarm                 Run on the nikhef farm
    -h --help                   Show this screen

"""

from docopt import docopt
args = docopt(__doc__)

import sys
import os
import math
import time
import re

#*****************************************************************

def execute_effmass_calc(args):
    """This function calls the function EffMhists

    args is an argument dictionary created by docopt.
    """

    #===========================================================
    # create the executable command for each selected file in summary dir
    #===========================================================

    summaryfiles = os.popen( "ls {}/{}".format(args['--summarydir'], args['--selstr']) ).read().split()

    gseagenfiles = []
    if (args['--gsgdir'] != ''):
        gseagenfiles = os.popen( "ls {}".format(args['--gsgdir']) ).read().split()

    triggerfiles = []
    if (args['--trigdir'] != ''):
        triggerfiles = os.popen( "ls {}".format(args['--trigdir']) ).read().split()

    flavors      = ['elec', 'muon', 'tau']
    interactions = ['NC','CC']

    jobcmds   = []

    for sf in summaryfiles:
        
        #--------------------------------------------------------------------------
        # extract flavor, NC/CC, energy range and file number from summary file name
        # and look for the corresponding gSeaGen file
        #--------------------------------------------------------------------------

        flav = inter = erange = fnr = ''

        flav  = [f for f in flavors if f in sf]
        inter = [i for i in interactions if i in sf]
        erange = sf[ sf[0:sf.index("GeV")].rfind('_')+1 : sf.index("GeV") ] + "GeV"
        fnr1 = sf[ sf.rfind('_') : sf.rfind('.')+1 ] + "evt" # file nr search str in format _123.evt
        fnr2 = "." + fnr1[1:]                                # file nr search str in format .123.evt

        if ( len(flav) != 1 or len(inter) != 1 or flav[0] == '' or inter[0] == '' or erange == '' or fnr1 == '' ):
            raise Exception("Flavor and interaction extraction failed!")

        # need to make an exception for tau's low-energy; summary file range 1-10, gSeaGen 3.4-10, jeez...
        if flav[0] == "tau" and float(erange.split('-')[0]) < 3.4:
            erange = erange[ erange.index('-'): ]

        #--------------------------------------------------------------------------
        # no trigger directory specified, look up a gSeaGen file for each summary file and
        # add an execution command to jobcmds list
        #--------------------------------------------------------------------------
        if len(triggerfiles) == 0:

            gsgfile = [g for g in gseagenfiles if ( flav[0] in g and inter[0] in g and erange in g 
                                                    and (fnr1 in g or fnr2 in g) ) ]
            
            if (len(gsgfile) != 1):
                print ( "WARNING! Could not find gSeaGen file for summary file {}".format(sf) )
                print ( "         search string: ;{};{};{};{};{};".format(flav[0], inter[0], erange, fnr1, fnr2) )
                print ( "         found files  : {}".format(gsgfile) )
                continue
            
            sumf = os.path.abspath(args['--summarydir']) + "/" + sf[ sf.rfind('/')+1 : ]
            gsgf = []
            gsgf.append( os.path.abspath(args['--gsgdir']) + "/" + gsgfile[0] )
            
            jobcmds.append( get_execution_cmd( sumf, gsgf, os.path.abspath(args['--odir']) ) )

        #--------------------------------------------------------------------------
        # trigger directory specified; look up trigger files corresponding to the summary file,
        # extract gSeaGen files and create commands
        #--------------------------------------------------------------------------

        else:
            
            fnr = fnr[:-1]+'-' # replace filenr search string to form _nr- instead of _nr.

            # fnr in searched in the part after GeV, otherwise the range _1-5 always matches for fnr=_1-
            trigfile = [ t for t in triggerfiles if ( flav[0] in t and inter[0] in t and erange in t and 
                                                      fnr in t[t.index("GeV"):] ) ]

            if (len(trigfile) != 1):
                print (trigfile)
                raise Exception("ERROR! Found the above list of trigger files, but should have 1 file exactly")

            sumf = os.path.abspath(args['--summarydir']) + "/" + sf[ sf.rfind('/')+1 : ]
            gsgf = get_gsg_filelist( os.path.abspath(args['--trigdir']) + "/" + trigfile[0] )

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
        elif args['--lyonfarm']:
            syscmd = "qsub -P P_km3net -l sps=1 -l ct=5:00:00 -o {0} -e {0} {1}".format(os.getcwd()+"/tmp", js)
        else:
            syscmd = "qsub -q short7 -o {0} -e {0} {1}".format(os.getcwd()+"/tmp", js)

        #print syscmd
        os.system(syscmd)

#*****************************************************************

def get_gsg_filelist(trigfile):
    """
    This function uses JPrintMeta to extract the list of gSeaGen files input to the trigger file
    """
    
    #-----------------------------------------------------------------------
    # read the meta-data
    #-----------------------------------------------------------------------
    metaout = os.popen( "JPrintMeta -f {}".format(trigfile) ).read()

    #-----------------------------------------------------------------------
    # file list starts with -f, ends whenever there is another argument (e.g. "-o "),
    # hence locate the next argument after the file-list and cut there
    #-----------------------------------------------------------------------

    flist = metaout[ metaout.index('-f ')+len('-f ') : ].split() # all text after '-f ', split at spaces
    regex = re.compile('-.')                                     # argument pattern "-*", where . stands for *
    flags = [str for str in flist if re.match(regex, str)]       # all arguments in the list
    
    # if no flags are found, raise error. This can actually happen if the -f filelist comes last,
    # but protect against that for now
    if len(flags) == 0:
        raise Exception("ERROR! get_gsg_filelist() could not find any more command line arguments")

    # finally, cut out everything after the first argument after the filelist
    flist = flist[ : flist.index(flags[0]) ]

    #-----------------------------------------------------------------------
    # check that all files exist
    #-----------------------------------------------------------------------
    for f in flist:
        if ( not os.path.isfile(f) ):
            raise Exception("ERROR! get_gsg_filelist() cannot find file {}".format(f))

    print ( "NOTICE found {} gSeaGen files input to trigger file {}".format(len(flist), trigfile) )

    return flist


#*****************************************************************
        
def get_execution_cmd(summaryfile, gseagenfiles, outputdir):
    """
    This function creates the execution command for the application.
    """
    
    syscmd  = os.getcwd() + "/./EffMhists "
    syscmd += "-d {} ".format(outputdir)
    syscmd += "-s {} ".format(summaryfile)

    for gsgf in gseagenfiles:
        syscmd += "-g {} ".format(gsgf)

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
    scriptf.write( "source {0}/setenv.sh\n\n".format( os.environ['MONADIR'] ) )
    
    for cmd in cmds:
        scriptf.write(cmd+"\n")

    scriptf.close()

    return scriptn

#*****************************************************************

if __name__ == "__main__":

    args = docopt(__doc__)
    execute_effmass_calc(args)
