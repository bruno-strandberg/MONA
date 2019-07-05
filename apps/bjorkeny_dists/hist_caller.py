#!/usr/bin/env python
"""
This script can be used to send Create_BY_hists application jobs to the farm.

Usage:
    hist_caller.py -g GSG_FILE (--local | --lyon | --nikhef) [-n FMIN]
    hist_caller.py -h

Option:
   -g GSG_FILE gSeaGen file(s); select several with -g '*some_files*'
   --local     Process locally
   --lyon      Process on the cc-lyon farm
   --nikhef    Process on the nikhef farm
   -n FMIN     If processed on the farm, minimum number of gsg files per job [default: 100]
   -h --help   Show this sreen

"""

import os
from docopt import docopt
args = docopt(__doc__)

#=======================================================================
#*******************FUNCTION DEFINITIONS********************************
#=======================================================================

def create_syscmd(gsg_file_list, output_name):
    """
    This function creates the command to execute Create_BY_hists application
    """
    syscmd  = "{}/./Create_BY_hists".format(cwd)
    syscmd += " -g {}".format(gsg_file_list)
    syscmd += " -o {}".format(output_name)

    return syscmd

#=======================================================================

def create_farm_script(syscmd, tmpdir, job_nr):
    """
    This script creates a farm script that can be executed on the farm
    """

    scriptn = "{0}/farm_job_{1}.sh".format(tmpdir, job_nr)
    scriptf = open(scriptn, "w")
    scriptf.write('#!/bin/bash\n\n')

    scriptf.write( "cd {0}\n".format( os.environ['JPP_DIR'] ) )
    scriptf.write( "source setenv.sh\n" )
    scriptf.write( "cd {0}\n".format(cwd) )
    scriptf.write( "source {0}/bin/thisroot.sh\n".format( os.environ['ROOTSYS'] ) )
    scriptf.write( 'export OSCPROBDIR="{}"\n'.format(os.environ['OSCPROBDIR']) )
    scriptf.write( "source {0}/setenv.sh\n\n".format( os.environ['MONADIR'] ) )
    scriptf.write(syscmd + "\n")

    scriptf.close()

    return os.path.abspath(scriptn)

#=======================================================================
#*******************START OF SCRIPT*************************************
#=======================================================================

#=======================================================================
# Read in the list of gSeaGen files to be processed, create tmp and output dir
#=======================================================================
gsg_flist = os.popen( "ls {}".format(args['-g']) ).read().split()

cwd = os.getcwd()
os.system("mkdir -p {}/tmp".format(cwd) )
os.system("mkdir -p {}/output".format(cwd) )

#=======================================================================
# loop over gsg files and split them to arguments for Create_BY_hists.C
#=======================================================================
if args['--local']:
    args['-n'] = len(gsg_flist) + 1

job_index     = 0
files_per_job = 0
files_total   = 0
arg_list = []

for gsgfile in gsg_flist:
    
    if (files_per_job == 0):
        job_gsglist_name = "{0}/tmp/gsglist_job_{1}.dat".format(cwd, job_index)
        job_outname      = "{0}/output/output_job_{1}.root".format(cwd, job_index)
        job_gsglist      = open(job_gsglist_name, 'w')
        arg_list.append( (job_gsglist_name, job_outname) )

    job_gsglist.write(os.path.abspath(gsgfile) + "\n")
    files_per_job += 1
    files_total   += 1
    
    if ( ( files_per_job == int(args['-n']) ) or ( files_total == len(gsg_flist) ) ):
        job_gsglist.close()
        job_index += 1
        files_per_job = 0

#=======================================================================
# use the arguments to execute Create_BY_hists, either on the farm or locally
#=======================================================================

if args['--local']:
    syscmd = create_syscmd(arg_list[0][0], arg_list[0][1])
    print syscmd
    os.system(syscmd)

else:
    for i, arg_pair in enumerate(arg_list):
        syscmd = create_syscmd(arg_pair[0], arg_pair[1])
        scriptname = create_farm_script(syscmd, cwd+"/tmp", i)

        if args['--lyon']:
            syscmd = "qsub -P P_km3net -l sps=1 -l vmem=2G -l ct=5:00:00 -o {0} -e {0} {1}".format(cwd+"/tmp", scriptname)
        elif args['--nikhef']:
            syscmd = "qsub -q short7 -o {0} -e {0} {1}".format(cwd+"/tmp", scriptname)
        else:
            raise Exception("ERROR! hist_caller.py: unknown farming resource requested.")

        os.system( syscmd )
