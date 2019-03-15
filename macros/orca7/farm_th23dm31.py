#!/usr/bin/python

"""
Script to farm theta23-dm31 pseudo-experiment fitting script.

Usage:
    farm_th23dm31.py [-n NFITS] [-j NJOBS] (--nikhef | --lyon)
    farm_th23dm31.py -h

Option:
    -n NFITS     Number of fits per job [default: 10]
    -j NJOBS     Number of farm jobs [default: 100]
    --nikhef     qsub command configured for running on nikhef farm
    --lyon       qsub command configured for running on lyon
    -h --help    Show this screen
"""

import os
from docopt import docopt
args = docopt(__doc__)

# get some necessary environment variables
jpp  = os.environ['JPP_DIR']
root = os.environ['ROOTSYS']
nmh  = os.environ['MONADIR']
oscp = os.environ['OSCPROBDIR']

# create temporary directory
curdir = os.getcwd()
os.system("mkdir -p {}/tmp".format(curdir))

# create farm jobs
for jnr in range ( 0, int(args['-j']) ):

    #----------------------------------------------
    # create root macro execution command
    #----------------------------------------------

    jobstr = '"job{}"'.format(jnr)
    expstr = '{}'.format(args['-n'])
    syscmd = "root -b -q '{}/th23dm31.C({},{})'".format(curdir,jobstr, expstr)

    #----------------------------------------------
    # create a bash script than can be sent to farm
    #----------------------------------------------

    script_name = "{0}/tmp/job_{1}.sh".format(curdir, jnr)
    script_file = open(script_name, 'w')

    script_file.write('#!/bin/bash\n\n')
    script_file.write('cd {0}\nsource setenv.sh\ncd {1}\n'.format(jpp, curdir)) #jpp
    script_file.write('source {}/bin/thisroot.sh\n'.format(root))               #root
    script_file.write("export OSCPROBDIR='{}'\n".format(oscp))                  #oscprob
    script_file.write('source {}/setenv.sh\n\n'.format(nmh))                    #nmh
    script_file.write(syscmd)
    
    script_file.close()

    #----------------------------------------------
    # create farm command
    #----------------------------------------------

    if args['--nikhef']:
        farm_cmd = "qsub -q short7 -o {0} -e {0} {1}".format(curdir+"/tmp", script_name)
    elif args['--lyon']:
        syscmd = "qsub -P P_km3net -l sps=1 -l ct=5:00:00 -o {0} -e {0} {1}".format(curdir+"/tmp", script_name)

    print ("Executing {}".format(farm_cmd))
    os.system(farm_cmd)
