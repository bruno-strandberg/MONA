#!/usr/bin/python

"""
Script to farm application `fit` for pseudo-experiment fitting.

Usage:
    farm_fit.py -t THETA23 [-s SEED] [--dm31prior] [--IO] [--WS] [-j NJOBS] [-d TMPDIR]
    farm_fit.py -h

Option:
    -t THETA23   Value for sinsqth23
    -s SEED      Seed for python script for generating random seeds for the `fit` application [default: 666]
    --dm31prior  Include a prior on dm31
    --IO         Inverted ordering data
    --WS         Include systematics
    -j NJOBS     Number of farm jobs [default: 100]
    -d TMPDIR    Temporary directory name where data is written [default: tmp-fit]
    -h --help    Show this screen
"""

import random
import os
from docopt import docopt
args = docopt(__doc__)

# create temporary directory
os.system( "mkdir -p {}/{}".format( os.getcwd(), args['-d'] ) )

# set seed for python random
random.seed( int(args['-s']) )

#========================================================
# create application execution commands
#========================================================

syscmds = []

for jnr in range ( 0, int(args['-j']) ):

    seed    = random.randint(1, 1e8)
    outfile = "{}/{}/fit_out_job_{}.root".format( os.getcwd(), args['-d'], jnr )
    syscmd  = "{}/./fit -t {} -S {} -o {}".format( os.getcwd(), args['-t'], seed, outfile )

    if args['--dm31prior']:
        syscmd += " -p"

    if args['--IO']:
        syscmd += " -i"

    if args['--WS']:
        syscmd += " -s"

    syscmds.append( [jnr, syscmd] )

#========================================================
# create job scripts
#========================================================

jobfiles = []

for syscmd in syscmds:

    jobfilename = "{}/{}/job_{}.sh".format(os.getcwd(), args['-d'], syscmd[0])
        
    jobfile = open(jobfilename, 'w')

    # copy lines from bashrc
    bashrcf = open("/user/bstrand/.bashrc", 'r')
    copylines = False
    for line in bashrcf:

        if "source" in line:
            copylines = True

        if copylines:
            jobfile.write(line)

    bashrcf.close()

    jobfile.write( syscmd[1] )
    jobfile.close()
    jobfiles.append( jobfilename )

#========================================================
# execute
#========================================================

for job in jobfiles:
    
    farmcmd = "qsub -q short7 -o {0} -e {0} {1}".format(os.getcwd()+"/" + args['-d'], job)
    os.system(farmcmd)
