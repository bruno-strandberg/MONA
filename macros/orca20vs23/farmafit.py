#!/usr/bin/python

"""
Script to farm asimovfit.C macro. At each theta-23 value the script calls the asimovfit command for ORCA20 and ORCA23 detectors.

Usage:
    farm_th23dm31.py --min TH23_MIN --max TH23_MAX --step TH23_STEP (--nikhef | --lyon)
    farm_th23dm31.py -h

Option:
    --min TH23_MIN     Start value of sin^2(th23)
    --max TH23_MAX     Stop value of sin^2(th23)
    --step TH23_STEP   Step size to go from min to max
    --nikhef           qsub command configured for running on nikhef farm
    --lyon             qsub command configured for running on lyon
    -h --help          Show this screen
"""

import os
from docopt import docopt
args = docopt(__doc__)

# get some necessary environment variables
jpp  = os.environ['JPP_DIR']
root = os.environ['ROOTSYS']
nmh  = os.environ['NMHDIR']
oscp = os.environ['OSCPROBDIR']

# create temporary directory
curdir = os.getcwd()
os.system("mkdir -p {}/tmp".format(curdir))

#==============================================================
# create commands to be executed on the farm
#==============================================================
jobcmds = []

th23 = float( args['--min'] )
while th23 <= float( args['--max'] ):

    outname_O20 = '"{}/rootfiles/out_ORCA20_th23_{}.root"'.format(curdir, th23)
    outname_O23 = '"{}/rootfiles/out_ORCA23_th23_{}.root"'.format(curdir, th23)

    jobcmd_O20 = "{}/./callfit -d ORCA20 -t {} -o {}".format(curdir, th23, outname_O20)
    jobcmd_O23 = "{}/./callfit -d ORCA23 -t {} -o {}".format(curdir, th23, outname_O23)

    jobcmds.append( jobcmd_O20 )
    jobcmds.append( jobcmd_O23 )
    th23 = th23 + float( args['--step'] )

#==============================================================
# create scripts to be sent to the farm
#==============================================================

for i, cmd in enumerate(jobcmds):

    # create script
    #--------------
    
    script_name = "{0}/tmp/job_{1}.sh".format(curdir, i)
    script_file = open(script_name, 'w')

    script_file.write('#!/bin/bash\n\n')
    script_file.write('cd {0}\nsource setenv.sh\ncd {1}\n'.format(jpp, curdir)) #jpp
    script_file.write('source {}/bin/thisroot.sh\n'.format(root))               #root
    script_file.write("export OSCPROBDIR='{}'\n".format(oscp))                  #oscprob
    script_file.write('source {}/setenv.sh\n\n'.format(nmh))                    #nmh
    script_file.write(cmd)
    
    script_file.close()

    # create farm command
    #--------------------

    if args['--nikhef']:
        farm_cmd = "qsub -q short7 -o {0} -e {0} {1}".format(curdir+"/tmp", script_name)
    elif args['--lyon']:
        syscmd = "qsub -P P_km3net -l sps=1 -l ct=5:00:00 -o {0} -e {0} {1}".format(curdir+"/tmp", script_name)

    print ("Executing {}".format(farm_cmd))
    os.system(farm_cmd)
