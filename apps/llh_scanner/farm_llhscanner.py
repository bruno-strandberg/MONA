#!/usr/bin/python

"""
Script to farm out LLH_scanner application. It creates a LLH scan job for each parameter.

Usage:
    farm_llhscanner [-d DATA_FILE] [-m EFFMASS_FILE] [-n NPOINTS] --evtresp -i ID_STR
    farm_llhscanner -h

Option:
    -d DATA_FILE      MC chain summary data file in MONA format
    -m EFFMASS_FILE   Effective mass file in MONA format
    -n NPOINTS        Number of points in LLH scan. If not specified, LLH_scanner default is used.
    --evtresp         Use event-by-event response
    -i ID_STR         Identifier string used to create job names and job outputs
    -h --help         Show this screen

"""

from docopt import docopt
args = docopt(__doc__)
import os

def create_script(script_name, envscript_name, cmd):
    
    script = open(script_name, 'w')
    
    script.write('#!/bin/bash\n\n')

    # write everything from the env script after the first source command
    envscript = open(envscript_name, 'r')
    ignoreall = True
    for line in envscript:
        if ( "source" in line ): ignoreall = False
        if (ignoreall): continue
        script.write(line)

    # add the execution command
    script.write("\n\n" + cmd + "\n\n")

    script.close()

    # return the script name
    return script_name
    

#***********************************************************************************

# parameters that are scanned in LLH
pars = ["SinsqTh12", "SinsqTh13", "SinsqTh23", "dcp", "Dm21", "Dm31", "E_tilt", "ct_tilt", "skew_mu_amu", "skew_e_ae", "skew_mu_e", "NC_norm", "Tau_norm", "E_scale"]

# create some directories
cwd = os.getcwd()
outdir = cwd + "/farm_out"
tmpdir = cwd + "/farm_tmp"
os.system("mkdir -p {}".format(outdir))
os.system("mkdir -p {}".format(tmpdir))

# create execution commands
app = "LLH_scanner"
scripts = []

for par in pars:

    jobname   = "{}_{}".format(args['-i'], par)
    jobscript = "{}/{}.sh".format(tmpdir, jobname) 
    output    = "{}/{}.root".format(outdir, jobname)

    cmd = "{0}/./{1} -p {2} -o {3}".format(cwd, app, par, output)

    if ( args['-d'] != None ):
        fpath = os.path.abspath(args['-D'])
        cmd += " -D {}".format(fpath)

    if ( args['-m'] != None ):
        fpath = os.path.abspath(args['-M'])
        cmd += " -M {}".format(fpath)

    if ( args['-n'] != None ):
        cmd += " -N {}".format( args['-n'] )

    if args['--evtresp']:
        cmd += ' -e'

    jobscript = create_script(jobscript, "/user/bstrand/.bashrc", cmd)
    scripts.append( jobscript )

# send all scripts to the farm
for script in scripts:

    syscmd = "qsub -q short7 -o {0} -e {0} {1}".format(tmpdir, script)
    os.system( syscmd )
