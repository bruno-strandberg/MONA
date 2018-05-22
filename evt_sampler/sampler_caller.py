#!/usr/bin/env python
"""
This script can be used to call the GSGSampler.C macro and execute it on the farm.
 
Usage:
    sampler_caller -f FLAVOR... -i INTERACTION... -a FLUX_LIST [-n SAMPLES]
    sampler_caller -h                                                                     

Option:
    -f FLAVOR         Neutrino flavor 0 - elec, 1 - muon, 2 - tau, may select several
    -i INTERACTION    Interaction 0 - nc, 1 - cc, may select several
    -a FLUX_LIST      Text file that lists the FluxChain.C output files input to GSGSampler.
                      flux_caller.py outputs such a file.
    -n SAMPLES        Samples per flux file [default: 5]
    -h --help         Show this screen

"""

import sys
import os
from docopt import docopt

#===========================================================================

def main(args):

    nmhdir  = os.environ['NMHDIR']
    gsg_dir = nmhdir + "/data/mc_start/data_atmnu/"  #gseagen files on sps dir
    flavs   = {0:"elec", 1:"muon", 2:"tau"}
    ints    = {0:"NC", 1:"CC"}

    for f in args['-f']:
        for i in args['-i']:

            inter = ints[ int(i) ] 
            flav  = flavs[ int(f) ]

            # for NC we only have elec-NC
            if (inter == "NC" and flav != "elec"): continue

            # make tmp directory for logfiles etc
            os.system("mkdir -p tmp/")

            # create the list of gsg files
            gsg_flist = "tmp/{0}_{1}_gsg_flist.dat".format( flav, inter )
            os.system( "ls {0}*{1}*{2}* > {3}".format(gsg_dir, flav, inter, gsg_flist) )

            # create the command to execute GSGSampler
            cmd  = "root -b -q 'GSGSampler.C+("
            cmd += '"' + args['-a'] + '", '  # flux file list
            cmd += '"' + gsg_flist + '", '   # gseagen file list
            cmd += str(int(f)) + ", "        # neutrino flavor
            cmd += str(int(i)) + ", "        # neutrino interaction
            cmd += str(int(args['-n']))      # number of samples per flux file
            cmd += ")'"

            # create a bash script than can be sent to farm
            cwd  = os.getcwd()
            jpp  = os.environ['JPP_DIR']
            root = os.environ['ROOTSYS']
            nmh  = os.environ['NMHDIR']

            script_name = "tmp/farm_job_{0}_{1}.sh".format(flav, inter)
            script_file = open(script_name, 'w')

            script_file.write('#!/bin/bash\n\n')
            script_file.write('cd {0}\nsource setenv.sh\ncd {1}\n'.format(jpp, cwd)) #jpp
            script_file.write('source {}/bin/thisroot.sh\n'.format(root))            #root
            script_file.write('source {}/setenv.sh -a\n\n'.format(nmh))              #nmh
            script_file.write(cmd)

            script_file.close()

            #os.system("qsub -P P_km3net -l sps=1 -o {0} -e {0} {1}".format(cwd+"/tmp", script_name))
#===========================================================================

if __name__ == "__main__":
    args = docopt(__doc__)
    main(args)