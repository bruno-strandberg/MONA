#!/usr/bin/python
"""
This macro can be used to farm out FitExps macro on the farm.

Usage:
    fitcaller.py [--nmin NMIN] [--nmax NMAX] [--odir OUTPUTDIR] [--infile INFILE]
    fitcaller.py -h --help

Option:
    --nmin   NMIN       Minimum number of samples [default: 1]
    --nmax   NMAX       Maximum number of samples [default: 10]
    --odir   OUTPUTDIR  Directory where farming outputs are written [default: tmp/]
    --infile INFILE     Input file to FitExps [default: OscResolution.root]
"""

from docopt import docopt
args = docopt(__doc__)
import os

os.system("mkdir -p {}/{}".format( os.getcwd(), args['--odir'] ) )
fin = "{}/{}".format(os.getcwd(), args['--infile'])

for n in range(int(args['--nmin']),int(args['--nmax'])+1):

    trkn = "trk_n{}".format(n)
    shwn = "shw_n{}".format(n)
    outn = os.getcwd() + "/{}/fit_n{}.root".format(args['--odir'],n)
    
    syscmd  = "root -b -q '{}/FitExps.C+(".format(os.getcwd())
    syscmd += '"{}",'.format(fin)
    syscmd += '"{}",'.format(trkn)
    syscmd += '"{}",'.format(shwn)
    syscmd += '"{}"'.format(outn)
    syscmd += ")'"
    
    #----------------------------------------------
    # create a bash script than can be sent to farm
    #----------------------------------------------

    script_name = "{0}/{1}/job_{2}.sh".format(os.getcwd(), args['--odir'], n)
    script_file = open(script_name, 'w')

    script_file.write('#!/bin/bash\n\n')
    script_file.write('cd {0}\nsource setenv.sh\ncd {1}\n'.format(os.environ["JPP_DIR"], os.getcwd())) #jpp
    script_file.write('source {}/bin/thisroot.sh\n'.format(os.environ["ROOTSYS"]))                     #root
    script_file.write("export OSCPROBDIR='{}'\n".format(os.environ["OSCPROBDIR"]))                     #oscprob
    script_file.write('source {}/setenv.sh\n\n'.format(os.environ["MONADIR"]))                         #nmh
    script_file.write(syscmd)
    
    script_file.close()

    farm_cmd = "qsub -q short7 -o {0} -e {0} {1}".format(os.getcwd()+"/" + args['--odir'], script_name)

    print ("Executing {}".format(farm_cmd))
    os.system(farm_cmd)


