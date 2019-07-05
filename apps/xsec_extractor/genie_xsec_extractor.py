#!/usr/bin/python
"""

This script extracts the neutrino-interaction cross-sections in sea water (oxygen and hydrogen) from gSeaGen/GENIE splines using GENIE's `gspl2root` tool. The output file is reduced by the `Reductor` application to contain only the graphs that are used by the `NuXsec.h/C` class of `MONA/common_software` 

Usage:
   genie_xsec_extractor.py [-g GSGENV] [-x SPLINE] [-e EMAX] [-o OUTFILE]
   genie_xsec_extractor.py -h

Option:
   -g GSGENV    gSeaGen enviroment script [default: /pbs/throng/km3net/src/gSeaGen/v5r1/setenv.sh]
   -x SPLINE    GENIE spline with merged muon, elec, tau data [default: /pbs/throng/km3net/src/gSeaGen/v5r1/dat/gxspl-seawater.xml]
   -e EMAX      Maximum energy for the xsec data [default: 1000]
   -o OUTFILE   Output file for input to NuXsec [default: xsec_gSeaGen_v5r1.root]
   -h --help    Show this screen

"""

from docopt import docopt
args = docopt(__doc__)

import os

#===================================================================================
# create a bash script to execute the gspl2root command
#===================================================================================

root = os.environ["ROOTSYS"]
scriptname = "tmp.sh"
tmpout     = "tmp.root"

scriptf = open(scriptname, 'w')
scriptf.write( '#!/bin/bash\n' )
scriptf.write( 'source {}\n'.format( args['-g'] ) )
scriptf.write( 'gspl2root -f {} -p 12,-12,14,-14,16,-16 -t 1000010010,1000080160 -e {} -o {}\n'.format( args['-x'], args['-e'], tmpout ) )
scriptf.write("source {}/bin/thisroot.sh".format(root))
scriptf.close()

#===================================================================================
# execute the script and execute reductor, remove temp files
#===================================================================================

os.system( "bash {}".format(scriptname) )
os.system( "./Reductor -e {} -f {} -g {} -o {} -x {}".format( args['-e'], tmpout, args['-g'], args['-o'], args['-x'] ) )
os.remove(scriptname)
os.remove(tmpout)

