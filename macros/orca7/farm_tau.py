#!/usr/bin/python

"""
Script to farm the application `apps/LLH_scanner` for tau sensitivity analysis in time.

Usage:
    farm_tau.py --tmin TMIN --tmax TMAX --step STEP [-d ODIR]
    farm_tau.py -h

Option:
    --tmin TMIN  Minimum ORCA7 run period in months [default: 1]
    --tmax TMAX  Maximum ORCA7 run period in months [default: 12]
    --step STEP  Step in time to scan from TMIN to TMAX [default: 1]
    -d ODIR      Output directory [default: tauscan/]
    -h --help    Show this screen
"""

import os, math
from docopt import docopt
args = docopt(__doc__)

# create output directory
odir = "{}/{}".format( os.getcwd(), args['-d'] )
os.system( "mkdir -p {}".format( odir ) )

#========================================================
# create application execution commands
#========================================================

appcmd  = "{}/apps/llh_scanner/./LLH_scanner ".format( os.environ['MONADIR'] )
appcmd += ' -M 0.05+0.01+0.05 '     # same muon cuts as in ORCA7 analysis
appcmd += ' -N 0.0+0.0+0.1 '        # same noise cuts as in ORCA7 analysis
appcmd += ' -P 0.0+0.3+0.7+1.0 '    # same PID ranges as in ORCA7 analysis
appcmd += ' -R shw -R shw -R comb ' # same reco types as in ORCA7 analysis
appcmd += ' -d {}/data/ORCA_MC_summary_ORCA7_23x9m_ECAP1018.root '.format( os.environ['MONADIR'] )
appcmd += ' -m {}/data/eff_mass/EffMass_ORCA7_23x9m_ECAP1018.root '.format( os.environ['MONADIR'] )
appcmd += ' -p Tau_norm '           # scan in tau normalisation
appcmd += ' -r 0.0+2.0 '            # in range 0 to 2
appcmd += ' -w '                    # request profile LLH scan
appcmd += ' -n 10 '                 # do the scan in 10 points

# set parameter values for data creation and fit start
appcmd += ' -x SinsqTh12 -x Dm21 -x SinsqTh13 -x SinsqTh23 -x Dm31 -x dcp '
appcmd += ' -y {} -y {} -y {} -y {} -y {} -y {} '.format( math.sin( math.radians(33.4) )**2,
                                                          7.53 * 1e-5,
                                                          math.sin( math.radians(8.42) )**2,
                                                          math.sin( math.radians(42.0) )**2,
                                                          2.44 * 1e-3 + 7.53 * 1e-5/2,
                                                          0. )

# fix a selection of parameters
appcmd += " -z SinsqTh12 -z SinsqTh13 -z Dm21 -z E_scale -z ct_tilt "


syscmds = []
runTime = float(args['--tmin'])
while runTime <= float(args['--tmax']):
    
    thiscmd  = appcmd
    thiscmd += ' -t {} '.format( runTime/12 )
    thiscmd += ' -o {}/tauProfile_t_{}_months.root '.format( odir, runTime )

    syscmds.append( thiscmd )

    runTime += float(args['--step'])

#========================================================
# create job scripts
#========================================================

jobfiles = []

for i,syscmd in enumerate(syscmds):

    jobfilename = "{}/job_{}.sh".format(odir, i)
        
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

    jobfile.write( syscmd )
    jobfile.close()
    jobfiles.append( jobfilename )

#========================================================
# execute
#========================================================

for job in jobfiles:
    
    farmcmd = "qsub -q short7 -o {0} -e {0} {1}".format(odir, job)
    os.system(farmcmd)
