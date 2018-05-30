#!/usr/bin/env python
"""
This script is to be used in conjuction with flux_caller.py and sampler_caller.py. If flux_caller.py has been run and the FLUX_LIST output by flux_caller.py has been used as input by sampler_caller.py to sample muon, elec, tau CC and all-flavor NC events, then this script can be used to merge the various outputs into files representing an experiment each.
 
Usage:
    merge_to_exps -f FLUX_LIST -l FLUX_LOG -n SAMPLES
    merge_to_exps -h                                                                     

Option:
    -f FLUX_LIST      List of flux files input to sampler_caller.py, output by flux_caller.py
    -l FLUX_LOG       Logfile output by flux_caller.py
    -n SAMPLES        Nr of samples input to sampler_caller.py
    -h --help         Show this screen

"""

import sys
import os
from random import shuffle
from docopt import docopt
args = docopt(__doc__)

#=======================================================================
# start of script
#=======================================================================

flavs = ['elec','muon','tau', 'allflavs']
flux_files    = open(args['-f'], 'r')
sampler_files = os.popen( "ls output/EvtSample_*" ).read().split()
shuffle(sampler_files) #this randomizes which samples (same flux file!) are combined into exps
            
# this file stores the info which files are merged together
merge_log = open('output/merge_log.dat', 'w')

for i, flux_f in enumerate(flux_files):
    
    # get GSG Sampler outputs associated with this flux file
    search_str = "_flux{}_".format( i )
    sf_subset = [sf for sf in sampler_files if search_str in sf ]

    # sort these by flavors
    sorted_subset = {key: list() for key in flavs}
    for flav in flavs:
        sorted_subset[flav] = [sf for sf in sf_subset if flav in sf]

    # combine the samples into experiments
    samples_left = 1
    for key, value in sorted_subset.iteritems():
        samples_left = samples_left * len(value)

    experiments  = []
    while (samples_left):
        experiment = []
        for key, value in sorted_subset.iteritems():
            experiment.append( value.pop() )
            samples_left = samples_left * len(value)
        experiments.append( experiment )

    # fetch the oscillation parameters associated with the experiment
    osc_pars = []
    log_file = open(args['-l'], 'r')
    for line in log_file:
        if (flux_f[:-1] not in line and len(osc_pars) == 0): continue      # lookup
        osc_pars.append(line)                                              # store osc parameters 
        if (line == '\n'): break                                           # next flux file, break
    log_file.close()

    # determine hierarchy
    mh = "NH"
    if ("NH0" in flux_f): mh = "IH"

    # create hadd commands, store info to log
    merge_log.write('#====================================================================\n')
    merge_log.write('# Experiments with same oscillation parameters\n')
    merge_log.write('#====================================================================\n')

    for n, exp in enumerate(experiments):

        syscmd = "hadd output/Experiment_oscpars{0}_sample_{1}_{2}.root".format(i, n, mh)
        for file in exp:
            syscmd += " {}".format(file)
            
        os.system(syscmd)
        merge_log.write(syscmd + "\n")

    merge_log.write("\n")
    merge_log.write("EvtSample flux{0} <-> flux file {1}\n".format(i, flux_f) )
    for op in osc_pars:
        merge_log.write(op)
            
flux_files.close()
merge_log.close()
