#!/usr/bin/env python
"""
This script can be used to call the FluxChain.C macro several times with different (random) oscillation parameter settings. If nothing but the IDSTR is specified, the macro will call FluxChain.C twice, once with NH, once with IH, with all OscPars at central values. For a single call of FluxChain.C use the root prompt.
 
Usage:
    flux_caller [-n NR_OF_CONFS] [-p PAR] [-s STEP] [-y YEARS] -i IDSTR
    flux_caller -h                                                                     

Option:
    -n NR_OF_CONFS    Number of randomly sampled oscillation parameter configurations. 
                      If -p specified, then the number of scan steps in the parameter
                      values from -3sigma to 3sigma. For each configuration normal and
                      inverted hierarchy is used.
    -p PAR            Name of the parameter to be scanned. Supported names are
                      sinsq_th12, sinsq_th23, sinsq_th13, dcp, dm21, dm32
    -s STEP           If -p specified, the step with which the par is scanned. If -p 
                      and -s are specified, then -n is ignored.
    -y YEARS          Number of years of data taking [default: 3]
    -i IDSTR          Identifier string added as prefix to files in output/

    -h --help         Show this screen
"""

import sys
import os
import numpy
from docopt import docopt

#=========================================================================

class osc_p():
    """
    Class to represent an oscillation parameter
    """
    
    def __init__(self, _val_NH, _min_NH, _max_NH, _val_IH=None, _min_IH=None, _max_IH=None):
        """
        Constructor: a value and 3sigma limits are to be provided. If NH and IH have
        different values, values for the inverted hierarchy can also be provided.
        """

        self.val_NH = _val_NH
        self.min_NH = _min_NH
        self.max_NH = _max_NH

        if (_val_IH != None):
            self.val_IH = _val_IH
        else:
            self.val_IH = _val_NH

        if (_min_IH != None):
            self.min_IH = _min_IH
        else:
            self.min_IH = _min_NH

        if (_max_IH != None):
            self.max_IH = _max_IH
        else:
            self.max_IH = _max_NH

    def get_val(self, NH):
        """
        Function to return value, depending on hierarchy.
        """
        if (NH):
            return self.val_NH
        else:
            return self.val_IH

    def get_random_gaus(self, NH):
        """
        Function to get a random value of the parameter following a gaussian density
        distribution.
        """

        mean  = (self.min_NH + self.max_NH)/2. # assume mean to be middle
        sigma = (self.max_NH - self.min_NH)/6. # min to max is 6 sigma

        if (not NH):
            mean  = (self.min_IH + self.max_IH)/2.
            sigma = (self.max_IH - self.min_IH)/6.

        return numpy.random.normal( mean, sigma )

    def get_random_flat(self, NH):
        """
        Function to get a random value of the parameter, flat distribution between
        minimum and maximum.
        """

        if (NH):
            return numpy.random.uniform(self.min_NH, self.max_NH)
        else:
            return numpy.random.uniform(self.min_IH, self.max_IH)

    def __str__(self):
        """
        Function to print the parameter.
        """
        return ( "{0} {1} {2} {3} {4} {5}".format(self.val_NH, self.min_NH, self.max_NH,
                                                  self.val_IH, self.min_IH, self.max_IH) )
        
#=========================================================================

def main(args):

    # define oscillation parameters
    op = { "sinsq_th12": osc_p(0.297, 0.250, 0.354), 
           "sinsq_th23": osc_p(0.425, 0.381, 0.615, 0.589, 0.384, 0.636), 
           "sinsq_th13": osc_p(0.0215, 0.0190, 0.0240, 0.0216, 0.0190, 0.0242), 
           "dcp": osc_p(1.38, 1., 1.9, 1.31, 0.92, 1.88), 
           "dm21": osc_p(7.37e-5, 6.93e-5, 7.96e-5), 
           "dm32": osc_p(2.56e-3 - 7.37e-5, 2.45e-3, 2.69e-3, 
                         -2.54e-3, -2.66e-3, -2.42e-3) }

    # command list to be executed; logfile where configurations are written
    cmds = []
    outputs = []
    logfile = open('output/{}_fluxcaller_log.dat'.format(args['-i']),'w')

    # no parameter scanning asked, generate random samples
    if (args['-p'] == None):

        # if -n specified create n samples        
        if args['-n'] != None:

            for n in range( 0, int(args['-n']) ):
                
                # each sample creates hists with normal and inverted hierarchy
                for NH in [False, True]:
                    
                    outname = "output/{0}_sample{1}_NH{2}.root".format(args['-i'], n, int(NH))
                    sinsq12 = op['sinsq_th12'].get_random_flat(NH)
                    sinsq23 = op['sinsq_th23'].get_random_flat(NH)
                    sinsq13 = op['sinsq_th13'].get_random_flat(NH)
                    dcp = op['dcp'].get_random_flat(NH)
                    dm21 = op['dm21'].get_random_flat(NH)
                    dm32 = op['dm32'].get_random_flat(NH)
            
                    cmd = get_fluxchain_cmd(args['-y'], outname, NH, sinsq12, sinsq23, sinsq13,
                                            dcp, dm21, dm32)
                    vars_to_log(logfile, cmd, sinsq12, sinsq23, sinsq13, dcp, dm21, dm32)
                    cmds.append(cmd)
                    outputs.append(outname)

        # -n not specified, create FluxChain outputs with one normal, one inverted hierarchy
        else:
            for NH in [False, True]:
                outname = "output/{0}_sample{1}_NH{2}.root".format(args['-i'], 0, int(NH))
                sinsq12 = op['sinsq_th12'].get_val(NH)
                sinsq23 = op['sinsq_th23'].get_val(NH)
                sinsq13 = op['sinsq_th13'].get_val(NH)
                dcp = op['dcp'].get_val(NH)
                dm21 = op['dm21'].get_val(NH)
                dm32 = op['dm32'].get_val(NH)
            
                cmd = get_fluxchain_cmd(args['-y'], outname, NH, sinsq12, sinsq23, sinsq13,
                                        dcp, dm21, dm32)
                vars_to_log(logfile, cmd, sinsq12, sinsq23, sinsq13, dcp, dm21, dm32)
                cmds.append(cmd)
                outputs.append(outname)

    # scan over one parameter
    else:
        
        # check for input
        if (args['-p'] not in op):
            print("Variable {} not recognised, supported variables are:".format(args['-p']))
            for key, value in op.iteritems(): print key
            raise Exception("Exiting.")

        p = op[ args['-p'] ]

        # if step defined use that, otherwise calculate step from n
        if (args['-s'] != None):
            step_NH = float(args['-s'])
            step_IH = step_NH
        else:
            step_NH = (p.max_NH - p.min_NH)/float(args['-n'])
            step_IH = (p.max_IH - p.min_IH)/float(args['-n'])

        # sample for NH
        n = 0
        p.val_NH = p.min_NH
        while (p.val_NH <= p.max_NH):
            outname = "output/{0}_step{1}_NH1.root".format(args['-i'], n)
            cmd = get_fluxchain_cmd(args['-y'], outname, 1, 
                                    op['sinsq_th12'].val_NH, op['sinsq_th23'].val_NH, 
                                    op['sinsq_th13'].val_NH, op['dcp'].val_NH, 
                                    op['dm21'].val_NH, op['dm32'].val_NH)
            vars_to_log(logfile, cmd, op['sinsq_th12'].val_NH, op['sinsq_th23'].val_NH, 
                        op['sinsq_th13'].val_NH, op['dcp'].val_NH, op['dm21'].val_NH, 
                        op['dm32'].val_NH)
            cmds.append(cmd)
            outputs.append(outname)
            p.val_NH += step_NH
            n += 1

        # sample for IH
        n = 0
        p.val_IH = p.min_IH
        while (p.val_IH <= p.max_IH):
            outname = "output/{0}_step{1}_NH0.root".format(args['-i'], n)
            cmd = get_fluxchain_cmd(args['-y'], outname, 0, 
                                    op['sinsq_th12'].val_IH, op['sinsq_th23'].val_IH, 
                                    op['sinsq_th13'].val_IH, op['dcp'].val_IH, 
                                    op['dm21'].val_IH, op['dm32'].val_IH)
            vars_to_log(logfile, cmd, op['sinsq_th12'].val_IH, op['sinsq_th23'].val_IH, 
                        op['sinsq_th13'].val_IH, op['dcp'].val_IH, op['dm21'].val_IH, 
                        op['dm32'].val_IH)
            cmds.append(cmd)
            outputs.append(outname)
            p.val_IH += step_IH
            n += 1
            
    # execute the commands; close logfile; write file list
    logfile.close()

    outlist = open('output/{}_output_list.dat'.format(args['-i']), 'w')
    for o in outputs:
        outlist.write(o + "\n")
    outlist.close()

    for cmd in cmds:
        os.system(cmd)


#=========================================================================

def get_fluxchain_cmd(years, outname, NH, sinsq12, sinsq23, sinsq13, dcp, dm21, dm32):
    """
    Function to return a command to execute the FluxChain.C macro.
    """

    syscmd  = "root -b -q 'FluxChain.C+("
    syscmd += str(years) + ", "      # 3 years running time
    syscmd += '"' + outname + '", '  # output name
    syscmd += str(20) + ","          # 20-fold oversampling in each bin
    syscmd += str(int(NH)) + ", "    # 0 - IH, 1 - NH
    syscmd += str(0) + ", "          # 0 - don't use Meff, 1 - use Meff
    syscmd += str(sinsq12) + ", "    # oscillation parameter values
    syscmd += str(sinsq23) + ", "
    syscmd += str(sinsq13) + ", "
    syscmd += str(dcp) + ", "
    syscmd += str(dm21) + ", "
    syscmd += str(dm32)
    syscmd += ")'"

    return syscmd

#=========================================================================

def vars_to_log(logfile, cmd, sinsq12, sinsq23, sinsq13, dcp, dm21, dm32):
    """
    Function that appends and command and the oscillation parameters to the logfile
    """
    logfile.write("{}\n".format(cmd))
    logfile.write("sin^2 th_12: {}\n".format(sinsq12) )
    logfile.write("sin^2 th_23: {}\n".format(sinsq23) )
    logfile.write("sin^2 th_13: {}\n".format(sinsq13) )
    logfile.write("delta cp: {}\n".format(dcp) )
    logfile.write("delta m_12: {}\n".format(dm21) )
    logfile.write("delta m_32: {}\n\n\n".format(dm32) )


if __name__ == "__main__":
    args = docopt(__doc__)
    main(args)
