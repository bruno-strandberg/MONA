#!/usr/bin/env python

# The print function for simplified errors is in pid_detres/simplified_error_calc/

import os

basedir = os.environ.get("NMHDIR")
basedir = basedir + "/fitter/"

def asym_print(filename):
    lt6 = False
    if "default" in filename:
        lt6 = True
    with open(filename) as f:
        for line in f:
            words = line.split()
            if "q<0.6" in line:
                lt6 = True
            elif "q<0.9" in line:
                lt6 = False
            if lt6:
                if "track" in line:
                    tr_num = float(words[5])
                    tr_err = float(words[7])
                elif "shower" in line:
                    sh_num = float(words[4])
                    sh_err = float(words[6])
                elif "combined" in line:
                    cm_num = float(words[4])
                    cm_err = float(words[6])
                else:
                    pass
    return [tr_num, tr_err, sh_num, sh_err, cm_num, cm_err]

bins_2 = asym_print(basedir + "default_detres/asym_output_correlations.txt")
bins_5 = asym_print(basedir + "pid_detres/pid_binning_5/asym_output_correlations.txt")
bins_10 = asym_print(basedir + "pid_detres/asym_output_correlations.txt")
bins_15 = asym_print(basedir + "pid_detres/pid_binning_15/asym_output_correlations.txt")
bins_20 = asym_print(basedir + "pid_detres/pid_binning_20/asym_output_correlations.txt")
bins_25 = asym_print(basedir + "pid_detres/pid_binning_25/asym_output_correlations.txt")
bins_30 = asym_print(basedir + "pid_detres/pid_binning_30/asym_output_correlations.txt")
bins_35 = asym_print(basedir + "pid_detres/pid_binning_35/asym_output_correlations.txt")
bins_40 = asym_print(basedir + "pid_detres/pid_binning_40/asym_output_correlations.txt")

bins = [bins_2, bins_5, bins_10, bins_15, bins_20, bins_25, bins_30, bins_35, bins_40]

tr_values = [b[0] for b in bins]
tr_errs = [b[1] for b in bins]
sh_values = [b[2] for b in bins]
sh_errs = [b[3] for b in bins]
cm_values = [b[4] for b in bins]
cm_errs = [b[5] for b in bins]

print "Track  values: ", tr_values
print "Track  errors: ", tr_errs
print "Shower values: ", sh_values
print "Shower errors: ", sh_errs
print "Comb   values: ", cm_values
print "Comb   errors: ", cm_errs
