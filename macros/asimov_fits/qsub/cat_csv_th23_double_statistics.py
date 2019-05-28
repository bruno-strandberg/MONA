#!/bin/env python

import os
import shutil

for order in ["no","io"]:
  for pid in [2,3,4,5,10]:
    path = "$MONADIR/macros/asimov_fits/output/csv/SensChi2Inf_20mCrossCheck/AsimovFit{0}Bins{1}Th23Range_PercentageOfMC/".format(pid, order.upper())
    file0_99 = "AsimovFit{0}Bins{1}Th23Range_PercentageOfMC_0-99.csv".format(pid, order.upper())
    file0_199 = "AsimovFit{0}Bins{1}Th23Range_PercentageOfMC_0-199.csv".format(pid, order.upper())
    os.system("mkdir " + path + "tmp/")
    os.system("mv " + path + file0_99 + " " + path + "tmp/" + file0_99)

    command = "cat " + path + "*.csv > " + path + file0_199

    #print(command.format(pid, order.upper()))
    os.system(command.format(pid, order.upper()))
