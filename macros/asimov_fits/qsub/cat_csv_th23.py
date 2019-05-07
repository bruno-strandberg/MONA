#!/bin/env python

import os


for order in ["no","io"]:
  for pid in [2,3,4,5,10]:
    path = "$MONADIR/macros/asimov_fits/output/csv/SensChi2Inf/AsimovFit{0}Bins{1}Th23Range_PercentageOfMC/"
    command = "cat " + path + "* > " + path + "AsimovFit{0}Bins{1}Th23Range_PercentageOfMC_0-99.csv"

    #print(command.format(pid, order.upper()))
    os.system(command.format(pid, order.upper()))
