#!/bin/env python

import os

job_file = "job_{0}.sh"

for job_nr in range(0,100):
  lines = ['#!/bin/bash', 'echo "input: {0}"'.format(job_nr),
           'command="root -q $NMHDIR/macros/asimov_fits/cross_check_macros/AsimovFitNO_PercentageOfMC.C\+\({0}\)"'.format(job_nr),
           'echo $command && eval $command']
  with open(job_file.format(job_nr), "w") as f:
    f.write("\n".join(lines))

  os.system("chmod +x {}".format(job_file.format(job_nr)))
