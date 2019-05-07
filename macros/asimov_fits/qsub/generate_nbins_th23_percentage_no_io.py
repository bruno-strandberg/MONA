#!/bin/env python

import os

job_file = "job_{0}_{1}_{2}.sh"

for order in ["no", "io"]:
  for pid in [2,3,4,5,10]:
    for job_nr in range(0, 100):
      lines = ['#!/bin/bash', 'echo "input jobnumber: {0} pidnumber: {1}"'.format(job_nr, pid),
               'command="$MONADIR/macros/asimov_fits/sensitivity_chi2_inf/AsimovQSub{2} {0} {1}"'.format(job_nr, pid, order.upper()),
               'echo $command && eval $command']
      with open(job_file.format(order, pid, job_nr), "w") as f:
        f.write("\n".join(lines))

      os.system("chmod +x {}".format(job_file.format(order, pid, job_nr)))
