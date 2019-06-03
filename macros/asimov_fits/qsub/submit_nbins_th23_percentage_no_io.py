#!/bin/env python
import os

for order in ["no", "io"]:
  for pid in [2,3,4,5,10,31,32,41]:
    for job_nr in range(0,200):

        qsub_submit = "qsub -V -q short\
                      -o $MONADIR/macros/asimov_fits/output/qsub/percentages_{0}_{1}_{2}.log\
                      -e $MONADIR/macros/asimov_fits/output/qsub/percentages_{0}_{1}_{2}.err\
                      $MONADIR/macros/asimov_fits/qsub/job_{0}_{1}_{2}.sh".format(order, pid, job_nr)

        os.system(qsub_submit)
