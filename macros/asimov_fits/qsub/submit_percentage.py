#!/bin/env python
import os

for job in range(0,100):

    qsub_submit = "qsub -V -q short\
                  -o $NMHDIR/macros/asimov_fits/output/qsub/percentages_{0}.log\
                  -e $NMHDIR/macros/asimov_fits/output/qsub/percentages_{0}.err\
                  $NMHDIR/macros/asimov_fits/qsub/job_{0}.sh".format(job)
 #   print(qsub_submit)
    os.system(qsub_submit)
