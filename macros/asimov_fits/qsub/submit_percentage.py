#!/bin/env python
import os

for job in range(0,100):

    qsub_submit = "qsub -V -q short\
                  -o $MONADIR/macros/asimov_fits/output/qsub/percentages_{0}.log\
                  -e $MONADIR/macros/asimov_fits/output/qsub/percentages_{0}.err\
                  $MONADIR/macros/asimov_fits/qsub/job_{0}.sh".format(job)
 #   print(qsub_submit)
    os.system(qsub_submit)
