#!/bin/env/python

import numpy as np
import pandas as pd
import sys

import matplotlib.pyplot as plt

filefolder = "./output/csv/"
filename = sys.argv[1] if (len(sys.argv) > 1) else "QRange.txt"

df = pd.read_csv(filefolder + filename)

y=["chi2_no", "chi2_io", "chi2_no_3.50", "chi2_io_3.50", "chi2_no_3.30", "chi2_io_3.30",
   "chi2_no_3.20","chi2_io_3.20"]
linestyle=['solid','dashed', 'solid', 'dashed', 'solid', 'dashed', 'solid', 'dashed']
marker=["o", "o", "o", "o", "o", "o", "o", "o"]
color=["darkred", "darkred", "royalblue", "royalblue", "forestgreen", "forestgreen", "indigo", "indigo"]

fig, ax=plt.subplots()
for y,s,m,c in zip(y,linestyle, marker, color):
  df.plot(x="q_cut", y=y, linestyle=s, marker=m, color=c, ax=ax,
          title="Asimov sensitivity at PDG values (3 years)", fontsize=12)
plt.xlabel(r"Edge of middle bin for (3 bins)", fontsize=14)
plt.ylabel(r"Sensitivity / $\sigma$", fontsize=14)
plt.ylim(2.5,4)
plt.grid()
plt.legend(["NO", "IO", "NO energy cut [3,50]", "IO energy cut [3,50]",
            "NO energy cut [3,30]", "IO energy cut [3,30]",
            "NO energy cut [3,20]", "IO energy cut [3,20]",], loc=2)
plt.show()

