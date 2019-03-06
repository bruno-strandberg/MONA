#!/bin/env/python

import numpy as np
import pandas as pd
import sys

import matplotlib.pyplot as plt

filefolder = "./output/csv/"
filename = sys.argv[1] if (len(sys.argv) > 1) else "BinRange.txt"

df = pd.read_csv(filefolder + filename)

df.plot(x="n_pid", y=["chi2_no", "chi2_io"], style=['o-','o--'])
plt.xlabel(r"Number of PID bins [default=2]", fontsize=12)
plt.ylabel(r"Sensitivity / $\sigma$ (3 years)", fontsize=12)
plt.grid()
plt.show()

