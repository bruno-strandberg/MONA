#!/bin/env/python

import numpy as np
import pandas as pd
import sys

import matplotlib.pyplot as plt

filefolder = "./"
filename = sys.argv[1] if (len(sys.argv) > 1) else "BinRange.csv"

df = pd.read_csv(filefolder + filename)

df.plot(x="n_pid", y=["chi2_no", "chi2_io", "chi2_no_3.60", "chi2_io_3.60", "chi2_no_5.15", 
        "chi2_io_5.15"], style=['o-','o--', 'o-', 'o--', 'o-', 'o--'])
plt.xlabel(r"Number of PID bins [default=2]", fontsize=12)
plt.ylabel(r"Sensitivity / $\sigma$ (3 years)", fontsize=12)
plt.grid()
plt.show()

