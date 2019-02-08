#!/bin/env/python

import numpy as np
import pandas as pd
import sys

import matplotlib.pyplot as plt

filefolder = "./"
filename1 = sys.argv[1]
filename2 = sys.argv[2]

df1 = pd.read_csv(filefolder + filename1)
df2 = pd.read_csv(filefolder + filename2)

df = pd.merge(df1, df2, left_on="th23", right_on="th23")

df['chi2_no'] = np.sqrt(df['n_chi2tr_no']**2 + df['n_chi2sh_no']**2)
df['chi2_io'] = np.sqrt(df['n_chi2tr_io']**2 + df['n_chi2sh_io']**2)

df.plot(x="th23", y=["n_chi2tr_no", "n_chi2sh_no", "chi2_no", "n_chi2tr_io", "n_chi2sh_io", "chi2_io"])
plt.xlabel(r"$\theta_{23}$")
plt.ylabel("Sensitivity")
plt.grid()
plt.show()
