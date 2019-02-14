#!/bin/env/python

import numpy as np
import pandas as pd
import sys

import matplotlib.pyplot as plt

filefolder = "./"
filename = sys.argv[1] if (len(sys.argv) > 1) else "AsimovFitTh23Range.txt"

df = pd.read_csv(filefolder + filename)

df['chi2_no'] = np.sqrt(df['n_chi2tr_no']**2 + df['n_chi2sh_no']**2)
df['chi2_io'] = np.sqrt(df['n_chi2tr_io']**2 + df['n_chi2sh_io']**2)

df.plot(x="th23", y=["n_chi2tr_no", "n_chi2sh_no", "chi2_no", "n_chi2tr_io", "n_chi2sh_io", "chi2_io"],
        style=['-','-','-','--','--','--'])
plt.xlabel(r"$\theta_{23}$", fontsize=12)
plt.ylabel(r"Sensitivity / $\sigma$ (3 years)", fontsize=12)
plt.grid()
plt.show()

