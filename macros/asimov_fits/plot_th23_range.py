#!/bin/env/python

import numpy as np
import pandas as pd
import sys

import matplotlib.pyplot as plt

filefolder = "./output/csv/"
savefolder = "./output/plots/"

df = pd.read_csv(filefolder + "AsimovFitTh23Range.csv")

try:
  df['n_chi2_no'] = np.sqrt(df['n_chi2tr_no']**2 + df['n_chi2sh_no']**2)
  df['n_chi2_io'] = np.sqrt(df['n_chi2tr_io']**2 + df['n_chi2sh_io']**2)
except:
  print("No track and shower data found, using combined data instead of calculating it")
  df['n_chi2_no'] = np.sqrt(df['n_chi2_no'])
  df['n_chi2_io'] = np.sqrt(df['n_chi2_io'])

df.plot(x="th23", y=["n_chi2_io", "n_chi2_no"], style=['-','-'], fontsize=12, 
        title="Asimov sensitivity (3 years)", color=["DarkRed", "DarkBlue"])
plt.xlabel(r"$\theta_{23}$ [deg]", fontsize=14)
plt.ylabel(r"Sensitivity / $\sigma$", fontsize=14)
plt.ylim(0,8)
plt.grid()
plt.legend(["NO nature","IO Nature"])
plt.savefig(savefolder + "AsimovSensitivity2Bins.pdf")
plt.show()
