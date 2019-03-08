#!/bin/env/python

import numpy as np
import pandas as pd
import sys

import matplotlib.pyplot as plt

"""
This script generates for these four scenarios simultaneously, which is NOT NECESSARY
if you only want the standard plot.

STANDARD PLOT: only have the output with Poisson errors ready in output/csv and run this
script with ADJUSTED folders and extensions

ALL PLOTS: have all four scenarios ready in the folders in
output/csv/{ExpectedError, PoissonError, Ebin24_Expected, Ebin24_Poisson}
and extensions {Expected, Poisson, ExpectedEbin24, PoissonEbin24}
to generate the four plots of the sensitivity comparison curves. These scenarios
mean the following:
Fit using Expected errors, Fit using Poisson errors, Use 24 Ebins instead of 40 for the fit
with Expected Errors, Use 24 Ebins for the fit with Poisson Errors.
"""

folders = ["ExpectedError", "PoissonError", "Ebin24_Expected", "Ebin24_Poisson"]
extensions = ["Expected", "Poisson", "ExpectedEbin24", "PoissonEbin24"]

for folder, fileext in zip(folders, extensions):
  filefolder = "./output/csv/{0}/".format(folder)
  savefolder = "./output/plots/"
  filename = sys.argv[1] if (len(sys.argv) > 1) else "AsimovFitTh23Range.csv"

  df = pd.read_csv(filefolder + filename)

  df['n_chi2_no'] = np.sqrt(df['n_chi2tr_no']**2 + df['n_chi2sh_no']**2)
  df['n_chi2_io'] = np.sqrt(df['n_chi2tr_io']**2 + df['n_chi2sh_io']**2)

  df.plot(x="th23", y=["n_chi2_io", "n_chi2_no"], style=['-','-'], fontsize=12,
          title="Asimov sensitivity (3 years)", color=["DarkRed", "DarkBlue"])
  ax1 = plt.gca()

  df_3bins = pd.read_csv(filefolder + "AsimovFit3BinsTh23Range.csv")
  df_3bins.plot(x="th23", y=["n_chi2_io", "n_chi2_no"], style=['--','--'], fontsize=12,
  # There is something wrong with my labeling...
                title="Asimov sensitivity (3 years)", color=["DarkRed", "DarkBlue"],
                ax=ax1)

  plt.xlabel(r"$\theta_{23}$ [deg]", fontsize=14)
  plt.ylabel(r"Sensitivity / $\sigma$", fontsize=14)
  plt.ylim(0,8)
  plt.grid()
  plt.legend(["NO nature [2 bins]", "IO nature [2 bins]", "NO nature [3 bins]", "IO nature [3 bins]"])
  plt.savefig(savefolder + "AsimovSensitivityComparison2vs3{0}.pdf".format(fileext))
  plt.show()
