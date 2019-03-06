#!/bin/env/python

import numpy as np
import pandas as pd
import sys

import matplotlib.pyplot as plt

filefolder = "./output/csv/CrossCheck/"
savefolder = "./output/plots/"
filenameNO = "paramNMHValuesNoSystematicsAllNO.csv"
filenameIO = "paramNMHValuesNoSystematicsAllIO.csv"

dfno = pd.read_csv(filefolder + filenameNO)
dfio = pd.read_csv(filefolder + filenameIO)

dfno = dfno.loc[dfno.groupby("theta23").apply(lambda x: x["sigma"].idxmin()).values]
dfio = dfio.loc[dfio.groupby("theta23").apply(lambda x: x["sigma"].idxmin()).values]

dfno.plot(x="theta23", y=["sigma"], style=['-'], fontsize=12, 
        title="Asimov sensitivity (3 years)", color=["DarkRed"])
ax1 = plt.gca()
dfio.plot(x="theta23", y=["sigma"], style=['-'], fontsize=12, 
        title="Asimov sensitivity (3 years)", color=["DarkBlue"], ax=ax1)

plt.xlabel(r"$\theta_{23}$ [deg]", fontsize=14)
plt.ylabel(r"Sensitivity / $\sigma$", fontsize=14)
plt.ylim(0,8)
plt.grid()
plt.legend(["NO nature [paramNMH]", "IO nature [paramNMH]"], loc=2)
plt.savefig(savefolder + "AsimovSensitivityComparisonParamNMH.pdf")
plt.show()
