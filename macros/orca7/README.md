orca7
======

Requirements
------------
To run these scripts, one requires:
1. Nikhef NMO package setup
2. Summary file and effective mass file for ORCA 7-line analysis in NNMO format

Scripts
-------
This directory contains some root macros for the analysis of the ORCA 7-line MC production. The macros are briefly documented in the source codes, short descriptions of the macros are presented below:

* `expectationplots.C` - this macro draws plots that illustrate noise suppression, depict expectation values and sensitivity to theta-23, dm31.

* `resolutions.C` - this histogram draws energy reco and direction reco plots for track and shower reconstructions.

* `thetafit.C` - this macro performs one fit to a pseudo-experiment corresponding to 1 year of data taking to extract theta-23, using only events with a high track score. It draws plots to depict the expectation value, the pseudo-experiment, a LLH scan, and a difference between model and data at different theta values.

* `contours.C` - this macro investigates the sensitivity to (theta-23, dm31) in various event selections by performing likelihood scans and drawing contour plots.

* `th23dm31.C` - this macro performs a specified number of fits to pseudo-experiments corresponding to 1 year of data talking to extract theta-23 and dm-31.

* `ana_th23dm31.C` - this macro analyses the output of `th23dm31.C` to create a plot of the true th23, dm31 point and the results from pseudo-experiments.