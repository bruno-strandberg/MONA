orca7
======

This directory contains some root macros for the analysis of the ORCA 7-line MC production. Short descriptions of the macros below:

* `expectationplots.C` - this macro draws expectation value histograms for tracks and showers, separated at track score 0.6. The histograms depict the total number of detected neutrionos in 1 year, projections to 1D and display atmospheric muon and noise contributions.

* `resolutions.C` - this histogram draws energy reco and direction reco plots for track and shower reconstructions.

* `thetafit.C` - this macro performs one fit to a pseudo-experiment corresponding to 1 year of data taking to extract theta-23, using only  events with a high track score.

* `th23dm31.C` - this macro performs a specified number of fits to pseudo-experiments corresponding to 1 year of data talking to extract theta-23 and dm-31.

* `ana_th23dm31.C` - this macro analyses the output(s) of `th23dm31.C`.