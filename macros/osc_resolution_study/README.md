Oscillation resolution study
=============================

Macros in this directory were used to investigate the effect of insufficient sampling in oscillation probability calculation. The results were presented at the ORCA phone meeting on 5.04.19 by B. Strandberg.

Quick quide
===========

1) `root OscResolution.C+`, this creates an output `OscResolution.root` with expectation value plots at different oscillation sampling density values.

2) `root`, `.x FitExps.C(...)`. This fits the expectation value histograms in the input file, which should be `OscResolution.root`. In practice, fitting can be distributed to nikhef farm by executing `fitcaller.py`.

3) `root`, `.x AnalyseFits.C(...)`. This macro makes plots with true and fitted parameter values.

More info/documentation in the macro codes.