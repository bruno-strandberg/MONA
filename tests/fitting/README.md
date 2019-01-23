Fitting examples
================

The three applications in this directory demonstrate the usage `common_software` and `fitter_software` libraries for NMO fits. The applications are documented in their source codes and summarised briefly here:

* `SampleFitter` is an example application to fit an experiment sample created with `apps/evt_sampler`
* `FitPseudoExperiment` is an example application that creates a pseudo-experiment and fits it using ROOT and RooFit for comparison.
* `Asymmetry` is an application that calculates the combined asymmetry in the track and shower channel between normal and inverted mass ordering at a specified theta value.

The applications were developed while analysing ORCA 23x9m production from 2015 and uses event selections specific to that production. For that reason the applications are not suitable for generic use with all productions, but they are hopefully reasonable examples for future developments.