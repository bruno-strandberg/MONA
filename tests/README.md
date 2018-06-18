Test scripts
============

* flux_check_inpol.py - script that compares the honda flux as read from the flux table file and estimated by the AtmFlux interpolator.

* flux_comp_mj.py - script that compares the atm. flux as estimated by the class AtmFlux and as stored in NMH/data/testing/mjongen_plots.root by Martijn's code.

* gsgparser.C - script to compare GSGParser event parsing with different methods.

* sanity_checks.C - script that takes the NMH/data/pid_result_XXX.root as input and draws some plots for sanity checks.

* xsec.py - script that draws xsec/E plots for various interactions.

* data_sorting.py - script to check that data_sorting/ scripts worked properly.

* evtid_functionality.C+ - script to check the find and sort algorithms used in evt_sampler/GSGSampler.C

* fileheader.C - a script to check/demostrate FileHeader functionality.