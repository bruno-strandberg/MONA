Main page {#mainpage}
===================

This software provides some tools for estimating the ORCA sensitivity to neutrino mass ordering. It relies on Monte-Carlo data from ORCA simulation chain (see `data/README.md` for a brief description of the simulation chain).

Prerequisities
==============
* ROOT 6: tested at Lyon with /usr/local/root/6.10.02
* Jpp with aanet, compiled against ROOT 6: tested with v10.1.11007
* [OscProb](https://github.com/joaoabcoelho/OscProb) package needs be available and compiled and an environment variable `OSCPROBDIR` needs to be set to point to the `OscProb` directory.
* Python scripts will require `docopt` package

More info
=========
* Regarding the data format(s) used in the analysis, consult `data/README.md` and `data_sorting/README.md`.
* Each subdirectory holds a `README.md` file that describes what the code can be used for.

Setup in lyon and elsewhere
===========================
TBD

How to use
==========

Sort data
----------
The starting point is the file `NMH/data/pid_result_XXX.root`. This needs to be converted to analysis format. Data sorting/conversion is handled by apps in `NMH/apps/data_sorting/`. Consult `NMH/apps/data_sorting/README.md`.

Effective mass
--------------
3D effective mass histograms can be created with the apps in `NMH/apps/effective_mass`. Consult `NMH/apps/effective_mass/README.md`. The combined output of the effective mass apps is stored in `data/eff_mass`.

Bjorken-y distributions
------------------------
TBD

Available applications
----------------------