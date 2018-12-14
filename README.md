Main page {#mainpage}
===================

This software provides some tools for estimating the ORCA sensitivity to neutrino mass ordering. It relies on Monte-Carlo data from ORCA simulation chain (see `data/README.md` for a brief description of the simulation chain).

Prerequisities
==============
* ROOT 6: tested at Lyon with /usr/local/root/6.10.02
* Jpp with aanet, compiled against ROOT 6: tested with v10.1.11007
* [OscProb](https://github.com/joaoabcoelho/OscProb) package needs be available and compiled and an environment variable `OSCPROBDIR` needs to be set to point to the `OscProb` directory.
* Python scripts will require `docopt` package

Documentation
=============
* Create with `doxygen doxyconf`, open `doxygen/html/index.html` in your favorite browser.
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

ORCA Monte-Carlo chain and data format
======================================

ORCA Monte-Carlo chain
----------------------

This directory stores the data necessary for the NMH sensitivity analysis. The raw starting point of the analysis is the file pid_result_XXX.root. The typical ORCA MC chain is as follows (see the following line in plain text in README.md):
```
gSeaGen->KM3Sim                JGandalf for tracks                                           
               \             /                      \                                       
		       JTE->         -  Recolns for tracks  -    -> merge, PID training -->summary  
               /             \                      /                                       
mupage->KM3                    Dusj reco for showers                                        
```
pid_result_XXX.root is the summary file, which contains all the information from the reco chains, plus the PID info and MC truth.

There are numerous gSeaGen and mupage files. We have *something* like this gSeaGen_<flavor>_<interaction>_<erange>_<runnr>.root, where flavor is muon/elec/tau, interaction is CC or NC, erange is 1--5 or 3--100, runnr is 1--1800.

This scheme persists until PID, however in the pid_result_XXX.root all flavours, interactions, energy ranges and run numbers are merged together, with special variables that help to re-trace the origin of each event. The merging is a necessary step for data input to ECAP Random Decision Forest PID.

For data sorting and quality purposes, having everything in one file is not optimal. Also, pid_result_XXX.root has about 300 branches, that occasionally change, depending on the version. Writing an analysis code that takes a changing and incomprehensible data format as input is not optimal. For that purpose the data in pid_result_XXX.root is converted to 'analysis' format (a smaller set of variables, organised tree, see below) and split up to match the file scheme used throughout the MC chain.

Analysis format
---------------

The analysis format is defined by the class ```NMH/common_software/SummaryEvent.h/C```. The class ```NMH/common_software/SummaryParser``` is set to read or write data in the analysis format by using the SummaryEvent class. The variables of the analysis format are described in SummaryEvent documentation.

Variables in tree
-----------------
* pos and dir refer to vertex position and neutrino direction, given for MC truth and 3 reco's.
* energy_nu and bjorkeny refer to the neutrion energy and Bjorken y, given for MC truth a 3 reco's (expect no bjorkeny for gandalf).
* Other variables are discussed in SummaryParser.h
* The variables <reco>_ql0(12) are described in detail below.

Quality levels
--------------

For each reconstruction there are branches called '<reco>_ql0, <reco>_ql1, ...' .  'ql' stands for 'quality level'.

* Quality level 0 is the lowest quality level, meaning somewhere in the reco chain it has been tested that the reconstruction worked in the most minimal way, e.g. that it did not return a NaN.

* Quality level 1 is the next lowest quality level. In this case the reconstruction has been tested against quality cuts set by the reconstruction specialists, as described on the wiki and the Moritz' script, linked in NMH/data/README.md. Additionally, in this case the event containment cut is loosened, such that the vertex can be up to 30 m outside the instrumented volume.

* Quality level 2 is the highest quality level. It is the same as above, but the vertex needs to be inside the instrumented volume.

We are usually (March 2018) advised to use ql1, as it selects more events and improves
statistics. For dusj there is no level 2 defined.