Main page {#mainpage}
===================

What does it do?
================

This software provides some tools for estimating the ORCA sensitivity to neutrino mass ordering and other ORCA analyses. It relies on Monte-Carlo data from ORCA simulation chain (discussed separately in section **ORCA Monte-Carlo chain and data format**). Broadly put, it provides tools to filter experiment data (either simulated or sea-data) and create a detector response to fit the filtered experiment data. Additionally, it provides a few small classes for the calculation of the number of expected atmospheric neutrino interactions at ORCA site.

Directory structure
===================
* `common_software/` - this directory holds the software that provides the base functionality of the package. It has the classes for the data format, atmospheric neutrino count calculation and data filtering. The software is documented carefully in `common_software/README.md`
* `fitter_software/` - software in this directory is built upon the software in `common_software/` and `RooFit` to provide tools for NMO fits. Documentation in `fitter_software/README.md`
* `apps/` - this directory has sub-directories for various applications that are built on top of `common_software` and `fitter_software`. The applications should be considered relatively dynamic and expected to develop/change/grow as the analyses become more complicated. Each application directory comes with a README-file that describes the purposes of the application. If the README-file is missing, it means the applications are not ready/meant for wider use.
* `macros/` - this directory has sub-directories for various analyses that are built using the `common_software/` and `fitter_software/` libraries. The `macros/` directory is the most dynamic working directory of the repository and typically documentation is sparse.
* `data/` - this directory holds data that is necessary for cross-section and atmospheric flux calculations. Additionally, it has sub-directories for Monte-Carlo data storage. More info in `data/README.md`.
* `doxygen/` - directory that is auto-populated by `doxygen` for documentation.
* `tests/` - directory that holds some example scripts/applications that also act as tests.

Prerequisities
==============
* ROOT 6: tested at Lyon with /usr/local/root/6.10.02
* Jpp with aanet, compiled against ROOT 6: tested with v10.1.11007 at Lyon
* [OscProb](https://github.com/joaoabcoelho/OscProb) package needs be available and compiled and an environment variable `OSCPROBDIR` needs to be set to point to the `OscProb` directory.
* Python scripts will require `docopt` package

Documentation
=============
* Create with `doxygen doxyconf`, open `doxygen/html/index.html` in your favorite browser (e.g. `firefox doxygen/html/index.html`)
* Each subdirectory holds a `README.md` file that describes what the code can be used for.

Versioning
==========
The numbers in the version string `vA.B.C-D`, where A to D are numbers, have the following meaning:
* A - version number for the `SummaryEvent` class only. This defines the data-format and should change rarely (but should be changed, if necessary). As data-format changes are likely to affect `apps/` and `macros/`, a specific number is held just for the data format.
* B - version number for the library in `common_software/`. Whenever significant changes occur, this should be incremented by 1. No upper limit, i.e. this can become larger that 9.
* C - version number for the library in `fitter_software/`, similar guidance as for B.
* D - a small version number for any small changes/improvements in B or C that did not require B or C to be incremented. This can grow indefinitely until A, B or C is increased, then reset to 0. For example, this number is useful for releasing certrain applications or macros as part of publishing plots.

Setup in lyon and elsewhere
===========================
* Add the following to your shell environment
```
export OSCPROBDIR=/my/path/to/oscprob/
source /my/path/to/NMH/setenv.sh
```
* Build the package by typing `make` in the master directory

Setup for usage
===============
This really depends on what is the targeted use of the software. For example, for the calculation of the atmosperic neutrino flux, the class `common_software/AtmFlux.h/C` can be used on the `root` prompt after compilation and nothing else is required. However, for something more advanced that requires ORCA MC data, the following tasks are usually required.

Data format conversion
----------------------
When Monte-Carlo data is in the equation, it needs to be in the analysis format, defined by `common_software/SummaryEvent.h/C` (reasons for this are discussed in detail in section **ORCA Monte-Carlo chain and data format** below). For example, consider ORCA MC summary file from ECAP from April 18 `/sps/km3net/users/shallman/ORCA_PID_OUTPUT/04_18/pid_result_shiftedVertexEventSelection.root`. The applications in `apps/data_sorting/` can be used to perform the conversion, see `apps/data_sorting/README.md` for more info.

Effective mass
--------------
*Effective mass* is required to predict the number of *selected* neutrino events in the detector in some time period. Typically, *selected* is defined as the events that would pass the trigger and simple atmospheric muon rejection cuts and end up in the ECAP PID output tree. The class `common_software/EffMass.h/C` provides a class that performs effective mass calculations, given an input file with histograms for *selected* (ECAP PID events) and *generated* (gSeaGen) events. The applications in `apps/effective_mass` can be used to create such a file, see `apps/effective_mass/README.md` for more info.

Bjorken-y distributions
------------------------
To distribute expectation values for *selected* events in bjorken-y (in addition to the conventional energy and cos-theta), knowledge of a 2D neutrino energy vs bjorken-y distribution is required. Such distributions can be generated from gSeaGen data, the applications in `apps/bjorkeny_dists` create such distributions, see `apps/bjorkeny_dists/README.md` for more info. Note that one such distribution file was generated from gSeaGen v4r1 data and comes with the repo. If significant updates are expected to the cross-section calculation in gSeaGen (this depends on the underlying GENIE version), the `apps/bjorkeny_dists` programs should be run again to update the distribution file `data/cross_sections_gSeaGen_v4r1/by_dists.root`. Otherwise, the existing file can be used. The `by_dists.root` file is used by `common_software/NuXsec.h/C` class.

Available applications
----------------------
Other available applications are in `apps/evt_sampler` and `apps/fitter`. The former directory contains applications to create data samples from gSeaGen data that are distributed as experiment data after a selected number of years. The latter directory contains applications to perform various NMO fits and analyses, using both `common_software` and `fitter_software`.

ORCA Monte-Carlo chain and *analysis format*
============================================

ORCA Monte-Carlo chain
----------------------

The typical ORCA MC chain consists of the applications and steps as illustrated below:
```
gSeaGen->KM3Sim                JGandalf for tracks                                           
               \             /                      \                                       
		       JTE->         -  Recolns for tracks  -    -> merge, PID training -->summary  
               /             \                      /                                       
mupage->KM3                    Dusj reco for showers                                        
```
The summary files come from ECAP and are described in more detail at https://wiki.km3net.de/index.php/ORCA_PID. Basically it is one huge root `TTree` with all of the simulated events that contains all the information from the reconstruction algorithms, plus PID info and MC truth.

There are numerous gSeaGen and mupage files. We have *something* like `gSeaGen_<flavor>_<interaction>_<erange>_<runnr>.root`, where flavor is `muon/elec/tau`, interaction is `CC` or `NC`, erange is `1-5` or `3-100`, runnr is `1-600`.

This scheme persists until PID. However, in the PID summary file all flavours, interactions, energy ranges and run numbers are merged together, with special variables that help to re-trace the origin of each event. The merging is a necessary step for data input to machine-learning PID. On the other hand, for data sorting and quality purposes, having everything in one file is not always optimal (e.g. for effective mass calculations). Also, the ECAP PID output tree has about 300 branches, which occasionally change, depending on the version. Writing an analysis code that takes a changing and incomprehensible data format as input is not optimal. Because of this, the ECAP PID summary data is converted to *analysis format* (a smaller set of variables, organised tree, see below) and split up to match the file scheme used throughout the MC chain. Importantly, event filtering and the creation of corresponding detector responses, which is one of the base functionalities of the package, could not have been implemented without a well-defined format for the data.

Analysis format
---------------
The analysis format is defined by the class `NMH/common_software/SummaryEvent.h/C`. The class `NMH/common_software/SummaryParser.h/C` is set to read or write data in the analysis format by using the `SummaryEvent` class. The variables of the analysis format are described in `SummaryEvent` documentation.

### Variables in tree
* pos and dir refer to vertex position and neutrino direction, respectively, and are provided for MC truth and 2 reco's (track and shower).
* energy and bjorkeny refer to the neutrion energy and Bjorken y, respectively, given for MC truth and 2 reco's (bjorkeny for tracks is currently empty).
* there are variables for various PID scores.
* For further info, see the documentation for `SummaryParser.h`
* The variables `<reco>_ql0(12)` define *quality levels* and are described in detail below.

### Quality levels

For each reconstruction there are branches called `<reco>_ql0, <reco>_ql1, ...`, where  `ql` stands for *quality level*.

* Quality level 0 is the lowest quality level, meaning somewhere in the reco chain it has been tested that the reconstruction worked in the most minimal way, e.g. that it did not return a NaN.

* Quality level 1 is the next lowest quality level. In this case the reconstruction has been tested against quality cuts set by the reconstruction specialists, as described on the wiki and the Moritz' script. Additionally, in this case the event containment cut is loosened, such that the vertex can be up to 30 m outside the instrumented volume.

* Quality level 2 is the highest quality level. It is the same as above, but the vertex needs to be inside the instrumented volume.

For ECAP April 2018 MC data We are usually advised to use ql1, as it selects more events and improves statistics. For shower there is no level 2 defined.

The quality levels are populated by the user in applications `apps/data_sorting/Alpha(Beta...)ToSummary`. In the end, it is up to the user to choose what sort of qualities these flags are used for. There is freedom to create custom quality cuts when PID tree is converted to summary formate, e.g. the user can create a new variable in the summary event 
```
Double_t fTrack_ql3 = (nhits > 30)
```
and use it when filtering events for `EventSelection` and `DetResponse`.
