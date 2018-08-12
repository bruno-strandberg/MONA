Main page {#mainpage}
===================
This software provides some tools for estimating the ORCA sensitivity to neutrino mass ordering. It relies on Monte-Carlo data from ORCA simulation chain (see `data/README.md` for a brief description of the simulation chain). The directory `NMH/common_software` holds classes that are used by various applications for recurring tasks. Other directories hold code for some specific analysis.

Prerequisities
==============
* ROOT 6
* To use `GSGParser` to parse `gSeaGen` files in `.evt` format, one also requires `Jpp` and `aanet` compiled against `ROOT 6`
* Python scripts will require `docopt` package

More info
=========
* Regarding the data format(s) used in the analysis, consult NMH/data/README.md and NMH/data_sorting/README.md.
* Each subdirectory holds a `README.md` file that describes what the code can be used for.

How to use
==========

Setup in lyon and elsewhere
-------------
* TBD
* ```source setenv.sh -a -o /oscprob/dir```. `-a` sets to compile against `aanet` (only required when reading `gSeaGen` files in `.evt` format). `-o` sets the `OscProb` directory that is required for oscillating neutrinos.
* ```cd common_software/ ; make ```

Sort data
----------
The starting point is the file `NMH/data/pid_result_XXX.root`. This needs to be converted to analysis format. Data sorting/conversion is handled by macros in `NMH/data_sorting/`. Consult `NMH/data_sorting/README.md`

Effective mass
--------------
2D effective mass histograms can be created with the scripts in `NMH/effective_mass`. Consult `NMH/data_sorting/README.md`. The combined outputs of the effective mass scripts are stored in `data/eff_mass`. There should be a file per each neutrino flavor, CC and electron NC.

Create experiment samples
-------------------------
The code in directory `NMH/evt_sampler/` can be used to create samples of Monte-Carlo data that look like event recorded by the detector over a certain amount of time. See `NMH/evt_sampler/README.md` for further information.

Asymmetry
---------
The code in directory `NMH/asymmetry` can be used to estimate the asymmetry between normal and inverted hierarchy. Documentation is pending.


Developments
============

For improvements
----------------

* Instead of mc_end, mc_start, you should have some discriminating name, eg mc_end_ECAP_PID_XXX, mc_start_LOI

* AtmFlux option enumerator should be part of the class, as in EventSelection.

* NuXsec uses data in a specific format. Would be nice to have something that can use GENIE output directly (not urgent).

* AtmFlux::ReadHondaFlux() makes a small approximation for the ranges abs(cosz) > 0.95 (not sure this can be easily improved, though).

* For both NuXsec and AtmFlux I foresee an option to add a scaler graph. The scaler graph can be initiated from a function or from a graph. It should provide a simple way to to scale/skew the crossection/flux.

* GSGSampler.C should output trees with GSG event IDs. These trees should be fed to a different application, that reads in the summary data and GSG samples and identifies events that made it to the end of the chain. In such a way you could disentangle the problem associated with the changing summary file format.

* You should make a separate application that cache's the gSeaGen data.

To-do, ideas
------------

* You could try to estimate the effective mass from the effective area, in which case you won't need to parse the gSeaGen files.

* Interactive python setup script that asks for OscProb, Jpp with root 6, etc, paths and creates the setenv.sh script.

Done
----

* Need aanet to parse tau gSeaGen files and calculate the effective mass --> compiling against aanet

* allow compilation without aanet/jpp --> controlled by setenv.sh script

* the test tests/flux_check_inpol.py illustrates that the interpolation works nicely in cosz and energy directions, no reason to put energy on log scale.

* SummaryParser should have a TChain instead of TTree to allow attaching several files --> done

* The can size, muon cut etc should be stored in the output of effective mass --> handled by FileHeader

* The above point applied also more generally: the input parameters should be saved to output for different applications, etc FluxChain, GSGSampler, ... --> handled by FileHeader

* pid_result conversion to summary format is still painful, especially when summary data format changes. One option would be to define summary data structure in a SummaryEvent() class. This would, however, break compatibility with most things I have written so far and, more importantly, will make accessing/cloning/etc of the data more difficult. Maybe a better option would be to have a script that, after DataReducer.C has run its script, it checks for branches in SummaryParser.h and advises on branches to be added? --> Improved by creating the SummaryEvent class, which provides better IO.

* It may be an idea that SummaryParser uses some event class instead of a flat tree -> uses SummaryEvent.
