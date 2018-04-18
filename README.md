Main page {#mainpage}
===================
Work in progress.

Prerequisities
==============
* ROOT 6
* To use GSGParser to parse gSeaGen files in .evt format, one also requires Jpp and aanet root6 branch.

How to use
==========

Setup
-----
* ```source setenv.sh```. Optionally, do ```source setenv.sh -a``` to compile against aanet.

Sort data
----------
The starting point is the file NMH/data/pid_result_XXX.root. This needs to be converted to analysis format. Data sorting/conversion is handled by macros in NMH/data_sorting/. Consult NMH/data_sorting/README.md

Effective mass
--------------
2D effective mass histograms can be created with the scripts in NMH/effective_mass. Consult NMH/data_sorting/README.md. The combined outputs of the effective mass scripts are stored in data/eff_mass. There should be a file per each neutrino flavor, CC and electron NC.

Developments
============

For improvements
----------------

* NuXsec uses data in a specific format. Would be nice to have something that can use GENIE output directly (not urgent).

* AtmFlux::ReadHondaFlux() makes a small approximation for the ranges abs(cosz) > 0.95 (not sure this can be easily improved, though).

* For both NuXsec and AtmFlux I foresee an option to add a scaler graph. The scaler graph can be initiated from a function or from a graph. It should provide a simple way to to scale/skew the crossection/flux.

To-do, ideas
------------

* You could try to estimate the effective mass from the effective area, in which case you won't need to parse the gSeaGen files.

* Once the event weight calculation is in place, it might be the time to give another thought to the data format. It may be smart to have a way to attach the weight to the summary files? Needs thought, data sorting-->effective mass-->event weight-->data sorting may not be ideal.

* SummaryParser should have a TChain instead of TTree to allow attaching several files.

* It may be an idea that SummaryParser uses some event class instead of a flat tree. Needs thought, though.

* Interactive python setup script that asks for OscProb, Jpp with root 6, etc, paths and creates the setenv.sh script.

Done
----

* Need aanet to parse tau gSeaGen files and calculate the effective mass

* allow compilation without aanet/jpp

* the test tests/flux_check_inpol.py illustrates that the interpolation works nicely in cosz and energy directions, no reason to put energy on log scale.