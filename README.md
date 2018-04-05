NMH code start page
===================
Work in progress.

Prerequisities
==============
* ROOT 6
* Jpp ROOT 6 branch with aanet




Developments
============

For improvements
----------------

* CScalculator uses data in a specific format. Would be nice to have something that can use GENIE output directly (not urgent).

* AtmFlux::ReadHondaFlux() makes a small approximation for the ranges abs(cosz) > 0.95 (not sure this can be easily improved, though).

* For both CScalculator and AtmFlux I foresee an option to add a scaler graph. The scaler graph can be initiated from a function or from a graph. It should provide a simple way to to scale/skew the crossection/flux.

* Maybe allow compilation without aanet/jpp

To-do, ideas
------------

* You could try to estimate the effective mass from the effective area, in which case you won't need to parse the gSeaGen files.

* Once the event weight calculation is in place, it might be the time to give another thought to the data format. It may be smart to have a way to attach the weight to the summary files? Needs thought, data sorting-->effective mass-->event weight-->data sorting may not be ideal.

* SummaryParser should have a TChain instead of TTree to allow attaching several files.

Done
----

* Need aanet to parse tau gSeaGen files and calculate the effective mass