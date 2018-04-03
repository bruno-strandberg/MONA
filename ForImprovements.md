For improvements
================

* CScalculator uses data in a specific format. Would be nice to have something that can use GENIE output directly (not urgent).

* AtmFlux::ReadHondaFlux() makes a small approximation for the ranges abs(cosz) > 0.95 (not sure this can be easily improved, though).

* For both CScalculator and AtmFlux I foresee an option to add a scaler graph. The scaler graph can be initiated from a function or from a graph. It should provide a simple way to to scale/skew the crossection/flux.

To-do, ideas
=============

* Need aanet to parse tau gSeaGen files and calculate the effective mass

* Alternatively, you could try to estimate the effective mass from the effective area, in which case you won't need to parse the gSeaGen files. At this stage you should also be ready to compare the oscillated event rates calculated with effective mass vs effective area. You will still need to get the header info from the .evt file, though...

* There may be a reason to have the AtmFlux data on a logarithmic energy scale. If the flux fall-off mimics exponential decay, then on a log scale this looks linear. Hence, the interpolation function is more accurate at higher energies, where there are fewer points. Should be trivial to implement. It may be a good idea to have the function still take the energy on a linear scale and take the log in the function to avoid mix-up between log2, log10

* Once the event weight calculation is in place, it might be the time to give another thought to the data format. It may be smart to have a way to attach the weight to the summary files? Needs thought, data sorting-->effective mass-->event weight-->data sorting may not be ideal.

* SummaryParser should have a TChain instead of TTree to allow attaching several files.