NMH analysis data directory
===========================


Directories and files
=====================

The directory structure (mcsummary/, gseagen/, ...) comes with the git repo. For most stuff to
work, one requires sorted data in mcsummary/... and gseagen/ data.

* pid_result_XXX.root          - summary of the mc events and reconstructed/pid'd mc events, ORCA detector.

* ORCA_MC_summary_XXX.root     - same events as pid_result.root, but converted to analysis format.

* mcsummary/                   - same data os ORCA_MC_summary_XXX.root, separated to files as in the rest of the MC chain.

* gseagen/                    - gSeaGen files, should be 1-1 correspondence with mcsummary/

* cross_sections_gSeaGen_v4r1/ - neutrino cc/nc cross-sections on n, p, O16 from GENIE. Gift from M. Jongen. 

* honda_flux/                  - honda flux tables for the atmospheric flux class.

* eff_mass/                    - files with effective mass histograms from the applications, NMH/effective_mass

* testing/                     - data used by/for NMH/tests scripts.

pid_result_XXX.root file versions
=================================

See https://wiki.km3net.de/mediawiki/index.php/ORCA_PID for info regarding various versions.

Info regarding the tree in pid_result_XXX.root
==============================================

The info in the tree is as follows (according to Liam, thanks!):

MC info
--------
* dir_x, dir_y, dir_z     -  for MC direction
* pos_x, pos_y, pos_z     -  for MC vertex
* Erange_min, Erange_max  -  tell you the lower and upper limit of the MC neutrinos
* bjorkeny                - mc bjorken y
* energy                  - mc energy
* is_neutrino             - separate neutrino/atm. muon events
* is_cc                   - separate cc/nc events
* type                    - neutrino type by PDG

Reco info
----------
Direction and position:
* reconame_pos(dir)_x(yz)

Neutrino energy:
* gandalf(dusj)_energy_corrected
* recolns_energy_neutrino

Bjorken y:
* dusj_best_DusjOrcaUsingProbabilitiesFinalFit_BjorkenY
* recolns_bjorken_y
* Not output by gandalf

Selection:
* reconame_is_good     - whether the reco converged (i.e. did not give NaN)
* reconame_is_selected - whether the reco passes selection cuts, as discussed in [selection cuts](https://wiki.km3net.de/mediawiki/index.php/Simulations/ORCA_productions#Default_Event_Selection_Cuts) and [Moritz' script](http://git.km3net.de/moritz/beluga/blob/master/beluga/cut_sets.py).
* reconame_loose_is_selected - a looser containment cut for event selection.