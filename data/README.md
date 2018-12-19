NMH data directory
==================

This is the data directory distributed with the repo. The sub-directories are described below.

* `cross_sections_gSeaGen_v4r1/` - neutrino cc/nc cross-sections on n, p, O16 from GENIE. Gift from M. Jongen. Bjorken-y distributions, created with `NMH/apps/bjorkeny_dists`, to help distribute events into bjorken-y bins. Files in this directory are used by `NMH/common_software/NuXsec.h/C` class.

* `honda_flux/` - honda flux tables for the class `NMH/common_software/AtmFlux.h/C`.

* `eff_mass/` - directory to store effective mass files, which are required for `NMH/common_software/EffMass.h/C`. The directory comes with the repo, but needs to be populated by running the applications in `NMH/apps/effective_mass/`.

* `mcsummary/` - Monte-Carlo chain summary data in `NMH/common_software/SummaryEvent.h/C` format. The directory comes with the repo, but needs to be populated by running the applications in `NMH/apps/data_sorting/`.

* `gseagen/` - gSeaGen files corresponding to the Monte-Carlo summary data in `mcsummary/`. The directory comes with the repo, but it is the users responsibility to copy `gSeaGen` files to this directory (note that, of course, there is no obligation to use this directory if the files are already centrally available somewhere).


ECAP PID output format
==============================================

File versions
--------------
See https://wiki.km3net.de/mediawiki/index.php/ORCA_PID for info regarding various versions.

Info in the ECAP PID `TTree`

###MC info
* `dir_x, dir_y, dir_z`     -  for MC direction
* `pos_x, pos_y, pos_z`     -  for MC vertex
* `Erange_min, Erange_max`  -  tell you the lower and upper limit of the MC neutrinos
* `bjorkeny`                - mc bjorken y
* `energy`                  - mc energy
* `is_neutrino`             - separate neutrino/atm. muon events
* `is_cc`                   - separate cc/nc events
* `type`                    - neutrino type by PDG

###Reco info
* `reconame_pos(dir)_x(yz)` -  for direction and position

* `gandalf(dusj)_energy_corrected, recolns_energy_neutrino` - for neutrino energy

* `dusj_best_DusjOrcaUsingProbabilitiesFinalFit_BjorkenY, recolns_bjorken_y` - for bjorken-y (nothing from gandalf, yet)

Selections:
* `reconame_is_good`     - whether the reco converged (i.e. did not give NaN)
* `reconame_is_selected` - whether the reco passes selection cuts, as discussed in [selection cuts](https://wiki.km3net.de/mediawiki/index.php/Simulations/ORCA_productions#Default_Event_Selection_Cuts) and [Moritz' script](http://git.km3net.de/moritz/beluga/blob/master/beluga/cut_sets.py).
* `reconame_loose_is_selected` - a looser containment cut for event selection.