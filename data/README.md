NMH analysis data directory
===========================

This directory stores the data necessary for the NMH sensitivity analysis. The raw starting point
of the analysis is the file pid_result_XXX.root. The typical ORCA MC chain is as follows (see the
following line in plain text in README.md):

gSeaGen->KM3Sim             JGandalf for tracks                                           
                 \        /                      \                                       
		   JTE-> -  Recolns for tracks      - -> merge, PID training -->summary  
                 /        \                      /                                       
mupage->KM3                 Dusj reco for showers                                        

pid_result_XXX.root is the summary file, which contains all the information from the reco chains,
plus the PID info and MC truth.

There are numerous gSeaGen and mupage files. We have *something* like this
gSeaGen_<flavor>_<interaction>_<erange>_<runnr>.root, where flavor is muon/elec/tau, interaction
is CC or NC, erange is 1--5 or 3--100, runnr is 1--1800.

This scheme persists until PID, however in the pid_result_XXX.root all flavours, interactions,
energy ranges and run numbers are merged together, with special variables that help to re-trace
the origin of each event.

For data sorting and quality purposes, having everything in one file is not optimal. Also,
pid_result_XXX.root has about 300 branches, that occasionally change, depending on the version.
Writing an analysis code that takes a changing and incomprehensible data format as input is not
optimal. For that purpose the data in pid_result_XXX.root is converted to 'analysis' format
(a smaller set of variables, organised tree) and split up to match the file scheme used throughout
the MC chain. Data sorting and the analysis format is discussed further in
NMH/data_sorting/README.md .

NB! The directory structure (mc_end/, mc_start/, ...) comes with the git repo. For most stuff to
work, one requires sorted data in mc_end/... and mc_start.


Directories and files
=====================

* pid_result_XXX.root          - summary of the mc events and reconstructed/pid'd mc events,
			         ORCA detector.

* ORCA_MC_summary_all_XXX.root - same events as pid_result.root, but converted to analysis format

* mc_end/                      - same data os ORCA_MC_summary_all.root, separated to files as in the
			         rest of the MC chain

* mc_start/                    - gSeaGen files, should be 1-1 correspondence with mc_end

* cross_sections_gSeaGen_v4r1/ - neutrino cc/nc cross-sections on n, p, O16 from GENIE. Gift from
                                 M. Jongen. 

* honda_flux/                  - honda flux tables for the atmospheric flux class.

* eff_mass/                    - files with effective mass histograms from the applications
  			       	 NMH/effective_mass

* testing/                     - data used by/for NMH/tests scripts.

pid_result_XXX.root file versions (from Steffen, thanks!)
=================================

pid_result_8Feb2018.root taken from
/sps/km3net/users/shallman/ORCA_PID/nemowater_with_mx_01_18/shiftedVertexEventSelection/pid_output/
See Steffen's email titled 'Fwd: ORCA PID update'.

pid_result_14Mar2018.root copied from
/sps/km3net/users/shallman/ORCA_PID_OUTPUT/pid_result_shiftedVertexEventSelection.root
More info in, 14.03.2018 version
https://wiki.km3net.de/mediawiki/index.php/ORCA_PID

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

PID
---

### PID 8Feb2018:
   * muon_probability  - separate atm. muon and neutrinos
   * track_probability - separate tracks and showers

### PID 14Mar2018:
   * muon_score  - separate atm. muon and neutrinos
   * track_score - separate tracks and showers
