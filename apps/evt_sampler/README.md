Evt sampler
==============

Applications in this directory can be used to generate experiment samples from gSeaGen events. Each experiment sample will be stored in a root file in the *analysis format* (see `NMH/README.md`). Each sample will contain a sample of muon, elec, tau CC and NC events, as would be recorded by running the detector for a certain amount of time. The experiment samples are to be subjected to subsequent analysis to estimate sensitivy to neutrino mass hierarchy.

Prerequisities
==============
* The applications require summary data in *analysis format*, usually in `NMH/data/mcsummary/TAG/data_atmnu/`, created with the applications in `apps/data_sorting`.
* The scripts require access to gSeaGen simulation data.
* The `FluxChain` application will require an effective mass file that can be read by the `NMH/common_software/EffMass.h/C` class. Such a file is created by the applications in `apps/effective_mass`.

Directories
===========
`output/` - where most of the outputs are stored
`cache/`  - `GSGSampler` reads in all simulated gSeaGen data, which is time consuming. It caches the data in this directory to speed things up for the next time it is called. Deleting everything from the cache directory is completely fine, but will slow down running `GSGSampler` the next time.
`tmp/`    - created automatically by `sampler_caller.py` to store scripts and info related to computing on the farm.

How to run
==========

Sampling experiments (use python scripts)
-----------------------------------------
* First, knowledge of the interacted neutrino flux at the detector is required. This is calculated by the app `FluxChain`. Use the python script `flux_caller.py` (do `./flux_caller.py -h` for usage) to create a desired number of flux samples with different oscillation parameter values.

* As a result of running the `flux_caller.py`, several flux output files corresponding to different oscillation parameters are created in `output/FluxChain/`, alongside an `output/FluxChain/{IDSTR}_output_list.dat` file that lists the created files and a log file `output/FluxChain/{IDSTR}_log.dat`. 

The flux data can now be used as an input to `GSGSampler` app, which 1) reads in gSeaGen data 2) creates samples from gSeaGen data depending on the input from a flux file 3) searches the summary data to determine which of the sampled gSeaGen events end up reconstructed and identified.

`GSGSampler` can be called with the script `sampler_caller.py` (do `./sampler_caller.py -h` for usage). The file `output/FluxChain/{IDSTR}_output_list.dat` will be one of the required inputs. This script will run `GSGSampler` on the computing farm at Lyon. Run `GSGSampler` for all flavors and interactions. Once this finishes, there will be root files `output/GSGSampler/EvtSample_{flavor}-{CC/NC}_flux{F}_sample{N}.root`. F corresponds to the sequence index of the input flux file in `output/FluxChain/{IDSTR}_output_list.dat`, whereas N corresponds to the sample number created with this flux file (for each flux file several samples can/should be created, to study statistical fluctuations related to sampling). For NC events the flavor will be allflav.

* After the previous two steps the script `merge_to_exps.py` (do `./merge_to_exps.py -h` for usage) has to be used to create root files with samples representing experiments. It will require `output/FluxChain/{IDSTR}_output_list.dat` and `output/FluxChain/{IDSTR}_log.dat` as input.

After running `merge_to_exps.py`, there are files `output/Experiments/Experiment_oscpars{i}_sample{j}.root` available. The index i stands for a combination of oscillation parameters and j for a sample with these oscillation parameters. The oscillation parameters associated with different files are readily available in `output/Experiments/merge_log.dat`, and are also stored in the headers (read with `common_software/FileHeader`) of the experiment files. No further weighting of the events are required, the events in these files can be treated as been recorded by the detector over a certain running time.

Standalone (use ROOT macros directly)
-------------------------------------
* `FluxChain` can be run standalone, if one wishes to create one set of flux histograms, do `./FluxChain -h!` for usage. See app documentation for more details.

* `GSGSampler` can be run standalone, if one wishes to create an event sample for one flavor, do `GSGSampler -h!` for usage. See app documentation for more details.

Outputs
==========

* `FluxChain` - a root file with E vs ct histograms. The atmospheric flux (atmflux) histograms will show the atmospheric neutrino count in operation time in units 1/m2. The oscillated flux (oscflux) histograms will show the oscillated neutrino count in operation time in units 1/m2. The iteracted flux (intflux) histograms will show the interacted neutrino count in operation time in units 1/MTon (hence, this needs to be multiplied with Vcan*rho to get interacted neutrino count inside the detector). Detected flux histograms (detflux/) will show the detected neutrino count in operation time (unitless).

* GSGSampler` - a root file with events in summary format in a distribution as expected from running the detector for a certain amount of time.
