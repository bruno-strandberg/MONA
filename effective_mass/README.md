Effective mass
==============

Scripts in this directory can be used to create effective mass histograms.

Prerequisities
==============
* The scripts use the code in NMH/common_software/.
* The scripts require ORCA simulation chain summary files in a specific format in NMH/data/mc_end
* The scripts require ORCA simulation chain start (gSeaGen files) in NMH/data/mc_start/

How to run
==========

Both EffMhists.C and EffMass.C macros can be run stand-alone. Run them in compiled mode
(do root, EffMhists.C+(...) ). Interfaces are documented in the code.

Assuming sorted data in NMH/data/... directory, a python script EMH_caller.py can be used to run
over several run numbers, flavors, energies, interactions. Do ```EMH_caller.py -h``` for help. It can also send jobs to the farm to save time. Combinations that do not exist (e.g. muon-NC) are ignored.

After EMH_caller.py has been run to produce effective mass histograms in output/ for each flavor, energy, interaction, the script EM_caller.py can be called. This will identify missing files, summarize the data and create effective mass histograms to combined_output.

Outputs
==========

EffMhists.C outputs histograms for 'Detected' events (events in the summary file) and 'Generated'
events (corresponding events in gSeaGen file).

EffMass.C takes the EffMhists.C output as input and divides Detected/Generated to write out
effective mass histograms. The division is done in a separate macro to falicitate adding several
EffMhists.C outputs together with 'hadd'.