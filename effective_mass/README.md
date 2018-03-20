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
(do root, EffMhists.C+(...) ).

Assuming sorted data in NMH/data/... directory, a python script caller.py can be used to run
over several run numbers, flavors etc. Do caller.py -h for usage.

Outputs
==========

EffMhists.C outputs histograms for 'Detected' events (events in the summary file) and 'Generated'
events (corresponding events in gSeaGen file).

EffMass.C takes the EffMhists.C output as input and divides Detected/Generated to write out
effective mass histograms. The division is done in a separate macro to falicitate adding several
EffMhists.C outputs together with 'hadd'.