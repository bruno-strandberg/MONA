Data sorting
============
This directory has the code to sort the data in NMH/data/pid_result_XXX.root to analysis format.

Prerequisities
==============
* The file NMH/data/pid_result_XXX.root.
* The compiled library NMH/common_software/libcommonsoft.so.

Note that if the analysis format, defined in DataReducer.C, changes, one needs to update
NMH/common_software/SummaryParser.h/C before running RestoreParity.C.

How to run
==========
1. Execute the script reduce_data.py. This creates the file NMH/data/ORCA_MC_summary_all.root, which
   will contain all the same events as NMH/data/pid_result_XXX.root, but will be in the analysis
   format (reduced and organised tree).
   
2. Run root, .x RestoreParity.C+("../data/data/ORCA_MC_summary_all.root"). This will split the data
   up to several files in NMH/data/mc_end/.. to match the file naming and numbering scheme used
   throughout the MC chain. This step will take about 3 hours to run if done on single machine
   (It is currently very ineffective).