Data sorting
============
This directory has the code to sort the data in NMH/data/pid_result_XXX.root to "analysis format" (see NMH/data/README.md for more info).

Prerequisities
==============
* The file NMH/data/pid_result_XXX.root.
* The compiled library NMH/common_software/libnmhsoft.so.

Maintenance
===========

If the format of the tree pid_result_XXX.root changes, ```RDFPIDReader``` class needs to be updated. If the PID tree in pid_result_XXX.root has changed, ```RDFPIDReader``` will complain about missing/wrong branch/variable names. The simplest way to update is to do ```root new_pid_result_XXX.root, PID->MakeClass(RDFPIDReader)``` and copy over the documenting header in ```RDFPIDReader.h```.

The RDFPID_to_Summary.C macro is a simple example how to map the RDFPID data format to summary data format. If at some point summary data takes a different form (currently I don't know what the ECAP deep learning PID will use), this macro serves as an example how to map the data to the analysis format.

How to run
==========
1. Execute ```root, RDFPID_to_Summary.C+("../data/pid_result_XXX.root","../data/ORCA_MC_summary_all_XXX.root")```. This creates the file NMH/data/ORCA_MC_summary_all_XXX.root, which will contain all the same events as NMH/data/pid_result_XXX.root, but will be in the analysis format (reduced and organised tree).
   
2. Run root, .x RestoreParity.C+("../data/data/ORCA_MC_summary_all_XXX.root"). This will split the data up to several files in NMH/data/mc_end/.. to match the file naming and numbering scheme used throughout the MC chain. This will take ~10 minutes.