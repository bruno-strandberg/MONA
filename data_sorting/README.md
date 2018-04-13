Data sorting
============
This directory has the code to sort the data in NMH/data/pid_result_XXX.root to "analysis format". The analysis format is defined in the function DataReducer::InitOutputTree(), which maps the tree in pid_result_XXX.root to a smaller tree with better variable naming convention. The class NMH/common_software/SummaryParser is set to parse data in the analysis format, the variables of the analysis format are described in SummaryParser documentation.

Prerequisities
==============
* The file NMH/data/pid_result_XXX.root.
* The compiled library NMH/common_software/libnmhsoft.so.

Maintenance
===========

If the analysis format, defined in DataReducer::InitOutputTree(), changes, one needs to update NMH/common_software/SummaryParser.h/C before running RestoreParity.C. This will also affect other applications (e.g. effective mass calculation) that use SummaryParser.

If the format of the tree pid_result_XXX.root changes, DataReducer class needs to be updated. If the PID tree in pid_result_XXX.root has changed, DataReducer will complain about missing/wrong branch/variable names. The "simplest" way to update is to do ```root new_pid_result_XXX.root, PID->MakeClass(NewDataReducer)``` and copy over the parts added by bstr. This is unpleasant, but caused by lack of agreement in the format of the summary data.

How to run
==========
1. Execute ```root, ReduceData.C+("../data/pid_result_XXX.root","../data/ORCA_MC_summary_all_XXX.root")```. This creates the file NMH/data/ORCA_MC_summary_all_XXX.root, which will contain all the same events as NMH/data/pid_result_XXX.root, but will be in the analysis format (reduced and organised tree).
   
2. Run root, .x RestoreParity.C+("../data/data/ORCA_MC_summary_all_XXX.root"). This will split the data up to several files in NMH/data/mc_end/.. to match the file naming and numbering scheme used throughout the MC chain. This will take ~10 minutes.