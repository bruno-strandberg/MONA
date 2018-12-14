Data sorting
============
This directory has the code to sort the summary data from the Monte-Carlo chain (*ECAP PID*, see below) to a *analysis format* that is used by the software in `NMH/common_software`. 

* *ECAP PID* is a usually a flat `TTree` that contains all of the variables from reconstruction algorithms + the PID scores for all simulated Monte-Carlo events that pass the initial atmospheric-muon rejection cuts. More info (not always up-to-date) can be found at [https://wiki.km3net.de/index.php/ORCA_PID]. 

* *analysis format* is defined by the class `common_software/SummaryEvent.h/C`. More info in the master `NMH/README.md` file.

Prerequisities
==============
* The file `pid_result_XXX.root`, which refers to a PID summary file from ECAP. Despite my several requests these files have not been uploaded to iRODS. Such files can/could be found at `/sps/km3net/users/shallman/ORCA_PID_OUTPUT/04_18`. Typically, to know exactly which file to use, one must contact someone active in the ORCA working group.

* The compiled library `NMH/common_software/`.

How to run
==========

1. Converter can be run as `./RDFPID_to_Summary -d ../../data/ -f ../../data/pid_result_XXX.root -t MY_TAG`, use `-h!` for more info. This creates the file `NMH/data/ORCA_MC_summary_MY_TAG.root`, which will contain all the same events as `NMH/data/pid_result_XXX.root`, but in the *analysis format*.
   
2. Run `RestoreParity -f ../../data/ORCA_MC_summary_MY_TAG.root -l 1+5 -u 3+100`, use `-h!` for more info. This will split the data up to several files in `NMH/data/mcsumary/MY_TAG/...` to match the file naming and numbering scheme used throughout the MC chain. This will take ~10 minutes. This step is required for the applications in `apps/effective_mass` to work (and some other applications) and allows for more flexibility in choosing how much MC data to use for various studies.

Maintenance
===========

If the format of the tree `pid_result_XXX.root` changes, `RDFPIDReader` class needs to be updated. If the PID tree in `pid_result_XXX.root` has changed, `RDFPIDReader` will complain about missing/wrong branch/variable names. The simplest way to update is to do `root new_pid_result_XXX.root, PID->MakeClass("RDFPIDReader")` and copy over the documenting header in `RDFPIDReader.h`.

The `RDFPID_to_Summary` application is a simple example how to map *ECAP PID* format to *analysis format*. If at some point *ECAP PID* takes a different form (currently I don't know what the ECAP deep learning PID will use), this application exemplifies how to map the data to *analysis format*.

