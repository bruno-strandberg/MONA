Data sorting
============
This directory has the code to sort the summary data from the Monte-Carlo chain (*ECAP PID*, see below) to an *analysis format* that is used by the software in `MONA/common_software`.

* *ECAP PID* is usually a flat `TTree` that contains all of the variables from reconstruction algorithms + the PID scores for all simulated Monte-Carlo events that pass the initial atmospheric-muon rejection cuts. More info (not always up-to-date) can be found at [https://wiki.km3net.de/index.php/ORCA_PID]. 

* *analysis format* is defined by the class `common_software/SummaryEvent.h/C`. More info is available in the master `MONA/README.md` file.

Prerequisities
==============
* The file `pid_result_XXX.root`, which refers to a PID summary file from ECAP. Such files can/could be found at `/sps/km3net/users/shallman/ORCA_PID_OUTPUT/04_18`, it is becoming a standard that these files are also available in iRODS. Typically, to know exactly which file to use, one must contact someone active in the ORCA working group.

* The configured and compiled library `MONA/common_software/`.

How to run
==========

This is a two-step process. First, the data from the PID summary file has to be converted to analysis format. This typically requires the creation of a new class to read the PID summary file (can be auto-generated by ROOT) and the creation of a new application that uses the reader and the `common_software/SummaryParser.h/C` class to map the data to *analysis format*. See section **Maintenance** below for instructions how to do so. Secondly, the data is to be split to several files by the application `RestoreParity` to match the file naming and numbering scheme used throughout the MC chain.

For example, the directory `parsers` contains a reader class `ECAP180401.h/C` and there is a corresponding application `ECAP180401_to_MONA.C` to convert the file `/in2p3/km3net/mc/atm_neutrino/KM3NeT_ORCA_115_23m_9m/v1.1.1/postprocessing/pid_180401/pid_output_atm_neutrino_atm_muon_pure_noise_shiftedVertexEventSelection_180401.root` to *analysis format*. The applications can be run as follows:

1. Converters can be run as `./ECAP180401_to_MONA -d ../../data/ -f ../../data/pid_result_XXX.root -t MY_TAG`, use `-h!` for more info. This creates the file `MONA/data/ORCA_MCsummary_SEvX_MY_TAG.root`, which will contain all the same events as `MONA/data/pid_result_XXX.root`, but in the *analysis format*. `X` is the version of the `SummaryEvent`.
   
2. Run `./RestoreParity -f ../../data/ORCA_MCsummary_SEvX_MY_TAG.root -l 1+5 -u 3+100`, use `-h!` for more info. This will split the data up to several files in `MONA/data/mcsumary/MY_TAG/...` to match the file naming and numbering scheme used throughout the MC chain. This will take ~10 minutes. This step is required for the applications in `apps/effective_mass` to work (and some other applications) and allows for more flexibility in choosing how much MC data to use for various studies.

Maintenance
===========

Typically, the format of the PID output tree from ECAP is not fixed and changes with every update. If the format of the tree `pid_result_XXX.root` changes, the `parsers/ECAP180401.h/C` class will no longer recongnize the `TTree` in the file. The application that uses the reader (`ECAP180401_to_MONA` in this case) will complain about missing/wrong branch names. Typically, a new reader and an application need to be created to convert a new PID summary file to *analysis format*. This can be achieved with the helper script `createconverter.py` or manually, both options are described below.

Helper script (recommended)
-------------
The python script `createconverter.py` can be used to create a `TTree` reader class and a skeleton application to convert the ECAP PID to *analysis format*. For example
~~~
./createconverter.py -s some_ecap_pid_file.root -n SomeReader -t PID -i iRODS_directory -c "comment"
~~~ 
will create files `SomeReader.h/C` and `SomeReaderToSummary.C`. The user has to modify the latter to define exactly how to map variables from the input root file to output file in *analysis format*. The branches that are available in the input file for transforming into summary format can be found in the header that is generated by `createconverter.py` and stored in `/parsers`. Once the user has edited `SomeReaderToSummary.C`, she/he should type `make` and execute `./SomeReaderToSummary -h!` for running guidance. The options `-i iRODS_directory` and `-c "comment"` are optional (recommended) and are used to automatically add the new application to the list of applications in this document.

Manual creation (not recommended)
----------------
Before the script `createconverter.py` was created, the following steps had to be done manually. Create a new reader by doing, e.g,
~~~
root new_pid_file.root
PID->MakeClass("NewReader")
~~~
Then, a new application `NewReader_to_MONA.C` can be created that uses the class `NewReader.h/C` to convert the data to analysis format. Practically, this can be achieved easily and quickly by doing
~~~
cp Other_to_MONA.C NewReader_to_MONA.C
~~~
and modifying the parts concerning the reader, data-mapping and documentation. Note that the new reader classes should be moved to the directory `parsers` for the make procedure to work.

The `ECAP180401_to_MONA` application is a simple example how to map *ECAP PID* format to *analysis format*. If at some point *ECAP PID* takes a different form (currently I don't know what the ECAP deep learning PID will use), this application exemplifies how to map the data to *analysis format*.

The `RestoreParity` application operates on the *analysis format* and hence should work on the output of `...to_MONA.C`, as long as the latter application has mapped the variables correctly.

Applications
=============
Application `ECAP180401_to_MONA.C` and parsers `parsers/ECAP180401.h/C` - comment: This application operates on the ECAP PID output from the indicated iRODS location to convert it into SummaryEvent format ; iRODS location of the ECAP PID output that the application operates on: /in2p3/km3net/mc/atm_neutrino/KM3NeT_ORCA_115_23m_9m/v1.1.1/postprocessing/pid_180401/pid_output_atm_neutrino_atm_muon_pure_noise_shiftedVertexEventSelection_180401.root

Application `ECAP181013_to_MONA.C` and parsers `parsers/ECAP181013.h/C` - comment: This application operates on the ECAP PID output from the indicated iRODS location to convert it into SummaryEvent format ; iRODS location of the ECAP PID output that the application operates on: /in2p3/km3net/mc/atm_neutrino/KM3NeT_ORCA_7_23m_9m/v0.0/postprocessing/pid_181013/pid_output_atm_neutrino_atm_muon_pure_noise_anyRecoUpgoing_181013.root

Application `ECAP190222_20m_to_MONA.C` and parsers `parsers/ECAP190222_20m.h/C` - comment: This application operates on the ECAP PID output from the indicated iRODS location to convert it into SummaryEvent format ; iRODS location of the ECAP PID output that the application operates on: /in2p3/km3net/mc/atm_neutrino/KM3NeT_ORCA_115_20m_9m/v5.0/postprocessing/pid_190222_v1.0

Application `ECAP190222_23m_to_MONA.C` and parsers `parsers/ECAP190222_23m.h/C` - comment: This application operates on the ECAP PID output from the indicated iRODS location to convert it into SummaryEvent format ; iRODS location of the ECAP PID output that the application operates on: /in2p3/km3net/mc/atm_neutrino/KM3NeT_ORCA_115_23m_9m/v5.0/postprocessing/pid_190222_v1.0
