Effective mass
==============

Programs in this directory can be used to create effective mass histograms. The output file of `Combine` application can be loaded by the `common_software/EffMass` class to calculate the effective masses.

Prerequisities
==============
* The scripts use the library in`common_software/` and Jpp headers.
* The scripts require ORCA simulation chain summary files in analysis format in `MONA/data/mcsummary/TAG/data_atmnu/`. These can be created with the `RestoreParity` applications in `apps/data_sorting/`.
* The scripts require ORCA simulation chain start (gSeaGen files) in `MONA/data/gseagen/TAG/data_atmnu`

How to run
==========

The application `EffMhists` creates effective mass histograms for one flavor for a given gSeaGen and summary file. Do `./EffMhists -h!` and read the doxygen doxumentation for more information.  NB: This application is applied to *one* summaryfile that is split up after `RestoreParity`. It can take multiple gSeaGen files as input, but only one summaryfile.

The script `EMH_caller.py` can be called to run `EffMhists` for all of the summary files in `MONA/data/mcsummary/TAG/data_atmnu/` directory, given that corresponding gSeaGen files are available.

The application `Combine` uses the `common_software/EffMass` class and combines all of the outputs from `EffMhists` to a single output file. That file can be loaded to another `EffMass` instance to calculate effective masses.

Optionally, the application `DataQuality` can be executed to create some data quality plots from the outputs of `EffMhists`.

Outputs
==========

`EffMhists` outputs 'generated' (gSeaGen events) and 'selected' (summary events) histograms for the given input files.

`Combine` created four combined outputs (`elec-CC, muon-CC, tau-CC, elec-NC`) from the outputs of `EffMhists` and additionally one file than is to be used with the `EffMass` class.
