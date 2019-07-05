Cross-section extractor
=======================

The scripts/applications in this directory can be used to extract neutrino cross-section data from GENIE/gSeaGen into a ROOT file for usage by the `common_software/NuXsec` class. The GENIE tool `gspl2root` is used to extract the ROOT file from xml spline data file. The scripts rely on access to a working gSeaGen installation, typically available on cc-lyon.

The created root file should be stored in `data/xsec_data`, the `NuXsec` constructor can be pointed to use a specific file. If significant cross-section updates occur in GENIE, a new xsec file should be created and the default constructor should be changed to point to the updated file by updating `common_software/NuXsec.h`.

Note that it is not foreseen that this step needs to be repeated often, the tools are made available for self-consistency of the MONA package.

How to run
==========

* Do `make` to build the application `Reductor`.

* Run `genie_xsec_extractor.py` to create a new xsec data file (do `genie_xsec_extractor.py -h` for usage).

* The created file may be copied to `MONA/data/xsec_data`, if necessary the default constructor of `NuXsec` can be modified to default to the new file by editing `common_software/NuXsec.h`.

Outputs
==========

* A root file with xsec data from GENIE/gSeaGen.