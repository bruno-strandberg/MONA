Cross-section data
==================

This directory contains root files with neutrino cross-section data that is used by the class `common_software/NuXsec`. There are the following files:

* `xsec_gSeaGen_XXX.root` holds `TGraph`'s with neutrino interaction cross-section data. The data is created from GENIE/gSeaGen cross-section splines with `gslp2root` tool. The application in `MONA/xsec_extractor` can be used to re-create such a file, see `apps/xsec_extractor/README.md` for more info.

* `by_dists_gSeaGen_XXX.root` holds neutrino energy vs bjorken-y distributions. This file was created by the applications in `MONA/apps/bjorkeny-dists`, see `apps/bjorkeny-dists/README.md` for more info.

Legacy
=======

The file `xsec_gSeaGen_v4r1.root` is a legacy from M. Jongen's NMH package. That file can, in principle, also be re-created, but in practice the script `apps/xsec_extractor/genie_xsec_extractor.py` only works from gSeaGen v5r1 onwards because of issues with the `setenv.sh` script of gSeaGen (it uses csh for setting environment variables, whereas `genie_xsec_extractor.py` assumes bash). It has been checked that `xsec_gSeaGen_v4r1.root` and `xsec_gSeaGen_v5r1.root` yield identical results.