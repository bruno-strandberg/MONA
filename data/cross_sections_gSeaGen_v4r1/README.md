Cross-section data
==================

This directory contains root files with neutrino cross-section data that is used by the class `NMHDIR/common_software/NuXsec`. There are two files

* `xsec.root` holds `TGraph`'s with neutrino interaction cross-section data. This file was extracted by Dr. M. Jongen (see below) and is pointed to by the default constructor of `NuXsec`.
* `by_dists.root` holds neutrino energy vs bjorken-y distributions. This file was created by the scripts in `NMHDIR/bjorkeny-dists`, see `NMHDIR/bjorkeny-dists/README.md` for more info. The default constructor of `NuXsec` points to this file.


Info from Dr. M. Jongen regarding `xsec.root`
============================================

Comment from B. Strandberg: this info and the file `xsec.root` was found from [http://git.km3net.de/mjongen/NikhefOrca].

These root files contain TGraphs representing the cross-sections used in gSeaGen v4r1.

gSeaGen v4r1 is used for the 2016 9m ORCA production (see http://wiki.km3net.de/index.php/Simulations/ORCA_productions#ORCA115_9m_Layout_2016).

These files were obtained by Martijn Jongen on 12 October 2017 in the following way.

The gSeaGen installation used for the MC production can be found in Lyon:

/afs/in2p3.fr/home/throng/km3net/src/gSeaGen/v4r1

This installation comes with the following cross-section files (found in the 'dat' directory):

> ls -ltrh dat/gxspl-seawater*

12M May 23  2016 dat/gxspl-seawater-tau-5000-genie2.10.2.xml
12M May 23  2016 dat/gxspl-seawater-muon-5000-genie2.10.2.xml
12M May 23  2016 dat/gxspl-seawater-elec-5000-genie2.10.2.xml
34M May 23  2016 dat/gxspl-seawater-5000-genie2.10.2.xml
93 Jul 26  2016 dat/gxspl-seawater-tau.xml -> /afs/in2p3.fr/home/throng/km3net/src/gSeaGen/v4r1/dat/gxspl-seawater-tau-5000-genie2.10.2.xml
94 Jul 26  2016 dat/gxspl-seawater-muon.xml -> /afs/in2p3.fr/home/throng/km3net/src/gSeaGen/v4r1/dat/gxspl-seawater-muon-5000-genie2.10.2.xml
94 Jul 26  2016 dat/gxspl-seawater-elec.xml -> /afs/in2p3.fr/home/throng/km3net/src/gSeaGen/v4r1/dat/gxspl-seawater-elec-5000-genie2.10.2.xml
89 Jul 26  2016 dat/gxspl-seawater.xml -> /afs/in2p3.fr/home/throng/km3net/src/gSeaGen/v4r1/dat/gxspl-seawater-5000-genie2.10.2.xml

These files are Genie cross-section splines in xml format. For more info on these, see chapter 5 of

"The GENIE Neutrino Monte Carlo Generator: Physics and User Manual"
(arXiv:1510.05494v1 [hep-ph]).

For the extraction, the xml files are converted to root format using the GENIE utility program gspl2root (see chapter 5.2.5). The spline data written-out have the energies given in GeV and the cross sections given in 10e−38 cm^2.

The scripts with which this was done are included in this folder. They can be run in Lyon:
- first do 'source setenv_gSeaGen.sh' to set up the environment for this gSeaGen version
- then do './make_cs.sh' to make a root file called 'xsec_full.root' containing the information from the gSeaGen xml files. It contains TGraphs of the cross-section for all individual subprocesses.
- then do 'root -b -q -l extract_total_cs.C' to extract only the total CC and NC cross-section TGraphs from 'xsec_full.root'. These are saved to the final output file called 'xsec.root' (this last step is just to save disk space)