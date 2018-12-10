Bjorken-y distributions
==============

Scripts in this directory can be used to create histograms that hold normalised bjorken-y distributions at given energies. The distributions are created from gSeaGen data. The histograms have neutrino energy on the x-axis and bjorken-y values on the y-axis. The normalisation is such that at a fixed energy \f[\sum_{bybins} = 1\f].

The output of these scripts is used by the class `NMHDIR/common_software/NuXsec`. The default constructor of `NuXsec` looks for a file in `NMHDIR/data/cross_sections_gSeaGen_v4r1/by_dists.root`. If files from a new gSeaGen version are expected to yield significantly different bjorken-y distributions, the scripts in this directory should be run and the `by_dists.root` file replaced. Alternatively, the default constructor of `NuXsec` should be changed to point to the new file.

Prerequisities
==============
* The scripts use the code in NMH/common_software/.
* To parse gSeaGen files in .evt format, common_software needs Jpp and aanet root6 branch.

How to run
==========

* First run `Create_BY_hists.C+(...)` from `ROOT` prompt (see next point). It takes as an input a list of gSeaGen files (create via e.g. `ls /path/to/gsg_files > list.dat`) and an output name. This application loops over gSeaGen events and fills the data to energy vs bjorken-y histograms. See also the documentation int the macro.

* Practically, running `Create_BY_hists.C` is better performed on the computing farm at Lyon. the python script `hist_caller.py` helps to do that. For example, run it as `./hist_caller.py -g 'my/dir/*gSeaGen*'`. See also documentation in the script.

* The output of `Create_BY_hists.C` needs to be normalised - this task is performed by `Normalised_BY_dists.C`. It takes as input a list of `Create_BY_hists.C` outputs and an output name. In practice, `hist_caller.py` creates an output per each farm job. `Create_BY_hists.C` can be run as
```
ls output/output_job* > output_list.dat
root
.x Normalised_BY_dists.C('output_list.dat', 'by_dists.root')
```
The file `by_dists.root` is the distribution file that should be copied to `NMHDIR/data/cross_sections_gSeaGen_v4r1` with some user-defined name and the `NuXsec` constructor should be updated to point to the new file.

Outputs
==========

* `Create_BY_hists.C` outputs a `root` file with events filled to {flavor}_{nu/nubar}_{nc/cc} histograms with energy on the X-axis and bjoren-y on the Y-axis.
* `Normalised_BY_dists.C` outputs a `root` file with {flavor}_{nu/nubar}_{nc/cc} histograms that represent energy vs bjorken-y distributions. The distributions have been normalised to 1 at each energy.