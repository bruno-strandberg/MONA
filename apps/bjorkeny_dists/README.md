Bjorken-y distributions
==============

Applications in this directory can be used to create histograms that hold normalised bjorken-y distributions at given energies. The distributions are created from gSeaGen data. The histograms have neutrino energy on the x-axis and bjorken-y values on the y-axis. The normalisation is such that at a fixed energy \f[\sum_{bybins} = 1\f].

The output of these scripts is used by the class `common_software/NuXsec`. The default constructor of `NuXsec` looks for a file in `data/cross_sections_gSeaGen_v4r1/by_dists.root`. If files from a new gSeaGen version are expected to yield significantly different bjorken-y distributions, the scripts in this directory should be run and the `by_dists.root` file replaced. Alternatively, the default constructor of `NuXsec` should be changed to point to the new file.

How to run
==========

* First run `Create_BY_hists` (do -h! for help). It takes as an input a list of gSeaGen files (create via e.g. `ls /path/to/gsg_files > list.dat`), an output name, the MC energy range and the binning configuration. This application loops over gSeaGen events and fills the data to energy vs bjorken-y histograms. See also the documentation int the app.

* Practically, running `Create_BY_hists` is better performed on the computing farm at Lyon. The python script `hist_caller.py` helps to do that. For example, run it as `./hist_caller.py -g 'my/dir/*gSeaGen*' --farm`. See also documentation in the script, use -h for more info.

* The output of `Create_BY_hists` needs to be normalised - this task is performed by `Normalised_BY_dists`. It takes as input a list of `Create_BY_hists` outputs and an output name. In practice, `hist_caller.py` creates an output per each farm job. `Create_BY_hists` can be run as
~~~
ls output/output_job* > output_list.dat
./Normalised_BY_dists -f output_list.dat -o by_dists.root
~~~
The file `by_dists.root` is the distribution file that should be copied to `data/cross_sections_gSeaGen_v4r1` with some user-defined name and the `NuXsec` constructor should be updated to point to the new file.

Outputs
==========

* `Create_BY_hists` outputs a `root` file with events filled to {flavor}_{nu/nubar}_{nc/cc} histograms with energy on the X-axis and bjoren-y on the Y-axis.
* `Normalised_BY_dists` outputs a `root` file with {flavor}_{nu/nubar}_{nc/cc} histograms that represent energy vs bjorken-y distributions. The distributions have been normalised to 1 at each energy.