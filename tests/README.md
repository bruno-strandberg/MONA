Test scripts
============

This directory holds ROOT macros and source code for applications that serve as examples, some of the applications also double as unit tests. The macros and directories are described briefly below.

## Directories

* `tests/common_software/` directory: these applications demonstrate the use of some of the classes in `NMH/common_software`. Each application in this directory is compiled by make and run as a test for continuous integration. Add applications here to increase the test suite of `NMH/common_software`, each application should return 0 at success and 1 at failure.

* `tests/fitter_software/` directory: these applications demonstrate the use of some of the classes in `NMH/fitter_software`. Each application in this directory is compiled by make and run as a test for continuous integration. Add applications here to increase the test suite of `NMH/fitter_software`, each application should return 0 at success and 1 at failure.

* `tests/fitting/` directory: these applications demonstrate the use of `NMH/common_softare` and `NMH/fitter_software` for NMO fits and analysis. These are somewhat more elaborate as the tests in `tests/fitter_software` and `tests/common_software`, rely on actual Monte-Carlo data files and are not run as part of a test suite.

* `tests/macros/` directory: hold small ROOT macros that demonstrate the usage of the classes of the base libraries in un-compiled macros.

  * `summaryparser.C` - a macro that demonstrates the usage of `SummaryParser` class for reading and writing data in `SummaryEvent` format.
  * `fileheader.C` - a macro that demonstrates the usage of `FileHeader` class for attaching simple headers to ROOT files.