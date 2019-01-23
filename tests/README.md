#Test scripts

This directory holds ROOT macros and source code for applications that serve as examples, some of the applications also double as unit tests. The macros and directories are described briefly below.

## Directories

* `event_filtering` directory: these applications demonstrate the use of `EventFilter`, `EventSelection` and `DetResponse`.

* `fitting` directory: these applications demonstrate the use of `common_softare` and `fitter_software` for NMO fits and analysis.

## Macros

* `summaryparser.C` - a macro that demonstrates the usage of `SummaryParser` class for reading and writing data in `SummaryEvent` format.
* `fileheader.C` - a macro that demonstrates the usage of `FileHeader` class for attaching simple headers to ROOT files.

# Suitable for unit tests

* All applications in `event_filtering`.