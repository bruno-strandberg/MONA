Common software
===============

This directory has some generic classes that are used by other scripts in the NMH directory. The following classes can be used stand-alone:

* SummaryParser - to parse summary files the analysis format as described below,
* GSGParser     - to parse gSeaGen files, either in .root or .evt format,
* CScalculator  - a class to get the neutrino interaction cross-section in water,
* AtmFlux       - a class to get the atmospheric neutrino flux.

Prerequisities
==============
* To compile, the script ../setenv.sh needs to be sourced before.
* Jpp with aanet, using root6, need to be configured.

How to run
==============

The classes can be used in root by doing:
1. > root
2. > .L libnmhsoft.so
3. > SummaryParser sp("/path/to/summary_file/")

The classes can also be used in compiled root macros, compiled programs etc., just like any other ROOT class.

The usage of data parsers is illustrated in NMH/examples/data_parsers.C. If there is gSeaGen data and summary data in mc_end and mc_start, one can, in the examples directory, do '''root data_parsers.C+ '''.

The analysis data format
========================

### Variables in tree
* pos and dir refer to vertex position and neutrino direction, given for MC truth and 3 reco's.
* energy_nu and bjorkeny refer to the neutrion energy and Bjorken y, given for MC truth a 3 reco's (expect no bjorkeny for gandalf).
* Other variables are discussed in SummaryParser.h
* The variables <reco>_ql0(12) are described in detail below.

### Quality levels

For each reconstruction there are branches called '<reco>_ql0, <reco>_ql1, ...' .  'ql' stands for 'quality level'.

* Quality level 0 is the lowest quality level, meaning somewhere in the reco chain it has been tested that the reconstruction worked in the most minimal way, e.g. that it did not return a NaN.

* Quality level 1 is the next lowest quality level. In this case the reconstruction has been tested against quality cuts set by the reconstruction specialists, as described on the wiki and the Moritz' script, linked in NMH/data/README.md. Additionally, in this case the event containment cut is loosened, such that the vertex can be up to 30 m outside the instrumented volume.

* Quality level 2 is the highest quality level. It is the same as above, but the vertex needs to be inside the instrumented volume.

We are usually (March 2018) advised to use ql1, as it selects more events and improves
statistics. For dusj there is no level 2 defined.