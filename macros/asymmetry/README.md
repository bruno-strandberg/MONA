Asymmetry Analysis
==================

This directory contains the root scripts that have been used to do the asymmetry analysis
into using multiple PID bins in the track/shower classification scheme. The basic idea is:
Instead of using track and shower as two bins to classify events into, use multiple, such
as 0-20% track, 20-40% track, .., 80-100% track, etc.

Prerequisites
=============
Many of the scripts, such as the detector responses and asymmetry calculations require 
output folders to be available to write into, so use `$ mkdir default_detres` or 
`$ mkdir pid_detres` for the standard 2 bins or respectively multibins scripts.

The scripts
===========
The scripts use the following naming convention: first it states what it does, and then
it states in which situation it does the thing.
 - `DetecectorResponse\*` generates a detector response, for example in the default 2 bin
   case or in the N pid bins case.
 - `Asymmetry\*` calculates the asymmetry of the already generated and saved detector 
   responses. It also saves some plots.
 - `Plot\*` generates and saves some plots.
 - `Print\*` prints out some statistics on the data in the root files to cross-check that
   for example the detector response is working as expected.

Running
=======
To run the scripts one must do the following:
`$ root ScriptName.C+(a, b, c)` to execute it using root with some input values `a, b, c`.
`$ root ScriptName.C++` to recompile the script after altering it.


