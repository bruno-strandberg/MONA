Common software
===============

Classes in this directory are compiled into a shared library object, which is automatically loaded by `ROOT` and provide tools for neutrino mass ordering analysis. Each class is carefully documented in the header and source file and can be explored in a web browser, once the doxygen documentation has been created (see `NMH/README.md`). Very brief descriptions of the classes is provided below.

Data format and data filtering
------------------------------
* `SummaryEvent`    - a class that defines the data format for this analysis framework
* `EventFilter`     - a class that defines an event filter that acts on `SummaryEvent`s; parent class of `EventSelection` and `DetResponse`
* `EventSelection`  - given a `TTree` of data in `SummaryEvent` format (either MC data or sea data), this class allows to filter events to a certain event class (e.g. filter events into track-like)
* `DetResponse`     - given a `TTree` of **monte-carlo** data in `SummaryEvent` format, this class allows to create a detector response for a certain event class (e.g. for track-like events)

These classes are the main work-horses for data manipulation. Given some *experiment data* (either simulated MC or sea data) and the corresponding Monte-Carlo data (both in `SummaryEvent` format), these classes provide tools to create selections of the data and corresponding detector responses. Given that the user has a model that predicts how many events are expected in a certain true (energy, cos-theta, bjorken-y) bin, the user can use the detector response to get the predicted number of events in a reco (energy, cos-theta, bjorken-y) bin. The model can be anything (e.g. normal ordering, inverted ordering, normal ordering with sterile) and it is up to the user to define. The point is that the `DetResponse` class can help to map events from true -> reco, and that the detector response can be set up easily for any event selection imposed on a `SummaryEvent`. The applications in `NMH/apps/fitter` demonstrate some use cases.

Calculator classes
------------------
* `AtmFlux`         - a class to get the atmospheric neutrino flux
* `NuXsec`          - a class to get the neutrino interaction cross-section in water
* `EffMass`         - a class to read/write effective mass histograms in a specific format and perform effective mass calculations

These classes act as calculators to get easy access to atmospheric flux (based on a honda table) or cross section (based on a ROOT file with xsec data from GENIE) predictions. The `EffMass` class is slightly more complicated, because effective mass needs to be calculated from the simulated data. Applications in `apps/effective_mass` can be used to create the required effective mass data file and the class `EffMass` can be subsequently used for effective mass calculations with the created data file as input.

Utility classes
---------------
* `SummaryParser`   - a class to read or write files in *analysis format* (`SummaryEvent.h/C`)
* `FileHeader`      - a class to create simple headers for `ROOT` files
* `NMHUtils`        - a misc collection of useful functions, e.g. for asymmetry calculation and for calculating logarithmic bins

These classes do not play an important role in terms of analysis concepts, but they do provide some utilities that make life substantially easier.

How to use?
===========

* The classes can be used readily from the ROOT prompt, e.g:
```
root
AtmFlux f
f.Flux_dE_dcosz(1, 0, 15, -0.8) 
```

* The classes can be used in un-compiled or compiled ROOT macros exactly like any other ROOT class.

* The classes can be used in compiled applications. In this case some experience with `Makefile` is required. For examples, consult the makefiles in `NMH/apps/../`.


Maintenance & outlook to the future
===================================

Changes to `SummaryEvent` class
-------------------------------
It is very likely that in the future, someone will discover that the variables in `SummaryEvent` are not enough and more/different variables are required. From the view-point of the library, this can be addressed trivially by changing the member variables of the `SummaryEvent` class and providing the corresponding getter and setter functions. After this, both the `EventSelection` and `DetResponse` can use the new/changed variable and the developments are finished. 

The applications (`NMH/apps/../`) that are based on the library are affected if some **existing** variables or interfaces (setters/getters) are changed. Addition of new variables will not affect them.

As the start of development of this package started in the early days of ORCA when we had 1 line in water, it was very difficult to forecast which variables will be required in the long-term. Thus, changes to `SummaryEvent` are expected and even changes that break all the applications should be considered lightly. The applications are easy to modify or delete and should be considered temporary anyway. The main value of the package is in this library and the `NMH/fitter_software` library, which do not make assumptions about the variable content of `SummaryEvent`.

MC chain data/sea data in native format
---------------------------------------
Sometime in the (far?) future we will agree that all data (sea data and monte-carlo summary data) will be in XXX format (when this package was started, no such agreement existed - some said `aanet` is the best, some said `hdf5` is the best, some said plain text is the best). Let it be assumed that, for example, `aanet` emerges victorious. How can one now use this package, which assumes `SummaryEvent` format? Should everything be deleted and started from scratch? Or should I convert each file I wish to use to `SummaryEvent` format? This sounds painful...

**Easy and sustainable solution**: TBD set up pseudo-code in `SummaryParser` to exemplify how `aanet` files could be read.