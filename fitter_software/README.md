Fitter software
===============

The functionality in this directory is built on top of the tools provided by the `common_software/` classes. The main aim of the code in this directory is to provide a fitter for the NMO analysis, using all the power available through `RooFit`.

How does this work?
-------------------

### The class `FitUtil`

`FitUtil` is a class that uses the machinery available through `common_software` to essentially provide a function that can predict the number of expected neutrino events in a given energy, cos-theta and bjorken-y bin. Such a function F is always necessary to fit a dataset D that represents an experiment (either MC or from the sea). Once one has such a function, one can:

   1. Take D and F and build a likelihood function oneself, pass the LLH function to a minimizer for fitting.
   2. Wrap F in a TF3 and fit it to a TH3 that holds D using ROOT. In this case ROOT automatically creates a chi2 or LLH function and passes it to a minimizer for fitting.
   3. Wrap F into a `RooFit` pdf class and fit the pdf to `RooDataSet/Hist` that holds D using `RooFit`. In this case `RooFit` creates the LLH (or chi2) function and passes it to a minimizer for fitting.

Here, approach 3. with `RooFit` is taken. The function that returns the expected number of neutrino events in a bin is implemented in `FitUtil::RecoEvts`. Other, mainly private or protected methods of `FitUtil` are mostly related to caching to speed-up the fitting process and book-keeping. The public functions of `FitUtil` provide an interface to allow manipulation of the fit parameters and functions that are to be called from inside `FitPDF` (discussed below).

### The class `FitPDF`

`FitPDF` is a skeleton code that was auto-generated by `RooFit` and modified to accommodate the specifics of the NMO analysis. `FitPDF` is a necessary wrapper to get access to the niceties offered by `RooFit`. The idea is to keep `FitPDF` a minimal and generic wrapper around `FitUtil`. Most importantly, the function `FitPDF::evaluate()` is the function `RooFit` uses to get the number of expected events in a bin when fitting. It can be used similarly to other `RooFit` pdf's.

The trickiest part in the code is the fit parameter map `proxymap_t`. Basically, a pdf in `RooFit` needs to know of the variables (parameters+observables) it depends on through proxies (which are `RooFit` internal pointers). In this case the list of proxies is created automatically from all of the parameters defined in `FitUtil`. This is an important bonus, because it means the user can create new parameters in `FitUtil` and inheriting classes and all she/he needs to do is add it to the parameter set `FitUtil::fParSet`. For every parameter in the parameter set of `FitUtil`, `FitPDF` will automatically create a proxy. Importantly, this means that modifications to `FitPDF` are not required to add parameters. 

### Why have `FitUtil` and `FitPDF` separately?

The idea is that we typically wish to perform several fits in parallel (e.g. fit a histogram for tracks and a histogram for showers simultaneously). Practically this means one has a histogram for tracks and a pdf to fit the tracks histogram and, correspondigly, a histogram for showers and a pdf to fit the showers histogram. However, both pdf's (tracks pdf and showers pdf) need to share the fit parameters, as the fit is performed simultaneously. For this reason the fit parameters are defined centrally in the `FitUtil` class, which is shared between the multiple pdf's. Another way of thinking about it is that if one uses a simple `RooGaussian` in `RooFit`, in this case also the parameters can only be defined externally (i.e. outside) of the pdf class `RooGaussian`.

Secondly, say one wishes to fit 10 histograms in parallel, instead of 2. Calculation of the oscillation probabilities is a CPU resource drain (approx 1/3 of the time). `FitUtil` has a central cache of the oscillation probabilities, shared between all 10 fits, which allows to save some time.

All of the fit parameters are defined in `FitUtil` and `FitPDF` is mainly just a wrapper around `FitUtil` functionalities. I hope this central store makes future development somewhat easier.

Finally - this is not a perfect implementation, because it is the first ever effort to use `RooFit` for the NMO analysis inside KM3NeT. It may be that all of the `FitUtil` functionality could be part of `FitPDF` and simultaneous fits, etc, would still work. However, I am hopeful that the implementation here can serve as some guidance to how `RooFit` can be used with quite a complex fit model.

How can I use it?
-----------------

Example applications that make use of the software in this directory are available in `tests/fitting` and `tests/fitter_software`. The software can also be used in ROOT macros like any other ROOT classes; some examples can be found in the analysis macros in `macros/` directory.

I need to modify the fit model. How to do that?
-----------------------------------------------
Extension/modification of the fit model for a specific analysis can be achieved by creating a new class that inherits from `FitUtil` and overloads (one or both of) the virtual functions `FitUtil::RecoEvts` and `FitUtil::TrueEvts`. After this, the new class can be used just like the original `FitUtil` class to initialise `FitPDF`'s and fit data. This is discussed in the documentation of `FitUtil`.