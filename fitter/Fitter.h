#ifndef Fitter_h
#define Fitter_h

#include "DetResponse.h"
#include "EventSelection.h"
#include "FitPDF.h"

#include <vector>
#include <tuple>

/** A namespace that collects functions and variables for the fitting */
namespace Fitter {

  // structure with parameters expected from the command line
  struct cmdpars {
    TString simdata_file;
    TString expdata_file;
    Bool_t  refill_response;
  };

  // functions
  cmdpars CommandLineParser(int argc, char **argv);
  void InitRespsAndSels();
  void FillRespsAndSels(TString simdata_file, TString expdata_file, Bool_t refill_response);
  void FillSelections();
  void FillExpectationValues(FitPDF *track_pdf, FitPDF *shower_pdf);
  void Cleanup();

  // binning variables
  static const Int_t    fENB    =  40;
  static const Double_t fEmin   =   1;
  static const Double_t fEmax   = 100;
  static const Int_t    fCtNB   =  40;
  static const Double_t fCtmin  =  -1;
  static const Double_t fCtmax  =   1;
  static const Int_t    fByNB   =   1;
  static const Double_t fBymin  =   0;
  static const Double_t fBymax  =   1;

  // member selections and responses
  DetResponse    *fTRres;
  DetResponse    *fSHres;
  EventSelection *fTRsel;
  EventSelection *fSHsel;

};

#endif
