#ifndef FitUtil_h
#define FitUtil_h

#include "RooRealVar.h"
#include "RooArgSet.h"
#include "RooRealProxy.h"
#include <map>

/**
   This is a class meant to be used fit `FitPDF`.

   This class defines all of the variables (observables+parameters) and defines argument sets to distinguish which variables are observables (energy, cos-theta, bjorken-y) and which variables are parameters.

 */
class FitUtil {

 public:
  
  FitUtil();
  ~FitUtil();
  Double_t GetValue(Double_t Ereco, Double_t ctreco, Double_t a, Double_t b);
  Double_t GetValue(const std::map<TString, RooRealProxy*> &parmap);

  RooArgSet   GetSet()                    { return fParSet; }
  RooArgSet   GetObs()                    { return fObsSet; }

 private:

  RooRealVar *Ereco;
  RooRealVar *Ctreco;
  RooRealVar *a;
  RooRealVar *b;
  RooArgSet   fParSet;
  RooArgSet   fObsSet;

};

#endif
