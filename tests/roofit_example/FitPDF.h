#ifndef FITPDF
#define FITPDF

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
#include "RooArgSet.h"

#include "TH2.h"

#include "FitUtil.h"

#include <map>
 
/**
   This class is a modified version of a RooFit auto-generated class for custom PDF.

   The template was generated by `RooClassFactory::makePdf("FitPDF", "a, b, ...")`. To understand the modifications to the template, the end-goal should first be described.

   Such a pdf will be used to perform several simultaneous fits. For simplest neutrino mass hierarchy fit, there are 8 parameters (energy, cos-theta, bjorken-y, th12, th13, th23, delta-cp, m12, m13). Recalculating everything each time I move to a new energy-costheta-bjorken-y bin will be expensive. For this one could implement a cache, such that the oscillation probabilities are re-calculated only when (th12, th13, th23, delta-cp, m12 or m23) change. The cache can be shared between all simultaneous fits, which hopefully saves time.

   For these purposes, the class `FitUtil` is created. All parameters and observables are declared in `FitUtil`, the purpose of this class is to create `RooRealProxy's for all variables and call a proper fuction from `FitUtil` in the `evaluate` method.

 */


class FitPDF : public RooAbsPdf {

public:

  /** Default constructor */
  FitPDF() {} ; 

  FitPDF(const char *name, const char *title, FitUtil *futil, TH2D *h=NULL);
  FitPDF(const FitPDF& other, const char* name=0);

  /** Clone function */
  virtual TObject* clone(const char* newname) const { return new FitPDF(*this,newname); }
  inline virtual ~FitPDF() { }

  double operator() (double *x, double *p);

  Int_t    getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const ;
  
protected:

  /// pointer to the fitted histogram for integration
  TH2 *fh;
  
  /// pointer to the fit utility that can be shared between several FitPDF instances
  FitUtil *fFitUtil;
  /// map of proxies to the RooRealVar's that the PDF depends on
  std::map<TString, RooRealProxy*> fProxies;

  Double_t evaluate() const ;

private:

  ClassDef(FitPDF,1)
};
 
#endif
