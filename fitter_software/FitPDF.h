#ifndef FITPDF
#define FITPDF

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"

#include "FitUtil.h"
#include "DetResponse.h"

#include "TH3.h"

#include <map>

/** This class is a modification of a code autogenerated with `RooClassFactory::makePdf` and provides access to fitting tools with `RooFit`.
    
    This class is meant to fit `EventSelection`'s that represent data from ORCA detector. The data is basically a 3D histogram with reconstructed energy on the x-axis, cos-theta on the y-axis and bjorken-y on the z-axis. In principle a fit to unbinned data could also be performed, but as the `DetResponse` is inherently binned, the validity of such an approach is questionable.

    For each `EventSelection`, a `DetResponse` can be defined (see the documentation of these classes) and a pointer to a `DetResponse` that corresponds to an `EventSelection` has to be provided to this class. Secondly, the class requires a pointer to a `FitUtil` class, where all of the calculations are performed (see documentation of `FitUtil`). Ideally the `FitPDF` class does not need to be modified, but acts as a wrapper around `FitUtil`.

    The class has also methods to get expectation value histograms (as used for asymmetry analysis) and to create simplistic pseudo-experiments (expectation value histograms that have been Poisson-smeared). These histograms can also be fitted with the PDF for toy studies.

*/
class FitPDF : public RooAbsPdf {

public:

  /** Default constructor */
  FitPDF() {} ;
  
  /** Destructor */
  inline virtual ~FitPDF() { }
  
  /** Clone function */
  virtual TObject* clone(const char* newname) const { return new FitPDF(*this,newname); }
  FitPDF(const char *name, const char *title, FitUtil *futil, DetResponse *resp);
  FitPDF(const FitPDF& other, const char* name=0);

  // public functions

  /** Get pointer the `DetResponse` used with this instance
      \return Pointer to `DetResponse`
   */
  DetResponse* GetResponse() { return fResponse; }
  /** Get pointer to `FitUtil` used with this instance
      \return Pointer to `FitUtil`
   */
  FitUtil*     GetUtil() { return fFitUtil; }
  /** Get a 3D histogram with expectation values for this fit model
      \param name Name of a range as used in `RooFit`, dummy for now
   */
  TH3D*        GetExpValHist(const char* name=0) {
    return fFitUtil->PdfExpectation(fProxies, fResponse, name);
  }


  TH3D*  SimplePseudoExp();
  double operator() (double *x, double *p);
  Int_t    getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const ;
  
protected:

  Double_t evaluate() const;

private:

  enum integrations {I_NUMERIC = 0, I_E_CT_BY}; //!< enumerator for integration types
  
  FitUtil     *fFitUtil;                     //!< pointer to the fit utility that can be shared between several `FitPDF` instances
  DetResponse *fResponse;                    //!< pointer to a specific fit response that describes the data to be fitted
  std::map<TString, RooRealProxy*> fProxies; //!< map of proxies to the RooRealVar's defined in `FitUtil`

  ClassDef(FitPDF,1)
};
 
#endif
