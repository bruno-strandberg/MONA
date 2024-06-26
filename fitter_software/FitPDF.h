#ifndef FITPDF
#define FITPDF

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"

#include "FitUtil.h"
#include "AbsResponse.h"

#include "TH3.h"
#include "TRandom3.h"

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
  FitPDF(const char *name, const char *title, FitUtil *futil, AbsResponse *resp);
  FitPDF(const FitPDF& other, const char* name=0);

  // public functions

  /** Get pointer the detector response used with this instance
      \return Pointer to the response
   */
  AbsResponse* GetResponse() { return fResponse; }

  /** Get pointer to `FitUtil` used with this instance
      \return Pointer to `FitUtil`
   */
  FitUtil*     GetUtil() { return fFitUtil; }

  /**Get reference to member `fProxies`, which stores the pointers to the variables in `RooFit`
     \return Reference to member `fProxies`
   */
  proxymap_t&    GetProxyMap() { return fProxies; }
  
  /** Get a 3D histogram with expectation values for this fit model.
      The histogram which is pointed to is created on the heap and it is the user's responsibility to delete the object.
      \param name Range string as used in `RooFit`, it will apply the defined rangeName to the E/Ct/By variables. If the range does not exists, the default ranges are used.
   */
  TH3D*        GetExpValHist(const char* name=0) const {
    return fFitUtil->Expectation(fResponse, fProxies, name);
  }

  TH3D*        GetExpValErrHist(const char* name=0); 

  /** Set the seed of the random generator (`TRandom3`) that is used for generating pseudo-experiments
      \param seed  Seed for the `TRandom3` generator `fRand`
   */
  void SetSeed(ULong_t seed = 0) { fRand.SetSeed(seed); }

  TH3D*  SimplePseudoExp(TString nametitle="pseudoexp", Bool_t IncludeStatErr=kFALSE, const char* rangeName=0);
  double operator() (double *x, double *p);
  Int_t    getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const ;

  /** This virtual function is required for extended LLH fits and needs to return false for RooFit to deal with normalisation correctly 
      \return False (pdf is not self-normalized)
   */
  virtual Bool_t     selfNormalized() const { return kFALSE ; }

  /** This function is used by RooFit internally to allow extended likelihood fits.

      By default, RooFit is a shape fitter, i.e. the overall normalisation is free. However, one can specify an option `Extended(kTRUE)` when fitting, in which case RooFit adds another term to the likelihood \f$ -\lambda + N log\lambda \f$, where \f$ \lambda \f$ is the the total number of events expected in the model and \f$ N \f$ is the total number of events in the fitted histogram. I.e, the total number of events is in this case also fitted.

      If this function returns `CanBeExtended`, RooFit will automatically add the above term. That is not ideal, because by default it is more conservative to keep the normalisation free. Also, this functionality was added later, meaning that existing applications will change default behaviour. To that end, the functions `FitPDF::IncludeNorm()` and `FitPDF::ExcludeNorm` can be used to control whether the overall normalisation is included in the fit or not. 
   */
  virtual ExtendMode extendMode()     const { return fExtMode ; }

  /** Function that returns the number of expected events in the model at given parameter values */
  virtual Double_t   expectedEvents(const RooArgSet* nset) const { return analyticalIntegral(I_E_CT_BY); }

  /** Function that returns the number of expected events in the model at given parameter values */
  virtual Double_t   expectedEvents(const RooArgSet& nset) const { return expectedEvents(&nset); }

  /** Tell RooFit to include the overall normalisation term in the likelihood.
      Note that if normalisation is made possible, it can be switched or or off during fitting by doing `pdf.fitTo(data, Extended(kFALSE/kTRUE))`.
   */
  void IncludeNorm() { fExtMode = RooAbsPdf::CanBeExtended;    }

  /** Tell RooFit to exclude the overall normalisation term in the likelihood.
      Note that doing `pdf.fitTo(data, Extended(kTRUE))` will have no effect.
   */
  void ExcludeNorm() { fExtMode = RooAbsPdf::CanNotBeExtended; }

protected:

  Double_t evaluate() const;

private:

  enum integrations {I_NUMERIC = 0, I_E_CT_BY}; //!< enumerator for integration types
  
  FitUtil     *fFitUtil;    //!< pointer to the fit utility that can be shared between several `FitPDF` instances
  AbsResponse *fResponse;   //!< pointer to a specific fit response that describes the data to be fitted
  proxymap_t   fProxies;    //!< map of proxies to the RooRealVar's defined in `FitUtil`
  TRandom3     fRand;       //!< random generator for pseudoexperiments
  RooAbsPdf::ExtendMode fExtMode; //!< this flag controls the behaviour whether RooFit peforms likelihood (default) or extended likelihood fits (see `FitPDF::ExtendMode`)

  ClassDef(FitPDF,1)
};
 
#endif
