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

  double operator() (double *x, double *p);
  DetResponse* GetResponse() { return fResponse; }

protected:

  Double_t evaluate() const;

private:

  void CheckBinning(TH3 *h1, TH3 *h2);

  FitUtil     *fFitUtil;  //!< pointer to the fit utility that can be shared between several `FitPDF` instances
  DetResponse *fResponse; //!< pointer to a specific fit response that describes the data to be fitted
  std::map<TString, RooRealProxy*> fProxies; //!< map of proxies to the RooRealVar's defined in `FitUtil`

  ClassDef(FitPDF,1)
};
 
#endif
