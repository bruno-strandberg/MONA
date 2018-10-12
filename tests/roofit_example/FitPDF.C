#include "Riostream.h" 

#include "FitPDF.h" 
#include "RooAbsReal.h" 
#include "RooAbsCategory.h" 
#include <math.h> 
#include "TMath.h" 

ClassImp(FitPDF); 

/** Copy constructor.

    \param other   Other instance of FitPDF
    \param name    Name of the copy
 */
FitPDF::FitPDF(const FitPDF& other, const char* name) : RooAbsPdf(other, name) { 
  
  // get the pointer to the fit utility
  fFitUtil = other.fFitUtil;
  
  // re-create the proxy map
  for (auto &p: other.fProxies) {    
    RooRealProxy *po = p.second;
    TString name     = p.first;
    fProxies.insert ( std::make_pair( name, new RooRealProxy(name, this, *po) ) );
  }
        
}

//*******************************************************************************

/** Constructor.

    \param name   Name of the pdf
    \param title  Title of the pdf
    \futil        Pointer to the FitUtil class

 */
FitPDF::FitPDF(const char *name, const char *title, FitUtil *futil) : RooAbsPdf(name, title) {

  // set the pointer and create a list of the parameters for iteration
  fFitUtil = futil;
  RooArgList pars( fFitUtil->GetSet() );

  // loop over parameters and create proxies
  for (Int_t i = 0; i < pars.getSize(); i++) {
    RooRealVar* v = (RooRealVar*)pars.at(i);
    TString name = v->GetName();
    fProxies.insert( std::make_pair( name, new RooRealProxy(name, name, this, *v) ) );
  }
  
}

//*******************************************************************************

/** This method is called by the minimiser to get the expected number of events in a bin.
 
    In this approach I just call a method of the fit utility.
 */
Double_t FitPDF::evaluate() const { 
  
  return fFitUtil->GetValue(fProxies);

} 

//*******************************************************************************

/** Implemented so that I can also use this class with ROOT TF1/2/3 fitting*/
double FitPDF::operator()(double *x, double *p) {

  Double_t E  = x[0];
  Double_t ct = x[1];
  Double_t a  = p[0];
  Double_t b  = p[1];

  return fFitUtil->GetValue(E, ct, a, b);
  
}
