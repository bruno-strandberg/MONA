#include "Riostream.h" 

#include "FitPDF.h" 
#include "RooAbsReal.h" 
#include "RooAbsCategory.h" 
#include <math.h> 
#include "TMath.h" 

#include <iostream>
using namespace std;

ClassImp(FitPDF); 

/** Copy constructor.

    \param other   Other instance of FitPDF
    \param name    Name of the copy
 */
FitPDF::FitPDF(const FitPDF& other, const char* name) : RooAbsPdf(other, name) { 
  
  // get the pointer to the fit utility
  fFitUtil = other.fFitUtil;
  fh = other.fh;
  
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
FitPDF::FitPDF(const char *name, const char *title, FitUtil *futil, TH2D* h) : RooAbsPdf(name, title) {

  // set the pointer and create a list of the parameters for iteration
  fFitUtil = futil;
  fh = h;
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

//*******************************************************************************

Int_t FitPDF::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const { 
  // LIST HERE OVER WHICH VARIABLES ANALYTICAL INTEGRATION IS SUPPORTED, 
  // ASSIGN A NUMERIC CODE FOR EACH SUPPORTED (SET OF) PARAMETERS 
  // THE EXAMPLE BELOW ASSIGNS CODE 1 TO INTEGRATION OVER VARIABLE X
  // YOU CAN ALSO IMPLEMENT MORE THAN ONE ANALYTICAL INTEGRAL BY REPEATING THE matchArgs 
  // EXPRESSION MULTIPLE TIMES

  if (fh == NULL) return 0;
  
  //if ( matchArgs( allVars, analVars, *( fProxies.at("Ereco") ), *( fProxies.at("Ctreco") ) ) ) { return 1; }
  if ( matchArgs( allVars, analVars, fFitUtil->GetObs() ) ) { return 1; }
  else if ( matchArgs( allVars, analVars, *(fProxies.at("Ereco")) ) ) { return 2; }
  else if ( matchArgs( allVars, analVars, *(fProxies.at("Ctreco")) ) ) { return 3; }
  else return 0;
  
  // if (matchArgs(allVars,analVars,x)) return 1 ; 
  //return 0 ; 
} 

//*******************************************************************************

Double_t FitPDF::analyticalIntegral(Int_t code, const char* rangeName) const { 
  // RETURN ANALYTICAL INTEGRAL DEFINED BY RETURN CODE ASSIGNED BY getAnalyticalIntegral
  // THE MEMBER FUNCTION x.min(rangeName) AND x.max(rangeName) WILL RETURN THE INTEGRATION
  // BOUNDARIES FOR EACH OBSERVABLE x
  
  if (code == 1) {
    cout << "Executing my analytical integral over E,ct" << endl;
    return fh->Integral();
  }
  else if (code == 2) {
    cout << "Executing my analytical integral over E" << endl;
    Int_t ctbin = fh->GetYaxis()->FindBin( *( fProxies.at("Ctreco") ) );
    Double_t integ = 0.;
    for (Int_t ebin = 1; ebin <= fh->GetXaxis()->GetNbins(); ebin++) {
      integ += fh->GetBinContent(ebin, ctbin);
    }
    return integ;
  }
  else if (code == 3) {
    cout << "Executing my analytical integral over Ct" << endl;
    Int_t ebin = fh->GetXaxis()->FindBin( *( fProxies.at("Ereco") ) );
    Double_t integ = 0.;
    for (Int_t ctbin = 1; ctbin <= fh->GetYaxis()->GetNbins(); ctbin++) {
      integ += fh->GetBinContent(ebin, ctbin);
    }
    return integ;
  }
  else {
    cout << "Guessing I am about to integrate numerically" << endl;
    return 0;
  }
  // assert(code==1) ; 
  // return (x.max(rangeName)-x.min(rangeName)) ; 
  //return 0 ; 
} 
