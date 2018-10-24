#include "Riostream.h" 

#include "FitPDF.h" 
#include "RooAbsReal.h" 
#include "RooAbsCategory.h" 
#include <math.h> 
#include "TMath.h" 

#include <stdexcept>

ClassImp(FitPDF); 

/** Copy constructor.

    \param other   Other instance of FitPDF
    \param name    Name of the copy
 */
FitPDF::FitPDF(const FitPDF& other, const char* name) : RooAbsPdf(other, name) { 
  
  // get the pointer to the fit utility and the response
  fFitUtil  = other.fFitUtil;
  fResponse = other.fResponse;

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
    \futil        Pointer to the `FitUtil` class
    \resp         Pointer to a `DetResponse`

 */
FitPDF::FitPDF(const char *name, const char *title, FitUtil *futil, DetResponse *resp) : RooAbsPdf(name, title) {

  // set the pointers
  fFitUtil  = futil;
  fResponse = resp;

  CheckBinning( fFitUtil->GetBinningHist(), fResponse->GetHist3D() );

  // create a list of the parameters for iteration and create proxies
  RooArgList pars( fFitUtil->GetSet() );

  for (Int_t i = 0; i < pars.getSize(); i++) {
    RooRealVar* v = (RooRealVar*)pars.at(i);
    TString name = v->GetName();
    fProxies.insert( std::make_pair( name, new RooRealProxy(name, name, this, *v) ) );
  }
  
}

//*******************************************************************************

/** This method is called by the minimiser to get the expected number of events in a bin */
Double_t FitPDF::evaluate() const { 
  
  return fFitUtil->PdfEvaluate(fProxies, fResponse);

} 

//*******************************************************************************

/** This funcion helps RooFit to decide which integration technique to use*/
Int_t FitPDF::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const {

  if ( matchArgs( allVars, analVars, fFitUtil->GetObs() ) ) { return I_E_CT_BY; }
  else                                                      { return (Int_t)I_NUMERIC; }

}

//*******************************************************************************

/** depending on the code returned by analyticalIntegral return the integral*/
Double_t FitPDF::analyticalIntegral(Int_t code, const char* rangeName) const {

  Double_t integral = 0.;
  
  if ( code == I_E_CT_BY ) {
    TH3D *hexp = fFitUtil->PdfGetExpValHist(fProxies, fResponse, rangeName);
    integral = hexp->Integral();
    delete hexp;
  }
  else {
    integral = 0.;
  }

  return integral;
  
}

//*******************************************************************************

/** Implementation of this function allows the class to be used with TF1/TF2/TF3

    \param x  Observable array (reconstructed E, cos-theta, bjorken-y)
    \param p  Parameter array (sinsqth12, 13, 23, delta-cp, dm21, dm31)
    \return   Number of events in the reco bin
    
 */
double FitPDF::operator()(double *x, double *p) {

  Double_t E  = x[0];
  Double_t ct = x[1];
  Double_t by = x[2];

  Double_t sinsqth12  = p[0];
  Double_t sinsqth13  = p[1];
  Double_t sinsqth23  = p[2];
  Double_t dcp        = p[3];
  Double_t dm21       = p[4];
  Double_t dm31       = p[5];

  return fFitUtil->RecoEvts(fResponse, E, ct, by, sinsqth12, sinsqth13, sinsqth23, dcp, dm21, dm31).first;
  
}

//*******************************************************************************

void FitPDF::CheckBinning(TH3 *h1, TH3 *h2) {

  TString xtitle =  h1->GetXaxis()->GetTitle();
  TString ytitle =  h1->GetYaxis()->GetTitle();
  TString ztitle =  h1->GetZaxis()->GetTitle();

  h1->GetXaxis()->SetTitle("x");
  h1->GetYaxis()->SetTitle("y");
  h1->GetZaxis()->SetTitle("z");

  std::vector< std::pair<TAxis*, TAxis*> > axps = { std::make_pair( h1->GetXaxis(), h2->GetXaxis() ), 
						    std::make_pair( h1->GetYaxis(), h2->GetYaxis() ),
						    std::make_pair( h1->GetZaxis(), h2->GetZaxis() ) };
  
  for (auto &axp: axps) {

    auto ax1 = axp.first;
    auto ax2 = axp.second;

    if ( ax1->GetNbins() != ax2->GetNbins() ) {
      throw std::invalid_argument("ERROR! FitPDF::CheckBinning() different number of bins on axis " + 
				  (string)ax1->GetTitle() );
    }

    Int_t nbins = ax1->GetNbins();

    for (Int_t bin = 1; bin <= nbins; bin++) {

      if ( ax1->GetBinLowEdge(bin) != ax2->GetBinLowEdge(bin) ) {
	throw std::invalid_argument("ERROR! FitPDF::CheckBinning() different bin low edges on axis " + 
				    (string)ax1->GetTitle() + ", bin number " + to_string(bin) );
      }

    }

  }

  h1->GetXaxis()->SetTitle(xtitle);
  h1->GetYaxis()->SetTitle(ytitle);
  h1->GetZaxis()->SetTitle(ztitle);

}
