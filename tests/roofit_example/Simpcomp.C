#include "TH2.h"
#include "TFile.h"
#include "TF2.h"
#include "TRandom3.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TMath.h"

#include "RooDataHist.h"
#include "RooFitResult.h"
#include "RooFormulaVar.h"
#include "RooGenericPdf.h"

#include <iostream>

using namespace RooFit;
using namespace std;

int main() {

  //-----------------------------------------------------
  // set up the roofit part
  //-----------------------------------------------------
  RooRealVar x("x","x",-1e3,1e3);
  RooRealVar y("y","y",-1e3,1e3);
  RooRealVar a("a","a",-1e3,1e3);
  RooRealVar b("b","b",-1e3,1e3);
  RooGenericPdf roofitpdf ("roofitpdf" , "(a+b)*x+(a-b)*y", RooArgList(x,y,a,b));
  RooFormulaVar roofitfunc("roofitfunc", "(a+b)*x+(a-b)*y", RooArgList(x,y,a,b));

  //-----------------------------------------------------
  // set up the roof fit function
  //-----------------------------------------------------
  TF2 rootfunc("rootfunc", "([0]+[1])*x+([0]-[1])*y", -100,100,-100,100);
  rootfunc.SetParNames("a","b");
  
  //-----------------------------------------------------
  // create data
  //-----------------------------------------------------

  //TH2D h1("h1","h1", 16, 0.5, 16.5, 16, 0.5, 16.5); // <== with this roofitfunc result matches ROOT result
  TH2D h1("h1","h1", 32, 0.5, 16.5, 16, 0.5, 16.5); // <== with this roofitfunc result overestimates x2
  h1.Sumw2();
  Double_t apar = 8.2;
  Double_t bpar = 4.3;
  
  for (Int_t xbin = 1; xbin <= h1.GetXaxis()->GetNbins(); xbin++) {
    for (Int_t ybin = 1; ybin <= h1.GetYaxis()->GetNbins(); ybin++) {

      Double_t E  = h1.GetXaxis()->GetBinCenter(xbin);
      Double_t ct = h1.GetYaxis()->GetBinCenter(ybin);

      Double_t bc1 = (apar+bpar)*E + (apar-bpar)*ct;

      h1.SetBinContent(xbin, ybin, bc1);
      h1.SetBinError(xbin, ybin, TMath::Sqrt(bc1));
    }
  }

  //-----------------------------------------------------
  // set guess values to parameters, both root and roofit
  //-----------------------------------------------------

  TRandom3 fRand(0);
  Double_t e = fRand.Uniform(-0.3, 0.3); // to add some error to the fitted pars

  rootfunc.SetParameters(apar + e*apar, bpar + e*bpar);
  a.setVal(apar + e*apar);
  b.setVal(bpar + e*bpar);
  
  //-----------------------------------------------------
  // import data to roofit and fit with both RooFormulaVar and GenericPdf
  //-----------------------------------------------------
  RooDataHist rf_h1("rf_h1", "rf_h1", RooArgList(x,y), Import(h1) );
  rf_h1.Print("v");
  
  RooFitResult *roofit_result1 = roofitfunc.chi2FitTo(rf_h1, Save(kTRUE), Verbose(kTRUE));
  RooArgSet final1 (roofit_result1->floatParsFinal() );

  cout << "******************************************************************" << endl;
  cout << "Finished fitting with RooFit RooRealVar" << endl;
  cout << "******************************************************************" << endl;

  // reset guess values
  a.setVal(apar + e*apar);
  b.setVal(bpar + e*bpar);

  RooFitResult *roofit_result2 = roofitpdf.chi2FitTo(rf_h1, Save(kTRUE), Verbose(kTRUE), PrintLevel(1000));
  RooArgSet final2 (roofit_result2->floatParsFinal() );
  
  cout << "******************************************************************" << endl;
  cout << "Finished fitting with RooFit RooGenericPdf" << endl;
  cout << "******************************************************************" << endl;

  //-----------------------------------------------------
  // do the fit in ROOT
  //-----------------------------------------------------
  TFitResultPtr root_result = h1.Fit("rootfunc", "VNS");

  cout << "******************************************************************" << endl;
  cout << "Finished fitting with ROOT" << endl;
  cout << "******************************************************************" << endl;
  
  //-----------------------------------------------------
  // print comparison
  //-----------------------------------------------------

  cout << "******************************************************************" << endl;
  cout << "ROOFIT RooRealVar Parameter a initial value, fitted value: "
       << apar << "\t" << ((RooRealVar*)final1.find("a"))->getVal() << endl;
  cout << "ROOFIT RooRealVar Parameter b initial value, fitted value: " 
       << bpar << "\t" << ((RooRealVar*)final1.find("b"))->getVal() << endl;
  cout << "******************************************************************" << endl;

  cout << "******************************************************************" << endl;
  cout << "ROOFIT RooGenericPdf Parameter a initial value, fitted value: "
       << apar << "\t" << ((RooRealVar*)final2.find("a"))->getVal() << endl;
  cout << "ROOFIT RooGenericPdf Parameter b initial value, fitted value: " 
       << bpar << "\t" << ((RooRealVar*)final2.find("b"))->getVal() << endl;
  cout << "******************************************************************" << endl;

  cout << "******************************************************************" << endl;
  cout << "ROOT Parameter a initial value, fitted value: "
       << apar << "\t" << root_result->Parameter(0) << endl;
  cout << "ROOT Parameter b initial value, fitted value: " 
       << bpar << "\t" << root_result->Parameter(1) << endl;
  cout << "******************************************************************" << endl;

  TFile fout("fitted-hist.root", "RECREATE");
  h1.Write();
  fout.Close();

}

void Simpcomp() {
  main();
}
