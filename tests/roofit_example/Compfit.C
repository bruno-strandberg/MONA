#include "FitUtil.h"
#include "FitPDF.h"

#include "TH2.h"
#include "TFile.h"
#include "TF2.h"
#include "TRandom3.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"

#include "RooDataHist.h"
#include "RooFitResult.h"
#include "RooGenericPdf.h"

#include <iostream>

using namespace RooFit;
using namespace std;

/**
   This program performs a fit to a simple 2D histogram with ROOT and with RooFit.

   It utilises the classes `FitUtil` and `FitPDF` to define a probability density function for RooFit and a TF2 for fitting with ROOT.
 */
int main() {
  
  //-----------------------------------------------------
  // create data
  //-----------------------------------------------------

  TH2D h1("h1","h1", 16, 0.5, 2.5, 16, 0.5, 2.5);

  FitUtil futil;

  Double_t apar = 8.2;
  Double_t bpar = 4.3;

  for (Int_t xbin = 1; xbin <= h1.GetXaxis()->GetNbins(); xbin++) {
    for (Int_t ybin = 1; ybin <= h1.GetYaxis()->GetNbins(); ybin++) {

      Double_t E  = h1.GetXaxis()->GetBinCenter(xbin);
      Double_t ct = h1.GetYaxis()->GetBinCenter(ybin);

      Double_t bc1 = futil.GetValue(E, ct, apar, bpar);

      h1.SetBinContent(xbin, ybin, bc1);

    }
  }

  //-----------------------------------------------------
  // set up the fitting with roofit
  //-----------------------------------------------------

  TRandom3 fRand(0);
  Double_t e = fRand.Uniform(-0.3, 0.3); // to add some error to the fitted pars

  // set guess values for the parameters
  ( (RooRealVar*)futil.GetSet().find("a") )->setVal(apar + e*apar);
  ( (RooRealVar*)futil.GetSet().find("b") )->setVal(bpar + e*bpar);

  // init the pdf
  FitPDF fpdf("mypdf","mypdf", &futil, &h1);

  // import data to RooFit - at this point RooFit will infer which variables defined in FitUtil
  // are observables and which are parameters
  RooDataHist rf_h1("rf_h1", "rf_h1", futil.GetObs(), Import(h1) );

  // perform fit in RooFit and save the result for printing
  RooFitResult *roofit_result = fpdf.chi2FitTo(rf_h1, Save(kTRUE));
  RooArgSet final (roofit_result->floatParsFinal() );

  cout << "******************************************************************" << endl;
  cout << "Finished fitting with RooFit method 1" << endl;
  cout << "******************************************************************" << endl;

  //-----------------------------------------------------
  // set up the fitting with roofit using the basics for comparison
  //-----------------------------------------------------
  RooGenericPdf genPDF("genPDF", "(a+b)*Ereco+a*b*Ctreco", futil.GetSet());
  RooFitResult *roofit_result2 = genPDF.chi2FitTo(rf_h1, Save(kTRUE));
  RooArgSet final2 (roofit_result2->floatParsFinal() );

  cout << "******************************************************************" << endl;
  cout << "Finished fitting with RooFit method 2" << endl;
  cout << "******************************************************************" << endl;
    
  //-----------------------------------------------------
  // do the fit in ROOT
  //-----------------------------------------------------

  TF2 *func = new TF2("fitfunc", fpdf, 0, 10, 0, 10, 2);
  func->SetParameters(apar + e*apar, bpar + e*bpar);
  func->SetParNames("a","b");

  TFitResultPtr root_result = h1.Fit("fitfunc", "VNS");

  //-----------------------------------------------------
  // print comparison
  //-----------------------------------------------------

  cout << "******************************************************************" << endl;
  cout << "ROOFIT FitPDF Parameter a initial value, fitted value: "
       << apar << "\t" << ((RooRealVar*)final.find("a"))->getVal() << endl;
  cout << "ROOFIT FitPDF Parameter b initial value, fitted value: " 
       << bpar << "\t" << ((RooRealVar*)final.find("b"))->getVal() << endl;
  cout << "******************************************************************" << endl;

  cout << "******************************************************************" << endl;
  cout << "ROOFIT genPDF Parameter a initial value, fitted value: "
       << apar << "\t" << ((RooRealVar*)final2.find("a"))->getVal() << endl;
  cout << "ROOFIT genPDF Parameter b initial value, fitted value: " 
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
