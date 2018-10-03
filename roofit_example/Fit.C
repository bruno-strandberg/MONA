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
#include "RooCategory.h"
#include "RooSimultaneous.h"

#include <iostream>

using namespace RooFit;
using namespace std;

int main() {
  
  //-----------------------------------------------------
  // create data
  //-----------------------------------------------------

  TH2D h1("h1","h1", 16, 0.5, 2.5, 16, 0.5, 2.5);
  TH2D h2("h2","h2", 16, 0.5, 2.5, 16, 0.5, 2.5);

  FitUtil futil;

  Double_t apar1 = 8.2;
  Double_t bpar1 = 4.3;
  Double_t apar2 = apar1-1; // apar should then be between two values in case of a simultaneous fit
  Double_t bpar2 = bpar1;

  for (Int_t xbin = 1; xbin <= h1.GetXaxis()->GetNbins(); xbin++) {
    for (Int_t ybin = 1; ybin <= h1.GetYaxis()->GetNbins(); ybin++) {

      Double_t E  = h1.GetXaxis()->GetBinCenter(xbin);
      Double_t ct = h1.GetYaxis()->GetBinCenter(ybin);

      Double_t bc1 = futil.GetValue(E, ct, apar1, bpar1);
      Double_t bc2 = futil.GetValue(E, ct, apar2, bpar2);

      h1.SetBinContent(xbin, ybin, bc1);
      h2.SetBinContent(xbin, ybin, bc2);

    }
  }

  //-----------------------------------------------------
  // set up the fitting with roofit
  //-----------------------------------------------------

  TRandom3 fRand(0);
  Double_t e = fRand.Uniform(-0.3, 0.3); // to add some error to the fitted pars

  // set guess values for the parameters
  ( (RooRealVar*)futil.GetSet().find("a") )->setVal(apar1 + e*apar1);
  ( (RooRealVar*)futil.GetSet().find("b") )->setVal(bpar1 + e*bpar1);

  // init the pdf's
  FitPDF fpdf1("mypdf1","mypdf1", &futil);
  FitPDF fpdf2("mypdf2","mypdf2", &futil);

  // import data to RooFit
  RooDataHist rf_h1("rf_h1", "rf_h1", futil.GetObs(), Import(h1) );
  RooDataHist rf_h2("rf_h2", "rf_h2", futil.GetObs(), Import(h2) );

  // create categories
  RooCategory cats("toydata", "my toy data");
  cats.defineType("hist1");
  cats.defineType("hist2");

  // combine data
  RooDataHist combData("combData", "combined data", futil.GetObs(), Index(cats), 
		      Import("hist1", rf_h1), Import("hist2", rf_h2) );
  
  // create a simultaneous pdf
  RooSimultaneous simPdf("simPdf", "simultaneous pdf", cats);
  simPdf.addPdf(fpdf1, "hist1");
  simPdf.addPdf(fpdf2, "hist2");

  // fit
  RooFitResult *roofit_result = simPdf.fitTo(combData, Save(kTRUE));
  RooArgSet final (roofit_result->floatParsFinal() );

  //-----------------------------------------------------
  // print comparison
  //-----------------------------------------------------

  cout << "******************************************************************" << endl;
  cout << "ROOFIT Parameter a initial values, fitted value: "
       << apar1 << "\t" << apar2 << "\t" << ((RooRealVar*)final.find("a"))->getVal() << endl;
  cout << "ROOFIT Parameter b initial value, fitted value: " 
       << bpar1 << "\t" << ((RooRealVar*)final.find("b"))->getVal() << endl;
  cout << "******************************************************************" << endl;

  TFile fout("fitted-hist.root", "RECREATE");
  h1.Write();
  fout.Close();

}
