
// ROOT headers
#include "TIterator.h"

// RooFit
#include "RooArgSet.h"
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooFitResult.h"
#include "RooRandom.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"

// NMH headers
#include "DetResponse.h"
#include "FitUtil.h"
#include "EffMass.h"
#include "FitPDF.h"
#include "AtmFlux.h"
#include "NuXsec.h"
#include "SummaryEvent.h"

// cpp headers
#include <iostream>

using namespace std;
using namespace RooFit;

/**
   This application tests simultaneous fitting. It is a deterministic application with some fixed seeds to check that the fit always gives the same result.
 */
int main(const int argc, const char **argv) {

  int failed = 0;
  
  //-------------------------------------------------------------
  // fill response with pseudo-data
  //-------------------------------------------------------------
  
  // create a responses and fill them with pseudo-data
  DetResponse trkR(DetResponse::track, "trackresp", 40, 1, 100, 40, -1, 1, 1, 0, 1);
  trkR.AddCut( &SummaryEvent::Get_RDF_track_score, std::greater<double>(), 0.6, true );
  trkR.AddCut( &SummaryEvent::Get_track_ql0      , std::greater<double>(), 0.5, true );

  DetResponse shwR(DetResponse::shower, "showerresp", 40, 1, 100, 40, -1, 1, 1, 0, 1);
  trkR.AddCut( &SummaryEvent::Get_RDF_track_score, std::less_equal<double>(), 0.6, true );
  trkR.AddCut( &SummaryEvent::Get_shower_ql0     , std::greater<double>()   , 0.5, true );
  
  SummaryEvent evt;
  evt.SetSeed(6928);
  
  for (Int_t i = 0; i < 1e7; i++) {
    evt.FillPseudoData();
    trkR.Fill(&evt);
    shwR.Fill(&evt);
  }

  cout << "NOTICE TestSimFit response ready" << endl;

  //-------------------------------------------------------------
  // init fitutil and fitpdf
  //-------------------------------------------------------------
  
  // init fitutil with dummy effective mass
  FitUtil futil(3, trkR.GetHist3D(), 3, 80, -1, 0, 0, 1, EffMass::DUMMYFILE);

  // init the pdfs
  FitPDF trkpdf("trkpdf","trkpdf", &futil, &trkR);
  FitPDF shwpdf("shwpdf","shwpdf", &futil, &shwR);

  //-------------------------------------------------------------
  // create expectation-value data at some oscillation parameter values
  //-------------------------------------------------------------

  RooRandom::randomGenerator()->SetSeed(99032);
  
  futil.SetNOlims();  
  futil.SetNOcentvals();
  futil.GetVar("SinsqTh13")->randomize();
  futil.GetVar("SinsqTh23")->randomize();
  futil.GetVar("Dm31")->randomize();

  // create expectation value histogram to be fitted
  TH3D* trkdata = trkpdf.GetExpValHist();
  TH3D* shwdata = shwpdf.GetExpValHist();

  // save true values
  RooArgSet *startvals = (RooArgSet*)futil.GetSet().snapshot(kTRUE);
  
  //-------------------------------------------------------------
  // set free parameter limits and start at NO central values; fix dcp, th12 and dm21; fit
  //-------------------------------------------------------------

  futil.FreeParLims();
  futil.SetNOcentvals();
  futil.GetVar("Dm21")->setConstant();
  futil.GetVar("SinsqTh12")->setConstant();
  futil.GetVar("dcp")->setConstant();  

  // prep for data import to RooFit
  TString trkcateg = "tracks";
  TString shwcateg = "showers";

  std::map<string, TH1* > hist_map = { { (string)trkcateg, trkdata }, 
				       { (string)shwcateg, shwdata } };
  
  RooCategory categs("categs","data categories");
  categs.defineType( trkcateg );
  categs.defineType( shwcateg );
  
  // create combined data set
  RooDataHist combData("combData", "combined data", futil.GetObs(), categs, hist_map);

  // create simultaneous pdf
  RooSimultaneous simPdf("simPdf", "simultaneous Pdf", categs);
  simPdf.addPdf(trkpdf, trkcateg );
  simPdf.addPdf(shwpdf, shwcateg );

  // create ranges for theta-23
  futil.GetVar("SinsqTh23")->setRange("firstq" , 0. , 0.5);
  futil.GetVar("SinsqTh23")->setRange("secondq", 0.5, 1);
  
  futil.GetVar("SinsqTh23")->setVal(0.4);
  RooFitResult *result1 = simPdf.fitTo( combData, Save(kTRUE), Range("firstq") );
  futil.SetNOcentvals();
  futil.GetVar("SinsqTh23")->setVal(0.6);
  RooFitResult *result2 = simPdf.fitTo( combData, Save(kTRUE), Range("secondq") );

  RooFitResult *result = ( result1->minNll() < result2->minNll() ) ? result1 : result2 ;
  
  //-------------------------------------------------------------
  // compare fitted and true values
  //-------------------------------------------------------------
  RooArgList fitvals = result->floatParsFinal();

  Double_t expected_th23 = 0.496957;
  Double_t expected_th13 = 0.0193736;
  Double_t expected_dm31 = 0.00257888;

  Double_t start_th23 = ((RooRealVar*)startvals->find("SinsqTh23"))->getVal();
  Double_t start_th13 = ((RooRealVar*)startvals->find("SinsqTh13"))->getVal();
  Double_t start_dm31 = ((RooRealVar*)startvals->find("Dm31")     )->getVal();
  
  Double_t fitted_th23 = ((RooRealVar*)fitvals.find("SinsqTh23"))->getVal();
  Double_t fitted_th13 = ((RooRealVar*)fitvals.find("SinsqTh13"))->getVal();
  Double_t fitted_dm31 = ((RooRealVar*)fitvals.find("Dm31")     )->getVal();

  cout << "Theta-23  start and fitted: " << start_th23 << "\t" << fitted_th23 << endl;
  cout << "Theta-13  start and fitted: " << start_th13 << "\t" << fitted_th13 << endl;
  cout << "DeltaM-31 start and fitted: " << start_dm31 << "\t" << fitted_dm31 << endl;

  // allow 0.5% change due to numerical accuracy
  if ( TMath::Abs(expected_th23 - fitted_th23)/expected_th23*100 > 0.5 ||
       TMath::Abs(expected_th13 - fitted_th13)/expected_th13*100 > 0.5 ||
       TMath::Abs(expected_dm31 - fitted_dm31)/expected_dm31*100 > 0.5 ) {
    cout << "Theta-23  expected and fitted: " << expected_th23 << "\t" << fitted_th23 << endl;
    cout << "Theta-13  expected and fitted: " << expected_th13 << "\t" << fitted_th13 << endl;
    cout << "DeltaM-31 expected and fitted: " << expected_dm31 << "\t" << fitted_dm31 << endl;
    failed = 1;
  }
  
  if (failed) cout << "NOTICE TestFit failed" << endl;
  else cout << "NOTICE TestFit passed" << endl;

  return failed;
  
}
