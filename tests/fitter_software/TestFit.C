
// ROOT headers
#include "TIterator.h"

// RooFit
#include "RooArgSet.h"
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooFitResult.h"
#include "RooRandom.h"

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

int main(const int argc, const char **argv) {

  int failed = 0;
  
  //-------------------------------------------------------------
  // fill response with pseudo-data
  //-------------------------------------------------------------
  
  // create a response and fill it with pseudo-data of muons
  DetResponse resp(DetResponse::mc_truth, "truthresp", 40, 1, 100, 40, -1, 1, 1, 0, 1);
  resp.AddCut( &SummaryEvent::Get_MC_type  , std::equal_to<double>() ,  14, false );
  resp.AddCut( &SummaryEvent::Get_MC_type  , std::equal_to<double>() , -14, false );
  resp.AddCut( &SummaryEvent::Get_MC_is_CC , std::greater<double>()  , 0.5, true );
  
  SummaryEvent evt;

  for (Int_t i = 0; i < 1e7; i++) {
    evt.FillPseudoData();
    resp.Fill(&evt);
  }

  cout << "NOTICE TestFit response ready" << endl;

  //-------------------------------------------------------------
  // init fitutil and fitpdf
  //-------------------------------------------------------------
  
  // init fitutil with dummy effective mass
  FitUtil futil(3, resp.GetHist3D(), 3, 80, -1, 0, 0, 1, EffMass::DUMMYFILE);

  // init the pdf
  FitPDF pdf("pdf","pdf", &futil, &resp);

  //-------------------------------------------------------------
  // create expectation-value data at some oscillation parameter values
  //-------------------------------------------------------------

  RooRandom::randomGenerator()->SetSeed(0);
  
  futil.SetNOlims();  
  futil.SetNOcentvals();
  futil.GetVar("SinsqTh13")->randomize();
  futil.GetVar("SinsqTh23")->randomize();
  futil.GetVar("Dm31")->randomize();

  // create expectation value histogram to be fitted
  TH3D* datah = pdf.GetExpValHist();

  // save true values
  RooArgSet *startvals = (RooArgSet*)futil.GetSet().snapshot(kTRUE);
  
  //-------------------------------------------------------------
  // set free parameter limits and start at NO central values and fix dcp; fit
  //-------------------------------------------------------------

  futil.FreeParLims();
  futil.SetNOcentvals();
  futil.GetVar("Dm21")->setConstant();
  futil.GetVar("SinsqTh12")->setConstant();
  futil.GetVar("dcp")->setConstant();

  futil.GetVar("SinsqTh23")->setRange("firstq" , 0. , 0.5);
  futil.GetVar("SinsqTh23")->setRange("secondq", 0.5, 1);
  
  // import data to RooFit and fit
  RooDataHist rfdatah("rfdata","rfdata", futil.GetObs(), Import(*datah) );

  futil.GetVar("SinsqTh23")->setVal(0.4);
  RooFitResult *result1 = pdf.fitTo( rfdatah, Save(kTRUE), Range("firstq") );
  futil.SetNOcentvals();
  futil.GetVar("SinsqTh23")->setVal(0.6);
  RooFitResult *result2 = pdf.fitTo( rfdatah, Save(kTRUE), Range("secondq") );

  RooFitResult *result = ( result1->minNll() < result2->minNll() ) ? result1 : result2 ;
  
  //-------------------------------------------------------------
  // compare fitted and true values
  //-------------------------------------------------------------
  RooArgList fitvals = result->floatParsFinal();

  TIterator *iter = fitvals.createIterator();
  RooRealVar *fitvar = NULL;
  while ( ( fitvar = (RooRealVar*)iter->Next() ) ) {
    RooRealVar *startvar = (RooRealVar*)startvals->find( fitvar->GetName() );

    Double_t start_value  = startvar->getVal();
    Double_t fitted_value = fitvar->getVal();

    // allow a 10% deviation from the true value
    if ( TMath::Abs(start_value - fitted_value)/start_value > 0.1 ) {
      cout << "NOTICE TestFit failed on variable: "
	   << fitvar->GetName() << " true: " << startvar->getVal() << " fitted: " << fitvar->getVal() << endl;
      failed = 1;
    }
    else {
      cout << "NOTICE TestFit passed on variable: "
	   << fitvar->GetName() << " true: " << startvar->getVal() << " fitted: " << fitvar->getVal() << endl;
    }
    
  }

  if (failed) cout << "NOTICE TestFit failed" << endl;
  else cout << "NOTICE TestFit passed" << endl;

  return failed;
  
}
