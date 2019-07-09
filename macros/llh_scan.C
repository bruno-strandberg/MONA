#include "SummaryParser.h"
#include "DetResponse.h"
#include "NMHUtils.h"
#include "FitUtil.h"
#include "FitPDF.h"
#include "RooDataHist.h"
#include "RooNLLVar.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include <iostream>

using namespace std;
using namespace RooFit;

/** This macro demonstrates the creation of a log-likelihood scan in the track channel. See README.md to see what is required to run this macro.*/

void llh_scan(TString dataf = (TString)getenv("MONADIR") + (TString)"/data/ORCA_MCsummary_SEv2_ORCA115_20x9m_ECAP190222.root",
	      TString effmf = (TString)getenv("MONADIR") + (TString)"/data/eff_mass/EffMass_ORCA115_20x9m_ECAP190222.root") {
  
  // initialise a detector response, the numbers configure the binning
  DetResponse DR(DetResponse::track, "DR", 40, 1, 100, 40, -1, 1, 1, 0, 1);

  // add selection cuts
  DR.AddCut( &SummaryEvent::Get_RDF_track_score, std::greater<double>(), 0.7, kTRUE ); // ask PID track score to be larger than > 0.7
  DR.AddCut( &SummaryEvent::Get_track_ql2      , std::greater<double>(), 0.5, kTRUE ); // ask for ql2 (gandalf_loose_is_selected), see `MONA/apps/data_sorting/ECAP190222_20m_to_MONA.C`

  // fill the response with MC events
  SummaryParser sp( dataf );
  for (Int_t i = 0; i < sp.GetTree()->GetEntries(); i++) {
    if (i%500000 == 0) cout << "Event: " << i << "/" << sp.GetTree()->GetEntries() << endl;
    DR.Fill( sp.GetEvt(i) );
  }
    
  // initialise fitter utility and a PDF with the above response; the first number is operation time in years,
  // the rest define the fit range. FitUtil is the worker class where calculations are done and parameters are defined
  FitUtil fu(3, DR.GetHist3DTrue(), DR.GetHist3DReco(), 3, 80, -1, -1e-5, 0, 1, effmf);

  // init the probability density function; this is just a wrapper class that uses FitUtil
  // and DetResponse to create a RooFit probability density function
  FitPDF  pdf("pdf","pdf", &fu, &DR);

  // create a pseudo-experiment, 3 years of data taking (run time determined at FitUtil construction)
  TH3D* D = pdf.SimplePseudoExp("Data");

  // Import data to RooFit
  RooDataHist RFD("RFD","RFD", fu.GetObs(), Import(*D) );

  // create a likelihood function between the pdf and the data
  RooAbsReal *nll = pdf.createNLL( RFD );

  // change th23 and evaluate the likelihood
  fu.GetVar("SinsqTh23")->setVal(0.6);
  cout << "NOTICE llh_scan() log-likelihood: " << nll->getVal() << endl;

  // let RooFit to scan the likelihood
  RooPlot *frame = fu.GetVar("SinsqTh23")->frame();
  nll->plotOn( frame, ShiftToZero() );

  TCanvas *c1 = new TCanvas("c1","c1",1);
  D->Project3D("yx")->Draw("colz");
  
  TCanvas *c2 = new TCanvas("c2","c2",1);
  frame->Draw();
  
}
