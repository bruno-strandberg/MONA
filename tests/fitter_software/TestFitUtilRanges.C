
// RooFit
#include "RooArgSet.h"
#include "RooCategory.h"
#include "RooDataHist.h"
#include "RooFitResult.h"
#include "RooMsgService.h"
#include "RooRandom.h"
#include "RooRealVar.h"
#include "RooSimultaneous.h"

// NMH headers
#include "DetResponse.h"
#include "FitUtil.h"
#include "EffMass.h"
#include "FitPDF.h"
#include "AtmFlux.h"
#include "NuXsec.h"
#include "SummaryEvent.h"
#include "SummaryParser.h"

// cpp headers
#include <iostream>

using namespace std;
using namespace RooFit;

int main(const int argc, const char **argv) {
  
  int failed = 0;

  // Turn off INFO:Eval output, it gets crowded!
  RooMsgService::instance().getStream(1).removeTopic(Eval) ;

  Int_t N_RANGES = 1E3;

  gRandom->SetSeed(0); 

  // Fitter ranges
  Int_t fitEMin  = 2;
  Int_t fitEMax  = 80;
  Int_t fitctMin = -1;
  Int_t fitctMax = 0;

  //----------------------------------------------------------
  // detector response for tracks 
  //----------------------------------------------------------
  std::vector<DetResponse> track_response_vector;

  DetResponse track_response(DetResponse::track, "track_response", 24, 1, 100, 40, -1, 1, 1, 0, 1);
  track_response.AddCut( &SummaryEvent::Get_track_ql0       , std::greater<double>()   , 0.5 , true );
  track_response.AddCut( &SummaryEvent::Get_track_ql1       , std::greater<double>()   , 0.5 , true );
  track_response.AddCut( &SummaryEvent::Get_RDF_track_score , std::greater<double>()   , 0.6 , true );
  track_response.AddCut( &SummaryEvent::Get_RDF_muon_score  , std::less_equal<double>(), 0.05, true );
  track_response.AddCut( &SummaryEvent::Get_RDF_noise_score , std::less_equal<double>(), 0.18, true );

  //-----------------------------------------------------
  // fill the detector response and event selection
  //-----------------------------------------------------

  auto summary_file = (TString)getenv("MONADIR") + "/data/ORCA_MC_summary_all_10Apr2018.root";
  SummaryParser sp(summary_file);
  for (Int_t i = 0; i < sp.GetTree()->GetEntries(); i++) {
    if (i % (Int_t)1e6 == 0) cout << "Event: " << i << endl;
    SummaryEvent *evt = sp.GetEvt(i);
    track_response.Fill(evt);
  }

  cout << "NOTICE detector response ready" << endl;

  auto meff_file = (TString)getenv("MONADIR") + "/data/eff_mass/EffMass_ORCA115_23x9m_ECAP0418.root";
  FitUtil *fitutil = new FitUtil(3, track_response.GetHist3D(), fitEMin, fitEMax, fitctMin, fitctMax, 0, 1, meff_file);
  FitPDF pdf_test("pdf_test", "pdf_test", fitutil, &track_response);

  // Set central values and make expectation histograms
  fitutil->SetNOcentvals();

  TH3D* test_exp_val = (TH3D*)pdf_test.GetExpValHist();
  TString test_name = "test_expval_norange";
  test_exp_val->SetName(test_name);

  cout << "NOTICE full histogram generated" << endl;

  //----------------------------------------------------------
  // set up data for simultaneous fitting and fit
  //----------------------------------------------------------

  cout << "NOTICE generating ranges and applying to full histogram" << endl;

  for (Int_t i = 0; i < N_RANGES; i++) {

    Double_t E_lo_range = gRandom->Uniform(fitEMin, fitEMax);
    Double_t E_hi_range = gRandom->Uniform(E_lo_range, fitEMax);
    Double_t ct_lo_range = gRandom->Uniform(fitctMin, fitctMax);
    Double_t ct_hi_range = gRandom->Uniform(ct_lo_range, fitctMax);

    TString rangeName = Form("range_%i", i);
    fitutil->GetVar("E_reco")->setRange(  rangeName, E_lo_range, E_hi_range);
    fitutil->GetVar("ct_reco")->setRange( rangeName, ct_lo_range, ct_hi_range);

    TH3D* test_exp_val_range = (TH3D*)pdf_test.GetExpValHist(rangeName) ;
    test_exp_val_range->SetName("test_expval_" + rangeName);

    Int_t E_bin_min  = test_exp_val->GetXaxis()->FindBin( E_lo_range );
    Int_t E_bin_max  = test_exp_val->GetXaxis()->FindBin( E_hi_range );
    Int_t ct_bin_min = test_exp_val->GetYaxis()->FindBin( ct_lo_range );
    Int_t ct_bin_max = test_exp_val->GetYaxis()->FindBin( ct_hi_range );

    // The Bjorken-y Z-axis contains only 1 bin, so the integral goes from 1 to 1.
    Double_t integral_full_range = test_exp_val->Integral();
    Double_t integral_rangeName = test_exp_val_range->Integral();
    Double_t integral_integralRange = test_exp_val_range->Integral( E_bin_min, E_bin_max, ct_bin_min, ct_bin_max, 1, 1 );
    //cout << "Integral of histo with no range is:       " << integral_full_range << endl;
    //cout << "Integral of histo with rangeName is:      " << integral_rangeName << endl;
    //cout << "Integral of histo with integral range is: " << integral_integralRange << endl;
 
    // If the integral of the full histogram is less than an integral of a histogram with a range, 
    // there is a problem. If the integral of the hist w/ range through the ExpVal method is 
    // different from the integral w/ range directly used in Integral, there is a problem.
    if ((integral_full_range < integral_rangeName) or (integral_rangeName != integral_integralRange)) {
      failed = 1;
    }

  }
  
  //-------------------------------------------------------------
  // evaluate if the test failed
  //-------------------------------------------------------------

  if (failed) cout << "NOTICE TestFit failed" << endl;
  else cout << "NOTICE TestFit passed" << endl;

  return failed;
  
}
