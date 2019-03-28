
// RooFit
#include "RooDataHist.h"
#include "RooMsgService.h"
#include "RooRandom.h"
#include "RooRealVar.h"

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

  Int_t N_RANGES = 1E3; // Number of ranges tested against.

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

  //-----------------------------------------------------
  // create fitutil and make pdfs for ExpValHists
  //-----------------------------------------------------

  auto meff_file = (TString)getenv("MONADIR") + "/data/eff_mass/EffMass_ORCA115_23x9m_ECAP0418.root";
  FitUtil *fitutil = new FitUtil(3, track_response.GetHist3D(), fitEMin, fitEMax, fitctMin, fitctMax, 0, 1, meff_file);
  FitPDF pdf_test("pdf_test", "pdf_test", fitutil, &track_response);

  // Set central values and make expectation histograms
  fitutil->SetNOcentvals();

  TH3D* test_exp_val = (TH3D*)pdf_test.GetExpValHist();
  TString test_name = "test_expval_norange";
  test_exp_val->SetName(test_name);

  cout << "NOTICE full histogram generated" << endl;

  // Evaluate range for RooDataHist
  RooRealVar E_reco("E_reco", "E_reco", 1, 100);
  RooRealVar ct_reco("ct_reco", "ct_reco", -1, 1);
  RooRealVar by_reco("by_reco", "by_reco", 0, 1);

  RooDataHist roo_test_hist("roo_test_hist", "roo_test_hist", RooArgSet(E_reco, ct_reco, by_reco), Import(*test_exp_val));
  // Check https://root.cern.ch/doc/master/classRooDataHist.html for sum() documentation
  Double_t rdh_integral = roo_test_hist.sum(0); 
  Double_t integral_full_range = test_exp_val->Integral();
  // If the difference is less than 1%, they are the same. (Normal values are 1E-13)
  if (TMath::Abs(integral_full_range - rdh_integral) > 1E-2) { 
    cout << "NOTICE Test failed on RooDataHist full integration vs. TH3D full integration" << endl;
    failed = 1;
  }

  //----------------------------------------------------------
  // set up data for range checking
  //----------------------------------------------------------

  cout << "NOTICE generating ranges and applying to full histogram" << endl;

  for (Int_t i = 0; i < N_RANGES; i++) {

    //-------------------------------------------------------------
    // Test how TH1 and GetExpValHist(rangeName) evaluate ranges
    //-------------------------------------------------------------

    // I took max / 2. as the upper bound to avoid getting empty intervals
    Double_t E_lo_range = gRandom->Uniform(fitEMin, fitEMax / 2.);
    Double_t E_hi_range = gRandom->Uniform(E_lo_range, fitEMax);
    Double_t ct_lo_range = gRandom->Uniform(fitctMin, fitctMax / 2.);
    Double_t ct_hi_range = gRandom->Uniform(ct_lo_range, fitctMax);

    // RooFit handles ranges on different RooRealVars with the same rangeName properly: they do not collide.
    TString rangeName = Form("range_%i", i);
    fitutil->GetVar("E_reco")->setRange(  rangeName, E_lo_range, E_hi_range);
    fitutil->GetVar("ct_reco")->setRange( rangeName, ct_lo_range, ct_hi_range);

    TH3D* test_exp_val_range = (TH3D*)pdf_test.GetExpValHist(rangeName) ;
    test_exp_val_range->SetName("test_expval_" + rangeName);

    // This takes the bin number in which the value falls. Similar to FitUtil::GetRange
    Int_t E_bin_min  = test_exp_val->GetXaxis()->FindBin( E_lo_range );
    Int_t E_bin_max  = test_exp_val->GetXaxis()->FindBin( E_hi_range );
    Int_t ct_bin_min = test_exp_val->GetYaxis()->FindBin( ct_lo_range );
    Int_t ct_bin_max = test_exp_val->GetYaxis()->FindBin( ct_hi_range );

    // The Bjorken-y Z-axis contains only 1 bin, so the integral goes from 1 to 1.
    Double_t integral_rangeName = test_exp_val_range->Integral();
    Double_t integral_integralRange = test_exp_val->Integral( E_bin_min, E_bin_max, ct_bin_min, ct_bin_max, 1, 1 );
 
    // If the integral of the full histogram is less than an integral of a histogram with a range, 
    // there is a problem. If the integral of the hist w/ range through the ExpVal method is 
    // different from the integral w/ range directly used in Integral, there is a problem.
    if ((integral_full_range < integral_rangeName) or (integral_rangeName != integral_integralRange)) {
      failed = 1;
    }

    //-------------------------------------------------------------
    // Test how RooDataHist handles ranges when integrating
    //-------------------------------------------------------------

    // RooDataHist integrates (with sum()) over continuous ranges. This behaviour can only be reconciled with 
    // TH1::Integral() behaviour by passing the RDH ranges that are exact bin edges.
    // Code snippet from https://root.cern.ch/doc/master/RooHistPdf_8cxx_source.html#l00341
    std::map<const RooAbsArg*, std::pair<Double_t, Double_t> > ranges;
    ranges[&E_reco]  = std::make_pair( test_exp_val->GetXaxis()->GetBinLowEdge( E_bin_min ), test_exp_val->GetXaxis()->GetBinUpEdge( E_bin_max ) );
    ranges[&ct_reco] = std::make_pair( test_exp_val->GetYaxis()->GetBinLowEdge( ct_bin_min ), test_exp_val->GetYaxis()->GetBinUpEdge( ct_bin_max ) );
    ranges[&by_reco] = std::make_pair( test_exp_val->GetZaxis()->GetBinLowEdge( 1 ), test_exp_val->GetZaxis()->GetBinUpEdge( 1 ) ); // By goes from bin 1 to bin 1

    Double_t rdh_integralRange = roo_test_hist.sum(RooArgSet(E_reco, ct_reco, by_reco), RooArgSet(E_reco, ct_reco, by_reco), kFALSE, kFALSE, ranges);
    // If the difference is more than 1%, they are not the same. (Normal values are 1E-13)
    if (TMath::Abs(integral_integralRange - rdh_integralRange) > 1E-2) {
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
