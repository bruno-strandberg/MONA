
// RooFit
#include "RooArgSet.h"
#include "RooCategory.h"
#include "RooDataHist.h"
#include "RooFitResult.h"
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

  Int_t N_RANGES = 10;

  std::vector< std::tuple<Double_t, Double_t> > fitRanges;
  fitRanges.push_back( std::make_tuple(3, 50) );
  fitRanges.push_back( std::make_tuple(3, 50) );
  fitRanges.push_back( std::make_tuple(3, 50) );

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
    if (i == (Int_t)1e6) break;
    SummaryEvent *evt = sp.GetEvt(i);
    track_response.Fill(evt);
  }

  cout << "NOTICE TestFit response ready" << endl;

  auto meff_file = (TString)getenv("MONADIR") + "/data/eff_mass/EffMass_ORCA115_23x9m_ECAP0418.root";
  FitUtil *fitutil = new FitUtil(3, track_response.GetHist3D(), fitEMin, fitEMax, fitctMin, fitctMax, 0, 1, meff_file);

  FitPDF pdf_test( "pdf_test", "pdf_test", fitutil, &track_response );
  fitutil->GetVar("E_reco")->setRange( "test_range", 4.56, 49.9);
  fitutil->GetVar("ct_reco")->setRange( "test_range", -0.8, -0.1);

  // Set IO values and make expectation histograms
  fitutil->SetNOcentvals();

  track_exp_val = (TH3D*)pdf_tracks.GetExpValHist() ;
  TString track_name = "tracks_expval";
  track_exp_val->SetName(track_name);

  //----------------------------------------------------------
  // set up data for simultaneous fitting and fit
  //----------------------------------------------------------

//  for (Int_t i = 0; i < N_RANGES; i++) {
//    fitutil->GetVar("E_reco")->setRange( Form("range_%i", i), std::get<0>(fitRanges[i]), std::get<1>(fitRanges[i]));
//  }

  Double_t fitrange_test_lower = fitutil->GetVar("E_reco")->getMin("fitRangeE_showers_expval_true_0.00");
  Double_t fitrange_test_upper = fitutil->GetVar("E_reco")->getMax("fitRangeE_showers_expval_true_0.00");

  cout << "Range test lower " << fitrange_test_lower << endl;
  cout << "Range test upper " << fitrange_test_upper << endl;







  TH3D* testhist_range    = (TH3D*)pdf_test.GetExpValHist("test_range") ;
  testhist_range->SetName("test_range");
  TH3D* testhist_no_range = (TH3D*)pdf_test.GetExpValHist() ;
  testhist_no_range->SetName("test_no_range");

  cout << "Integral of histo with fitRange is:    " << testhist_range->Integral() << endl;
  cout << "Integral of histo with no fitRange is: " << testhist_no_range->Integral() << endl;

  TH2D* testhist_range_2d    = (TH2D*)testhist_range->Project3D("yx");
  TH2D* testhist_no_range_2d = (TH2D*)testhist_no_range->Project3D("yx");

  cout << "Integral of 2d with fitRange    " <<  testhist_range_2d->Integral() << endl;
  cout << "Integral of 2d with no fitRange " <<  testhist_no_range_2d->Integral() << endl;

  TFile fout("test.root", "RECREATE");
  testhist_range_2d->Write("testhist_range");
  testhist_no_range_2d->Write("testhist_no_range");
  fout.Close();





  
  //-------------------------------------------------------------
  // compare fitted and true values
  //-------------------------------------------------------------
//  RooArgList fitvals = result->floatParsFinal();


  if (failed) cout << "NOTICE TestFit failed" << endl;
  else cout << "NOTICE TestFit passed" << endl;

  return failed;
  
}
