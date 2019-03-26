
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

  const int N_PID_CLASSES = 3;
  const Double_t PID_CUT = 0.6;

  std::map<Int_t, Double_t> pid_map;
  cout << "NOTICE: Set PID case " << N_PID_CLASSES << endl;
  pid_map.insert(std::make_pair(0, 0.0)); // shower
  pid_map.insert(std::make_pair(1, 0.4)); // middle group: shower
  pid_map.insert(std::make_pair(2, 0.6)); // track
  pid_map.insert(std::make_pair(3, 1.0)); // upper limit 

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
  // detector response for tracks and showers
  //----------------------------------------------------------
  std::vector<DetResponse> track_response_vector;
  std::vector<DetResponse> shower_response_vector;

  for (Int_t i = 0; i < N_PID_CLASSES; i++) {
    std::function<bool(double, double)> comparison_operator;
    if (i == 0) { comparison_operator = std::greater_equal<double>(); // The first bin needs to include the lower limit.
    } else { comparison_operator = std::greater<double>(); }

    DetResponse track_response(DetResponse::track, Form("track_response_%.2f", pid_map[i]), 24, 1, 100, 40, -1, 1, 1, 0, 1);
    track_response.AddCut( &SummaryEvent::Get_track_ql0       , std::greater<double>()   , 0.5         , true );
    track_response.AddCut( &SummaryEvent::Get_track_ql1       , std::greater<double>()   , 0.5         , true );
    track_response.AddCut( &SummaryEvent::Get_RDF_track_score , comparison_operator      , pid_map[i]  , true );
    track_response.AddCut( &SummaryEvent::Get_RDF_track_score , std::less_equal<double>(), pid_map[i+1], true );
    track_response.AddCut( &SummaryEvent::Get_RDF_muon_score  , std::less_equal<double>(), 0.05        , true );
    track_response.AddCut( &SummaryEvent::Get_RDF_noise_score , std::less_equal<double>(), 0.18        , true );
    track_response_vector.push_back(track_response);

    DetResponse shower_response(DetResponse::shower, Form("shower_response_%.2f", pid_map[i]), 24, 1, 100, 40, -1, 1, 1, 0, 1);
    shower_response.AddCut( &SummaryEvent::Get_shower_ql0     , std::greater<double>()   , 0.5         , true );
    shower_response.AddCut( &SummaryEvent::Get_shower_ql1     , std::greater<double>()   , 0.5         , true );
    shower_response.AddCut( &SummaryEvent::Get_RDF_track_score, comparison_operator      , pid_map[i]  , true );
    shower_response.AddCut( &SummaryEvent::Get_RDF_track_score, std::less_equal<double>(), pid_map[i+1], true );
    shower_response.AddCut( &SummaryEvent::Get_RDF_muon_score , std::less_equal<double>(), 0.05        , true );
    shower_response.AddCut( &SummaryEvent::Get_RDF_noise_score, std::less_equal<double>(), 0.5         , true );
    shower_response_vector.push_back(shower_response);
  }

  //-----------------------------------------------------
  // fill the detector response and event selection
  //-----------------------------------------------------

  auto summary_file = (TString)getenv("MONADIR") + "/data/ORCA_MC_summary_all_10Apr2018.root";
  SummaryParser sp(summary_file);
  for (Int_t i = 0; i < sp.GetTree()->GetEntries(); i++) {
    if (i % (Int_t)1e6 == 0) cout << "Event: " << i << endl;
    if (i == (Int_t)1e6) break;
    SummaryEvent *evt = sp.GetEvt(i);
    for (Int_t i = 0; i < N_PID_CLASSES; i++) {
      track_response_vector[i].Fill(evt);
      shower_response_vector[i].Fill(evt);
    }
  }

  cout << "NOTICE TestFit response ready" << endl;

  auto meff_file = (TString)getenv("MONADIR") + "/data/eff_mass/EffMass_ORCA115_23x9m_ECAP0418.root";
  FitUtil *fitutil = new FitUtil(3, track_response_vector[0].GetHist3D(), fitEMin, fitEMax, fitctMin, fitctMax, 0, 1, meff_file);

  std::vector<TH3D*> track_vector_true;
  std::vector<TH3D*> shower_vector_true;

  std::vector<FitPDF> pdf_tracks_vector;
  std::vector<FitPDF> pdf_showers_vector;

  for (int i = 0; i < N_PID_CLASSES; i++) {
    FitPDF pdf_tracks(  Form("pdf_tracks_%.2f", pid_map[i]),  "pdf_tracks",  fitutil, &track_response_vector[i]);
    FitPDF pdf_showers( Form("pdf_showers_%.2f", pid_map[i]), "pdf_showers", fitutil, &shower_response_vector[i]);
    pdf_tracks_vector.push_back( pdf_tracks );
    pdf_showers_vector.push_back( pdf_showers );

    // Set IO values and make expectation histograms
    fitutil->SetIOcentvals();

    track_vector_true.push_back(  (TH3D*)pdf_tracks.GetExpValHist("track") );
    shower_vector_true.push_back( (TH3D*)pdf_showers.GetExpValHist("shower") );
    TString track_name_true  = Form("tracks_expval_true_%.2f", pid_map[i]);
    TString shower_name_true = Form("showers_expval_true_%.2f", pid_map[i]);
    track_vector_true[i] ->SetName(track_name_true);
    shower_vector_true[i]->SetName(shower_name_true);

    // Set parameters to constant
    fitutil->GetVar("SinsqTh12")->setConstant(kTRUE);
    fitutil->GetVar("SinsqTh13")->setConstant(kTRUE);
    fitutil->GetVar("dcp")->setConstant(kTRUE);
    fitutil->GetVar("Dm21")->setConstant(kTRUE);
  }

  // create quantiles for theta23
  fitutil->GetVar("SinsqTh23")->setRange("firstq" , 0., 0.5);
  fitutil->GetVar("SinsqTh23")->setRange("secondq", 0.5, 1.);

  //----------------------------------------------------------
  // set up data for simultaneous fitting and fit
  //----------------------------------------------------------
  cout << "NOTICE Fitter started fitting" << endl;
  
  TStopwatch timer;

  // Fit under NO model, IO data
  fitutil->SetNOlims();
  fitutil->GetVar("SinsqTh23")->setMin(0.);
  fitutil->GetVar("SinsqTh23")->setMax(1.);
  fitutil->SetNOcentvals();

  std::map<string, TH1*> hist_map;
  std::vector<TString> fitRangeCategories;
  RooCategory cats("categories", "data categories");
  for (Int_t i = 0; i < N_PID_CLASSES; i++) {
    if (pid_map[i] < PID_CUT) {
      hist_map.insert( {(string)shower_vector_true[i]->GetName(), shower_vector_true[i] } );
      cats.defineType( shower_vector_true[i]->GetName() );
      fitRangeCategories.push_back( shower_vector_true[i]->GetName() );

      cout << "NOTICE: Added hist and cat to shower" << endl;
    }
    else {
      hist_map.insert( {(string)track_vector_true[i]->GetName(), track_vector_true[i] } );
      cats.defineType( track_vector_true[i]->GetName() );
      fitRangeCategories.push_back( track_vector_true[i]->GetName() );

      cout << "NOTICE: Added hist and cat to track" << endl;
    }
  }

  RooSimultaneous simPdf("simPdf", "simultaneous Pdf for IO", cats);
  for (Int_t i = 0; i < N_PID_CLASSES; i++) {
    if (pid_map[i] < PID_CUT) {
      simPdf.addPdf(pdf_showers_vector[i], shower_vector_true[i]->GetName() );
      cout << "NOTICE: Added simpdf to shower " << shower_vector_true[i]->GetName() << endl;
    }
    else {
      simPdf.addPdf(pdf_tracks_vector[i],  track_vector_true[i]->GetName() );
      cout << "NOTICE: Added simpdf to track " << track_vector_true[i]->GetName() <<  endl;
    }
  }

  RooDataHist data_hists("data_hists", "track and shower data", fitutil->GetObs(), cats, hist_map);

  // create ranges for fitter
  //fitutil->GetVar("E_reco")->setRange( TString("fitRangeE"), fitEMin, fitEMax);
  for (Int_t i = 0; i < N_PID_CLASSES; i++) {
    fitutil->GetVar("E_reco")->setRange( TString("fitRangeE_") + fitRangeCategories[i],  // ITS NOT THE TSTRING
                                        std::get<0>(fitRanges[i]), std::get<1>(fitRanges[i]));
  }

  Double_t fitrange_test_lower = fitutil->GetVar("E_reco")->getMin("fitRangeE_showers_expval_true_0.00");
  Double_t fitrange_test_upper = fitutil->GetVar("E_reco")->getMax("fitRangeE_showers_expval_true_0.00");

  cout << "Range test lower " << fitrange_test_lower << endl;
  cout << "Range test upper " << fitrange_test_upper << endl;

  FitPDF pdf_test( "pdf_test", "pdf_test", fitutil, &shower_response_vector[0] );
  fitutil->GetVar("E_reco")->setRange( "test_range", 4.56, 49.9);
  fitutil->GetVar("ct_reco")->setRange( "test_range", -0.8, -0.1);
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
