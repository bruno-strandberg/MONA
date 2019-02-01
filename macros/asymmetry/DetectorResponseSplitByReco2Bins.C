#include "TSystem.h"
#include "TROOT.h"
#include "TDatime.h"
#include "TH3.h"
#include "TFile.h"
#include "TStopwatch.h"
#include "TString.h"

#include "DetResponse.h"
#include "SummaryParser.h"
#include "SummaryEvent.h"

#include "FitUtil.h"
#include "FitPDF.h"

#include "RooRealVar.h"


#include <iostream>
using namespace std;
using namespace RooFit;

void DetectorResponseSplitByReco2Bins() {

  const int N_PID_CLASSES = 2;
  TString filefolder = "./quality_detres/";

  gSystem->Load("$OSCPROBDIR/libOscProb.so");

  std::map<int, float> cut_map;
  cut_map.insert(std::make_pair(0, 0.0));
  cut_map.insert(std::make_pair(1, 0.6));
  cut_map.insert(std::make_pair(2, 1.0));

  //----------------------------------------------------------
  // detector response for tracks and showers
  //----------------------------------------------------------
  std::vector<DetResponse> track_response_vector;
  std::vector<DetResponse> shower_response_vector;
  std::vector<DetResponse> mc_response_vector;
  for (int i = 0; i < N_PID_CLASSES; i++) {
    std::function<bool(double, double)> comparison_operator;
    if (i == 0) { comparison_operator = std::greater_equal<double>(); // The first bin needs to include the lower limit.
    } else { comparison_operator = std::greater<double>(); }

    DetResponse track_response(DetResponse::track, TString::Format("track_response_%.2f", cut_map[i]), 40, 1, 100, 40, -1, 1, 1, 0, 1);
    track_response.AddCut( &SummaryEvent::Get_track_ql0       , std::greater<double>()   , 0.5         , true );
    track_response.AddCut( &SummaryEvent::Get_track_ql1       , std::greater<double>()   , 0.5         , true );
    track_response.AddCut( &SummaryEvent::Get_RDF_track_score , comparison_operator      , cut_map[i]  , true );
    track_response.AddCut( &SummaryEvent::Get_RDF_track_score , std::less_equal<double>(), cut_map[i+1], true );
    track_response.AddCut( &SummaryEvent::Get_RDF_muon_score  , std::less_equal<double>(), 0.05        , true );
    track_response.AddCut( &SummaryEvent::Get_RDF_noise_score , std::less_equal<double>(), 0.18        , true );
    track_response_vector.push_back(track_response);

    DetResponse shower_response(DetResponse::shower, TString::Format("shower_response_%.2f", cut_map[i]), 40, 1, 100, 40, -1, 1, 1, 0, 1);
    shower_response.AddCut( &SummaryEvent::Get_shower_ql0     , std::greater<double>()   , 0.5         , true );
    shower_response.AddCut( &SummaryEvent::Get_shower_ql1     , std::greater<double>()   , 0.5         , true );
    shower_response.AddCut( &SummaryEvent::Get_RDF_track_score, comparison_operator      , cut_map[i]  , true );
    shower_response.AddCut( &SummaryEvent::Get_RDF_track_score, std::less_equal<double>(), cut_map[i+1], true );
    shower_response.AddCut( &SummaryEvent::Get_RDF_muon_score , std::less_equal<double>(), 0.05        , true );
    shower_response.AddCut( &SummaryEvent::Get_RDF_noise_score, std::less_equal<double>(), 0.5         , true );
    shower_response_vector.push_back(shower_response);

    DetResponse mc_response(DetResponse::mc_truth, TString::Format("mc_response_%.2f", cut_map[i]), 40, 1, 100, 40, -1, 1, 1, 0, 1);
    mc_response.AddCut( &SummaryEvent::Get_RDF_track_score, comparison_operator      , cut_map[i]  , true );
    mc_response.AddCut( &SummaryEvent::Get_RDF_track_score, std::less_equal<double>(), cut_map[i+1], true );
    mc_response.AddCut( &SummaryEvent::Get_RDF_muon_score , std::less_equal<double>(), 0.05        , true );
    mc_response.AddCut( &SummaryEvent::Get_RDF_noise_score, std::less_equal<double>(), 0.5         , true );
    mc_response_vector.push_back(mc_response);
  }

  auto summary_file = (TString)getenv("NMHDIR") + "/data/ORCA_MC_summary_all_10Apr2018.root";
  SummaryParser sp(summary_file);
  for (Int_t i = 0; i < sp.GetTree()->GetEntries(); i++) {
    if (i % (Int_t)1e6 == 0) cout << "Event: " << i << endl;
    SummaryEvent *evt = sp.GetEvt(i);
    for (int i = 0; i < N_PID_CLASSES; i++) {
      track_response_vector[i].Fill(evt);
      shower_response_vector[i].Fill(evt);
      mc_response_vector[i].Fill(evt);
    }
  }

  for (int i = 0; i < N_PID_CLASSES; i++) {
    track_response_vector[i].WriteToFile(filefolder + TString::Format("track_response_%.2f.root" , cut_map[i]));
    shower_response_vector[i].WriteToFile(filefolder + TString::Format("shower_response_%.2f.root", cut_map[i]));
    mc_response_vector[i].WriteToFile(filefolder + TString::Format("mc_response_%.2f.root", cut_map[i]));
  }

  cout << "NOTICE: Finished filling response" << endl;

  //----------------------------------------------------------
  // set up the PDFs and static oscillation parameters
  //----------------------------------------------------------

  auto meff_file = (TString)getenv("NMHDIR") + "/data/eff_mass/EffMass_ORCA115_23x9m_ECAP0418.root";
  FitUtil *fitutil = new FitUtil(3, track_response_vector[0].GetHist3D(),
                 1, 100, -1, 1, 0, 1, meff_file);

  Double_t sinsqth12 = TMath::Power( TMath::Sin( 33.4 * TMath::Pi()/180. ), 2 );
  Double_t sinsqth13 = TMath::Power( TMath::Sin( 8.42 * TMath::Pi()/180. ), 2 );
  Double_t sinsqth23 = TMath::Power( TMath::Sin( 45   * TMath::Pi()/180. ), 2 );
  Double_t dcp       = 0.;
  Double_t dm32      = 2.44e-3;
  Double_t dm21      = 7.53e-5;
  Double_t DM        = dm32 + 0.5*dm21;

  // deconstrain th23 and dm31, when fitting you want constraints, otherwise you dont.
  fitutil->GetVar("Dm31")->setMin( -1 );
  fitutil->GetVar("Dm31")->setMax(  1 );
  fitutil->GetVar("SinsqTh23")->setMin( -1 );
  fitutil->GetVar("SinsqTh23")->setMax(  1 );

  // set parameter values 
  fitutil->GetVar("SinsqTh12")->setVal( sinsqth12 );
  fitutil->GetVar("SinsqTh13")->setVal( sinsqth13 );
  fitutil->GetVar("dcp")->setVal( dcp );
  fitutil->GetVar("Dm21")->setVal( dm21 );
  fitutil->GetVar("SinsqTh23")->setVal( sinsqth23 );

  for (int i = 0; i < N_PID_CLASSES; i++){
    FitPDF pdf_tracks("pdf_tracks", "pdf_tracks"   , fitutil, &track_response_vector[i]);
    FitPDF pdf_showers("pdf_showers", "pdf_showers", fitutil, &shower_response_vector[i]);
    FitPDF pdf_mc("pdf_mc", "pdf_mc", fitutil, &mc_response_vector[i]);
  
    //----------------------------------------------------------
    // set normal hierarchy
    //----------------------------------------------------------
    Double_t dm31 = DM + 0.5*dm21;
    fitutil->GetVar("Dm31")->setVal( dm31 );

    TH2D *tracks_NO  = (TH2D*)pdf_tracks.GetExpValHist()->Project3D("yx");
    TH2D *showers_NO = (TH2D*)pdf_showers.GetExpValHist()->Project3D("yx");
    TH2D *mc_NO      = (TH2D*)pdf_mc.GetExpValHist()->Project3D("yx");

    tracks_NO->SetNameTitle("detected_tracks", "detected_tracks");
    showers_NO->SetNameTitle("detected_showers", "detected_showers");
    mc_NO->SetNameTitle("detected_mc", "detected_mc");

    TH2D *tracks_NO_err  = (TH2D*)pdf_tracks.GetExpValErrHist()->Project3D("yx");
    TH2D *showers_NO_err = (TH2D*)pdf_showers.GetExpValErrHist()->Project3D("yx");
    TH2D *mc_NO_err      = (TH2D*)pdf_mc.GetExpValErrHist()->Project3D("yx");
    //----------------------------------------------------------
    // set inverted hierarchy
    //----------------------------------------------------------
    dm31 = -DM + 0.5*dm21;
    fitutil->GetVar("Dm31")->setVal( dm31 );

    TH2D *tracks_IO  = (TH2D*)pdf_tracks.GetExpValHist()->Project3D("yx");
    TH2D *showers_IO = (TH2D*)pdf_showers.GetExpValHist()->Project3D("yx");
    TH2D *mc_IO      = (TH2D*)pdf_mc.GetExpValHist()->Project3D("yx");

    tracks_IO->SetNameTitle("detected_tracks", "detected_tracks");
    showers_IO->SetNameTitle("detected_showers", "detected_showers");
    mc_IO->SetNameTitle("detected_mc", "detected_mc");

    TH2D *tracks_IO_err  = (TH2D*)pdf_tracks.GetExpValErrHist()->Project3D("yx");
    TH2D *showers_IO_err = (TH2D*)pdf_showers.GetExpValErrHist()->Project3D("yx");
    TH2D *mc_IO_err      = (TH2D*)pdf_mc.GetExpValErrHist()->Project3D("yx");
    //----------------------------------------------------------
    // save output
    //----------------------------------------------------------
    TString output_NO = TString::Format("split_expected_evts_NO_%.2f.root", cut_map[i]);
    TFile fout_NO(filefolder + output_NO,"RECREATE");
    auto hists_NO = {tracks_NO, showers_NO, mc_NO,
                     tracks_NO_err, showers_NO_err, mc_NO_err};
    for (auto hist: hists_NO) {
      hist->Write();
    }
    fout_NO.Close();

    TString output_IO = TString::Format("split_expected_evts_IO_%.2f.root", cut_map[i]);
    TFile fout_IO(filefolder + output_IO,"RECREATE");
    auto hists_IO = {tracks_IO, showers_IO, mc_IO,
                     tracks_IO_err, showers_IO_err, mc_IO_err};
    for (auto hist: hists_IO) {
      hist->Write();
    }
    fout_IO.Close();
  }
}
