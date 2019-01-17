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

#include "CustomEventFilters.h"

#include <iostream>
using namespace std;
using namespace RooFit;

void DetectorResponseSplitByRecoPID() {

  const int N_PID_CLASSES = 10;
  Double_t PID_step = 1 / float(N_PID_CLASSES);
  TString filefolder = TString::Format("./energy_detres/pid_bins_%i/", N_PID_CLASSES);

  gSystem->Load("$OSCPROBDIR/libOscProb.so");

  //----------------------------------------------------------
  // detector response for tracks and showers
  //----------------------------------------------------------
  // Good tracks only, no overlap with good showers
  // gt = good track, gs = good shower, ge = good event
  // The track response has been split into these 3 types, shower has not.
  std::vector<DetResponse> response_good_track_vector; // good tracks in track channel
  std::vector<DetResponse> response_good_shower_vector; // good showers in track channel
  std::vector<DetResponse> response_good_event_vector; // good events in track channel
  std::vector<DetResponse> response_shower_vector; // showers
  std::vector<DetResponse> response_mc_vector; // mc

  for (int i = 0; i < N_PID_CLASSES; i++) {
    std::function<bool(double, double)> comparison_operator;
    if (i == 0) { comparison_operator = std::greater_equal<double>(); // The first bin needs to include the lower limit.
    } else { comparison_operator = std::greater<double>(); }
    DetResponse response_good_track(DetResponse::customreco, TString::Format("track_response_gt_%.2f", PID_step * i), 40, 1, 100, 40, -1, 1, 1, 0, 1);
    response_good_track.SetObsFuncPtrs( &CUSTOMEF::HybridEnergy, &CUSTOMEF::TrackDir, &CUSTOMEF::TrackPos, &CUSTOMEF::TrackBY );

    response_good_track.AddCut( &SummaryEvent::Get_track_ql0       , std::greater<double>()   ,  0.5, true );
    response_good_track.AddCut( &SummaryEvent::Get_track_ql1       , std::greater<double>()   ,  0.5, true );
    response_good_track.AddCut( &SummaryEvent::Get_shower_ql0      , std::less<double>()      ,  0.5, true );
    response_good_track.AddCut( &SummaryEvent::Get_shower_ql1      , std::less<double>()      ,  0.5, true );
    response_good_track.AddCut( &SummaryEvent::Get_RDF_track_score , comparison_operator      , PID_step * i    , true );
    response_good_track.AddCut( &SummaryEvent::Get_RDF_track_score , std::less_equal<double>(), PID_step * (i+1), true );
    response_good_track.AddCut( &SummaryEvent::Get_RDF_muon_score  , std::less_equal<double>(), 0.05, true );
    response_good_track.AddCut( &SummaryEvent::Get_RDF_noise_score , std::less_equal<double>(), 0.18, true );
    response_good_track_vector.push_back(response_good_track);

    // Good showers only, no overlap with good tracks
    DetResponse response_good_shower(DetResponse::customreco, TString::Format("track_response_gs_%.2f", PID_step * i), 40, 1, 100, 40, -1, 1, 1, 0, 1);
    response_good_shower.SetObsFuncPtrs( &CUSTOMEF::HybridEnergy, &CUSTOMEF::TrackDir, &CUSTOMEF::TrackPos, &CUSTOMEF::TrackBY );

    response_good_shower.AddCut( &SummaryEvent::Get_track_ql0       , std::greater<double>()   ,  0.5, true ); // This would cut almost all events, since they all have fTrack_ql0 > 0.5
    response_good_shower.AddCut( &SummaryEvent::Get_track_ql1       , std::less<double>()      ,  0.5, true ); // Moreover: the numbers are consistent with resolution_plot_flav_complementary_events.C
    response_good_shower.AddCut( &SummaryEvent::Get_shower_ql0      , std::greater<double>()   ,  0.5, true );
    response_good_shower.AddCut( &SummaryEvent::Get_shower_ql1      , std::greater<double>()   ,  0.5, true );
    response_good_shower.AddCut( &SummaryEvent::Get_RDF_track_score , comparison_operator      , PID_step * i    , true );
    response_good_shower.AddCut( &SummaryEvent::Get_RDF_track_score , std::less_equal<double>(), PID_step * (i+1), true );
    response_good_shower.AddCut( &SummaryEvent::Get_RDF_muon_score  , std::less_equal<double>(), 0.05, true );
    response_good_shower.AddCut( &SummaryEvent::Get_RDF_noise_score , std::less_equal<double>(), 0.18, true );
    response_good_shower_vector.push_back(response_good_shower);

    // Good events only, no overlap with bad tracks or bad showers
    DetResponse response_good_event(DetResponse::customreco, TString::Format("track_response_ge_%.2f", PID_step * i), 40, 1, 100, 40, -1, 1, 1, 0, 1);
    response_good_event.SetObsFuncPtrs( &CUSTOMEF::HybridEnergy, &CUSTOMEF::TrackDir, &CUSTOMEF::TrackPos, &CUSTOMEF::TrackBY );

    response_good_event.AddCut( &SummaryEvent::Get_track_ql0       , std::greater<double>()   ,  0.5, true );
    response_good_event.AddCut( &SummaryEvent::Get_track_ql1       , std::greater<double>()   ,  0.5, true );
    response_good_event.AddCut( &SummaryEvent::Get_shower_ql0      , std::greater<double>()   ,  0.5, true );
    response_good_event.AddCut( &SummaryEvent::Get_shower_ql1      , std::greater<double>()   ,  0.5, true );
    response_good_event.AddCut( &SummaryEvent::Get_RDF_track_score , comparison_operator      , PID_step * i    , true );
    response_good_event.AddCut( &SummaryEvent::Get_RDF_track_score , std::less_equal<double>(), PID_step * (i+1), true );
    response_good_event.AddCut( &SummaryEvent::Get_RDF_muon_score  , std::less_equal<double>(), 0.05, true );
    response_good_event.AddCut( &SummaryEvent::Get_RDF_noise_score , std::less_equal<double>(), 0.18, true );
    response_good_event_vector.push_back(response_good_event);

    DetResponse response_shower(DetResponse::customreco, TString::Format("shower_response_%.2f", PID_step * i), 40, 1, 100, 40, -1, 1, 1, 0, 1);
    response_shower.SetObsFuncPtrs( &CUSTOMEF::HybridEnergy, &CUSTOMEF::TrackDir, &CUSTOMEF::TrackPos, &CUSTOMEF::TrackBY );

    response_shower.AddCut( &SummaryEvent::Get_shower_ql0     , std::greater<double>()   ,  0.5, true );
    response_shower.AddCut( &SummaryEvent::Get_shower_ql1     , std::greater<double>()   ,  0.5, true );
    response_shower.AddCut( &SummaryEvent::Get_RDF_track_score , comparison_operator      , PID_step * i    , true );
    response_shower.AddCut( &SummaryEvent::Get_RDF_track_score , std::less_equal<double>(), PID_step * (i+1), true );
    response_shower.AddCut( &SummaryEvent::Get_RDF_muon_score , std::less_equal<double>(), 0.05, true );
    response_shower.AddCut( &SummaryEvent::Get_RDF_noise_score, std::less_equal<double>(),  0.5, true );
    response_shower_vector.push_back(response_shower);

    DetResponse response_mc(DetResponse::customreco, TString::Format("mc_response_%.2f", PID_step * i), 40, 1, 100, 40, -1, 1, 1, 0, 1);
    response_mc.SetObsFuncPtrs( &CUSTOMEF::HybridEnergy, &CUSTOMEF::TrackDir, &CUSTOMEF::TrackPos, &CUSTOMEF::TrackBY );

    response_mc.AddCut( &SummaryEvent::Get_RDF_track_score , comparison_operator      , PID_step * i    , true );
    response_mc.AddCut( &SummaryEvent::Get_RDF_track_score , std::less_equal<double>(), PID_step * (i+1), true );
    response_mc.AddCut( &SummaryEvent::Get_RDF_muon_score , std::less_equal<double>(), 0.05, true );
    response_mc.AddCut( &SummaryEvent::Get_RDF_noise_score, std::less_equal<double>(),  0.5, true );
    response_mc_vector.push_back(response_mc);
  }

  auto summary_file = (TString)getenv("NMHDIR") + "/data/ORCA_MC_summary_all_10Apr2018.root";
  SummaryParser sp(summary_file);
  for (Int_t i = 0; i < sp.GetTree()->GetEntries(); i++) {
    if (i % (Int_t)1e6 == 0) cout << "Event: " << i << endl;
    SummaryEvent *evt = sp.GetEvt(i);
    for (int i = 0; i < N_PID_CLASSES; i++) {
      response_good_track_vector[i].Fill(evt);
      response_good_shower_vector[i].Fill(evt);
      response_good_event_vector[i].Fill(evt);
      response_shower_vector[i].Fill(evt);
      response_mc_vector[i].Fill(evt);
    }
  }

  for (int i = 0; i < N_PID_CLASSES; i++) { 
    response_good_track_vector[i].WriteToFile(filefolder + TString::Format("track_response_gt_%.2f.root", PID_step * i));
    response_good_shower_vector[i].WriteToFile(filefolder + TString::Format("track_response_gs_%.2f.root", PID_step * i));
    response_good_event_vector[i].WriteToFile(filefolder + TString::Format("track_response_ge_%.2f.root", PID_step * i));
    response_shower_vector[i].WriteToFile(filefolder  + TString::Format("shower_response_%.2f.root", PID_step * i));
    response_mc_vector[i].WriteToFile(filefolder  + TString::Format("mc_response_%.2f.root", PID_step * i));
  }

  cout << "NOTICE: Finished filling response" << endl;

  //----------------------------------------------------------
  // set up the PDFs and static oscillation parameters
  //----------------------------------------------------------

  auto meff_file = (TString)getenv("NMHDIR") + "/data/eff_mass/EffMass_ORCA115_23x9m_ECAP0418.root";
  FitUtil *fitutil = new FitUtil(3, response_good_track.GetHist3D(),
                 1, 100, -1, 1, 0, 1, meff_file);

  Double_t sinsqth12 = TMath::Power( TMath::Sin( 33.4 * TMath::Pi()/180. ), 2 );
  Double_t sinsqth13 = TMath::Power( TMath::Sin( 8.42 * TMath::Pi()/180. ), 2 );
  Double_t sinsqth23 = TMath::Power( TMath::Sin( 45   * TMath::Pi()/180. ), 2 );
  Double_t dcp       = 0.;
  Double_t dm32      = 2.44e-3;
  Double_t dm21      = 7.53e-5;
  Double_t DM        = dm32 + 0.5*dm21;

  // deconstrain th23 and dm31, when fitting you want constraints, otherwise you dont.
  ( (RooRealVar*)fitutil->GetSet().find("Dm31") )->setMin( -1 );
  ( (RooRealVar*)fitutil->GetSet().find("Dm31") )->setMax(  1 );
  ( (RooRealVar*)fitutil->GetSet().find("SinsqTh23") )->setMin( -1 );
  ( (RooRealVar*)fitutil->GetSet().find("SinsqTh23") )->setMax(  1 );

  // set parameter values 
  ( (RooRealVar*)fitutil->GetSet().find("SinsqTh12") )->setVal( sinsqth12 );
  ( (RooRealVar*)fitutil->GetSet().find("SinsqTh13") )->setVal( sinsqth13 );
  ( (RooRealVar*)fitutil->GetSet().find("dcp") )->setVal( dcp );
  ( (RooRealVar*)fitutil->GetSet().find("Dm21") )->setVal( dm21 );
  ( (RooRealVar*)fitutil->GetSet().find("SinsqTh23") )->setVal( sinsqth23 );


  //----------------------------------------------------------
  // set normal hierarchy
  //----------------------------------------------------------
  Double_t dm31 = DM + 0.5*dm21;
  ( (RooRealVar*)fitutil->GetSet().find("Dm31") )->setVal( dm31 );

  for (int i = 0; i < N_PID_CLASSES; i++) {
    FitPDF pdf_track_gt("pdf_tracks_gt", "pdf_tracks_gt", fitutil, &track_response_gt[i]);
    FitPDF pdf_tracks_gs("pdf_tracks_gs", "pdf_tracks_gs", fitutil, &track_response_gs[i]);
    FitPDF pdf_tracks_ge("pdf_tracks_ge", "pdf_tracks_ge", fitutil, &track_response_ge[i]);
    FitPDF pdf_showers("pdf_showers", "pdf_showers", fitutil, &shower_response[i]);
    FitPDF pdf_mc("pdf_mc", "pdf_mc", fitutil, &mc_response[i]);

    TH2D *tracks_gt_NO = (TH2D*)pdf_tracks_gt.GetExpValHist()->Project3D("yx");
    TH2D *tracks_gs_NO = (TH2D*)pdf_tracks_gs.GetExpValHist()->Project3D("yx");
    TH2D *tracks_ge_NO = (TH2D*)pdf_tracks_ge.GetExpValHist()->Project3D("yx");
    TH2D *showers_NO   = (TH2D*)pdf_showers.GetExpValHist()->Project3D("yx");
    TH2D *mc_NO        = (TH2D*)pdf_mc.GetExpValHist()->Project3D("yx");

    tracks_gt_NO->SetNameTitle("detected_tracks_gt", "detected_tracks_gt");
    tracks_gs_NO->SetNameTitle("detected_tracks_gs", "detected_tracks_gs");
    tracks_ge_NO->SetNameTitle("detected_tracks_ge", "detected_tracks_ge");
    showers_NO->SetNameTitle("detected_showers", "detected_showers");
    mc_NO->SetNameTitle("detected_mc", "detected_mc");
 
    //----------------------------------------------------------
    // set inverted hierarchy
    //----------------------------------------------------------
    Double_t dm31 = -DM + 0.5*dm21;
    TH2D *tracks_gt_IO = (TH2D*)pdf_tracks_gt.GetExpValHist()->Project3D("yx");
    TH2D *tracks_gs_IO = (TH2D*)pdf_tracks_gs.GetExpValHist()->Project3D("yx");
    TH2D *tracks_ge_IO = (TH2D*)pdf_tracks_ge.GetExpValHist()->Project3D("yx");
    TH2D *showers_IO   = (TH2D*)pdf_showers.GetExpValHist()->Project3D("yx");
    TH2D *mc_IO        = (TH2D*)pdf_mc.GetExpValHist()->Project3D("yx");

    tracks_gt_IO->SetNameTitle("detected_tracks_gt", "detected_tracks_gt");
    tracks_gs_IO->SetNameTitle("detected_tracks_gs", "detected_tracks_gs");
    tracks_ge_IO->SetNameTitle("detected_tracks_ge", "detected_tracks_ge");
    showers_IO->SetNameTitle("detected_showers", "detected_showers");
    mc_IO->SetNameTitle("detected_mc", "detected_mc");

    //----------------------------------------------------------
    // save output
    //----------------------------------------------------------
    TString output_NO = TString::Format("split_expected_evts_NO_%.2f.root", N_PID_CLASSES * i);
    TFile fout_NO(filefolder + output_NO,"RECREATE");
    tracks_gt_NO->Write();
    tracks_gs_NO->Write();
    tracks_ge_NO->Write();
    showers_NO->Write();
    mc_NO->Write();
    fout_NO.Close();

    TString output_IO = TString::Format("split_expected_evts_IO_%.2f.root", N_PID_CLASSES * i);
    TFile fout_IO(filefolder + output_IO,"RECREATE");
    tracks_gt_IO->Write();
    tracks_gs_IO->Write();
    tracks_ge_IO->Write();
    showers_IO->Write();
    mc_IO->Write();
  }
}
