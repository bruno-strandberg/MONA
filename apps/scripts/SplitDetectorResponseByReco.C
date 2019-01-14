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

void SplitDetectorResponseByReco() {

  TString filefolder = "./energy_detres/";

  gSystem->Load("$OSCPROBDIR/libOscProb.so");

  //----------------------------------------------------------
  // detector response for tracks and showers
  //----------------------------------------------------------
  // Good tracks only, no overlap with good showers
  // gt = good track, gs = good shower, ge = good event
  // The track response has been split into these 3 types, shower has not.
  DetResponse track_response_gt(DetResponse::hybridE, "track_track_response_gt", 40, 1, 100, 40, -1, 1, 1, 0, 1);
  track_response_gt.AddCut( &SummaryEvent::Get_track_ql0       , std::greater<double>()   ,  0.5, true );
  track_response_gt.AddCut( &SummaryEvent::Get_track_ql1       , std::greater<double>()   ,  0.5, true );
  track_response_gt.AddCut( &SummaryEvent::Get_shower_ql0      , std::less<double>()      ,  0.5, true );
  track_response_gt.AddCut( &SummaryEvent::Get_shower_ql1      , std::less<double>()      ,  0.5, true );
  track_response_gt.AddCut( &SummaryEvent::Get_RDF_track_score , std::greater<double>()   ,  0.6, true );
  track_response_gt.AddCut( &SummaryEvent::Get_RDF_muon_score  , std::less_equal<double>(), 0.05, true );
  track_response_gt.AddCut( &SummaryEvent::Get_RDF_noise_score , std::less_equal<double>(), 0.18, true );

  // Good showers only, no overlap with good tracks
  DetResponse track_response_gs(DetResponse::hybridE, "track_track_response_gs", 40, 1, 100, 40, -1, 1, 1, 0, 1);
  track_response_gs.AddCut( &SummaryEvent::Get_track_ql0       , std::greater<double>()   ,  0.5, true ); // This would cut almost all events, since they all have fTrack_ql0 > 0.5
  track_response_gs.AddCut( &SummaryEvent::Get_track_ql1       , std::less<double>()      ,  0.5, true ); // Moreover: the numbers are consistent with resolution_plot_flav_complementary_events.C
  track_response_gs.AddCut( &SummaryEvent::Get_shower_ql0      , std::greater<double>()   ,  0.5, true );
  track_response_gs.AddCut( &SummaryEvent::Get_shower_ql1      , std::greater<double>()   ,  0.5, true );
  track_response_gs.AddCut( &SummaryEvent::Get_RDF_track_score , std::greater<double>()   ,  0.6, true );
  track_response_gs.AddCut( &SummaryEvent::Get_RDF_muon_score  , std::less_equal<double>(), 0.05, true );
  track_response_gs.AddCut( &SummaryEvent::Get_RDF_noise_score , std::less_equal<double>(), 0.18, true );

  // Good events only, no overlap with bad tracks or bad showers
  DetResponse track_response_ge(DetResponse::hybridE, "track_track_response_ge", 40, 1, 100, 40, -1, 1, 1, 0, 1);
  track_response_ge.AddCut( &SummaryEvent::Get_track_ql0       , std::greater<double>()   ,  0.5, true );
  track_response_ge.AddCut( &SummaryEvent::Get_track_ql1       , std::greater<double>()   ,  0.5, true );
  track_response_ge.AddCut( &SummaryEvent::Get_shower_ql0      , std::greater<double>()   ,  0.5, true );
  track_response_ge.AddCut( &SummaryEvent::Get_shower_ql1      , std::greater<double>()   ,  0.5, true );
  track_response_ge.AddCut( &SummaryEvent::Get_RDF_track_score , std::greater<double>()   ,  0.6, true );
  track_response_ge.AddCut( &SummaryEvent::Get_RDF_muon_score  , std::less_equal<double>(), 0.05, true );
  track_response_ge.AddCut( &SummaryEvent::Get_RDF_noise_score , std::less_equal<double>(), 0.18, true );

  DetResponse shower_response(DetResponse::hybridE, "shower_response", 40, 1, 100, 40, -1, 1, 1, 0, 1);
  shower_response.AddCut( &SummaryEvent::Get_shower_ql0     , std::greater<double>()   ,  0.5, true );
  shower_response.AddCut( &SummaryEvent::Get_shower_ql1     , std::greater<double>()   ,  0.5, true );
  shower_response.AddCut( &SummaryEvent::Get_RDF_track_score, std::less_equal<double>(),  0.6, true );
  shower_response.AddCut( &SummaryEvent::Get_RDF_muon_score , std::less_equal<double>(), 0.05, true );
  shower_response.AddCut( &SummaryEvent::Get_RDF_noise_score, std::less_equal<double>(),  0.5, true );

  DetResponse mc_response(DetResponse::mc_truth, "mc_response", 40, 1, 100, 40, -1, 1, 1, 0, 1);
  mc_response.AddCut( &SummaryEvent::Get_RDF_muon_score , std::less_equal<double>(), 0.05, true );
  mc_response.AddCut( &SummaryEvent::Get_RDF_noise_score, std::less_equal<double>(),  0.5, true );

  auto summary_file = (TString)getenv("NMHDIR") + "/data/ORCA_MC_summary_all_10Apr2018.root";
  SummaryParser sp(summary_file);
  for (Int_t i = 0; i < sp.GetTree()->GetEntries(); i++) {
    if (i % (Int_t)1e6 == 0) cout << "Event: " << i << endl;
    SummaryEvent *evt = sp.GetEvt(i);
    track_response_gt.Fill(evt);
    track_response_gs.Fill(evt);
    track_response_ge.Fill(evt);
    shower_response.Fill(evt);
    response_mc.Fill(evt);
  }

  track_response_gt.WriteToFile(filefolder + "track_response_gt.root");
  track_response_gs.WriteToFile(filefolder + "track_response_gs.root");
  track_response_ge.WriteToFile(filefolder + "track_response_ge.root");
  shower_response.WriteToFile(filefolder  + "shower_response.root");
  mc_response.WriteToFile(filefolder  + "mc_response.root");

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

  FitPDF pdf_track_gt("pdf_tracks_gt", "pdf_tracks_gt", fitutil, &track_response_gt);
  FitPDF pdf_tracks_gs("pdf_tracks_gs", "pdf_tracks_gs", fitutil, &track_response_gs);
  FitPDF pdf_tracks_ge("pdf_tracks_ge", "pdf_tracks_ge", fitutil, &track_response_ge);
  FitPDF pdf_showers("pdf_showers", "pdf_showers", fitutil, &shower_response);
  FitPDF pdf_mc("pdf_mc", "pdf_mc", fitutil, &mc_response);

  //----------------------------------------------------------
  // set normal hierarchy
  //----------------------------------------------------------
  Double_t dm31 = DM + 0.5*dm21;
  ( (RooRealVar*)fitutil->GetSet().find("Dm31") )->setVal( dm31 );

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
  ( (RooRealVar*)fitutil->GetSet().find("Dm31") )->setVal( dm31 );

  TH2D *tracks_gt_IO = (TH2D*)pdf_tracks_gt.GetExpValHist()->Project3D("yx");
  TH2D *tracks_gs_IO = (TH2D*)pdf_tracks_gs.GetExpValHist()->Project3D("yx");
  TH2D *tracks_ge_IO = (TH2D*)pdf_tracks_ge.GetExpValHist()->Project3D("yx");
  TH2D *showers_IO   = (TH2D*)pdf_showers.GetExpValHist()->Project3D("yx");
  TH2D *mc_IO        = (TH2D*)pdf_mc.GetExpValHist()->Project3D("yx");

  tracks_gt_NO->SetNameTitle("detected_tracks_gt", "detected_tracks_gt");
  tracks_gs_NO->SetNameTitle("detected_tracks_gs", "detected_tracks_gs");
  tracks_ge_NO->SetNameTitle("detected_tracks_ge", "detected_tracks_ge");
  showers_NO->SetNameTitle("detected_showers", "detected_showers");
  mc_NO->SetNameTitle("detected_mc", "detected_mc");

  //----------------------------------------------------------
  // save output
  //----------------------------------------------------------
  TString output_NO = "split_expected_evts_NO.root";
  TFile fout_NO(filefolder + output_NO,"RECREATE");
  tracks_gt_NO->Write();
  tracks_gs_NO->Write();
  tracks_ge_NO->Write();
  showers_NO->Write();
  mc_NO->Write();
  fout_NO.Close();

  TString output_IO = "split_expected_evts_IO.root";
  TFile fout_IO(filefolder + output_IO,"RECREATE");
  tracks_gt_IO->Write();
  tracks_gs_IO->Write();
  tracks_ge_IO->Write();
  showers_IO->Write();
  mc_IO->Write();
}
