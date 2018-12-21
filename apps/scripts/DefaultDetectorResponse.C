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

void DefaultDetectorResponse() {

  TString filefolder = "./default_detres/RooFit/";

  gROOT->ProcessLine(".L FitFunction.C+");
  gSystem->Load("$OSCPROBDIR/libOscProb.so");

  //----------------------------------------------------------
  // detector response for tracks and showers
  //----------------------------------------------------------
  DetResponse track_response(DetResponse::track, "track_response", 40, 1, 100, 40, -1, 1, 1, 0, 1);
  track_response.AddCut( &SummaryEvent::Get_track_ql0       , std::greater<double>()   ,  0.5, true );
  track_response.AddCut( &SummaryEvent::Get_track_ql1       , std::greater<double>()   ,  0.5, true );
  track_response.AddCut( &SummaryEvent::Get_RDF_track_score , std::greater<double>()   ,  0.6, true );
  track_response.AddCut( &SummaryEvent::Get_RDF_muon_score  , std::less_equal<double>(), 0.05, true );
  track_response.AddCut( &SummaryEvent::Get_RDF_noise_score , std::less_equal<double>(), 0.18, true );

  DetResponse shower_response(DetResponse::shower, "shower_response", 40, 1, 100, 40, -1, 1, 1, 0, 1);
  shower_response.AddCut( &SummaryEvent::Get_shower_ql0     , std::greater<double>()   ,  0.5, true );
  shower_response.AddCut( &SummaryEvent::Get_shower_ql1     , std::greater<double>()   ,  0.5, true );
  shower_response.AddCut( &SummaryEvent::Get_RDF_track_score, std::less_equal<double>(),  0.6, true );
  shower_response.AddCut( &SummaryEvent::Get_RDF_muon_score , std::less_equal<double>(), 0.05, true );
  shower_response.AddCut( &SummaryEvent::Get_RDF_noise_score, std::less_equal<double>(),  0.5, true );

  DetResponse mc_response(DetResponse::mc_truth, "mc_response", 40, 1, 100, 40, -1, 1, 1, 0, 1); // Are these cuts necessary in MC?
  mc_response.AddCut( &SummaryEvent::Get_RDF_muon_score , std::less_equal<double>(), 0.05, true );
  mc_response.AddCut( &SummaryEvent::Get_RDF_noise_score, std::less_equal<double>(),  0.5, true );

  auto summary_file = (TString)getenv("NMHDIR") + "/data/ORCA_MC_summary_all_10Apr2018.root";
  SummaryParser sp(summary_file);
  for (Int_t i = 0; i < sp.GetTree()->GetEntries(); i++) {
    if (i % (Int_t)1e6 == 0) cout << "Event: " << i << endl;
    sp.GetTree()->GetEntry(i);
    SummaryEvent *evt = sp.GetEvt();
    if (evt->Get_MC_is_neutrino() < 0.5) continue; // Is this still relevant?
    track_response.Fill(evt);
    shower_response.Fill(evt);
    mc_response.Fill(evt);
  }

  track_response.WriteToFile(filefolder + "track_response.root");
  shower_response.WriteToFile(filefolder + "shower_response.root");
  mc_response.WriteToFile(filefolder + "mc_response.root");

  cout << "NOTICE: Finished filling response" << endl;

  //----------------------------------------------------------
  // set up the PDFs and static oscillation parameters
  //----------------------------------------------------------

  auto meff_file = (TString)getenv("NMHDIR") + "/data/eff_mass/EffMass_ORCA115_23x9m_ECAP0418.root";
  FitUtil *fitutil = new FitUtil(3, track_response.GetHist3D(),
                 1, 100, -1, 1, 0, 1, meff_file);

  FitPDF pdf_tracks("pdf_tracks", "pdf_tracks"   , fitutil, &track_response);
  FitPDF pdf_showers("pdf_showers", "pdf_showers", fitutil, &shower_response);
  FitPDF pdf_mc("pdf_mc", "pdf_mc", fitutil, &mc_response);

  Double_t sinsqth12 = TMath::Power( TMath::Sin( 33.4 * TMath::Pi()/180. ), 2 );
  Double_t sinsqth13 = TMath::Power( TMath::Sin( 8.42 * TMath::Pi()/180. ), 2 );
  Double_t sinsqth23 = TMath::Power( TMath::Sin( 45   * TMath::Pi()/180. ), 2 );
  Double_t dcp       = 0.;
  Double_t dm32      = 2.44e-3;
  Double_t dm21      = 7.53e-5;
  Double_t DM        = dm32 + 0.5*dm21;

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

  TH2D *tracks_NO  = (TH2D*)pdf_tracks.GetExpValHist()->Project3D("yx");
  TH2D *showers_NO = (TH2D*)pdf_showers.GetExpValHist()->Project3D("yx");
  TH2D *mc_NO      = (TH2D*)pdf_mc.GetExpValHist()->Project3D("yx");

  tracks_NO->SetNameTitle("detected_tracks", "detected_tracks");
  showers_NO->SetNameTitle("detected_showers", "detected_showers");
  mc_NO->SetNameTitle("detected_mc", "detected_mc");

  //----------------------------------------------------------
  // set inverted hierarchy
  //----------------------------------------------------------
  dm31 = -DM + 0.5*dm21; //IO
  ( (RooRealVar*)fitutil->GetSet().find("Dm31") )->setVal( dm31 );

  TH2D *tracks_IO  = (TH2D*)pdf_tracks.GetExpValHist()->Project3D("yx");
  TH2D *showers_IO = (TH2D*)pdf_showers.GetExpValHist()->Project3D("yx");
  TH2D *mc_IO      = (TH2D*)pdf_mc.GetExpValHist()->Project3D("yx");

  tracks_IO->SetNameTitle("detected_tracks", "detected_tracks");
  showers_IO->SetNameTitle("detected_showers", "detected_showers");
  mc_IO->SetNameTitle("detected_mc", "detected_mc");

  //----------------------------------------------------------
  // save output
  //----------------------------------------------------------
  TString output_NO = "default_expectated_evts_NO.root";
  TFile fout_NO(filefolder + output_NO,"RECREATE");
  tracks_NO->Write();
  showers_NO->Write();
  mc_NO->Write();
  fout_NO.Close();

  TString output_IO = "default_expectated_evts_IO.root";
  TFile fout_IO(filefolder + output_IO,"RECREATE");
  tracks_IO->Write();
  showers_IO->Write();
  mc_IO->Write();
  fout_IO.Close();
}
