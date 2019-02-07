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

/* Macro to generate the detector responses in the default scheme: a track response and a shower
 * response. The shower response are events with a quality factor of < 0.6, tracks are events with
 * a quality factor of > 0.6.
 *
 * The responses and all other root files are saved into `filefolder`.
 */

void DetectorResponseDefault() {

  TString filefolder = "./default_detres/";

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

  DetResponse mc_response(DetResponse::mc_truth, "mc_response", 40, 1, 100, 40, -1, 1, 1, 0, 1);
  mc_response.AddCut( &SummaryEvent::Get_RDF_muon_score , std::less_equal<double>(), 0.05, true );
  mc_response.AddCut( &SummaryEvent::Get_RDF_noise_score, std::less_equal<double>(),  0.5, true );

  auto summary_file = (TString)getenv("NMHDIR") + "/data/ORCA_MC_summary_all_10Apr2018.root";
  SummaryParser sp(summary_file);
  for (Int_t i = 0; i < sp.GetTree()->GetEntries(); i++) {
    if (i % (Int_t)1e6 == 0) cout << "Event: " << i << endl;
    SummaryEvent *evt = sp.GetEvt(i);
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

  TH2D* tracks_err_NO  = (TH2D*)pdf_tracks.GetExpValErrHist()->Project3D("yx");
  TH2D* showers_err_NO = (TH2D*)pdf_showers.GetExpValErrHist()->Project3D("yx");
  TH2D* mc_err_NO      = (TH2D*)pdf_mc.GetExpValErrHist()->Project3D("yx");
  //----------------------------------------------------------
  // set inverted hierarchy
  //----------------------------------------------------------
  dm31 = -DM + 0.5*dm21; //IO
  fitutil->GetVar("Dm31")->setVal( dm31 );

  TH2D *tracks_IO  = (TH2D*)pdf_tracks.GetExpValHist()->Project3D("yx");
  TH2D *showers_IO = (TH2D*)pdf_showers.GetExpValHist()->Project3D("yx");
  TH2D *mc_IO      = (TH2D*)pdf_mc.GetExpValHist()->Project3D("yx");

  tracks_IO->SetNameTitle("detected_tracks", "detected_tracks");
  showers_IO->SetNameTitle("detected_showers", "detected_showers");
  mc_IO->SetNameTitle("detected_mc", "detected_mc");

  TH2D* tracks_err_IO  = (TH2D*)pdf_tracks.GetExpValErrHist()->Project3D("yx");
  TH2D* showers_err_IO = (TH2D*)pdf_showers.GetExpValErrHist()->Project3D("yx");
  TH2D* mc_err_IO      = (TH2D*)pdf_mc.GetExpValErrHist()->Project3D("yx");
  //----------------------------------------------------------
  // save output
  //----------------------------------------------------------
  TString output_NO = "default_expected_evts_NO.root";
  TFile fout_NO(filefolder + output_NO,"RECREATE");
  auto hists_NO = {tracks_NO, showers_NO, mc_NO,
                   tracks_err_NO, showers_err_NO, mc_err_NO};
  for (auto hist: hists_NO) {
    hist->Write();
  }
  fout_NO.Close();

  TString output_IO = "default_expected_evts_IO.root";
  TFile fout_IO(filefolder + output_IO,"RECREATE");
  auto hists_IO = {tracks_IO, showers_IO, mc_IO,
                   tracks_err_IO, showers_err_IO, mc_err_IO};
  for (auto hist: hists_IO) {
    hist->Write();
  }
  fout_IO.Close();
}
