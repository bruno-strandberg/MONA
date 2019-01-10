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

void SplitDetectorResponse() {

  const int N_PID_CLASSES = 10;
  Double_t PID_step = 1 / float(N_PID_CLASSES);
  TString filefolder = TString::Format("./pid_detres/RooFit/pid_binning_%i/", N_PID_CLASSES);

  gSystem->Load("$OSCPROBDIR/libOscProb.so");

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

    DetResponse track_response(DetResponse::track, TString::Format("track_response_%.2f", PID_step * i), 40, 1, 100, 40, -1, 1, 1, 0, 1);
    track_response.AddCut( &SummaryEvent::Get_track_ql0       , std::greater<double>()   , 0.5             , true );
    track_response.AddCut( &SummaryEvent::Get_track_ql1       , std::greater<double>()   , 0.5             , true );
    track_response.AddCut( &SummaryEvent::Get_RDF_track_score , comparison_operator      , PID_step * i    , true );
    track_response.AddCut( &SummaryEvent::Get_RDF_track_score , std::less_equal<double>(), PID_step * (i+1), true );
    track_response.AddCut( &SummaryEvent::Get_RDF_muon_score  , std::less_equal<double>(), 0.05            , true );
    track_response.AddCut( &SummaryEvent::Get_RDF_noise_score , std::less_equal<double>(), 0.18            , true );
    track_response_vector.push_back(track_response);

    DetResponse shower_response(DetResponse::shower, TString::Format("shower_response_%.2f", PID_step * i), 40, 1, 100, 40, -1, 1, 1, 0, 1);
    shower_response.AddCut( &SummaryEvent::Get_shower_ql0     , std::greater<double>()   , 0.5             , true );
    shower_response.AddCut( &SummaryEvent::Get_shower_ql1     , std::greater<double>()   , 0.5             , true );
    shower_response.AddCut( &SummaryEvent::Get_RDF_track_score, comparison_operator      , PID_step * i    , true );
    shower_response.AddCut( &SummaryEvent::Get_RDF_track_score, std::less_equal<double>(), PID_step * (i+1), true );
    shower_response.AddCut( &SummaryEvent::Get_RDF_muon_score , std::less_equal<double>(), 0.05            , true );
    shower_response.AddCut( &SummaryEvent::Get_RDF_noise_score, std::less_equal<double>(), 0.5             , true );
    shower_response_vector.push_back(shower_response);

    DetResponse mc_response(DetResponse::mc_truth, TString::Format("mc_response_%.2f", PID_step * i), 40, 1, 100, 40, -1, 1, 1, 0, 1);
    mc_response.AddCut( &SummaryEvent::Get_RDF_track_score, comparison_operator      , PID_step * i    , true );
    mc_response.AddCut( &SummaryEvent::Get_RDF_track_score, std::less_equal<double>(), PID_step * (i+1), true );
    mc_response.AddCut( &SummaryEvent::Get_RDF_muon_score , std::less_equal<double>(), 0.05            , true );
    mc_response.AddCut( &SummaryEvent::Get_RDF_noise_score, std::less_equal<double>(), 0.5             , true ); // Do these cuts make sense?? Yes, we want to compare the DR with "real DR"
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
    track_response_vector[i].WriteToFile(filefolder + TString::Format("track_response_%.2f.root" , PID_step * i));       
    shower_response_vector[i].WriteToFile(filefolder + TString::Format("shower_response_%.2f.root", PID_step * i));       
    mc_response_vector[i].WriteToFile(filefolder + TString::Format("mc_response_%.2f.root", PID_step * i));     
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

  for (int i = 0; i < N_PID_CLASSES; i++){
    FitPDF pdf_tracks("pdf_tracks", "pdf_tracks"   , fitutil, &track_response_vector[i]);
    FitPDF pdf_showers("pdf_showers", "pdf_showers", fitutil, &shower_response_vector[i]);
    FitPDF pdf_mc("pdf_mc", "pdf_mc", fitutil, &mc_response_vector[i]);
  
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
    TString output_NO = TString::Format("split_expectated_evts_NO_%.2f.root", PID_step * i);
    TFile fout_NO(filefolder + output_NO,"RECREATE");
    tracks_NO->Write();
    showers_NO->Write();
    mc_NO->Write();
    fout_NO.Close();

    TString output_IO = TString::Format("split_expectated_evts_IO_%.2f.root", PID_step * i);
    TFile fout_IO(filefolder + output_IO,"RECREATE");
    tracks_IO->Write();
    showers_IO->Write();
    mc_IO->Write();
  }
}
