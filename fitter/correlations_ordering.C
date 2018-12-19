#include "TSystem.h"
#include "TROOT.h"
#include "TDatime.h"
#include "TH3.h"
#include "TFile.h"
#include "TStopwatch.h"

#include "DetResponse.h"
#include "FitFunction.h"
#include "SummaryParser.h"
#include "SummaryEvent.h"
#include "HelperFunctions.C"

#include <iostream>
using namespace std;

void correlations_ordering() {
  const int N_PID_CLASSES = 10;
  Double_t PID_step = 1 / float(N_PID_CLASSES);
  TString filefolder = "./pid_detres/";

  gROOT->ProcessLine(".L FitFunction.C+");
  gSystem->Load("$OSCPROBDIR/libOscProb.so");

  //----------------------------------------------------------
  // detector response for tracks and showers
  //----------------------------------------------------------
  std::vector<DetResponse> drts_NO;
  std::vector<DetResponse> drts_IO;
  std::vector<DetResponse> drss_NO;
  std::vector<DetResponse> drss_IO;
  std::vector<DetResponse> drms_NO;
  std::vector<DetResponse> drms_IO;

  for (int i = 0; i < N_PID_CLASSES; i++) {
    DetResponse drt_NO(DetResponse::track, TString::Format("track_response_NO_%.2f", PID_step * i), 40, 1, 100, 40, -1, 1, 1, 0, 1);
    drt_NO.AddCut( &SummaryEvent::Get_track_ql0       , std::greater<double>()   , 0.5             , true );
    drt_NO.AddCut( &SummaryEvent::Get_track_ql1       , std::greater<double>()   , 0.5             , true );
    drt_NO.AddCut( &SummaryEvent::Get_RDF_track_score , std::greater<double>()   , PID_step * i    , true );
    drt_NO.AddCut( &SummaryEvent::Get_RDF_track_score , std::less_equal<double>(), PID_step * (i+1), true );
    drt_NO.AddCut( &SummaryEvent::Get_RDF_muon_score  , std::less_equal<double>(), 0.05            , true );
    drt_NO.AddCut( &SummaryEvent::Get_RDF_noise_score , std::less_equal<double>(), 0.18            , true );
    drts_NO.push_back(drt_NO);

    DetResponse drs_NO(DetResponse::shower, TString::Format("shower_response_NO_%.2f", PID_step * i), 40, 1, 100, 40, -1, 1, 1, 0, 1);
    drs_NO.AddCut( &SummaryEvent::Get_shower_ql0     , std::greater<double>()   , 0.5             , true );
    drs_NO.AddCut( &SummaryEvent::Get_shower_ql1     , std::greater<double>()   , 0.5             , true );
    drs_NO.AddCut( &SummaryEvent::Get_RDF_track_score, std::greater<double>()   , PID_step * i    , true );
    drs_NO.AddCut( &SummaryEvent::Get_RDF_track_score, std::less_equal<double>(), PID_step * (i+1), true );
    drs_NO.AddCut( &SummaryEvent::Get_RDF_muon_score , std::less_equal<double>(), 0.05            , true );
    drs_NO.AddCut( &SummaryEvent::Get_RDF_noise_score, std::less_equal<double>(), 0.5             , true );
    drss_NO.push_back(drs_NO);

    DetResponse drm_NO(DetResponse::mc_truth, TString::Format("mc_response_NO_%.2f", PID_step * i), 40, 1, 100, 40, -1, 1, 1, 0, 1);
    drm_NO.AddCut( &SummaryEvent::Get_RDF_track_score, std::greater<double>()   , PID_step * i    , true );
    drm_NO.AddCut( &SummaryEvent::Get_RDF_track_score, std::less_equal<double>(), PID_step * (i+1), true );
    drm_NO.AddCut( &SummaryEvent::Get_RDF_muon_score , std::less_equal<double>(), 0.05            , true );
    drm_NO.AddCut( &SummaryEvent::Get_RDF_noise_score, std::less_equal<double>(), 0.5             , true );
    drms_NO.push_back(drm_NO);

    DetResponse drt_IO(DetResponse::track, TString::Format("track_response_IO_%.2f", PID_step * i), 40, 1, 100, 40, -1, 1, 1, 0, 1);
    drt_IO.AddCut( &SummaryEvent::Get_track_ql0       , std::greater<double>()   , 0.5             , true );
    drt_IO.AddCut( &SummaryEvent::Get_track_ql1       , std::greater<double>()   , 0.5             , true );
    drt_IO.AddCut( &SummaryEvent::Get_RDF_track_score , std::greater<double>()   , PID_step * i    , true );
    drt_IO.AddCut( &SummaryEvent::Get_RDF_track_score , std::less_equal<double>(), PID_step * (i+1), true );
    drt_IO.AddCut( &SummaryEvent::Get_RDF_muon_score  , std::less_equal<double>(), 0.05            , true );
    drt_IO.AddCut( &SummaryEvent::Get_RDF_noise_score , std::less_equal<double>(), 0.18            , true );
    drts_IO.push_back(drt_IO);

    DetResponse drs_IO(DetResponse::shower, TString::Format("shower_response_IO_%.2f", PID_step * i), 40, 1, 100, 40, -1, 1, 1, 0, 1);
    drs_IO.AddCut( &SummaryEvent::Get_shower_ql0     , std::greater<double>()   , 0.5             , true );
    drs_IO.AddCut( &SummaryEvent::Get_shower_ql1     , std::greater<double>()   , 0.5             , true );
    drs_IO.AddCut( &SummaryEvent::Get_RDF_track_score, std::greater<double>()   , PID_step * i    , true );
    drs_IO.AddCut( &SummaryEvent::Get_RDF_track_score, std::less_equal<double>(), PID_step * (i+1), true );
    drs_IO.AddCut( &SummaryEvent::Get_RDF_muon_score , std::less_equal<double>(), 0.05            , true );
    drs_IO.AddCut( &SummaryEvent::Get_RDF_noise_score, std::less_equal<double>(), 0.5             , true );
    drss_IO.push_back(drs_IO);
    
    DetResponse drm_IO(DetResponse::mc_truth, TString::Format("mc_response_IO_%.2f", PID_step * i), 40, 1, 100, 40, -1, 1, 1, 0, 1);
    drm_IO.AddCut( &SummaryEvent::Get_RDF_track_score, std::greater<double>()   , PID_step * i    , true );
    drm_IO.AddCut( &SummaryEvent::Get_RDF_track_score, std::less_equal<double>(), PID_step * (i+1), true );
    drm_IO.AddCut( &SummaryEvent::Get_RDF_muon_score , std::less_equal<double>(), 0.05            , true );
    drm_IO.AddCut( &SummaryEvent::Get_RDF_noise_score, std::less_equal<double>(), 0.5             , true );
    drms_IO.push_back(drm_IO);

  }

  for (int i = 0; i < N_PID_CLASSES; i++) {
    drts_NO[i].ReadFromFile(filefolder + TString::Format("track_response_NO_%.2f.root" , PID_step * i));
    drts_IO[i].ReadFromFile(filefolder + TString::Format("track_response_IO_%.2f.root" , PID_step * i));
    drss_NO[i].ReadFromFile(filefolder + TString::Format("shower_response_NO_%.2f.root", PID_step * i));
    drss_IO[i].ReadFromFile(filefolder + TString::Format("shower_response_IO_%.2f.root", PID_step * i));
    drms_NO[i].ReadFromFile(filefolder + TString::Format("mc_response_NO_%.2f.root", PID_step * i));
    drms_IO[i].ReadFromFile(filefolder + TString::Format("mc_response_IO_%.2f.root", PID_step * i));
  }  

  Double_t bybin = 0.5; // Theres only one by

  for (int i = 0; i < N_PID_CLASSES; i++){
    FitFunction tfitf_NO(&drts_NO[i], 3, "../data/eff_mass/EffMhists_elec_CC.root", "../data/eff_mass/EffMhists_muon_CC.root", "../data/eff_mass/EffMhists_tau_CC.root", "../data/eff_mass/EffMhists_elec_NC.root");
    FitFunction tfitf_IO(&drts_IO[i], 3, "../data/eff_mass/EffMhists_elec_CC.root", "../data/eff_mass/EffMhists_muon_CC.root", "../data/eff_mass/EffMhists_tau_CC.root", "../data/eff_mass/EffMhists_elec_NC.root");
    FitFunction sfitf_NO(&drss_NO[i], 3, "../data/eff_mass/EffMhists_elec_CC.root", "../data/eff_mass/EffMhists_muon_CC.root", "../data/eff_mass/EffMhists_tau_CC.root", "../data/eff_mass/EffMhists_elec_NC.root");
    FitFunction sfitf_IO(&drss_IO[i], 3, "../data/eff_mass/EffMhists_elec_CC.root", "../data/eff_mass/EffMhists_muon_CC.root", "../data/eff_mass/EffMhists_tau_CC.root", "../data/eff_mass/EffMhists_elec_NC.root");
    FitFunction mfitf_NO(&drms_NO[i], 3, "../data/eff_mass/EffMhists_elec_CC.root", "../data/eff_mass/EffMhists_muon_CC.root", "../data/eff_mass/EffMhists_tau_CC.root", "../data/eff_mass/EffMhists_elec_NC.root");
    FitFunction mfitf_IO(&drms_IO[i], 3, "../data/eff_mass/EffMhists_elec_CC.root", "../data/eff_mass/EffMhists_muon_CC.root", "../data/eff_mass/EffMhists_tau_CC.root", "../data/eff_mass/EffMhists_elec_NC.root");

    TH3D *h_cor_tracks  = (TH3D*)drts_NO[i].GetHist3D()->Clone("correlation_track_bins");
    TH3D *h_cor_showers = (TH3D*)drss_NO[i].GetHist3D()->Clone("correlation_shower_bins");
    TH3D *h_cor_mc      = (TH3D*)drms_NO[i].GetHist3D()->Clone("correlation_mc_bins");
    h_cor_tracks->Reset();
    h_cor_showers->Reset();
    h_cor_mc->Reset();

    for (Int_t ebin = 1; ebin <= h_cor_tracks->GetXaxis()->GetNbins(); ebin++) {
      for (Int_t ctbin = 1; ctbin <= h_cor_tracks->GetYaxis()->GetNbins(); ctbin++) {
        for (Int_t bybin = 1; bybin <= h_cor_tracks->GetZaxis()->GetNbins(); bybin++) {
  
          Double_t E  = h_cor_tracks->GetXaxis()->GetBinCenter( ebin );
          Double_t ct = h_cor_tracks->GetYaxis()->GetBinCenter( ctbin );
          Double_t by = h_cor_tracks->GetZaxis()->GetBinCenter( bybin );

          double x[] = {E, ct, by};
  
          std::vector<Double_t> track_weights_NO = tfitf_NO.Weights(x);
          std::vector<Double_t> track_weights_IO = tfitf_IO.Weights(x);
          std::vector<Double_t> shower_weights_NO = sfitf_NO.Weights(x);
          std::vector<Double_t> shower_weights_IO = sfitf_IO.Weights(x);
          std::vector<Double_t> mc_weights_NO = mfitf_NO.Weights(x);
          std::vector<Double_t> mc_weights_IO = mfitf_IO.Weights(x);

          Double_t correlation_track = normed_innerproduct(track_weights_NO, track_weights_IO);
          Double_t correlation_shower = normed_innerproduct(shower_weights_NO, shower_weights_IO);
          Double_t correlation_mc = normed_innerproduct(mc_weights_NO, mc_weights_IO);

          h_cor_tracks->SetBinContent(  ebin, ctbin, bybin, correlation_track );
          h_cor_showers->SetBinContent( ebin, ctbin, bybin, correlation_shower );
          h_cor_mc->SetBinContent( ebin, ctbin, bybin, correlation_shower );

        }
      }
    }
    TString output = TString::Format("correlations_%.2f.root", PID_step * i);
    TFile fout(filefolder + output, "RECREATE");
    h_cor_tracks->Write();
    h_cor_showers->Write();
    h_cor_mc->Write();
    fout.Close();
  }
}
