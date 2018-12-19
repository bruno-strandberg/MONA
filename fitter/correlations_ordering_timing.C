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

void correlations_ordering_timing() {
  TString filefolder = "./default_detres/shower_energy_in_tracks/";

  gROOT->ProcessLine(".L FitFunction.C+");
  gSystem->Load("$OSCPROBDIR/libOscProb.so");

  //----------------------------------------------------------
  // detector response for tracks and showers
  //----------------------------------------------------------

  DetResponse drt_NO(DetResponse::track, "track_response_NO", 40, 1, 100, 40, -1, 1, 1, 0, 1);
  drt_NO.AddCut( &SummaryEvent::Get_track_ql0       , std::greater<double>()   , 0.5             , true );
  drt_NO.AddCut( &SummaryEvent::Get_track_ql1       , std::greater<double>()   , 0.5             , true );
  drt_NO.AddCut( &SummaryEvent::Get_RDF_track_score , std::greater<double>()   , 0.6             , true );
  drt_NO.AddCut( &SummaryEvent::Get_RDF_muon_score  , std::less_equal<double>(), 0.05            , true );
  drt_NO.AddCut( &SummaryEvent::Get_RDF_noise_score , std::less_equal<double>(), 0.18            , true );

  DetResponse drs_NO(DetResponse::shower, "shower_response_NO", 40, 1, 100, 40, -1, 1, 1, 0, 1);
  drs_NO.AddCut( &SummaryEvent::Get_shower_ql0     , std::greater<double>()   , 0.5             , true );
  drs_NO.AddCut( &SummaryEvent::Get_shower_ql1     , std::greater<double>()   , 0.5             , true );
  drs_NO.AddCut( &SummaryEvent::Get_RDF_track_score, std::less_equal<double>(),  0.6            , true );
  drs_NO.AddCut( &SummaryEvent::Get_RDF_muon_score , std::less_equal<double>(), 0.05            , true );
  drs_NO.AddCut( &SummaryEvent::Get_RDF_noise_score, std::less_equal<double>(), 0.5             , true );

  DetResponse drm_NO(DetResponse::mc_truth, "mc_response_NO", 40, 1, 100, 40, -1, 1, 1, 0, 1);
  drm_NO.AddCut( &SummaryEvent::Get_RDF_muon_score , std::less_equal<double>(), 0.05            , true );
  drm_NO.AddCut( &SummaryEvent::Get_RDF_noise_score, std::less_equal<double>(), 0.5             , true );

  DetResponse drt_IO(DetResponse::track, "track_response_IO", 40, 1, 100, 40, -1, 1, 1, 0, 1);
  drt_IO.AddCut( &SummaryEvent::Get_track_ql0       , std::greater<double>()   , 0.5             , true );
  drt_IO.AddCut( &SummaryEvent::Get_track_ql1       , std::greater<double>()   , 0.5             , true );
  drt_NO.AddCut( &SummaryEvent::Get_RDF_track_score , std::greater<double>()   , 0.6             , true );
  drt_IO.AddCut( &SummaryEvent::Get_RDF_muon_score  , std::less_equal<double>(), 0.05            , true );
  drt_IO.AddCut( &SummaryEvent::Get_RDF_noise_score , std::less_equal<double>(), 0.18            , true );

  DetResponse drs_IO(DetResponse::shower, "shower_response_IO", 40, 1, 100, 40, -1, 1, 1, 0, 1);
  drs_IO.AddCut( &SummaryEvent::Get_shower_ql0     , std::greater<double>()   , 0.5             , true );
  drs_IO.AddCut( &SummaryEvent::Get_shower_ql1     , std::greater<double>()   , 0.5             , true );
  drs_IO.AddCut( &SummaryEvent::Get_RDF_track_score, std::less_equal<double>(),  0.6            , true );
  drs_IO.AddCut( &SummaryEvent::Get_RDF_muon_score , std::less_equal<double>(), 0.05            , true );
  drs_IO.AddCut( &SummaryEvent::Get_RDF_noise_score, std::less_equal<double>(), 0.5             , true );
  
  DetResponse drm_IO(DetResponse::mc_truth, "mc_response_IO", 40, 1, 100, 40, -1, 1, 1, 0, 1);
  drm_IO.AddCut( &SummaryEvent::Get_RDF_muon_score , std::less_equal<double>(), 0.05            , true );
  drm_IO.AddCut( &SummaryEvent::Get_RDF_noise_score, std::less_equal<double>(), 0.5             , true );


  drt_NO.ReadFromFile(filefolder + "track_response_timing_NO.root") ;
  drt_IO.ReadFromFile(filefolder + "track_response_timing_IO.root") ;
  drs_NO.ReadFromFile(filefolder + "shower_response_timing_NO.root");
  drs_IO.ReadFromFile(filefolder + "shower_response_timing_IO.root");
  drm_NO.ReadFromFile(filefolder + "mc_response_timing_NO.root");
  drm_IO.ReadFromFile(filefolder + "mc_response_timing_IO.root");

  Double_t bybin = 0.5; // Theres only one by

  FitFunction tfitf_NO(&drt_NO, 3, "../data/eff_mass/EffMhists_elec_CC.root", "../data/eff_mass/EffMhists_muon_CC.root", "../data/eff_mass/EffMhists_tau_CC.root", "../data/eff_mass/EffMhists_elec_NC.root");
  FitFunction tfitf_IO(&drt_IO, 3, "../data/eff_mass/EffMhists_elec_CC.root", "../data/eff_mass/EffMhists_muon_CC.root", "../data/eff_mass/EffMhists_tau_CC.root", "../data/eff_mass/EffMhists_elec_NC.root");
  FitFunction sfitf_NO(&drs_NO, 3, "../data/eff_mass/EffMhists_elec_CC.root", "../data/eff_mass/EffMhists_muon_CC.root", "../data/eff_mass/EffMhists_tau_CC.root", "../data/eff_mass/EffMhists_elec_NC.root");
  FitFunction sfitf_IO(&drs_IO, 3, "../data/eff_mass/EffMhists_elec_CC.root", "../data/eff_mass/EffMhists_muon_CC.root", "../data/eff_mass/EffMhists_tau_CC.root", "../data/eff_mass/EffMhists_elec_NC.root");
  FitFunction mfitf_NO(&drm_NO, 3, "../data/eff_mass/EffMhists_elec_CC.root", "../data/eff_mass/EffMhists_muon_CC.root", "../data/eff_mass/EffMhists_tau_CC.root", "../data/eff_mass/EffMhists_elec_NC.root");
  FitFunction mfitf_IO(&drm_IO, 3, "../data/eff_mass/EffMhists_elec_CC.root", "../data/eff_mass/EffMhists_muon_CC.root", "../data/eff_mass/EffMhists_tau_CC.root", "../data/eff_mass/EffMhists_elec_NC.root");

  TH3D *h_cor_tracks  = (TH3D*)drt_NO.GetHist3D()->Clone("correlation_track_bins");
  TH3D *h_cor_showers = (TH3D*)drs_NO.GetHist3D()->Clone("correlation_shower_bins");
  TH3D *h_cor_mc      = (TH3D*)drm_NO.GetHist3D()->Clone("correlation_mc_bins");
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

        Double_t correlation_track = normed_innerproduct(track_weights_NO, track_weights_IO); // This is imported from correlations_ordering.so
        Double_t correlation_shower = normed_innerproduct(shower_weights_NO, shower_weights_IO);
        Double_t correlation_mc = normed_innerproduct(mc_weights_NO, mc_weights_IO);

        h_cor_tracks->SetBinContent(  ebin, ctbin, bybin, correlation_track );
        h_cor_showers->SetBinContent( ebin, ctbin, bybin, correlation_shower );
        h_cor_mc->SetBinContent( ebin, ctbin, bybin, correlation_shower );

      }
    }
  }
  TString output = "correlations.root";
  TFile fout(filefolder + output, "RECREATE");
  h_cor_tracks->Write();
  h_cor_showers->Write();
  h_cor_mc->Write();
  fout.Close();
}
