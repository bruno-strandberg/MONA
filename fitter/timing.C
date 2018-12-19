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

#include <iostream>
using namespace std;

void timing() {


  bool use_half_data = false;
  TString filefolder = "./default_detres/";

  gROOT->ProcessLine(".L FitFunction.C+");
  gSystem->Load("$OSCPROBDIR/libOscProb.so");

  //----------------------------------------------------------
  // detector response for tracks and showers
  //----------------------------------------------------------
  DetResponse drt(DetResponse::track, "track_response", 40, 1, 100, 40, -1, 1, 1, 0, 1);
  drt.AddCut( &SummaryEvent::Get_track_ql0       , std::greater<double>()   ,  0.5, true );
  drt.AddCut( &SummaryEvent::Get_track_ql1       , std::greater<double>()   ,  0.5, true );
  drt.AddCut( &SummaryEvent::Get_RDF_track_score , std::greater<double>()   ,  0.6, true );
  drt.AddCut( &SummaryEvent::Get_RDF_muon_score  , std::less_equal<double>(), 0.05, true );
  drt.AddCut( &SummaryEvent::Get_RDF_noise_score , std::less_equal<double>(), 0.18, true );

  DetResponse drs(DetResponse::shower, "shower_response", 40, 1, 100, 40, -1, 1, 1, 0, 1);
  drs.AddCut( &SummaryEvent::Get_shower_ql0     , std::greater<double>()   ,  0.5, true );
  drs.AddCut( &SummaryEvent::Get_shower_ql1     , std::greater<double>()   ,  0.5, true );
  drs.AddCut( &SummaryEvent::Get_RDF_track_score, std::less_equal<double>(),  0.6, true );
  drs.AddCut( &SummaryEvent::Get_RDF_muon_score , std::less_equal<double>(), 0.05, true );
  drs.AddCut( &SummaryEvent::Get_RDF_noise_score, std::less_equal<double>(),  0.5, true );

  DetResponse drm(DetResponse::mc_truth, "mc_response", 40, 1, 100, 40, -1, 1, 1, 0, 1); // In multi pid bin, we make mc detres for all
  drm.AddCut( &SummaryEvent::Get_RDF_muon_score , std::less_equal<double>(), 0.05, true );
  drm.AddCut( &SummaryEvent::Get_RDF_noise_score, std::less_equal<double>(),  0.5, true );

  Double_t n_allowed = n_evts;
  SummaryParser sp("../data/ORCA_MC_summary_all_10Apr2018.root");
  for (Int_t i = 0; i < sp.GetTree()->GetEntries(); i++) {
    sp.GetTree()->GetEntry(i);
    SummaryEvent *evt = sp.GetEvt();
    if (evt->Get_MC_is_neutrino() < 0.5) continue;
    if (use_half_data) {
      if ((std::abs(sp.GetEvt()->Get_MC_type()) == 11) or (std::abs(sp.GetEvt()->Get_MC_type()) == 12) or 
          (std::abs(sp.GetEvt()->Get_MC_type()) == 13) or (std::abs(sp.GetEvt()->Get_MC_type()) == 14)) {
        if (evt->Get_MC_runID() <=300) continue ;
      }
      if ((std::abs(sp.GetEvt()->Get_MC_type()) == 15) or (std::abs(sp.GetEvt()->Get_MC_type()) == 16)) {
        if (evt->Get_MC_runID() <=900) continue ;
      }
    }
    drt.Fill(evt);
    drs.Fill(evt);
    drm.Fill(evt);
  }

  drt.WriteToFile(filefolder + "track_response_timing.root");
  drs.WriteToFile(filefolder + "shower_response_timing.root");
  drm.WriteToFile(filefolder + "mc_response_timing.root");

//  drt.ReadFromFile("track_response_timing.root");
//  drs.ReadFromFile("shower_response_timing.root");
//  drm.ReadFromFile("mc_response_timing.root");
  
  cout << "NOTICE: Finished filling response" << endl;

  //----------------------------------------------------------
  // fill a detected histogram using the fit function
  //----------------------------------------------------------

  FitFunction tfitf(&drt, 3, "../data/eff_mass/EffMhists_elec_CC.root", "../data/eff_mass/EffMhists_muon_CC.root", "../data/eff_mass/EffMhists_tau_CC.root", "../data/eff_mass/EffMhists_elec_NC.root");
  FitFunction sfitf(&drs, 3, "../data/eff_mass/EffMhists_elec_CC.root", "../data/eff_mass/EffMhists_muon_CC.root", "../data/eff_mass/EffMhists_tau_CC.root", "../data/eff_mass/EffMhists_elec_NC.root");
  FitFunction mfitf(&drm, 3, "../data/eff_mass/EffMhists_elec_CC.root", "../data/eff_mass/EffMhists_muon_CC.root", "../data/eff_mass/EffMhists_tau_CC.root", "../data/eff_mass/EffMhists_elec_NC.root");

  TH3D *hdet_tracks  = (TH3D*)drt.GetHist3D()->Clone("detected_tracks");
  TH3D *hdet_showers = (TH3D*)drs.GetHist3D()->Clone("detected_showers");
  TH3D *hdet_mc      = (TH3D*)drm.GetHist3D()->Clone("detected_mc");
  hdet_tracks->Reset();
  hdet_showers->Reset();
  hdet_mc->Reset();

  std::vector<std::tuple<string, Double_t>> orderings; // Value is dm32
  orderings.push_back(std::make_tuple("NO", 2.52e-3));
  orderings.push_back(std::make_tuple("IO", -2.44e-3));

  for (auto order: orderings) {
    double p[] = {0.303, 0.0214, 0.5, 0, 7.53e-5, std::get<1>(order)}; // dm32 is last parameter.
    string order_string = std::get<0>(order);  

    TStopwatch timer;

    for (Int_t ebin = 1; ebin <= hdet_tracks->GetXaxis()->GetNbins(); ebin++) {
      for (Int_t ctbin = 1; ctbin <= hdet_tracks->GetYaxis()->GetNbins(); ctbin++) {
        for (Int_t bybin = 1; bybin <= hdet_tracks->GetZaxis()->GetNbins(); bybin++) {

          Double_t E  = hdet_tracks->GetXaxis()->GetBinCenter( ebin );
          Double_t ct = hdet_tracks->GetYaxis()->GetBinCenter( ctbin );
          Double_t by = hdet_tracks->GetZaxis()->GetBinCenter( bybin );

          double x[] = {E, ct, by};

          // With Errors
          std::tuple<double, double> track_response  = tfitf.operator()(x, p);
          std::tuple<double, double> shower_response = sfitf.operator()(x, p);
          std::tuple<double, double> mc_response     = mfitf.operator()(x, p);

          hdet_tracks->SetBinContent(  ebin, ctbin, bybin, std::get<0>(track_response) );
          hdet_tracks->SetBinError(    ebin, ctbin, bybin, std::get<1>(track_response) );
          hdet_showers->SetBinContent( ebin, ctbin, bybin, std::get<0>(shower_response) );
          hdet_showers->SetBinError(   ebin, ctbin, bybin, std::get<1>(shower_response) );
          hdet_mc->SetBinContent( ebin, ctbin, bybin, std::get<0>(mc_response) );
          hdet_mc->SetBinError(   ebin, ctbin, bybin, std::get<1>(mc_response) );

        }
      }
    }
    
    cout << "NOTICE: Finished filling histograms" << endl;
    cout << "NOTICE: Time for filling hists " << (Double_t)timer.RealTime() << endl;
    

    TString output = Form("timing_%s.root", order_string.c_str());
    TFile fout(filefolder + output,"RECREATE");
    hdet_tracks->Write();
    hdet_showers->Write();
    hdet_mc->Write();
    fout.Close();

  }
}
