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

void SplitQuality(bool mass_ordering=true) {
  const int N_PID_CLASSES = 2;
  Double_t PID_step = 1 / float(N_PID_CLASSES);
  TString filefolder = "./quality_detres/";
  string order_string;
  double p5;
  // GET RID OF THIS boolean STUFF
  if (mass_ordering) { // True = NO, +1; False = IO, -1
    order_string = "NO";
    p5 = 2.52e-3;
  } 
  else { 
    order_string = "IO";
    p5 = -2.44e-3;
  }

  gROOT->ProcessLine(".L FitFunction.C+");
  gSystem->Load("$OSCPROBDIR/libOscProb.so");

  std::map<int, float> cut_map;
  cut_map.insert(std::make_pair(0, 0.0));
  cut_map.insert(std::make_pair(1, 0.6));
  cut_map.insert(std::make_pair(2, 1.0));

  //----------------------------------------------------------
  // detector response for tracks and showers
  //----------------------------------------------------------
  std::vector<DetResponse> track_response_vector;
  for (int i = 0; i < N_PID_CLASSES; i++) {
    DetResponse track_response(DetResponse::track, TString::Format("track_response_%s_%.2f", order_string.c_str(), cut_map[i]), 40, 1, 100, 40, -1, 1, 1, 0, 1);
    track_response.AddCut( &SummaryEvent::Get_track_ql0       , std::greater<double>()   , 0.5         , true );
    track_response.AddCut( &SummaryEvent::Get_track_ql1       , std::greater<double>()   , 0.5         , true );
    track_response.AddCut( &SummaryEvent::Get_RDF_track_score , std::greater<double>()   , cut_map[i]  , true );
    track_response.AddCut( &SummaryEvent::Get_RDF_track_score , std::less_equal<double>(), cut_map[i+1], true );
    track_response.AddCut( &SummaryEvent::Get_RDF_muon_score  , std::less_equal<double>(), 0.05        , true );
    track_response.AddCut( &SummaryEvent::Get_RDF_noise_score , std::less_equal<double>(), 0.18        , true );
    track_response_vector.push_back(track_response);
  }

  std::vector<DetResponse> shower_response_vector;
  for (int i = 0; i < N_PID_CLASSES; i++) {
    DetResponse shower_response(DetResponse::shower, TString::Format("shower_response_%s_%.2f", order_string.c_str(), cut_map[i]), 40, 1, 100, 40, -1, 1, 1, 0, 1);
    shower_response.AddCut( &SummaryEvent::Get_shower_ql0     , std::greater<double>()   , 0.5         , true );
    shower_response.AddCut( &SummaryEvent::Get_shower_ql1     , std::greater<double>()   , 0.5         , true );
    shower_response.AddCut( &SummaryEvent::Get_RDF_track_score, std::greater<double>()   , cut_map[i]  , true );
    shower_response.AddCut( &SummaryEvent::Get_RDF_track_score, std::less_equal<double>(), cut_map[i+1], true );
    shower_response.AddCut( &SummaryEvent::Get_RDF_muon_score , std::less_equal<double>(), 0.05        , true );
    shower_response.AddCut( &SummaryEvent::Get_RDF_noise_score, std::less_equal<double>(), 0.5         , true );
    shower_response_vector.push_back(shower_response);
  }

  SummaryParser sp("../../data/ORCA_MC_summary_all_10Apr2018.root");
  bool writeFiles = true;
  for (Int_t i = 0; i < sp.GetTree()->GetEntries(); i++) {
    sp.GetTree()->GetEntry(i);
    SummaryEvent *evt = sp.GetEvt();
    if (evt->Get_MC_is_neutrino() < 0.5) continue;
    for (int i = 0; i < N_PID_CLASSES; i++) {
      track_response_vector[i].Fill(evt);
      shower_response_vector[i].Fill(evt);
    }
  }

  for (int i = 0; i < N_PID_CLASSES; i++) { 
    if (writeFiles) { // This flag is to not accidently overwrite with empty files
      track_response_vector[i].WriteToFile(filefolder + TString::Format("track_response_%s_%.2f.root" , order_string.c_str(), cut_map[i]));
      shower_response_vector[i].WriteToFile(filefolder + TString::Format("shower_response_%s_%.2f.root", order_string.c_str(), cut_map[i]));
    }
  }

  cout << "NOTICE: Finished filling response" << endl;

  //----------------------------------------------------------
  // fill a detected histogram using the fit function
  //----------------------------------------------------------


  auto effmass_folder = (TString)getenv("NMHDIR") + "/data/eff_mass/";
  for (int i = 0; i < N_PID_CLASSES; i++){
    FitFunction tfitf(&track_response_vector[i], 3, effmass_folder + "EffMhists_elec_CC.root", effmass_folder + "EffMhists_muon_CC.root", effmass_folder + "EffMhists_tau_CC.root", effmass_folder + "EffMhists_elec_NC.root");
    FitFunction sfitf(&shower_response_vector[i], 3, effmass_folder + "EffMhists_elec_CC.root", effmass_folder + "EffMhists_muon_CC.root", effmass_folder + "EffMhists_tau_CC.root", effmass_folder + "EffMhists_elec_NC.root");
  
    TH3D *hdet_tracks  = (TH3D*)track_response_vector[i].GetHist3D()->Clone("detected_tracks");
    TH3D *hdet_showers = (TH3D*)shower_response_vector[i].GetHist3D()->Clone("detected_showers");
    hdet_tracks->Reset();
    hdet_showers->Reset();
  
    // Mass ordering implied by the oscillation parameters
    double p[] = {0.303, 0.0214, 0.5, 0, 7.53e-5, p5}; // Taking p5 from the boolean, this is only difference between NO and IO.
  
    TStopwatch timer;
  
    for (Int_t ebin = 1; ebin <= hdet_tracks->GetXaxis()->GetNbins(); ebin++) {
      for (Int_t ctbin = 1; ctbin <= hdet_tracks->GetYaxis()->GetNbins(); ctbin++) {
        for (Int_t bybin = 1; bybin <= hdet_tracks->GetZaxis()->GetNbins(); bybin++) {
  
          Double_t E  = hdet_tracks->GetXaxis()->GetBinCenter( ebin );
          Double_t ct = hdet_tracks->GetYaxis()->GetBinCenter( ctbin );
          Double_t by = hdet_tracks->GetZaxis()->GetBinCenter( bybin );
  
          double x[] = {E, ct, by};
 
          std::tuple<double, double> track_response  = tfitf.operator()(x, p);
          std::tuple<double, double> shower_response = sfitf.operator()(x, p);

          hdet_tracks->SetBinContent(  ebin, ctbin, bybin, std::get<0>(track_response) );
          hdet_tracks->SetBinError(    ebin, ctbin, bybin, std::get<1>(track_response) );
          hdet_showers->SetBinContent( ebin, ctbin, bybin, std::get<0>(shower_response) );
          hdet_showers->SetBinError(   ebin, ctbin, bybin, std::get<1>(shower_response) );
        }
      }
    }
    
    cout << "NOTICE: Finished filling histograms" << endl;
    cout << "NOTICE: Time for filling hists " << (Double_t)timer.RealTime() << endl;
  
    TString output = TString::Format("split_pid_%s_%.2f.root", order_string.c_str(), cut_map[i]);
    TFile fout(filefolder + output, "RECREATE");
    hdet_tracks->Write();
    hdet_showers->Write();
    fout.Close();
    cout << "NOTICE: Written " << output << endl;
  }
}
