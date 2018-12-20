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

void weight_hist(Double_t ebin=10, Double_t ctbin=-0.8) {
  const int N_PID_CLASSES = 10;
  Double_t PID_step = 1 / float(N_PID_CLASSES);
  TString filefolder = TString::Format("./pid_detres/pid_binning_%i/", N_PID_CLASSES);

  gROOT->ProcessLine(".L FitFunction.C+");
  gSystem->Load("$OSCPROBDIR/libOscProb.so");

  //----------------------------------------------------------
  // detector response for tracks and showers
  //----------------------------------------------------------
  std::vector<DetResponse> track_response_vector;
  std::vector<DetResponse> shower_response_vector;

  for (int i = 0; i < N_PID_CLASSES; i++) {
    DetResponse track_response(DetResponse::track, TString::Format("track_response_%.2f", PID_step * i), 40, 1, 100, 40, -1, 1, 1, 0, 1);
    track_response.AddCut( &SummaryEvent::Get_track_ql0       , std::greater<double>()   , 0.5             , true );
    track_response.AddCut( &SummaryEvent::Get_track_ql1       , std::greater<double>()   , 0.5             , true );
    track_response.AddCut( &SummaryEvent::Get_RDF_track_score , std::greater<double>()   , PID_step * i    , true );
    track_response.AddCut( &SummaryEvent::Get_RDF_track_score , std::less_equal<double>(), PID_step * (i+1), true );
    track_response.AddCut( &SummaryEvent::Get_RDF_muon_score  , std::less_equal<double>(), 0.05            , true );
    track_response.AddCut( &SummaryEvent::Get_RDF_noise_score , std::less_equal<double>(), 0.18            , true );
    track_response_vector.push_back(track_response);
  }

  for (int i = 0; i < N_PID_CLASSES; i++) {
    DetResponse shower_response(DetResponse::shower, TString::Format("shower_response_%.2f", PID_step * i), 40, 1, 100, 40, -1, 1, 1, 0, 1);
    shower_response.AddCut( &SummaryEvent::Get_shower_ql0     , std::greater<double>()   , 0.5             , true );
    shower_response.AddCut( &SummaryEvent::Get_shower_ql1     , std::greater<double>()   , 0.5             , true );
    shower_response.AddCut( &SummaryEvent::Get_RDF_track_score, std::greater<double>()   , PID_step * i    , true );
    shower_response.AddCut( &SummaryEvent::Get_RDF_track_score, std::less_equal<double>(), PID_step * (i+1), true );
    shower_response.AddCut( &SummaryEvent::Get_RDF_muon_score , std::less_equal<double>(), 0.05            , true );
    shower_response.AddCut( &SummaryEvent::Get_RDF_noise_score, std::less_equal<double>(), 0.5             , true );
    shower_response_vector.push_back(shower_response);
  }

  for (int i = 0; i < N_PID_CLASSES; i++) {
    track_response_vector[i].ReadFromFile(filefolder + TString::Format("track_response_%.2f.root" , PID_step * i));
    shower_response_vector[i].ReadFromFile(filefolder + TString::Format("shower_response_%.2f.root", PID_step * i));
  }  

  cout << "NOTICE: Finished filling response" << endl;

  //----------------------------------------------------------
  // fill a detected histogram using the fit function
  //----------------------------------------------------------

  Double_t bybin = 0.5; // Theres only one by

  std::vector<TH1D> h_weights_track; 
  std::vector<TH1D> h_weights_shower; 
  auto effmass_folder = (TString)getenv("NMHDIR") + "/data/eff_mass/";
  for (int i = 0; i < N_PID_CLASSES; i++){
    FitFunction tfitf(&track_response_vector[i], 3, effmass_folder + "EffMhists_elec_CC.root", effmass_folder + "EffMhists_muon_CC.root", effmass_folder + "EffMhists_tau_CC.root", effmass_folder + "EffMhists_elec_NC.root");
  
    FitFunction sfitf(&shower_response_vector[i], 3, effmass_folder + "EffMhists_elec_CC.root", effmass_folder + "EffMhists_muon_CC.root", effmass_folder + "EffMhists_tau_CC.root", effmass_folder + "EffMhists_elec_NC.root");
  
    std::vector<TrueB> track_weights;
    std::vector<TrueB> shower_weights;
    
    
    track_weights  = track_response_vector[i].GetBinWeights(ebin, ctbin, bybin);
    shower_weights = shower_response_vector[i].GetBinWeights(ebin, ctbin, bybin);
 

    TH1D t = TH1D(Form("h_w_t_q%i", i), Form("Weights track q%i", i),  100, 0, 0.05); 
    TH1D s = TH1D(Form("h_w_s_q%i", i), Form("Weights shower q%i", i), 100, 0, 0.05);
    for (int j = 0; j < track_weights.size(); j++) {
      t.Fill(track_weights[j].fW);
    }
    for (int j = 0; j < shower_weights.size(); j++) {
      s.Fill(shower_weights[j].fW);
    }
    h_weights_track.push_back(t);
    h_weights_shower.push_back(s);

  
  }
  TString output = TString::Format("weights_E%.2f_ct%.2f.root", ebin, ctbin);
  TFile fout(filefolder + output, "RECREATE");
  for (int i = 0; i < N_PID_CLASSES; i++){
    h_weights_track[i].Write();
    h_weights_shower[i].Write();
  }
  fout.Close();
  cout << "NOTICE: Written " << output << endl;
}
