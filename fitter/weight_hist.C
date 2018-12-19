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
  TString filefolder = "./pid_detres/";

  gROOT->ProcessLine(".L FitFunction.C+");
  gSystem->Load("$OSCPROBDIR/libOscProb.so");

  //----------------------------------------------------------
  // detector response for tracks and showers
  //----------------------------------------------------------
  std::vector<DetResponse> drts;
  std::vector<DetResponse> drss;

  for (int i = 0; i < N_PID_CLASSES; i++) {
    DetResponse drt(DetResponse::track, TString::Format("track_response_%.1f", PID_step * i), 40, 1, 100, 40, -1, 1, 1, 0, 1);
    drt.AddCut( &SummaryEvent::Get_track_ql0       , std::greater<double>()   , 0.5             , true );
    drt.AddCut( &SummaryEvent::Get_track_ql1       , std::greater<double>()   , 0.5             , true );
    drt.AddCut( &SummaryEvent::Get_RDF_track_score , std::greater<double>()   , PID_step * i    , true );
    drt.AddCut( &SummaryEvent::Get_RDF_track_score , std::less_equal<double>(), PID_step * (i+1), true );
    drt.AddCut( &SummaryEvent::Get_RDF_muon_score  , std::less_equal<double>(), 0.05            , true );
    drt.AddCut( &SummaryEvent::Get_RDF_noise_score , std::less_equal<double>(), 0.18            , true );
    drts.push_back(drt);
  }

  for (int i = 0; i < N_PID_CLASSES; i++) {
    DetResponse drs(DetResponse::shower, TString::Format("shower_response_%.1f", PID_step * i), 40, 1, 100, 40, -1, 1, 1, 0, 1);
    drs.AddCut( &SummaryEvent::Get_shower_ql0     , std::greater<double>()   , 0.5             , true );
    drs.AddCut( &SummaryEvent::Get_shower_ql1     , std::greater<double>()   , 0.5             , true );
    drs.AddCut( &SummaryEvent::Get_RDF_track_score, std::greater<double>()   , PID_step * i    , true );
    drs.AddCut( &SummaryEvent::Get_RDF_track_score, std::less_equal<double>(), PID_step * (i+1), true );
    drs.AddCut( &SummaryEvent::Get_RDF_muon_score , std::less_equal<double>(), 0.05            , true );
    drs.AddCut( &SummaryEvent::Get_RDF_noise_score, std::less_equal<double>(), 0.5             , true );
    drss.push_back(drs);
  }

  for (int i = 0; i < N_PID_CLASSES; i++) {
    drts[i].ReadFromFile(filefolder + TString::Format("track_response_%.1f.root" , PID_step * i));
    drss[i].ReadFromFile(filefolder + TString::Format("shower_response_%.1f.root", PID_step * i));
  }  

  cout << "NOTICE: Finished filling response" << endl;

  //----------------------------------------------------------
  // fill a detected histogram using the fit function
  //----------------------------------------------------------

  Double_t bybin = 0.5; // Theres only one by

  std::vector<TH1D> h_weights_track; 
  std::vector<TH1D> h_weights_shower; 
  for (int i = 0; i < N_PID_CLASSES; i++){
    FitFunction tfitf(&drts[i], 3, "../data/eff_mass/EffMhists_elec_CC.root", "../data/eff_mass/EffMhists_muon_CC.root", "../data/eff_mass/EffMhists_tau_CC.root", "../data/eff_mass/EffMhists_elec_NC.root");
  
    FitFunction sfitf(&drss[i], 3, "../data/eff_mass/EffMhists_elec_CC.root", "../data/eff_mass/EffMhists_muon_CC.root", "../data/eff_mass/EffMhists_tau_CC.root", "../data/eff_mass/EffMhists_elec_NC.root");
  
    std::vector<TrueB> track_weights;
    std::vector<TrueB> shower_weights;
    
    
    track_weights  = drts[i].GetBinWeights(ebin, ctbin, bybin);
    shower_weights = drss[i].GetBinWeights(ebin, ctbin, bybin);
 

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
  TString output = TString::Format("weights_E%.1f_ct%.1f.root", ebin, ctbin);
  TFile fout(filefolder + output, "RECREATE");
  for (int i = 0; i < N_PID_CLASSES; i++){
    h_weights_track[i].Write();
    h_weights_shower[i].Write();
  }
  fout.Close();
  cout << "NOTICE: Written " << output << endl;
}
