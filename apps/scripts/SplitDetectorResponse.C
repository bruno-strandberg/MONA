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

void SplitDetectorResponse() {

  const int N_PID_CLASSES = 10;
  Double_t PID_step = 1 / float(N_PID_CLASSES);
  TString filefolder = TString::Format("./pid_detres/pid_binning_%i", N_PID_CLASSES);

  gROOT->ProcessLine(".L FitFunction.C+");
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

  //auto simdata_file = (TString)getenv("NMHDIR") + "/data/ORCA_MC_summary_ORCA115_23x9m_ECAP0418.root"; // To be used once cross-check of file below is done.
  SummaryParser sp("../../data/ORCA_MC_summary_all_10Apr2018.root");
  bool writeFiles = true;
  for (Int_t i = 0; i < sp.GetTree()->GetEntries(); i++) {
    if (i == 1000) break;
    if (i % (Int_t)1e6 == 0) cout << "Event: " << i << endl;
    sp.GetTree()->GetEntry(i);
    SummaryEvent *evt = sp.GetEvt();
    if (evt->Get_MC_is_neutrino() < 0.5) continue; // Is this still relevant?
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
  // fill a detected histogram using the fit function
  //----------------------------------------------------------


  auto effmass_folder = (TString)getenv("NMHDIR") + "/data/eff_mass/";
  for (int i = 0; i < N_PID_CLASSES; i++){
    FitFunction tfitf(&track_response_vector[i], 3, effmass_folder + "EffMhists_elec_CC.root", effmass_folder + "EffMhists_muon_CC.root", effmass_folder + "EffMhists_tau_CC.root", effmass_folder + "EffMhists_elec_NC.root");
    FitFunction sfitf(&shower_response_vector[i], 3, effmass_folder + "EffMhists_elec_CC.root", effmass_folder + "EffMhists_muon_CC.root", effmass_folder + "EffMhists_tau_CC.root", effmass_folder + "EffMhists_elec_NC.root");
    FitFunction mfitf(&mc_response_vector[i], 3, effmass_folder + "EffMhists_elec_CC.root", effmass_folder + "EffMhists_muon_CC.root", effmass_folder + "EffMhists_tau_CC.root", effmass_folder + "EffMhists_elec_NC.root");
  
    TH3D *hdet_tracks  = (TH3D*)track_response_vector[i].GetHist3D()->Clone("detected_tracks");
    TH3D *hdet_showers = (TH3D*)shower_response_vector[i].GetHist3D()->Clone("detected_showers");
    TH3D *hdet_mc      = (TH3D*)mc_response_vector[i].GetHist3D()->Clone("detected_mc");
    hdet_tracks->Reset();
    hdet_showers->Reset();
    hdet_mc->Reset();
  
    // The construction of TString is complaining about a conversion when inside a tuple...
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
 
            std::tuple<double, double> track_response  = tfitf.operator()(x, p);
            std::tuple<double, double> shower_response = sfitf.operator()(x, p);
            std::tuple<double, double> mc_response     = mfitf.operator()(x, p);

            hdet_tracks->SetBinContent(  ebin, ctbin, bybin, std::get<0>(track_response) );
            hdet_tracks->SetBinError(    ebin, ctbin, bybin, std::get<1>(track_response) );
            hdet_showers->SetBinContent( ebin, ctbin, bybin, std::get<0>(shower_response) );
            hdet_showers->SetBinError(   ebin, ctbin, bybin, std::get<1>(shower_response) );
            hdet_mc->SetBinContent(      ebin, ctbin, bybin, std::get<0>(mc_response) );
            hdet_mc->SetBinError(        ebin, ctbin, bybin, std::get<1>(mc_response) );
          }
        }
      }
    
      cout << "NOTICE: Finished filling histograms" << endl;
      cout << "NOTICE: Time for filling hists " << (Double_t)timer.RealTime() << endl;
    
      TString output = TString::Format("split_expected_evts_%s_%.2f.root", order_string.c_str(), PID_step * i);
      TFile fout(filefolder + output, "RECREATE");
      hdet_tracks->Write();
      hdet_showers->Write();
      hdet_mc->Write();
      fout.Close();
      cout << "NOTICE: Written " << output << endl;
    }
  }
}
