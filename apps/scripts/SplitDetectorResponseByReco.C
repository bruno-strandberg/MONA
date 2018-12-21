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

void SplitDetectorResponseByReco() {

  TString filefolder = "./energy_detres/";

  gROOT->ProcessLine(".L FitFunction.C+");
  gSystem->Load("$OSCPROBDIR/libOscProb.so");

  //----------------------------------------------------------
  // detector response for tracks and showers
  //----------------------------------------------------------
  // Good tracks only, no overlap with good showers
  // gt = good track, gs = good shower, ge = good event
  DetResponse response_good_track(DetResponse::hybridE, "track_response_gt", 40, 1, 100, 40, -1, 1, 1, 0, 1);
  response_good_track.AddCut( &SummaryEvent::Get_track_ql0       , std::greater<double>()   ,  0.5, true );
  response_good_track.AddCut( &SummaryEvent::Get_track_ql1       , std::greater<double>()   ,  0.5, true );
  response_good_track.AddCut( &SummaryEvent::Get_shower_ql0      , std::less<double>()      ,  0.5, true );
  response_good_track.AddCut( &SummaryEvent::Get_shower_ql1      , std::less<double>()      ,  0.5, true );
  response_good_track.AddCut( &SummaryEvent::Get_RDF_track_score , std::greater<double>()   ,  0.6, true );
  response_good_track.AddCut( &SummaryEvent::Get_RDF_muon_score  , std::less_equal<double>(), 0.05, true );
  response_good_track.AddCut( &SummaryEvent::Get_RDF_noise_score , std::less_equal<double>(), 0.18, true );

  // Good showers only, no overlap with good tracks
  DetResponse response_good_shower(DetResponse::hybridE, "track_response_gs", 40, 1, 100, 40, -1, 1, 1, 0, 1);
  response_good_shower.AddCut( &SummaryEvent::Get_track_ql0       , std::greater<double>()   ,  0.5, true ); // This would cut almost all events, since they all have fTrack_ql0 > 0.5
  response_good_shower.AddCut( &SummaryEvent::Get_track_ql1       , std::less<double>()      ,  0.5, true ); // Moreover: the numbers are consistent with resolution_plot_flav_complementary_events.C
  response_good_shower.AddCut( &SummaryEvent::Get_shower_ql0      , std::greater<double>()   ,  0.5, true );
  response_good_shower.AddCut( &SummaryEvent::Get_shower_ql1      , std::greater<double>()   ,  0.5, true );
  response_good_shower.AddCut( &SummaryEvent::Get_RDF_track_score , std::greater<double>()   ,  0.6, true );
  response_good_shower.AddCut( &SummaryEvent::Get_RDF_muon_score  , std::less_equal<double>(), 0.05, true );
  response_good_shower.AddCut( &SummaryEvent::Get_RDF_noise_score , std::less_equal<double>(), 0.18, true );

  // Good events only, no overlap with bad tracks or bad showers
  DetResponse response_good_event(DetResponse::hybridE, "track_response_ge", 40, 1, 100, 40, -1, 1, 1, 0, 1);
  response_good_event.AddCut( &SummaryEvent::Get_track_ql0       , std::greater<double>()   ,  0.5, true );
  response_good_event.AddCut( &SummaryEvent::Get_track_ql1       , std::greater<double>()   ,  0.5, true );
  response_good_event.AddCut( &SummaryEvent::Get_shower_ql0      , std::greater<double>()   ,  0.5, true );
  response_good_event.AddCut( &SummaryEvent::Get_shower_ql1      , std::greater<double>()   ,  0.5, true );
  response_good_event.AddCut( &SummaryEvent::Get_RDF_track_score , std::greater<double>()   ,  0.6, true );
  response_good_event.AddCut( &SummaryEvent::Get_RDF_muon_score  , std::less_equal<double>(), 0.05, true );
  response_good_event.AddCut( &SummaryEvent::Get_RDF_noise_score , std::less_equal<double>(), 0.18, true );

  DetResponse response_shower(DetResponse::hybridE, "shower_response", 40, 1, 100, 40, -1, 1, 1, 0, 1);
  response_shower.AddCut( &SummaryEvent::Get_shower_ql0     , std::greater<double>()   ,  0.5, true );
  response_shower.AddCut( &SummaryEvent::Get_shower_ql1     , std::greater<double>()   ,  0.5, true );
  response_shower.AddCut( &SummaryEvent::Get_RDF_track_score, std::less_equal<double>(),  0.6, true );
  response_shower.AddCut( &SummaryEvent::Get_RDF_muon_score , std::less_equal<double>(), 0.05, true );
  response_shower.AddCut( &SummaryEvent::Get_RDF_noise_score, std::less_equal<double>(),  0.5, true );

  DetResponse response_mc(DetResponse::mc_truth, "mc_response", 40, 1, 100, 40, -1, 1, 1, 0, 1);
  response_mc.AddCut( &SummaryEvent::Get_RDF_muon_score , std::less_equal<double>(), 0.05, true );
  response_mc.AddCut( &SummaryEvent::Get_RDF_noise_score, std::less_equal<double>(),  0.5, true );

  auto summary_file = (TString)getenv("NMHDIR") + "/data/ORCA_MC_summary_all_10Apr2018.root";
  SummaryParser sp(summary_file);
  for (Int_t i = 0; i < sp.GetTree()->GetEntries(); i++) {
    if (i % (Int_t)1e6 == 0) cout << "Event: " << i << endl;
    sp.GetTree()->GetEntry(i);
    SummaryEvent *evt = sp.GetEvt();
    if (evt->Get_MC_is_neutrino() < 0.5) continue;
    response_good_track.Fill(evt);
    response_good_shower.Fill(evt);
    response_good_event.Fill(evt);
    response_shower.Fill(evt);
    response_mc.Fill(evt);
  }

  response_good_track.WriteToFile(filefolder + "track_response_gt.root");
  response_good_shower.WriteToFile(filefolder + "track_response_gs.root");
  response_good_event.WriteToFile(filefolder + "track_response_ge.root");
  response_shower.WriteToFile(filefolder  + "shower_response.root");
  response_mc.WriteToFile(filefolder  + "mc_response.root");

  cout << "NOTICE: Finished filling response" << endl;

  //----------------------------------------------------------
  // fill a detected histogram using the fit function
  //----------------------------------------------------------

  auto effmass_folder = (TString)getenv("NMHDIR") + "/data/eff_mass/";
  FitFunction gtfitf(&response_good_track, 3, effmass_folder + "EffMhists_elec_CC.root", effmass_folder + "EffMhists_muon_CC.root", effmass_folder + "EffMhists_tau_CC.root", effmass_folder + "EffMhists_elec_NC.root");
  FitFunction gsfitf(&response_good_shower, 3, effmass_folder + "EffMhists_elec_CC.root", effmass_folder + "EffMhists_muon_CC.root", effmass_folder + "EffMhists_tau_CC.root", effmass_folder + "EffMhists_elec_NC.root");
  FitFunction gefitf(&response_good_event, 3, effmass_folder + "EffMhists_elec_CC.root", effmass_folder + "EffMhists_muon_CC.root", effmass_folder + "EffMhists_tau_CC.root", effmass_folder + "EffMhists_elec_NC.root");
  FitFunction sfitf(&response_shower, 3, effmass_folder + "EffMhists_elec_CC.root", effmass_folder + "EffMhists_muon_CC.root", effmass_folder + "EffMhists_tau_CC.root", effmass_folder + "EffMhists_elec_NC.root");
  FitFunction mfitf(&response_mc, 3, effmass_folder + "EffMhists_elec_CC.root", effmass_folder + "EffMhists_muon_CC.root", effmass_folder + "EffMhists_tau_CC.root", effmass_folder + "EffMhists_elec_NC.root");

  TH3D *hdet_good_t  = (TH3D*)response_good_track.GetHist3D()->Clone("detected_tracks_gt"); // Track channel, using good tracks, etc.
  TH3D *hdet_good_s  = (TH3D*)response_good_shower.GetHist3D()->Clone("detected_tracks_gs");
  TH3D *hdet_good_e  = (TH3D*)response_good_event.GetHist3D()->Clone("detected_tracks_ge");
  TH3D *hdet_showers = (TH3D*)response_shower.GetHist3D()->Clone("detected_showers");
  TH3D *hdet_mc      = (TH3D*)response_mc.GetHist3D()->Clone("detected_mc");
  hdet_good_t->Reset();
  hdet_good_s->Reset();
  hdet_good_e->Reset();
  hdet_showers->Reset();
  hdet_mc->Reset();

  // This is so rediculous, the construction of TString is complaining about a conversion when inside a tuple...
  std::vector<std::tuple<string, Double_t>> orderings; // Value is dm32
  orderings.push_back(std::make_tuple("NO", 2.52e-3));
  orderings.push_back(std::make_tuple("IO", -2.44e-3));

  for (auto order: orderings) {
    double p[] = {0.303, 0.0214, 0.5, 0, 7.53e-5, std::get<1>(order)}; // dm32 is last parameter.
    string order_string = std::get<0>(order);  

    TStopwatch timer;
  
    for (Int_t ebin = 1; ebin <= hdet_showers->GetXaxis()->GetNbins(); ebin++) {
      for (Int_t ctbin = 1; ctbin <= hdet_showers->GetYaxis()->GetNbins(); ctbin++) {
        for (Int_t bybin = 1; bybin <= hdet_showers->GetZaxis()->GetNbins(); bybin++) {
  
          Double_t E  = hdet_showers->GetXaxis()->GetBinCenter( ebin );
          Double_t ct = hdet_showers->GetYaxis()->GetBinCenter( ctbin );
          Double_t by = hdet_showers->GetZaxis()->GetBinCenter( bybin );
  
          double x[] = {E, ct, by};
  
          // With Errors
          std::tuple<double, double> track_response_gt = gtfitf.operator()(x, p);
          std::tuple<double, double> track_response_gs = gsfitf.operator()(x, p);
          std::tuple<double, double> track_response_ge = gefitf.operator()(x, p);
          std::tuple<double, double> shower_response   = sfitf.operator()(x, p);
          std::tuple<double, double> mc_response       = mfitf.operator()(x, p);
  
          hdet_good_t->SetBinContent(  ebin, ctbin, bybin, std::get<0>(track_response_gt) );
          hdet_good_t->SetBinError(    ebin, ctbin, bybin, std::get<1>(track_response_gt) );
          hdet_good_s->SetBinContent(  ebin, ctbin, bybin, std::get<0>(track_response_gs) );
          hdet_good_s->SetBinError(    ebin, ctbin, bybin, std::get<1>(track_response_gs) );
          hdet_good_e->SetBinContent(  ebin, ctbin, bybin, std::get<0>(track_response_ge) );
          hdet_good_e->SetBinError(    ebin, ctbin, bybin, std::get<1>(track_response_ge) );
          hdet_showers->SetBinContent( ebin, ctbin, bybin, std::get<0>(shower_response) );
          hdet_showers->SetBinError(   ebin, ctbin, bybin, std::get<1>(shower_response) );
          hdet_mc->SetBinContent( ebin, ctbin, bybin, std::get<0>(mc_response) );
          hdet_mc->SetBinError(   ebin, ctbin, bybin, std::get<1>(mc_response) );
  
        }
      }
    }
    
    cout << "NOTICE: Finished filling histograms" << endl;
    cout << "NOTICE: Time for filling hists " << (Double_t)timer.RealTime() << endl;
 
 
    TString output = Form("split_expected_evts_by_energy_%s.root", order_string.c_str());
    TFile fout(filefolder + output,"RECREATE");
    hdet_good_t->Write();
    hdet_good_s->Write();
    hdet_good_e->Write();
    hdet_showers->Write();
    hdet_mc->Write();
    fout.Close();
  
  }
}
