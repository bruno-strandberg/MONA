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

void SplitEnergyPID() {

  TString filefolder = "./energy_detres/pid_bins_10/";
  const int N_PID_CLASSES = 10;
  Double_t PID_step = 1 / float(N_PID_CLASSES);

  gROOT->ProcessLine(".L FitFunction.C+");
  gSystem->Load("$OSCPROBDIR/libOscProb.so");

  //----------------------------------------------------------
  // detector response for tracks and showers
  //----------------------------------------------------------
  // Good tracks only, no overlap with good showers
  // gt = good track, gs = good shower, ge = good event
  std::vector<DetResponse> drgts; // good tracks in track channel
  std::vector<DetResponse> drgss; // good showers in track channel
  std::vector<DetResponse> drges; // good events in track channel
  std::vector<DetResponse> drss; // showers
  std::vector<DetResponse> drms; // mc

  for (int i = 0; i < N_PID_CLASSES; i++) {
    std::function<bool(double, double)> comparison_operator;
    if (i == 0) { comparison_operator = std::greater_equal<double>(); // The first bin needs to include the lower limit.
    } else { comparison_operator = std::greater<double>(); }
    DetResponse drgt(DetResponse::hybridE, TString::Format("track_response_gt_%.2f", PID_step * i), 40, 1, 100, 40, -1, 1, 1, 0, 1);
    drgt.AddCut( &SummaryEvent::Get_track_ql0       , std::greater<double>()   ,  0.5, true );
    drgt.AddCut( &SummaryEvent::Get_track_ql1       , std::greater<double>()   ,  0.5, true );
    drgt.AddCut( &SummaryEvent::Get_shower_ql0      , std::less<double>()      ,  0.5, true );
    drgt.AddCut( &SummaryEvent::Get_shower_ql1      , std::less<double>()      ,  0.5, true );
    drgt.AddCut( &SummaryEvent::Get_RDF_track_score , comparison_operator      , PID_step * i    , true );
    drgt.AddCut( &SummaryEvent::Get_RDF_track_score , std::less_equal<double>(), PID_step * (i+1), true );
    drgt.AddCut( &SummaryEvent::Get_RDF_muon_score  , std::less_equal<double>(), 0.05, true );
    drgt.AddCut( &SummaryEvent::Get_RDF_noise_score , std::less_equal<double>(), 0.18, true );
    drgts.push_back(drgt);

    // Good showers only, no overlap with good tracks
    DetResponse drgs(DetResponse::hybridE, TString::Format("track_response_gs_%.2f", PID_step * i), 40, 1, 100, 40, -1, 1, 1, 0, 1);
    drgs.AddCut( &SummaryEvent::Get_track_ql0       , std::greater<double>()   ,  0.5, true ); // This would cut almost all events, since they all have fTrack_ql0 > 0.5
    drgs.AddCut( &SummaryEvent::Get_track_ql1       , std::less<double>()      ,  0.5, true ); // Moreover: the numbers are consistent with resolution_plot_flav_complementary_events.C
    drgs.AddCut( &SummaryEvent::Get_shower_ql0      , std::greater<double>()   ,  0.5, true );
    drgs.AddCut( &SummaryEvent::Get_shower_ql1      , std::greater<double>()   ,  0.5, true );
    drgs.AddCut( &SummaryEvent::Get_RDF_track_score , comparison_operator      , PID_step * i    , true );
    drgs.AddCut( &SummaryEvent::Get_RDF_track_score , std::less_equal<double>(), PID_step * (i+1), true );
    drgs.AddCut( &SummaryEvent::Get_RDF_muon_score  , std::less_equal<double>(), 0.05, true );
    drgs.AddCut( &SummaryEvent::Get_RDF_noise_score , std::less_equal<double>(), 0.18, true );
    drgss.push_back(drgs);

    // Good events only, no overlap with bad tracks or bad showers
    DetResponse drge(DetResponse::hybridE, TString::Format("track_response_ge_%.2f", PID_step * i), 40, 1, 100, 40, -1, 1, 1, 0, 1);
    drge.AddCut( &SummaryEvent::Get_track_ql0       , std::greater<double>()   ,  0.5, true );
    drge.AddCut( &SummaryEvent::Get_track_ql1       , std::greater<double>()   ,  0.5, true );
    drge.AddCut( &SummaryEvent::Get_shower_ql0      , std::greater<double>()   ,  0.5, true );
    drge.AddCut( &SummaryEvent::Get_shower_ql1      , std::greater<double>()   ,  0.5, true );
    drge.AddCut( &SummaryEvent::Get_RDF_track_score , comparison_operator      , PID_step * i    , true );
    drge.AddCut( &SummaryEvent::Get_RDF_track_score , std::less_equal<double>(), PID_step * (i+1), true );
    drge.AddCut( &SummaryEvent::Get_RDF_muon_score  , std::less_equal<double>(), 0.05, true );
    drge.AddCut( &SummaryEvent::Get_RDF_noise_score , std::less_equal<double>(), 0.18, true );
    drges.push_back(drge);

    DetResponse drs(DetResponse::hybridE, TString::Format("shower_response_%.2f", PID_step * i), 40, 1, 100, 40, -1, 1, 1, 0, 1);
    drs.AddCut( &SummaryEvent::Get_shower_ql0     , std::greater<double>()   ,  0.5, true );
    drs.AddCut( &SummaryEvent::Get_shower_ql1     , std::greater<double>()   ,  0.5, true );
    drs.AddCut( &SummaryEvent::Get_RDF_track_score , comparison_operator      , PID_step * i    , true );
    drs.AddCut( &SummaryEvent::Get_RDF_track_score , std::less_equal<double>(), PID_step * (i+1), true );
    drs.AddCut( &SummaryEvent::Get_RDF_muon_score , std::less_equal<double>(), 0.05, true );
    drs.AddCut( &SummaryEvent::Get_RDF_noise_score, std::less_equal<double>(),  0.5, true );
    drss.push_back(drs);

    DetResponse drm(DetResponse::hybridE, TString::Format("mc_response_%.2f", PID_step * i), 40, 1, 100, 40, -1, 1, 1, 0, 1);
    drm.AddCut( &SummaryEvent::Get_RDF_track_score , comparison_operator      , PID_step * i    , true );
    drm.AddCut( &SummaryEvent::Get_RDF_track_score , std::less_equal<double>(), PID_step * (i+1), true );
    drm.AddCut( &SummaryEvent::Get_RDF_muon_score , std::less_equal<double>(), 0.05, true );
    drm.AddCut( &SummaryEvent::Get_RDF_noise_score, std::less_equal<double>(),  0.5, true );
    drms.push_back(drm);
  }

  SummaryParser sp("../data/ORCA_MC_summary_all_10Apr2018.root");
  bool writeFiles = true;
  for (Int_t i = 0; i < sp.GetTree()->GetEntries(); i++) {
    sp.GetTree()->GetEntry(i);
    SummaryEvent *evt = sp.GetEvt();
    if (evt->Get_MC_is_neutrino() < 0.5) continue;
    for (int i = 0; i < N_PID_CLASSES; i++) {
      drgts[i].Fill(evt);
      drgss[i].Fill(evt);
      drges[i].Fill(evt);
      drss[i].Fill(evt);
      drms[i].Fill(evt);
    }
  }

  for (int i = 0; i < N_PID_CLASSES; i++) { 
    if (writeFiles) { // This flag is to not accidently overwrite with empty files
      drgts[i].WriteToFile(filefolder + TString::Format("track_response_gt_%.2f.root", PID_step * i));
      drgss[i].WriteToFile(filefolder + TString::Format("track_response_gs_%.2f.root", PID_step * i));
      drges[i].WriteToFile(filefolder + TString::Format("track_response_ge_%.2f.root", PID_step * i));
      drss[i].WriteToFile(filefolder  + TString::Format("shower_response_%.2f.root", PID_step * i));
      drms[i].WriteToFile(filefolder  + TString::Format("mc_response_%.2f.root", PID_step * i));
    }
  }

//  for (int i = 0; i < N_PID_CLASSES; i++) {
//    drgts[i].ReadFromFile(filefolder + TString::Format("track_response_gt_%.2f.root", PID_step * i));
//    drgss[i].ReadFromFile(filefolder + TString::Format("track_response_gs_%.2f.root", PID_step * i));
//    drges[i].ReadFromFile(filefolder + TString::Format("track_response_ge_%.2f.root", PID_step * i));
//    drss[i].ReadFromFile(filefolder  + TString::Format("shower_response_%.2f.root", PID_step * i));
//    drms[i].ReadFromFile(filefolder  + TString::Format("mc_response_%.2f.root", PID_step * i));
//  }

  cout << "NOTICE: Finished filling response" << endl;

  //----------------------------------------------------------
  // fill a detected histogram using the fit function
  //----------------------------------------------------------

  for (int i = 0; i < N_PID_CLASSES; i++){
    FitFunction gtfitf(&drgts[i], 3, "../data/eff_mass/EffMhists_elec_CC.root", "../data/eff_mass/EffMhists_muon_CC.root", "../data/eff_mass/EffMhists_tau_CC.root", "../data/eff_mass/EffMhists_elec_NC.root");
    FitFunction gsfitf(&drgss[i], 3, "../data/eff_mass/EffMhists_elec_CC.root", "../data/eff_mass/EffMhists_muon_CC.root", "../data/eff_mass/EffMhists_tau_CC.root", "../data/eff_mass/EffMhists_elec_NC.root");
    FitFunction gefitf(&drges[i], 3, "../data/eff_mass/EffMhists_elec_CC.root", "../data/eff_mass/EffMhists_muon_CC.root", "../data/eff_mass/EffMhists_tau_CC.root", "../data/eff_mass/EffMhists_elec_NC.root");
    FitFunction sfitf(&drss[i], 3, "../data/eff_mass/EffMhists_elec_CC.root", "../data/eff_mass/EffMhists_muon_CC.root", "../data/eff_mass/EffMhists_tau_CC.root", "../data/eff_mass/EffMhists_elec_NC.root");
    FitFunction mfitf(&drms[i], 3, "../data/eff_mass/EffMhists_elec_CC.root", "../data/eff_mass/EffMhists_muon_CC.root", "../data/eff_mass/EffMhists_tau_CC.root", "../data/eff_mass/EffMhists_elec_NC.root");

    TH3D *hdet_good_t  = (TH3D*)drgts[i].GetHist3D()->Clone("detected_tracks_gt"); // Track channel, using good tracks, etc.
    TH3D *hdet_good_s  = (TH3D*)drgss[i].GetHist3D()->Clone("detected_tracks_gs");
    TH3D *hdet_good_e  = (TH3D*)drges[i].GetHist3D()->Clone("detected_tracks_ge");
    TH3D *hdet_showers = (TH3D*)drss[i].GetHist3D()->Clone("detected_showers");
    TH3D *hdet_mc      = (TH3D*)drms[i].GetHist3D()->Clone("detected_mc");
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
 
 
      TString output = Form("split_energy_%s_%.2f.root", order_string.c_str(), PID_step * i);
      TFile fout(filefolder + output,"RECREATE");
      hdet_good_t->Write();
      hdet_good_s->Write();
      hdet_good_e->Write();
      hdet_showers->Write();
      hdet_mc->Write();
      fout.Close();
    }
  }
}
