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

void expectation() {

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

  SummaryParser sp("../data/ORCA_MC_summary_all_10Apr2018.root");
  for (Int_t i = 0; i < sp.GetTree()->GetEntries(); i++) {
    sp.GetTree()->GetEntry(i);
    SummaryEvent *evt = sp.GetEvt();
    drt.Fill(evt);
    drs.Fill(evt);
  }

  drt.WriteToFile("track_response_timing.root");
  drs.WriteToFile("shower_response_timing.root");

  // drt.ReadFromFile("track_response_timing.root");
  // drs.ReadFromFile("shower_response_timing.root");
  
  cout << "NOTICE: Finished filling response" << endl;

  //----------------------------------------------------------
  // fill a detected histogram using the fit function
  //----------------------------------------------------------

  FitFunction tfitf(&drt, 3, "../data/eff_mass/EffMhists_elec_CC.root", "../data/eff_mass/EffMhists_muon_CC.root", "../data/eff_mass/EffMhists_tau_CC.root", "../data/eff_mass/EffMhists_elec_NC.root");

  FitFunction sfitf(&drs, 3, "../data/eff_mass/EffMhists_elec_CC.root", "../data/eff_mass/EffMhists_muon_CC.root", "../data/eff_mass/EffMhists_tau_CC.root", "../data/eff_mass/EffMhists_elec_NC.root");

  TH3D *hdet_tracks  = (TH3D*)drt.GetHist3D()->Clone("detected_tracks");
  TH3D *hdet_showers = (TH3D*)drs.GetHist3D()->Clone("detected_showers");
  hdet_tracks->Reset();
  hdet_showers->Reset();

  double p[] = {0.297, 0.0215, 0.425, 1.38, 7.37e-5, 2.56e-3};

  TStopwatch timer;

  for (Int_t ebin = 1; ebin <= hdet_tracks->GetXaxis()->GetNbins(); ebin++) {
    for (Int_t ctbin = 1; ctbin <= hdet_tracks->GetYaxis()->GetNbins(); ctbin++) {
      for (Int_t bybin = 1; bybin <= hdet_tracks->GetZaxis()->GetNbins(); bybin++) {

	Double_t E  = hdet_tracks->GetXaxis()->GetBinCenter( ebin );
	Double_t ct = hdet_tracks->GetYaxis()->GetBinCenter( ctbin );
	Double_t by = hdet_tracks->GetZaxis()->GetBinCenter( bybin );

	double x[] = {E, ct, by};

	hdet_tracks->SetBinContent( ebin, ctbin, bybin, tfitf.operator()(x, p) );
	hdet_showers->SetBinContent( ebin, ctbin, bybin, sfitf.operator()(x, p) );
      }
    }
  }
  
  cout << "NOTICE: Finished filling histograms" << endl;
  cout << "NOTICE: Time for filling hists " << (Double_t)timer.RealTime() << endl;
  

  TFile fout("expectation_values.root","RECREATE");
  hdet_tracks->Write();
  hdet_showers->Write();
  fout.Close();

}
