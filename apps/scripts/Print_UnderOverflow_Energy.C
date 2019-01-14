#include "TSystem.h"
#include "TROOT.h"
#include "TDatime.h"
#include "TH2.h"
#include "TH3.h"
#include "TFile.h"
#include "TStopwatch.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"

#include "DetResponse.h"
#include "FitFunction.h"
#include "NMHUtils.h"
#include "SummaryParser.h"
#include "SummaryEvent.h"
#include "HelperFunctions.C"

#include <iostream>
using namespace std;

void Print_UnderOverflow_Energy(TString summary_file=(TString)getenv("NMHDIR") + "/data/ORCA_MC_summary_all_10Apr2018.root") { 
  SummaryParser sp(sum_file);
  
  bool print = false;
  int n_bins = 40;
  std::vector<Double_t> e_edges  = NMHUtils::GetLogBins(n_bins, 1, 100);
  std::vector<Double_t> ct_edges = NMHUtils::GetBins(n_bins, -1, 0);

  std::vector<TH2D*> h2elcc;
  std::vector<TH2D*> h2elnc;
  std::vector<TH2D*> h2mucc;
  std::vector<TH2D*> h2munc;
  std::vector<TH2D*> h2tacc;
  std::vector<TH2D*> h2tanc;
  for (int i = 0; i < 10; i++){
    h2elcc.push_back(new TH2D(Form("h2elcc_%i", i), Form("Energy resolution nu_e_CC_q_0.%i", i), n_bins, &e_edges[0], n_bins, &e_edges[0]));
    h2elnc.push_back(new TH2D(Form("h2elnc_%i", i), Form("Energy resolution nu_e_NC_q_0.%i", i), n_bins, &e_edges[0], n_bins, &e_edges[0]));
    h2mucc.push_back(new TH2D(Form("h2mucc_%i", i), Form("Energy resolution nu_m_CC_q_0.%i", i), n_bins, &e_edges[0], n_bins, &e_edges[0]));
    h2munc.push_back(new TH2D(Form("h2munc_%i", i), Form("Energy resolution nu_m_NC_q_0.%i", i), n_bins, &e_edges[0], n_bins, &e_edges[0]));
    h2tacc.push_back(new TH2D(Form("h2tacc_%i", i), Form("Energy resolution nu_t_CC_q_0.%i", i), n_bins, &e_edges[0], n_bins, &e_edges[0]));
    h2tanc.push_back(new TH2D(Form("h2tanc_%i", i), Form("Energy resolution nu_t_NC_q_0.%i", i), n_bins, &e_edges[0], n_bins, &e_edges[0]));
  } 

  Double_t q;
  Int_t good_track_and_shower = 0;
  Int_t bad_track = 0;
  Int_t bad_shower = 0;
  Int_t bad_track_and_shower = 0;
  Int_t total_events = 0;
  for (Int_t i = 0; i < sp.GetTree()->GetEntries(); i++) {
    total_events++;

    SummaryEvent *evt = sp.GetEvt(i);
    // Filters, bad implementation of EventFilters for this purpose.
    bool good_tr = false;
    bool good_sh = false;
    // Watch out! In showers the RDF noise is different than in tracks
    if ((evt->Get_RDF_muon_score() > 0.05) or (evt->Get_RDF_noise_score() > 0.5)) { continue; } // Filters for the events
    if ((evt->Get_track_ql0() > 0.5) and (evt->Get_track_ql1() > 0.5)) { good_tr = true; }
    if ((evt->Get_shower_ql0() > 0.5) and (evt->Get_shower_ql1() > 0.5)) { good_sh = true; }
    if ((not good_tr) and (not good_sh)) { bad_track_and_shower++; continue; }
    else if (not good_tr) { bad_track++; continue; }
    else if (not good_sh) { bad_shower++; continue; }
    good_track_and_shower++; // If it passes to here, the reconstruction is good for both track and shower

    // We order the histograms by quality: 0.0-0.1 --> 0, 0.1-0.2 --> 1, etc.
    q = evt->Get_RDF_track_score();
    int index = (int)(TMath::Floor(q * 10));
    if (index == 10) { index = 9; } // A perfect track will get an index 10, which does not exist

    // Fill events
    if ((std::abs(evt->Get_MC_type()) == 11) or (std::abs(evt->Get_MC_type()) == 12)) { // Abs for anti-particles
      if (evt->Get_MC_is_CC()) h2elcc[index]->Fill(evt->Get_MC_energy(), evt->Get_shower_energy());
      else                     h2elnc[index]->Fill(evt->Get_MC_energy(), evt->Get_shower_energy());
    }
    if ((std::abs(evt->Get_MC_type()) == 13) or (std::abs(evt->Get_MC_type()) == 14)) {
      if (evt->Get_MC_is_CC()) h2mucc[index]->Fill(evt->Get_MC_energy(), evt->Get_track_energy());
      else                     h2munc[index]->Fill(evt->Get_MC_energy(), evt->Get_shower_energy());
    }
    if ((std::abs(evt->Get_MC_type()) == 15) or (std::abs(evt->Get_MC_type()) == 16)) {
      if (evt->Get_MC_is_CC()) h2tacc[index]->Fill(evt->Get_MC_energy(), evt->Get_shower_energy());
      else                     h2tanc[index]->Fill(evt->Get_MC_energy(), evt->Get_shower_energy());
    }
  }
  if (print) {
    cout << "Total events: " << total_events << endl;
    cout << "Good events : " << good_track_and_shower << endl;
    cout << "Bad tracks  : " << bad_track << endl;
    cout << "Bad showers : " << bad_shower << endl;
    cout << "Bad both    : " << bad_track_and_shower << endl;
  }
  
  for (int i = 0; i < 10; i++) {
    std::tuple<Double_t, Double_t, Double_t, Double_t> under_over_flow = GetOverFlow2D(h2elcc[i]);
    cout << "Under- and overflow for elcc at q" << i << " are X: " << std::get<0>(under_over_flow) << " " << std::get<1>(under_over_flow) 
                                                         << " Y: " << std::get<2>(under_over_flow) << " " << std::get<3>(under_over_flow) << endl;
  }
  for (int i = 0; i < 10; i++) {
    std::tuple<Double_t, Double_t, Double_t, Double_t> under_over_flow = GetOverFlow2D(h2elnc[i]);
    cout << "Under- and overflow for elnc at q" << i << " are X: " << std::get<0>(under_over_flow) << " " << std::get<1>(under_over_flow) 
                                                         << " Y: " << std::get<2>(under_over_flow) << " " << std::get<3>(under_over_flow) << endl;
  }
  for (int i = 0; i < 10; i++) {
    std::tuple<Double_t, Double_t, Double_t, Double_t> under_over_flow = GetOverFlow2D(h2mucc[i]);
    cout << "Under- and overflow for mucc at q" << i << " are X: " << std::get<0>(under_over_flow) << " " << std::get<1>(under_over_flow) 
                                                         << " Y: " << std::get<2>(under_over_flow) << " " << std::get<3>(under_over_flow) << endl;
  }
  for (int i = 0; i < 10; i++) {
    std::tuple<Double_t, Double_t, Double_t, Double_t> under_over_flow = GetOverFlow2D(h2munc[i]);
    cout << "Under- and overflow for munc at q" << i << " are X: " << std::get<0>(under_over_flow) << " " << std::get<1>(under_over_flow) 
                                                         << " Y: " << std::get<2>(under_over_flow) << " " << std::get<3>(under_over_flow) << endl;
  }
  for (int i = 0; i < 10; i++) {
    std::tuple<Double_t, Double_t, Double_t, Double_t> under_over_flow = GetOverFlow2D(h2tacc[i]);
    cout << "Under- and overflow for tacc at q" << i << " are X: " << std::get<0>(under_over_flow) << " " << std::get<1>(under_over_flow) 
                                                         << " Y: " << std::get<2>(under_over_flow) << " " << std::get<3>(under_over_flow) << endl;
  }
  for (int i = 0; i < 10; i++) {
    std::tuple<Double_t, Double_t, Double_t, Double_t> under_over_flow = GetOverFlow2D(h2tanc[i]);
    cout << "Under- and overflow for tanc at q" << i << " are X: " << std::get<0>(under_over_flow) << " " << std::get<1>(under_over_flow) 
                                                         << " Y: " << std::get<2>(under_over_flow) << " " << std::get<3>(under_over_flow) << endl;
  }
}
