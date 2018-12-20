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

void under_over_flow_ct_flav(TString sum_file="../../data/ORCA_MC_summary_all_10Apr2018.root") { 
  SummaryParser sp(sum_file);
  
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
    h2elcc.push_back(new TH2D(Form("h2elcc_%i", i), Form("Angular resolution nu_e_CC_q_0.%i", i), n_bins, &ct_edges[0], n_bins, &ct_edges[0]));
    h2elnc.push_back(new TH2D(Form("h2elnc_%i", i), Form("Angular resolution nu_e_NC_q_0.%i", i), n_bins, &ct_edges[0], n_bins, &ct_edges[0]));
    h2mucc.push_back(new TH2D(Form("h2mucc_%i", i), Form("Angular resolution nu_m_CC_q_0.%i", i), n_bins, &ct_edges[0], n_bins, &ct_edges[0]));
    h2munc.push_back(new TH2D(Form("h2munc_%i", i), Form("Angular resolution nu_m_NC_q_0.%i", i), n_bins, &ct_edges[0], n_bins, &ct_edges[0]));
    h2tacc.push_back(new TH2D(Form("h2tacc_%i", i), Form("Angular resolution nu_t_CC_q_0.%i", i), n_bins, &ct_edges[0], n_bins, &ct_edges[0]));
    h2tanc.push_back(new TH2D(Form("h2tanc_%i", i), Form("Angular resolution nu_t_NC_q_0.%i", i), n_bins, &ct_edges[0], n_bins, &ct_edges[0]));
  } 

  Double_t q;
  Int_t good_track_and_shower = 0;
  Int_t bad_track = 0;
  Int_t bad_shower = 0;
  Int_t bad_track_and_shower = 0;
  Int_t total_events = 0;
  for (Int_t i = 0; i < sp.GetTree()->GetEntries(); i++) {
    total_events++;
    // Filter shyte
    bool good_tr = false;
    bool good_sh = false;
    sp.GetTree()->GetEntry(i);
    if ((sp.GetEvt()->Get_RDF_muon_score() > 0.05) or (sp.GetEvt()->Get_RDF_noise_score() > 0.5)) { continue; } // Filters for the events
    if ((sp.GetEvt()->Get_track_ql0() > 0.5) and (sp.GetEvt()->Get_track_ql1() > 0.5)) { good_tr = true; }
    if ((sp.GetEvt()->Get_shower_ql0() > 0.5) and (sp.GetEvt()->Get_shower_ql1() > 0.5)) { good_sh = true; }
    if ((not good_tr) and (not good_sh)) { bad_track_and_shower++; continue; }
    else if (not good_tr) { bad_track++; continue; }
    else if (not good_sh) { bad_shower++; continue; }
    good_track_and_shower++;

    // We order the histograms by quality: 0.0-0.1 --> 0, 0.1-0.2 --> 1, etc.
    q = sp.GetEvt()->Get_RDF_track_score();
    int index = (int)(TMath::Floor(q * 10));
    if (index == 10) { index = 9; } // A perfect track will get an index 10, which does not exist

    // Fill events
    if ((std::abs(sp.GetEvt()->Get_MC_type()) == 11) or (std::abs(sp.GetEvt()->Get_MC_type()) == 12)) { // Abs for anti-particles
      if (sp.GetEvt()->Get_MC_is_CC()) h2elcc[index]->Fill(-sp.GetEvt()->Get_MC_dir_z(), -sp.GetEvt()->Get_shower_dir_z());
      else                             h2elnc[index]->Fill(-sp.GetEvt()->Get_MC_dir_z(), -sp.GetEvt()->Get_shower_dir_z());
    }
    if ((std::abs(sp.GetEvt()->Get_MC_type()) == 13) or (std::abs(sp.GetEvt()->Get_MC_type()) == 14)) {
      if (sp.GetEvt()->Get_MC_is_CC()) h2mucc[index]->Fill(-sp.GetEvt()->Get_MC_dir_z(), -sp.GetEvt()->Get_track_dir_z());
      else                             h2munc[index]->Fill(-sp.GetEvt()->Get_MC_dir_z(), -sp.GetEvt()->Get_shower_dir_z());
    }
    if ((std::abs(sp.GetEvt()->Get_MC_type()) == 15) or (std::abs(sp.GetEvt()->Get_MC_type()) == 16)) {
      if (sp.GetEvt()->Get_MC_is_CC()) h2tacc[index]->Fill(-sp.GetEvt()->Get_MC_dir_z(), -sp.GetEvt()->Get_shower_dir_z());
      else                             h2tanc[index]->Fill(-sp.GetEvt()->Get_MC_dir_z(), -sp.GetEvt()->Get_shower_dir_z());
    }

  }
  cout << "Total events: " << total_events << endl;
  cout << "Good events : " << good_track_and_shower << endl;
  cout << "Bad tracks  : " << bad_track << endl;
  cout << "Bad showers : " << bad_shower << endl;
  cout << "Bad both    : " << bad_track_and_shower << endl;

  
  for (int i = 0; i < 10; i++) {
    std::tuple<Double_t, Double_t, Double_t, Double_t> over_under_elcc = GetOverFlow2D(h2elcc[i]);
    cout << "Under- and overflow for elcc at q" << i << " are X: " << std::get<0>(over_under_elcc) << " " << std::get<1>(over_under_elcc) 
                                                         << " Y: " << std::get<2>(over_under_elcc) << " " << std::get<3>(over_under_elcc) << endl;
  }
  for (int i = 0; i < 10; i++) {
    std::tuple<Double_t, Double_t, Double_t, Double_t> over_under_elcc = GetOverFlow2D(h2elnc[i]);
    cout << "Under- and overflow for elnc at q" << i << " are X: " << std::get<0>(over_under_elcc) << " " << std::get<1>(over_under_elcc) 
                                                         << " Y: " << std::get<2>(over_under_elcc) << " " << std::get<3>(over_under_elcc) << endl;
  }
  for (int i = 0; i < 10; i++) {
    std::tuple<Double_t, Double_t, Double_t, Double_t> over_under_elcc = GetOverFlow2D(h2mucc[i]);
    cout << "Under- and overflow for mucc at q" << i << " are X: " << std::get<0>(over_under_elcc) << " " << std::get<1>(over_under_elcc) 
                                                         << " Y: " << std::get<2>(over_under_elcc) << " " << std::get<3>(over_under_elcc) << endl;
  }
  for (int i = 0; i < 10; i++) {
    std::tuple<Double_t, Double_t, Double_t, Double_t> over_under_elcc = GetOverFlow2D(h2munc[i]);
    cout << "Under- and overflow for munc at q" << i << " are X: " << std::get<0>(over_under_elcc) << " " << std::get<1>(over_under_elcc) 
                                                         << " Y: " << std::get<2>(over_under_elcc) << " " << std::get<3>(over_under_elcc) << endl;
  }
  for (int i = 0; i < 10; i++) {
    std::tuple<Double_t, Double_t, Double_t, Double_t> over_under_elcc = GetOverFlow2D(h2tacc[i]);
    cout << "Under- and overflow for tacc at q" << i << " are X: " << std::get<0>(over_under_elcc) << " " << std::get<1>(over_under_elcc) 
                                                         << " Y: " << std::get<2>(over_under_elcc) << " " << std::get<3>(over_under_elcc) << endl;
  }
  for (int i = 0; i < 10; i++) {
    std::tuple<Double_t, Double_t, Double_t, Double_t> over_under_elcc = GetOverFlow2D(h2tanc[i]);
    cout << "Under- and overflow for tanc at q" << i << " are X: " << std::get<0>(over_under_elcc) << " " << std::get<1>(over_under_elcc) 
                                                         << " Y: " << std::get<2>(over_under_elcc) << " " << std::get<3>(over_under_elcc) << endl;
  }
}
