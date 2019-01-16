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
#include "NMHUtils.h"
#include "SummaryParser.h"
#include "SummaryEvent.h"
#include "HelperFunctions.C"

#include <iostream>
using namespace std;


void Plot_Resolution_Energy_Flav(TString summary_file=(TString)getenv("NMHDIR") + "/data/ORCA_MC_summary_all_10Apr2018.root") {
  SummaryParser sp(sum_file);
  
  bool plot = false;
  bool print = false;
  int n_bins = 40;
  std::vector<Double_t> e_edges  = NMHUtils::GetLogBins(n_bins, 1, 100);
  std::vector<Double_t> ct_edges = NMHUtils::GetBins(n_bins, -1, 1);

  // Split up scheme: DR are split into PID based bins
  std::vector<TH2D*> h2elcc;
  std::vector<TH2D*> h2elnc;
  std::vector<TH2D*> h2mucc;
  std::vector<TH2D*> h2munc;
  std::vector<TH2D*> h2tacc;
  std::vector<TH2D*> h2tanc;
  for (int i = 0; i < 10; i++){
    h2elcc.push_back(new TH2D(Form("h2elcc_%i", i), Form("Energy resolution nu_e_CC q_0.%i", i), n_bins, &e_edges[0], n_bins, &e_edges[0]));
    h2elnc.push_back(new TH2D(Form("h2elnc_%i", i), Form("Energy resolution nu_e_NC q_0.%i", i), n_bins, &e_edges[0], n_bins, &e_edges[0]));
    h2mucc.push_back(new TH2D(Form("h2mucc_%i", i), Form("Energy resolution nu_m_CC q_0.%i", i), n_bins, &e_edges[0], n_bins, &e_edges[0]));
    h2munc.push_back(new TH2D(Form("h2munc_%i", i), Form("Energy resolution nu_m_NC q_0.%i", i), n_bins, &e_edges[0], n_bins, &e_edges[0]));
    h2tacc.push_back(new TH2D(Form("h2tacc_%i", i), Form("Energy resolution nu_t_CC q_0.%i", i), n_bins, &e_edges[0], n_bins, &e_edges[0]));
    h2tanc.push_back(new TH2D(Form("h2tanc_%i", i), Form("Energy resolution nu_t_NC q_0.%i", i), n_bins, &e_edges[0], n_bins, &e_edges[0]));
  }

  Double_t q;
  Int_t good_track_and_shower = 0;
  Int_t bad_track = 0;
  Int_t bad_shower = 0;
  Int_t bad_track_and_shower = 0;
  Int_t rejected = 0;
  Int_t total_events = 0;
  for (Int_t i = 0; i < sp.GetTree()->GetEntries(); i++) {
    total_events++;
    // Filters
    bool good_tr = false;
    bool good_sh = false;

    SummaryEvent *evt = sp.GetEvt(i);
    if ((evt->Get_RDF_muon_score() > 0.05) or (evt->Get_RDF_noise_score() > 0.5)) { rejected++; continue; } // Filters for the events
    if ((evt->Get_track_ql0() > 0.5) and (evt->Get_track_ql1() > 0.5)) { good_tr = true; }
    if ((evt->Get_shower_ql0() > 0.5) and (evt->Get_shower_ql1() > 0.5)) { good_sh = true; }
    if ((not good_tr) and (not good_sh)) { bad_track_and_shower++; continue; }
    else if (not good_tr) { bad_track++; }
    else if (not good_sh) { bad_shower++; }
    good_track_and_shower++;

    // We order the histograms by quality: 0.0-0.1 --> 0, 0.1-0.2 --> 1, etc.
    q = evt->Get_RDF_track_score();
    int index = (int)(TMath::Floor(q * 10));
    if (index == 10) { index = 9; } // A perfect track will get an index 10, which does not exist

    // Fill events PID binned
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
    cout << "Rejected evt: " << rejected << endl;
  }

  if (plot) {
    TCanvas *c1 = new TCanvas("c1", "c1", 1800, 500);
    TCanvas *c2 = new TCanvas("c2", "c2", 1800, 500);
    TCanvas *c3 = new TCanvas("c3", "c3", 1800, 500); 
    TCanvas *c4 = new TCanvas("c4", "c4", 1800, 500); 
    TCanvas *c5 = new TCanvas("c5", "c5", 1800, 500); 
    TCanvas *c6 = new TCanvas("c6", "c6", 1800, 500); 
    gStyle->SetPalette(kBird);
    c1->Divide(5,2);
    c2->Divide(5,2);
    c3->Divide(5,2);
    c4->Divide(5,2);
    c5->Divide(5,2);
    c6->Divide(5,2);
    for (int i = 0; i < 10; i++) {
      c1->cd(i+1);
      GetNormalizedSlicesY(h2elcc[i]);
      h2elcc[i]->Draw("colz");
      h2elcc[i]->GetXaxis()->SetTitle("E_{true} [GeV]");
      h2elcc[i]->GetYaxis()->SetTitle("E_{reco} [GeV]");
      h2elcc[i]->GetXaxis()->SetRangeUser(3,100);
      h2elcc[i]->GetYaxis()->SetRangeUser(3,100);
      c1->cd(i+1)->SetLogx();
      c1->cd(i+1)->SetLogy();
      c1->cd(i+1)->SetLogz();

      c2->cd(i+1);
      GetNormalizedSlicesY(h2elnc[i]);
      h2elnc[i]->Draw("colz");
      h2elnc[i]->GetXaxis()->SetTitle("E_{true} [GeV]");
      h2elnc[i]->GetYaxis()->SetTitle("E_{reco} [GeV]");
      h2elnc[i]->GetXaxis()->SetRangeUser(3,100);
      h2elnc[i]->GetYaxis()->SetRangeUser(3,100);
      c2->cd(i+1)->SetLogx();
      c2->cd(i+1)->SetLogy();
      c2->cd(i+1)->SetLogz();

      c3->cd(i+1);
      GetNormalizedSlicesY(h2mucc[i]);
      h2mucc[i]->Draw("colz");
      h2mucc[i]->GetXaxis()->SetTitle("E_{true} [GeV]");
      h2mucc[i]->GetYaxis()->SetTitle("E_{reco} [GeV]");
      h2mucc[i]->GetXaxis()->SetRangeUser(3,100);
      h2mucc[i]->GetYaxis()->SetRangeUser(3,100);
      c3->cd(i+1)->SetLogx();
      c3->cd(i+1)->SetLogy();
      c3->cd(i+1)->SetLogz();

      c4->cd(i+1);
      GetNormalizedSlicesY(h2munc[i]);
      h2munc[i]->Draw("colz");
      h2munc[i]->GetXaxis()->SetTitle("E_{true} [GeV]");
      h2munc[i]->GetYaxis()->SetTitle("E_{reco} [GeV]");
      h2munc[i]->GetXaxis()->SetRangeUser(3,100);
      h2munc[i]->GetYaxis()->SetRangeUser(3,100);
      c4->cd(i+1)->SetLogx();
      c4->cd(i+1)->SetLogy();

      c5->cd(i+1);
      GetNormalizedSlicesY(h2tacc[i]);
      h2tacc[i]->Draw("colz");
      h2tacc[i]->GetXaxis()->SetTitle("E_{true} [GeV]");
      h2tacc[i]->GetYaxis()->SetTitle("E_{reco} [GeV]");
      h2tacc[i]->GetXaxis()->SetRangeUser(3,100);
      h2tacc[i]->GetYaxis()->SetRangeUser(3,100);
      c5->cd(i+1)->SetLogx();
      c5->cd(i+1)->SetLogy();
      c5->cd(i+1)->SetLogz();

      c6->cd(i+1);
      GetNormalizedSlicesY(h2tanc[i]);
      h2tanc[i]->Draw("colz");
      h2tanc[i]->GetXaxis()->SetTitle("E_{true} [GeV]");
      h2tanc[i]->GetYaxis()->SetTitle("E_{reco} [GeV]");
      h2tanc[i]->GetXaxis()->SetRangeUser(3,100);
      h2tanc[i]->GetYaxis()->SetRangeUser(3,100);
      c6->cd(i+1)->SetLogx();
      c6->cd(i+1)->SetLogy();
    }

    c1->SaveAs("./pid_detres/energy_resolution_plots_most_evts_used/energy_resolution_elcc_per_q.pdf");
    c2->SaveAs("./pid_detres/energy_resolution_plots_most_evts_used/energy_resolution_elnc_per_q.pdf");
    c3->SaveAs("./pid_detres/energy_resolution_plots_most_evts_used/energy_resolution_mucc_per_q.pdf");
    c4->SaveAs("./pid_detres/energy_resolution_plots_most_evts_used/energy_resolution_munc_per_q.pdf");
    c5->SaveAs("./pid_detres/energy_resolution_plots_most_evts_used/energy_resolution_tacc_per_q.pdf");
    c6->SaveAs("./pid_detres/energy_resolution_plots_most_evts_used/energy_resolution_tanc_per_q.pdf");
  }
}
