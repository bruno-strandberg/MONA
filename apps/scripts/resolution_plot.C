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

void resolution_plot(TString sum_file="../../data/ORCA_MC_summary_all_10Apr2018.root") { 

  SummaryParser sp(sum_file);
  
  int n_bins = 40;
  std::vector<Double_t> e_edges  = NMHUtils::GetLogBins(n_bins, 1, 100);
  std::vector<Double_t> ct_edges = NMHUtils::GetBins(n_bins, -1, 1);

  // Default scheme: all events
  TH2D* ht_def = new TH2D("ht_def", "Reconstructed events", n_bins, &e_edges[0], n_bins, &e_edges[0]);
  TH2D* hs_def = new TH2D("hs_def", "Reconstructed events", n_bins, &e_edges[0], n_bins, &e_edges[0]);

  std::vector<TH1D*> ht;
  std::vector<TH1D*> hs;
  std::vector<TH1D*> hm;
  std::vector<TH2D*> h2t;
  std::vector<TH2D*> h2s;
  std::vector<TH2D*> h2elcc;
  std::vector<TH2D*> h2elnc;
  std::vector<TH2D*> h2mucc;
  std::vector<TH2D*> h2munc;

  for (int i = 0; i < 10; i++){
    ht.push_back(new TH1D(Form("ht_%i", i), Form("Reconstructed events_q_0.%i", i), n_bins, &e_edges[0]));
    hs.push_back(new TH1D(Form("hs_%i", i), Form("Reconstructed events_q_0.%i", i), n_bins, &e_edges[0]));
    hm.push_back(new TH1D(Form("hm_%i", i), Form("MC events_q_0.%i", i), n_bins, &e_edges[0]));
    h2t.push_back(new TH2D(Form("h2t_%i", i), Form("Energy resolution tracks_q_0.%i", i), n_bins, &e_edges[0], n_bins, &e_edges[0]));
    h2s.push_back(new TH2D(Form("h2s_%i", i), Form("Energy resolution showers_q_0.%i", i), n_bins, &e_edges[0], n_bins, &e_edges[0]));
    h2elcc.push_back(new TH2D(Form("h2elcc_%i", i), Form("Energy resolution nu_e_CC_q_0.%i", i), n_bins, &e_edges[0], n_bins, &e_edges[0]));
    h2elnc.push_back(new TH2D(Form("h2elnc_%i", i), Form("Energy resolution nu_e_NC_q_0.%i", i), n_bins, &e_edges[0], n_bins, &e_edges[0]));
    h2mucc.push_back(new TH2D(Form("h2mucc_%i", i), Form("Energy resolution nu_m_CC_q_0.%i", i), n_bins, &e_edges[0], n_bins, &e_edges[0]));
    h2munc.push_back(new TH2D(Form("h2munc_%i", i), Form("Energy resolution nu_m_NC_q_0.%i", i), n_bins, &e_edges[0], n_bins, &e_edges[0]));
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
    if ((sp.GetEvt()->Get_RDF_muon_score() > 0.5) or (sp.GetEvt()->Get_RDF_noise_score() > 0.5)) { continue; } // Filters for the events
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
    ht[index]->Fill(sp.GetEvt()->Get_track_energy());
    hs[index]->Fill(sp.GetEvt()->Get_shower_energy());
    hm[index]->Fill(sp.GetEvt()->Get_MC_energy());
    if ((std::abs(sp.GetEvt()->Get_MC_type()) == 11) or (std::abs(sp.GetEvt()->Get_MC_type() == 12))) { // Abs for anti-particles
      if (sp.GetEvt()->Get_MC_is_CC()) h2elcc[index]->Fill(sp.GetEvt()->Get_MC_energy(), sp.GetEvt()->Get_track_energy());
      else                             h2elnc[index]->Fill(sp.GetEvt()->Get_MC_energy(), sp.GetEvt()->Get_track_energy());
    }
    if ((std::abs(sp.GetEvt()->Get_MC_type() == 13)) or (std::abs(sp.GetEvt()->Get_MC_type() == 14))) {
      if (sp.GetEvt()->Get_MC_is_CC()) h2mucc[index]->Fill(sp.GetEvt()->Get_MC_energy(), sp.GetEvt()->Get_track_energy());
      else                             h2munc[index]->Fill(sp.GetEvt()->Get_MC_energy(), sp.GetEvt()->Get_track_energy());
    }
    h2t[index]->Fill(sp.GetEvt()->Get_MC_energy(), sp.GetEvt()->Get_track_energy());
    h2s[index]->Fill(sp.GetEvt()->Get_MC_energy(), sp.GetEvt()->Get_shower_energy());

    if (q > 0.6) { ht_def->Fill(sp.GetEvt()->Get_MC_energy(), sp.GetEvt()->Get_track_energy()); }
    if (q < 0.6) { hs_def->Fill(sp.GetEvt()->Get_MC_energy(), sp.GetEvt()->Get_shower_energy()); }
  }
  cout << "Total events: " << total_events << endl;
  cout << "Good events : " << good_track_and_shower << endl;
  cout << "Bad tracks  : " << bad_track << endl;
  cout << "Bad showers : " << bad_shower << endl;
  cout << "Bad both    : " << bad_track_and_shower << endl;

  TCanvas *c1 = new TCanvas("c1", "c1", 1800, 500);
  TCanvas *c2 = new TCanvas("c2", "c2", 1800, 500);
  TCanvas *c3 = new TCanvas("c3", "c3", 1800, 500); 
  TCanvas *c4 = new TCanvas("c4", "c4", 600, 500); // track default
  TCanvas *c5 = new TCanvas("c5", "c5", 600, 500); // shower default
  gStyle->SetPalette(kLightTemperature);
  c1->Divide(5,2);
  c2->Divide(5,2);
  c3->Divide(5,2);
  for (int i = 0; i < 10; i++) {
    c1->cd(i+1);
    ht[i]->Draw("");
    hs[i]->Draw("same");
    hm[i]->Draw("same");
    hs[i]->SetLineColor(kRed);
    hm[i]->SetLineColor(kBlack);
    ht[i]->GetXaxis()->SetTitle("Energy [GeV]");
    c1->cd(i+1)->SetLogx();

    TLegend* leg = new TLegend(0.1, 0.75, 0.3, 0.9);
    leg->AddEntry(ht[i], "Track", "l");
    leg->AddEntry(hs[i], "Shower", "l");
    leg->AddEntry(hm[i], "MC", "l");
    leg->Draw();


    c2->cd(i+1);
    GetNormalizedSlicesY(h2t[i]);
    h2t[i]->Draw("colz");
    h2t[i]->GetXaxis()->SetTitle("E_{true} [GeV]");
    h2t[i]->GetYaxis()->SetTitle("E_{reco} [GeV]");
    h2t[i]->GetXaxis()->SetRangeUser(3,100);
    h2t[i]->GetYaxis()->SetRangeUser(3,100);
    c2->cd(i+1)->SetLogx();
    c2->cd(i+1)->SetLogy();
    c2->cd(i+1)->SetLogz();

    c3->cd(i+1);
    GetNormalizedSlicesY(h2s[i]);
    h2s[i]->Draw("colz");
    h2s[i]->GetXaxis()->SetTitle("E_{true} [GeV]");
    h2s[i]->GetYaxis()->SetTitle("E_{reco} [GeV]");

    h2s[i]->GetXaxis()->SetRangeUser(3,100);
    h2s[i]->GetYaxis()->SetRangeUser(3,100);
    c3->cd(i+1)->SetLogx();
    c3->cd(i+1)->SetLogy();
    c3->cd(i+1)->SetLogz();
  }

  c4->cd();
  GetNormalizedSlicesY(ht_def);
  ht_def->Draw("colz");
  ht_def->GetXaxis()->SetTitle("E_{true} [GeV]");
  ht_def->GetYaxis()->SetTitle("E_{reco} [GeV]");

  ht_def->GetXaxis()->SetRangeUser(3,100);
  ht_def->GetYaxis()->SetRangeUser(3,100);
  c4->SetLogx();
  c4->SetLogy();
  c4->SetLogz();
  
  c5->cd();
  GetNormalizedSlicesY(hs_def);
  hs_def->Draw("colz");
  hs_def->GetXaxis()->SetTitle("E_{true} [GeV]");
  hs_def->GetYaxis()->SetTitle("E_{reco} [GeV]");

  hs_def->GetXaxis()->SetRangeUser(3,100);
  hs_def->GetYaxis()->SetRangeUser(3,100);
  c5->SetLogx();
  c5->SetLogy();
  c5->SetLogz();

  c1->SaveAs("./pid_detres/reconstructed_evts_per_q.pdf");
  c2->SaveAs("./pid_detres/energy_resolution_track_per_q.pdf");
  c3->SaveAs("./pid_detres/energy_resolution_shower_per_q.pdf");
  c4->SaveAs("./default_detres/energy_resolution_track.pdf");
  c5->SaveAs("./default_detres/energy_resolution_shower.pdf");

}
