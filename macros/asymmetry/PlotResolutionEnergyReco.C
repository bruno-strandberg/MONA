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

void PlotResolutionEnergyReco(TString summary_file=(TString)getenv("NMHDIR") + "/data/ORCA_MC_summary_all_10Apr2018.root") {
  SummaryParser sp(sum_file);
  
  bool plot = false;

  int n_bins = 40;
  std::vector<Double_t> e_edges  = NMHUtils::GetLogBins(n_bins, 1, 100);
  std::vector<Double_t> ct_edges = NMHUtils::GetBins(n_bins, -1, 1);

  // Default scheme: all events
  TH2D* ht_def = new TH2D("ht_def", "Reconstructed events", n_bins, &e_edges[0], n_bins, &e_edges[0]);
  TH2D* hs_def = new TH2D("hs_def", "Reconstructed events", n_bins, &e_edges[0], n_bins, &e_edges[0]);

  std::vector<TH2D*> h2t;
  std::vector<TH2D*> h2s;

  for (int i = 0; i < 10; i++){
    h2t.push_back(new TH2D(Form("h2t_%i", i), Form("Energy resolution tracks_q_0.%i", i), n_bins, &e_edges[0], n_bins, &e_edges[0]));
    h2s.push_back(new TH2D(Form("h2s_%i", i), Form("Energy resolution showers_q_0.%i", i), n_bins, &e_edges[0], n_bins, &e_edges[0]));
  } 

  Double_t q;
  for (Int_t i = 0; i < sp.GetTree()->GetEntries(); i++) {
    total_events++;
    // Filters 
    bool good_tr = false;
    bool good_sh = false;

    SummaryEvent *evt = sp.GetEvt(i);
    if ((evt->Get_RDF_muon_score() > 0.5) or (evt->Get_RDF_noise_score() > 0.5)) { continue; } // Filters for the events
    if ((evt->Get_track_ql0() > 0.5) and (evt->Get_track_ql1() > 0.5)) { good_tr = true; }
    if ((evt->Get_shower_ql0() > 0.5) and (evt->Get_shower_ql1() > 0.5)) { good_sh = true; }
    if ((not good_tr) and (not good_sh)) { continue; }
    else if (not good_tr) { continue; }
    else if (not good_sh) { continue; }

    // We order the histograms by quality: 0.0-0.1 --> 0, 0.1-0.2 --> 1, etc.
    q = evt->Get_RDF_track_score();
    int index = (int)(TMath::Floor(q * 10));
    if (index == 10) { index = 9; } // A perfect track will get an index 10, which does not exist

    // Fill events
    h2t[index]->Fill(evt->Get_MC_energy(), evt->Get_track_energy());
    h2s[index]->Fill(evt->Get_MC_energy(), evt->Get_shower_energy());

    if (q > 0.6) { ht_def->Fill(evt->Get_MC_energy(), evt->Get_track_energy()); }
    if (q < 0.6) { hs_def->Fill(evt->Get_MC_energy(), evt->Get_shower_energy()); }
  }

  if (plot) {
    TCanvas *c1 = new TCanvas("c1", "c1", 1800, 500);
    TCanvas *c2 = new TCanvas("c2", "c2", 1800, 500); 
    TCanvas *c3 = new TCanvas("c3", "c3", 600, 500); // track default
    TCanvas *c4 = new TCanvas("c4", "c4", 600, 500); // shower default
    gStyle->SetPalette(kBird);
    c1->Divide(5,2);
    c2->Divide(5,2);
    for (int i = 0; i < 10; i++) {
      c1->cd(i+1);
      GetNormalizedSlicesY(h2t[i]);
      h2t[i]->Draw("colz");
      h2t[i]->GetXaxis()->SetTitle("E_{true} [GeV]");
      h2t[i]->GetYaxis()->SetTitle("E_{reco} [GeV]");
      h2t[i]->GetXaxis()->SetRangeUser(3,100);
      h2t[i]->GetYaxis()->SetRangeUser(3,100);
      c1->cd(i+1)->SetLogx();
      c1->cd(i+1)->SetLogy();
      c1->cd(i+1)->SetLogz();

      c2->cd(i+1);
      GetNormalizedSlicesY(h2s[i]);
      h2s[i]->Draw("colz");
      h2s[i]->GetXaxis()->SetTitle("E_{true} [GeV]");
      h2s[i]->GetYaxis()->SetTitle("E_{reco} [GeV]");

      h2s[i]->GetXaxis()->SetRangeUser(3,100);
      h2s[i]->GetYaxis()->SetRangeUser(3,100);
      c2->cd(i+1)->SetLogx();
      c2->cd(i+1)->SetLogy();
      c2->cd(i+1)->SetLogz();
    }

    c3->cd();
    GetNormalizedSlicesY(ht_def);
    ht_def->Draw("colz");
    ht_def->GetXaxis()->SetTitle("E_{true} [GeV]");
    ht_def->GetYaxis()->SetTitle("E_{reco} [GeV]");

    ht_def->GetXaxis()->SetRangeUser(3,100);
    ht_def->GetYaxis()->SetRangeUser(3,100);
    c3->SetLogx();
    c3->SetLogy();
    c3->SetLogz();
    
    c4->cd();
    GetNormalizedSlicesY(hs_def);
    hs_def->Draw("colz");
    hs_def->GetXaxis()->SetTitle("E_{true} [GeV]");
    hs_def->GetYaxis()->SetTitle("E_{reco} [GeV]");

    hs_def->GetXaxis()->SetRangeUser(3,100);
    hs_def->GetYaxis()->SetRangeUser(3,100);
    c4->SetLogx();
    c4->SetLogy();
    c4->SetLogz();

    c1->SaveAs("./pid_detres/energy_resolution_track_events.pdf");
    c2->SaveAs("./pid_detres/energy_resolution_shower_events.pdf");
    c3->SaveAs("./default_detres/energy_resolution_track_events.pdf");
    c4->SaveAs("./default_detres/energy_resolution_shower_events.pdf");
  }
}
