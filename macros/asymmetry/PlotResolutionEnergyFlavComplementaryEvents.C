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

void PlotResolutionEnergyFlavComplementaryEvents(TString summary_file=(TString)getenv("MONADIR") + "/data/ORCA_MC_summary_all_10Apr2018.root") {
  SummaryParser sp(sum_file);
  
  bool plot = false;
  int n_bins = 40;
  std::vector<Double_t> e_edges  = NMHUtils::GetLogBins(n_bins, 1, 100);
  std::vector<Double_t> ct_edges = NMHUtils::GetBins(n_bins, -1, 1);

  // Different selections: only track, only shower, good events (good tr, good sh)
  std::vector<TH2D*> h2_g_tr;
  std::vector<TH2D*> h2_g_sh;
  std::vector<TH2D*> h2_g_ev;
  for (int i = 0; i < 10; i++){
    h2_g_tr.push_back(new TH2D(Form("h2_g_tr_%i", i), Form("Energy resolution nu_m_CC q_0.%i [good tr]", i), n_bins, &e_edges[0], n_bins, &e_edges[0]));
    h2_g_sh.push_back(new TH2D(Form("h2_g_sh_%i", i), Form("Energy resolution nu_m_CC q_0.%i [good sh]", i), n_bins, &e_edges[0], n_bins, &e_edges[0]));
    h2_g_ev.push_back(new TH2D(Form("h2_g_ev_%i", i), Form("Energy resolution nu_m_CC q_0.%i [good ev]", i), n_bins, &e_edges[0], n_bins, &e_edges[0]));
  }

  Double_t q;
  // Good Tracks only, no good showers
  for (Int_t i = 0; i < sp.GetTree()->GetEntries(); i++) {
    // Filters 
    bool good_tr = false;
    bool good_sh = false;

    SummaryEvent *evt = sp.GetEvt(i);
    if ((evt->Get_RDF_muon_score() > 0.05) or (evt->Get_RDF_noise_score() > 0.5)) { continue; } // Filters for the events
    if ((evt->Get_track_ql0() > 0.5) and (evt->Get_track_ql1() > 0.5)) { good_tr = true; }
    if ((evt->Get_shower_ql0() > 0.5) and (evt->Get_shower_ql1() > 0.5)) { good_sh = true; }
    if ((good_tr) and (not good_sh)) { 
      q = evt->Get_RDF_track_score();
      int index = (int)(TMath::Floor(q * 10));
      if (index == 10) { index = 9; } // A perfect track will get an index 10, which does not exist
      if ((std::abs(evt->Get_MC_type()) == 13) or (std::abs(evt->Get_MC_type()) == 14)) {
        if (evt->Get_MC_is_CC()) h2_g_tr[index]->Fill(evt->Get_MC_energy(), evt->Get_track_energy());
      }
    }
    if ((good_sh) and (not good_tr)) {
      q = evt->Get_RDF_track_score();
      int index = (int)(TMath::Floor(q * 10));
      if (index == 10) { index = 9; }
      if ((std::abs(evt->Get_MC_type()) == 13) or (std::abs(evt->Get_MC_type()) == 14)) {
        if (evt->Get_MC_is_CC()) h2_g_sh[index]->Fill(evt->Get_MC_energy(), evt->Get_shower_energy());
      }
    }
    if ((good_sh) and (good_tr)) {
      q = evt->Get_RDF_track_score();
      int index = (int)(TMath::Floor(q * 10));
      if (index == 10) { index = 9; }
      if ((std::abs(evt->Get_MC_type()) == 13) or (std::abs(evt->Get_MC_type()) == 14)) {
        if (evt->Get_MC_is_CC()) h2_g_ev[index]->Fill(evt->Get_MC_energy(), evt->Get_shower_energy());
      }
    }
  }

  if (plot) {
    TCanvas *c1 = new TCanvas("c1", "c1", 1800, 500);
    TCanvas *c2 = new TCanvas("c2", "c2", 1800, 500);
    TCanvas *c3 = new TCanvas("c3", "c3", 1800, 500); 
    gStyle->SetPalette(kLightTemperature);
    c1->Divide(5,2);
    c2->Divide(5,2);
    c3->Divide(5,2);
    for (int i = 0; i < 10; i++) {
      c1->cd(i+1);
      GetNormalizedSlicesY(h2_g_tr[i]);
      h2_g_tr[i]->Draw("colz");
      h2_g_tr[i]->GetXaxis()->SetTitle("E_{true} [GeV]");
      h2_g_tr[i]->GetYaxis()->SetTitle("E_{reco} [GeV]");
      h2_g_tr[i]->GetXaxis()->SetRangeUser(3,100);
      h2_g_tr[i]->GetYaxis()->SetRangeUser(3,100);
      c1->cd(i+1)->SetLogx();
      c1->cd(i+1)->SetLogy();
      c1->cd(i+1)->SetLogz();

      c2->cd(i+1);
      GetNormalizedSlicesY(h2_g_sh[i]);
      h2_g_sh[i]->Draw("colz");
      h2_g_sh[i]->GetXaxis()->SetTitle("E_{true} [GeV]");
      h2_g_sh[i]->GetYaxis()->SetTitle("E_{reco} [GeV]");
      h2_g_sh[i]->GetXaxis()->SetRangeUser(3,100);
      h2_g_sh[i]->GetYaxis()->SetRangeUser(3,100);
      c2->cd(i+1)->SetLogx();
      c2->cd(i+1)->SetLogy();
      c2->cd(i+1)->SetLogz();

      c3->cd(i+1);
      GetNormalizedSlicesY(h2_g_ev[i]);
      h2_g_ev[i]->Draw("colz");
      h2_g_ev[i]->GetXaxis()->SetTitle("E_{true} [GeV]");
      h2_g_ev[i]->GetYaxis()->SetTitle("E_{reco} [GeV]");
      h2_g_ev[i]->GetXaxis()->SetRangeUser(3,100);
      h2_g_ev[i]->GetYaxis()->SetRangeUser(3,100);
      c3->cd(i+1)->SetLogx();
      c3->cd(i+1)->SetLogy();
      c3->cd(i+1)->SetLogz();
    }

    c1->SaveAs("./pid_detres/energy_resolution_complementing_events/energy_resolution_mucc_track_events.pdf");
    c2->SaveAs("./pid_detres/energy_resolution_complementing_events/energy_resolution_mucc_shower_events.pdf");
    c3->SaveAs("./pid_detres/energy_resolution_complementing_events/energy_resolution_mucc_good_events.pdf");
  }
}
