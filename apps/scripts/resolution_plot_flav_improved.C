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

void resolution_plot_flav_improved(TString sum_file="../../data/ORCA_MC_summary_all_10Apr2018.root") { 
  SummaryParser sp(sum_file);
  
  bool plot = true;
  int n_bins = 40;
  std::vector<Double_t> e_edges  = NMHUtils::GetLogBins(n_bins, 1, 100);
  std::vector<Double_t> ct_edges = NMHUtils::GetBins(n_bins, -1, 1);

  std::vector<TH2D*> h2mucc; // THIS IS SHOWER THEN TRACK FILLING
  std::vector<TH2D*> h2munc;
  std::vector<TH2D*> h2mucc_good_track;
  std::vector<TH2D*> h2mucc_tr_then_sh;
  for (int i = 0; i < 10; i++){
    h2mucc.push_back(new TH2D(Form("h2mucc_%i", i), Form("Energy resolution nu_m_CC sh then tr q_0.%i", i), n_bins, &e_edges[0], n_bins, &e_edges[0]));
    h2munc.push_back(new TH2D(Form("h2munc_%i", i), Form("Energy resolution nu_m_NC q_0.%i", i), n_bins, &e_edges[0], n_bins, &e_edges[0]));
    h2mucc_good_track.push_back(new TH2D(Form("h2mucc_g_t_%i", i), Form("Energy resolution nu_m_CC good track q_0.%i", i), n_bins, &e_edges[0], n_bins, &e_edges[0]));
    h2mucc_tr_then_sh.push_back(new TH2D(Form("h2mucc_t_s_%i", i), Form("Energy resolution nu_m_CC tr then sh q_0.%i", i), n_bins, &e_edges[0], n_bins, &e_edges[0]));
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

    // We order the histograms by quality: 0.0-0.1 --> 0, 0.1-0.2 --> 1, etc.
    q = sp.GetEvt()->Get_RDF_track_score();
    int index = (int)(TMath::Floor(q * 10));
    if (index == 10) { index = 9; } // A perfect track will get an index 10, which does not exist

    if ((std::abs(sp.GetEvt()->Get_MC_type()) == 13) or (std::abs(sp.GetEvt()->Get_MC_type()) == 14)) {
      if (good_sh) { 
        if (sp.GetEvt()->Get_MC_is_CC()) h2mucc[index]->Fill(sp.GetEvt()->Get_MC_energy(), sp.GetEvt()->Get_shower_energy());
        else                             h2munc[index]->Fill(sp.GetEvt()->Get_MC_energy(), sp.GetEvt()->Get_shower_energy());
      }
      else {
        if (good_tr) {
          if (sp.GetEvt()->Get_MC_is_CC()) { 
            h2mucc[index]->Fill(sp.GetEvt()->Get_MC_energy(), sp.GetEvt()->Get_track_energy());
            h2mucc_good_track[index]->Fill(sp.GetEvt()->Get_MC_energy(), sp.GetEvt()->Get_track_energy());
          }
          else h2munc[index]->Fill(sp.GetEvt()->Get_MC_energy(), sp.GetEvt()->Get_shower_energy());
        }
      }
    }
    if ((std::abs(sp.GetEvt()->Get_MC_type()) == 13) or (std::abs(sp.GetEvt()->Get_MC_type()) == 14)) {
      if (good_tr) { 
        if (sp.GetEvt()->Get_MC_is_CC()) h2mucc_tr_then_sh[index]->Fill(sp.GetEvt()->Get_MC_energy(), sp.GetEvt()->Get_track_energy());
      }
      else {
        if (good_sh) {
          if (sp.GetEvt()->Get_MC_is_CC()) h2mucc_tr_then_sh[index]->Fill(sp.GetEvt()->Get_MC_energy(), sp.GetEvt()->Get_track_energy());
        }
      }
    }
    if ((not good_tr) and (not good_sh)) { bad_track_and_shower++; continue; }
    else if (not good_tr) { bad_track++; continue; }
    else if (not good_sh) { bad_shower++; continue; }
    good_track_and_shower++;

  }
  cout << "Total events: " << total_events << endl;
  cout << "Good events : " << good_track_and_shower << endl;
  cout << "Bad tracks  : " << bad_track << endl;
  cout << "Bad showers : " << bad_shower << endl;
  cout << "Bad both    : " << bad_track_and_shower << endl;

  if (plot) {
    TCanvas *c3 = new TCanvas("c3", "c3", 1800, 500); 
    TCanvas *c4 = new TCanvas("c4", "c4", 1800, 500); 
    TCanvas *c5 = new TCanvas("c5", "c5", 1800, 500); 
    TCanvas *c6 = new TCanvas("c6", "c6", 1800, 500); 
    gStyle->SetPalette(kLightTemperature);
    c3->Divide(5,2);
    c4->Divide(5,2);
    c5->Divide(5,2);
    c6->Divide(5,2);
    for (int i = 0; i < 10; i++) {
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
      GetNormalizedSlicesY(h2mucc_good_track[i]);
      h2mucc_good_track[i]->Draw("colz");
      h2mucc_good_track[i]->GetXaxis()->SetTitle("E_{true} [GeV]");
      h2mucc_good_track[i]->GetYaxis()->SetTitle("E_{reco} [GeV]");
      h2mucc_good_track[i]->GetXaxis()->SetRangeUser(3,100);
      h2mucc_good_track[i]->GetYaxis()->SetRangeUser(3,100);
      c5->cd(i+1)->SetLogx();
      c5->cd(i+1)->SetLogy();

      c6->cd(i+1);
      GetNormalizedSlicesY(h2mucc_tr_then_sh[i]);
      h2mucc_tr_then_sh[i]->Draw("colz");
      h2mucc_tr_then_sh[i]->GetXaxis()->SetTitle("E_{true} [GeV]");
      h2mucc_tr_then_sh[i]->GetYaxis()->SetTitle("E_{reco} [GeV]");
      h2mucc_tr_then_sh[i]->GetXaxis()->SetRangeUser(3,100);
      h2mucc_tr_then_sh[i]->GetYaxis()->SetRangeUser(3,100);
      c6->cd(i+1)->SetLogx();
      c6->cd(i+1)->SetLogy();
    }

    c3->SaveAs("./pid_detres/resolution_different_energy_schemes/energy_resolution_mucc_per_q_sh_then_tr.pdf");
    c4->SaveAs("./pid_detres/resolution_different_energy_schemes/energy_resolution_munc_per_q.pdf");
    c5->SaveAs("./pid_detres/resolution_different_energy_schemes/energy_resolution_mucc_per_q_good_tracks.pdf");
    c6->SaveAs("./pid_detres/resolution_different_energy_schemes/energy_resolution_mucc_per_q_tr_then_sh.pdf");
  }
}
