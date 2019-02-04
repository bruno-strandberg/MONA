#include "TSystem.h"
#include "TROOT.h"
#include "TDatime.h"
#include "TH2.h"
#include "TH3.h"
#include "TFile.h"
#include "TStopwatch.h"
#include "TCanvas.h"
#include "TStyle.h"

#include "DetResponse.h"
#include "NMHUtils.h"
#include "SummaryParser.h"
#include "SummaryEvent.h"

#include <iostream>
using namespace std;

/* This script calculates the asymmetry for track energy binned data in good track, good shower and good event.
 * Showers are binned like showers as in the default detector response scheme.
 *
 * The output files are saved into `filefolder`
 */


void AsymmetrySplitByReco() {

  bool b_plot = false;
  TString filefolder = "./energy_detres/";

  Double_t asym_t_gt;
  Double_t asym_t_gs;
  Double_t asym_t_ge;
  Double_t asym_s;
  Double_t asym_m;

  Double_t asym_t_gt_err;
  Double_t asym_t_gs_err;
  Double_t asym_t_ge_err;
  Double_t asym_s_err;
  Double_t asym_m_err;

  TString file_NO = filefolder + "split_expected_evts_NO.root";
  TString file_IO = filefolder + "split_expected_evts_IO.root";
  TString output  = filefolder + "asymmetry_split_by_energy.root";
  
  TFile *f_NO  = TFile::Open(file_NO, "READ");
  TFile *f_IO  = TFile::Open(file_IO, "READ");

  TH2D *h_t_gt_NO = (TH2D*)f_NO->Get("detected_tracks_gt");
  TH2D *h_t_gs_NO = (TH2D*)f_NO->Get("detected_tracks_gs");
  TH2D *h_t_ge_NO = (TH2D*)f_NO->Get("detected_tracks_ge");
  TH2D *h_s_NO = (TH2D*)f_NO->Get("detected_showers");
  TH2D *h_m_NO = (TH2D*)f_NO->Get("detected_mc");

  h_t_gt_NO->SetName("detected_tracks_gt_NO");
  h_t_gs_NO->SetName("detected_tracks_gs_NO");
  h_t_ge_NO->SetName("detected_tracks_ge_NO");
  h_s_NO->SetName("detected_showers_NO");
  h_m_NO->SetName("detected_mc_NO");

  TH2D *h_t_gt_IO = (TH2D*)f_IO->Get("detected_tracks_gt");
  TH2D *h_t_gs_IO = (TH2D*)f_IO->Get("detected_tracks_gs");
  TH2D *h_t_ge_IO = (TH2D*)f_IO->Get("detected_tracks_ge");
  TH2D *h_s_IO = (TH2D*)f_IO->Get("detected_showers");
  TH2D *h_m_IO = (TH2D*)f_IO->Get("detected_mc");

  h_t_gt_IO->SetName("detected_tracks_gt_IO");
  h_t_gs_IO->SetName("detected_tracks_gs_IO");
  h_t_ge_IO->SetName("detected_tracks_ge_IO");
  h_s_IO->SetName("detected_showers_IO");
  h_m_IO->SetName("detected_mc_IO");

  auto asym_t_gt_tuple = NMHUtils::Asymmetry(h_t_gt_NO, h_t_gt_IO, "asymmetry_track_gt", 2, 80, -1, 0);
  auto asym_t_gs_tuple = NMHUtils::Asymmetry(h_t_gs_NO, h_t_gs_IO, "asymmetry_track_gs", 2, 80, -1, 0);
  auto asym_t_ge_tuple = NMHUtils::Asymmetry(h_t_ge_NO, h_t_ge_IO, "asymmetry_track_ge", 2, 80, -1, 0);
  auto asym_s_tuple    = NMHUtils::Asymmetry(h_s_NO, h_s_IO, "asymmetry_shower", 2, 80, -1, 0);
  auto asym_m_tuple    = NMHUtils::Asymmetry(h_m_NO, h_m_IO, "asymmetry_mc", 2, 80, -1, 0);

  TH2D* h_asym_t_gt = std::get<0>(asym_t_gt_tuple);
  TH2D* h_asym_t_gs = std::get<0>(asym_t_gs_tuple);
  TH2D* h_asym_t_ge = std::get<0>(asym_t_ge_tuple);
  TH2D* h_asym_s    = std::get<0>(asym_s_tuple);
  TH2D* h_asym_m    = std::get<0>(asym_m_tuple);

  TCanvas *c1 = new TCanvas("c1", "c1", 2000, 1200); // has to be outside of the if, so that the writing below works.
  gStyle->SetPalette(kBird);
  if (b_plot) { 
    c1->Divide(5,3);

    vector<TH2D*> plots = {h_t_gt_NO, h_t_gs_NO, h_t_ge_NO, h_s_NO, h_m_NO, 
                           h_t_gt_IO, h_t_gs_IO, h_t_ge_IO, h_s_IO, h_m_IO, 
                           h_asym_t_gt, h_asym_t_gs, h_asym_t_ge, h_asym_s, h_asym_m};
    int i = 1;
    for (auto plot: plots) { 
      c1->cd(i);
      plot->Draw("colz");
      plot->GetXaxis()->SetTitle("Energy [GeV]");
      plot->GetYaxis()->SetTitle("cos(#theta)");
      plot->GetXaxis()->SetRangeUser(3,100);
      plot->GetYaxis()->SetRangeUser(-1,0);
      if (i>10) { plot->GetZaxis()->SetRangeUser(-1,1); }
      c1->cd(i)->SetLogx();

      i++;
    }
  }

  auto asym_t_gt = std::make_pair(std::get<1>(asym_t_gt_tuple), std::get<2>(asym_t_gt_tuple));
  auto asym_t_gs = std::make_pair(std::get<1>(asym_t_gs_tuple), std::get<2>(asym_t_gs_tuple));
  auto asym_t_ge = std::make_pair(std::get<1>(asym_t_ge_tuple), std::get<2>(asym_t_ge_tuple));
  auto asym_s    = std::make_pair(std::get<1>(asym_s_tuple), std::get<2>(asym_s_tuple));
  auto asym_m    = std::make_pair(std::get<1>(asym_m_tuple), std::get<2>(asym_m_tuple));
  
  PrintAsymmetryWithErrors("tracks (good tracks)", asym_t_gt.first, asym_t_gt.second);
  PrintAsymmetryWithErrors("tracks (good showrs)", asym_t_gs.first, asym_t_gs.second);
  PrintAsymmetryWithErrors("tracks (good events)", asym_t_ge.first, asym_t_ge.second);
  PrintAsymmetryWithErrors("showers             ", asym_s.first, asym_s.second);
  PrintAsymmetryWithErrors("mc                  ", asym_m.first, asym_m.second);

  std::tuple<Double_t, Double_t> sq_sum = NMHUtils::SquaredSumErrorProp({asym_t_gt, asym_t_gs, asym_t_ge, asym_s, asym_m}); 
  double asym_tot_value = std::get<0>(sq_sum);
  double asym_tot_err   = std::get<1>(sq_sum);
  cout << "-----------------------------" << endl;
  PrintAsymmetryWithErrors("total combined      ", asym_tot_value, asym_tot_err);
  
  TFile fout(output, "RECREATE");
  c1->Write();
  h_t_gt_NO->Write();
  h_t_gt_IO->Write();
  h_t_gs_NO->Write();
  h_t_gs_IO->Write();
  h_t_ge_NO->Write();
  h_t_ge_IO->Write();
  h_s_NO->Write();
  h_s_IO->Write();
  h_m_NO->Write();
  h_m_IO->Write();
  h_asym_t_gt->Write();
  h_asym_t_gs->Write();
  h_asym_t_ge->Write();
  h_asym_s->Write();
  h_asym_m->Write();
  fout.Write();
  fout.Close();

  if (b_plot) { // To keep the memory low (esp. over tunnels) close the canvas.
    c1->Close(); 
    gSystem->ProcessEvents(); 
  }
}
