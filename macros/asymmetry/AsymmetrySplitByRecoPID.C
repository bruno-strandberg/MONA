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
#include "HelperFunctions.C"
#include "NMHUtils.h"
#include "SummaryParser.h"
#include "SummaryEvent.h"

#include <iostream>
using namespace std;

/* This script calculates the asymmetry for track energy binned data in good track, good shower and good event
 * and simultaneously using N (default=10) PID bins to split the data along the quality axis.
 * Showers are binned like showers as in the default detector response scheme, but also have the N PID bins.
 *
 * The output files are saved into `filefolder`
 */



void AsymmetrySplitByRecoPID() {

  bool b_plot = false;
  const int N_PID_CLASSES = 10;
  Double_t PID_step = 1 / float(N_PID_CLASSES);
  TString filefolder = Form("./energy_detres/pid_bins_%i/", N_PID_CLASSES);

  vector<Double_t> asym_t_gt(N_PID_CLASSES);
  vector<Double_t> asym_t_gs(N_PID_CLASSES);
  vector<Double_t> asym_t_ge(N_PID_CLASSES);
  vector<Double_t> asym_s(N_PID_CLASSES);
  vector<Double_t> asym_m(N_PID_CLASSES);

  vector<Double_t> asym_t_gt_err(N_PID_CLASSES);
  vector<Double_t> asym_t_gs_err(N_PID_CLASSES);
  vector<Double_t> asym_t_ge_err(N_PID_CLASSES);
  vector<Double_t> asym_s_err(N_PID_CLASSES);
  vector<Double_t> asym_m_err(N_PID_CLASSES);

  cout << "Asymmetries per single bin: " << endl;
  for (int i = 0; i < N_PID_CLASSES; i++) {
    TString file_NO = filefolder + TString::Format("split_energy_NO_%.2f.root", PID_step * i);
    TString file_IO = filefolder + TString::Format("split_energy_IO_%.2f.root", PID_step * i);
    TString output  = filefolder + TString::Format("asym_energy_%.2f.root", PID_step * i);
    
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

    auto asym_t_gt_tuple = NMHUtils::Asymmetry(h_t_gt_NO, h_t_gt_IO, TString::Format("asymmetry_track_gt_%.2f", i * PID_step), 2, 80, -1, 0);
    auto asym_t_gs_tuple = NMHUtils::Asymmetry(h_t_gs_NO, h_t_gs_IO, TString::Format("asymmetry_track_gs_%.2f", i * PID_step), 2, 80, -1, 0);
    auto asym_t_ge_tuple = NMHUtils::Asymmetry(h_t_ge_NO, h_t_ge_IO, TString::Format("asymmetry_track_ge_%.2f", i * PID_step), 2, 80, -1, 0);
    auto asym_s_tuple    = NMHUtils::Asymmetry(h_s_NO, h_s_IO, TString::Format("asymmetry_shower_%.2f", i * PID_step), 2, 80, -1, 0);
    auto asym_m_tuple    = NMHUtils::Asymmetry(h_m_NO, h_m_IO, TString::Format("asymmetry_mc_%.2f", i * PID_step), 2, 80, -1, 0);

    TH2D* h_asym_t_gt = std::get<0>(asym_t_gt_tuple);
    TH2D* h_asym_t_gs = std::get<0>(asym_t_gs_tuple);
    TH2D* h_asym_t_ge = std::get<0>(asym_t_ge_tuple);
    TH2D* h_asym_s    = std::get<0>(asym_s_tuple);
    TH2D* h_asym_m    = std::get<0>(asym_m_tuple);

    TCanvas *c1 = new TCanvas(TString::Format("c1_%.2f", i * PID_step), "c1", 2000, 1200); // has to be outside of the if, so that the writing below works.
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

    // Lets print the number of events per category
    Double_t num_NO_gt = ((TH2D*)f_NO->Get("detected_tracks_gt"))->Integral();
    Double_t num_NO_gs = ((TH2D*)f_NO->Get("detected_tracks_gs"))->Integral();
    Double_t num_NO_ge = ((TH2D*)f_NO->Get("detected_tracks_ge"))->Integral();
    Double_t num_NO_s  = ((TH2D*)f_NO->Get("detected_showers"))->Integral();
    Double_t num_NO_m  = ((TH2D*)f_NO->Get("detected_mc"))->Integral();

    Double_t num_IO_gt = ((TH2D*)f_IO->Get("detected_tracks_gt"))->Integral();
    Double_t num_IO_gs = ((TH2D*)f_IO->Get("detected_tracks_gs"))->Integral();
    Double_t num_IO_ge = ((TH2D*)f_IO->Get("detected_tracks_ge"))->Integral();
    Double_t num_IO_s  = ((TH2D*)f_IO->Get("detected_showers"))->Integral();
    Double_t num_IO_m  = ((TH2D*)f_IO->Get("detected_mc"))->Integral();

    asym_t_gt[i] = std::get<1>(asym_t_gt_tuple);
    asym_t_gs[i] = std::get<1>(asym_t_gs_tuple);
    asym_t_ge[i] = std::get<1>(asym_t_ge_tuple);
    asym_s[i]    = std::get<1>(asym_s_tuple);
    asym_m[i]    = std::get<1>(asym_m_tuple);
    asym_t_gt_err[i] = std::get<2>(asym_t_gt_tuple);
    asym_t_gs_err[i] = std::get<2>(asym_t_gs_tuple);
    asym_t_ge_err[i] = std::get<2>(asym_t_ge_tuple);
    asym_s_err[i]    = std::get<2>(asym_s_tuple);
    asym_m_err[i]    = std::get<2>(asym_m_tuple);

    PrintAsymmetryWithErrors("tracks (good tracks)", asym_t_gt[i], asym_t_gt_err[i], num_NO_gt, num_IO_gt);
    PrintAsymmetryWithErrors("tracks (good showrs)", asym_t_gs[i], asym_t_gs_err[i], num_NO_gs, num_IO_gs);
    PrintAsymmetryWithErrors("tracks (good events)", asym_t_ge[i], asym_t_ge_err[i], num_NO_ge, num_IO_ge);
    PrintAsymmetryWithErrors("showers             ", asym_s[i], asym_s_err[i], num_NO_s, num_IO_s);
    PrintAsymmetryWithErrors("mc                  ", asym_m[i], asym_m_err[i], num_NO_m, num_NO_m);

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

  cout << "Single counting (q<0.6): " << endl;
  int offset = N_PID_CLASSES * 0.6;

  std::pair<Double_t, Double_t> asym_t_gt_q6;
  std::pair<Double_t, Double_t> asym_t_gs_q6;
  std::pair<Double_t, Double_t> asym_t_ge_q6;
  std::pair<Double_t, Double_t> asym_s_q6;
  for (Int_t i = 0; i < N_PID_CLASSES; i++) {
    if (i < offset) {
      asym_s_q6.push_back(std::make_pair(asym_s[i], asym_s_err[i]));
    } 
    else {
      asym_t_gt_q6.push_back(std::make_pair(asym_t_gt[i], asym_t_gt_err[i]));
      asym_t_gs_q6.push_back(std::make_pair(asym_t_gs[i], asym_t_gs_err[i]));
      asym_t_ge_q6.push_back(std::make_pair(asym_t_ge[i], asym_t_ge_err[i]));
    }
  }

  std::tuple<Double_t, Double_t> track_gt_value_squared_q6 = 
    NMHUtils::SquaredSumErrorProp(asym_t_gt_q6);
  std::tuple<Double_t, Double_t> track_gs_value_squared_q6 = 
    NMHUtils::SquaredSumErrorProp(asym_t_gs_q6);
  std::tuple<Double_t, Double_t> track_ge_value_squared_q6 = 
    NMHUtils::SquaredSumErrorProp(asym_t_ge_q6);
  std::tuple<Double_t, Double_t> shower_value_squared_q6 = 
    NMHUtils::SquaredSumErrorProp(asym_s_q6);
  std::tuple<Double_t, Double_t> total_value_squared_q6  = 
    NMHUtils::SquaredSumErrorProp( { std::make_pair(std::get<0>(track_gt_value_squared_q6), std::get<1>(track_gt_value_squared_q6)),
                                     std::make_pair(std::get<0>(track_gs_value_squared_q6), std::get<1>(track_gs_value_squared_q6)), 
                                     std::make_pair(std::get<0>(track_ge_value_squared_q6), std::get<1>(track_ge_value_squared_q6)), 
                                     std::make_pair(std::get<0>(shower_value_squared_q6), std::get<1>(shower_value_squared_q6)) } );

  Double_t track_gt_value_q6  = std::get<0>(track_gt_value_squared_q6);
  Double_t track_gs_value_q6  = std::get<0>(track_gs_value_squared_q6);
  Double_t track_ge_value_q6  = std::get<0>(track_ge_value_squared_q6);
  Double_t shower_value_q6    = std::get<0>(shower_value_squared_q6);
  Double_t total_value_q6     = std::get<0>(total_value_squared_q6);

  Double_t track_gt_error_q6  = std::get<1>(track_gt_value_squared_q6);
  Double_t track_gs_error_q6  = std::get<1>(track_gs_value_squared_q6);
  Double_t track_ge_error_q6  = std::get<1>(track_ge_value_squared_q6);
  Double_t shower_error_q6    = std::get<1>(shower_value_squared_q6);
  Double_t total_error_q6     = std::get<1>(total_value_squared_q6);

  PrintAsymmetryWithErrors("tracks (good tracks)", track_gt_value_q6, track_gt_error_q6); 
  PrintAsymmetryWithErrors("tracks (good showrs)", track_gs_value_q6, track_gs_error_q6); 
  PrintAsymmetryWithErrors("tracks (good events)", track_ge_value_q6, track_ge_error_q6); 
  PrintAsymmetryWithErrors("total showers       ", shower_value_q6, shower_error_q6); 
  PrintAsymmetryWithErrors("total combined      ", total_value_q6, total_error_q6); 
}
