#include "TCanvas.h"
#include "TDatime.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH2.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"

#include "DetResponse.h"
#include "NMHUtils.h"
#include "SummaryParser.h"
#include "SummaryEvent.h"

#include <iostream>
using namespace std;

// This script plots the errors per E,CT bin into a graph as function of the PID number of bins

void plot_asymmetry_error_all_bins() {

  bool print = false;
  bool normalize = true;
  std::vector<Int_t> N_PID_CLASSES = {5, 10, 15, 20, 25, 30, 35, 40};
  std::vector<Double_t> N_PID_CLASSES_D(N_PID_CLASSES.begin(), N_PID_CLASSES.end())  ;
  Double_t Q_CUTOFF = 0.6;

  // Make two vectors that contain vectors with the average errors for tracks and showers:
  // { {average_5_bin_1, ..., average_5_bin_1600}, {average_10_bin_1, ... } , ..., { ..., average_40_bin_1600} }
  std::vector<std::vector<Double_t>> matrix_tr_errors;
  std::vector<std::vector<Double_t>> matrix_sh_errors;
  for (auto n_pid: N_PID_CLASSES) {

    Double_t PID_step = 1 / float(n_pid);
    TString basefolder = std::getenv("NMHDIR");
    TString filefolder = TString::Format(basefolder + "/fitter/pid_detres/simplified_error_calc/pid_binning_%i/", n_pid);
  
    std::vector<Double_t> asym_tr_err(1600); // I know we use 40x40 bins, however this is magic number, bad!
    std::vector<Double_t> asym_sh_err(1600);
    std::vector<Double_t> average_tr_err(1600, 0); // Total for errors
    std::vector<Double_t> average_sh_err(1600, 0);
  
    for (int i = 0; i < n_pid; i++) {
      TString input_file = filefolder + TString::Format("asym_pid_%.2f.root", PID_step * i);
      
      TFile *f = TFile::Open(input_file, "READ");
  
      TH2D* h_asym_tr = (TH2D*)f->Get(TString::Format("asymmetry_track_%.2f", PID_step * i));
      TH2D* h_asym_sh = (TH2D*)f->Get(TString::Format("asymmetry_shower_%.2f", PID_step * i));
  
      for (Int_t xb = 1; xb <= h_asym_tr->GetXaxis()->GetNbins(); xb++) {
        for (Int_t yb = 1; yb <= h_asym_tr->GetYaxis()->GetNbins(); yb++) {
          if (PID_step * i < Q_CUTOFF ) { // If below 0.6 take shower, else track.
            Double_t sh_error = h_asym_sh->GetBinError(xb, yb);
            asym_sh_err[(yb - 1) + (xb - 1) * h_asym_tr->GetYaxis()->GetNbins()] = sh_error;
          } else {
            Double_t tr_error = h_asym_tr->GetBinError(xb, yb);
            asym_tr_err[(yb - 1) + (xb - 1) * h_asym_tr->GetYaxis()->GetNbins()] = tr_error;
          }
        }
      }
      // Calculate the average error per (E, CT) BIN for a certain PID binning.
      for (unsigned int i = 0; i < average_tr_err.size(); i++) {
        average_tr_err[i] = average_tr_err[i] + asym_tr_err[i] / n_pid;
        average_sh_err[i] = average_sh_err[i] + asym_sh_err[i] / n_pid;
      }
    }

    matrix_tr_errors.push_back(average_tr_err);
    matrix_sh_errors.push_back(average_sh_err);
  }
  if (normalize) {
    for (int i = matrix_tr_errors.size() - 1; i >= 0; i--) {
      for (unsigned int j = 0; j < matrix_tr_errors[0].size(); j++) {
        if (matrix_tr_errors[0][j] != 0) {
          matrix_tr_errors[i][j] = matrix_tr_errors[i][j] / matrix_tr_errors[0][j];
          matrix_sh_errors[i][j] = matrix_sh_errors[i][j] / matrix_sh_errors[0][j];
        }
      }
    }
  }
  if (print) {
    for (unsigned int i = 0; i < N_PID_CLASSES.size(); i++) {
      for (unsigned int j = 0; j < matrix_tr_errors[i].size(); j++) {
        cout << matrix_tr_errors[i][j] << endl;
      }
    }
    for (unsigned int i = 0; i < N_PID_CLASSES.size(); i++) {
      for (unsigned int j = 0; j < matrix_tr_errors[i].size(); j++) {
        cout << matrix_sh_errors[i][j] << endl;
      }
    }
  }

  TH2D* h_track = new TH2D("h_tr", "h_tr", N_PID_CLASSES_D.size() + 1, 0, 45, 100, 0, 5);
  TH2D* h_shwer = new TH2D("h_sh", "h_sh", N_PID_CLASSES_D.size() + 1, 0, 45, 100, 0, 5);
  
  //for (int i = matrix_tr_errors.size() - 1; i >= 0; i--) {
  //  for (unsigned int j = 0; j < matrix_tr_errors[0].size(); j++) {
  //    if (matrix_tr_errors[i][j] == 0) continue;
  //    if (matrix_tr_errors[i][j] == 1) continue;
  //    h_track->Fill(N_PID_CLASSES_D[i], matrix_tr_errors[i][j]);
  //  }
  //}
  //h_track->Draw("colz");

  for (int i = matrix_sh_errors.size() - 1; i >= 0; i--) {
    for (unsigned int j = 0; j < matrix_sh_errors[0].size(); j++) {
      if (matrix_sh_errors[i][j] < 0.05) continue;
      if (matrix_sh_errors[i][j] == 1) continue;
      h_shwer->Fill(N_PID_CLASSES_D[i], matrix_sh_errors[i][j]);
    }
  }
  h_shwer->Draw("colz");

  // Plotting
//  TCanvas *c1 = new TCanvas("c1", "c1", 200, 10, 700, 500);
//  c1->SetGrid();
//
//  std::vector<TGraph*> graph_vector;
//  const int n = N_PID_CLASSES_D.size();
//  // track
//  for (unsigned int i = 0; i < matrix_sh_errors[0].size(); i++) {
//    TGraph* gtr = new TGraph(n);
//    if (matrix_tr_errors[0][i] == 0) continue;
//    for (int j = 0; j < n; j++) {
//      gtr->SetPoint(j, N_PID_CLASSES_D[j], matrix_tr_errors[j][i]);
//    }
//    graph_vector.push_back(gtr);
//  }
//  cout << graph_vector.size() << endl;
//  
//  graph_vector[0]->Draw("AL");
//  for (auto gtr: graph_vector) {
//    gtr->Draw("L");
//    gtr->SetTitle("Error on bins as function of number of PID bins");
//    gtr->GetXaxis()->SetTitle("Number of bins used");
//    gtr->GetYaxis()->SetTitle("Asymmetry error");
//    gtr->SetLineColorAlpha(kBlue+2, 0.01);
//    gtr->SetMarkerColor(kYellow+2);
//    gtr->SetMarkerStyle(4);
//    gtr->SetMarkerSize(1);
//    //gtr->GetYaxis()->SetRangeUser(0, 5e-3);
//    gtr->GetYaxis()->SetRangeUser(0, 5);
//    c1->Update();
//  }
//  graph_vector[0]->SaveAs("./all_errors_track.pdf");

//  // shower
//  TGraph *gsh = new TGraph(n);
//  for (int i=0; i < n; ++i) {
//    gsh->SetPoint(i, N_PID_CLASSES_D.at(i), average_sh_errors.at(i));
//  }
//  gsh->Draw("L");
//  gsh->SetLineColor(kBlue+2);
//  gsh->SetMarkerColor(kBlue+2);
//  gsh->SetMarkerStyle(4);
//  c1->Update();
//
//  auto legend = new TLegend(0.1,0.7,0.48,0.9);
//  legend->SetHeader("Channel","C"); // option "C" allows to center the header
//  legend->AddEntry(gtr, "Track", "p");
//  legend->AddEntry(gsh, "Shower", "p");
//  legend->Draw();

}
