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

// This script calculates the asymmetry for 10 binned data in pid_detres.

void plot_asymmetry_error_one_bin(Int_t bin=700) {

  bool print = false;
  bool normalize = true;
  std::vector<Int_t> N_PID_CLASSES = {5, 10, 15, 20, 25, 30, 35, 40};
  std::vector<Double_t> N_PID_CLASSES_D(N_PID_CLASSES.begin(), N_PID_CLASSES.end())  ;
  Double_t Q_CUTOFF = 0.6;

  std::vector<Double_t> average_tr_errors;
  std::vector<Double_t> average_sh_errors;
  for (auto n_pid: N_PID_CLASSES) {

    Double_t PID_step = 1 / float(n_pid);
    TString basefolder = std::getenv("NMHDIR");
    TString filefolder = TString::Format(basefolder + "/fitter/pid_detres/simplified_error_calc/pid_binning_%i/", n_pid);
  
    std::vector<Double_t> asym_tr_err(N_PID_CLASSES.size());
    std::vector<Double_t> asym_sh_err(N_PID_CLASSES.size());
    Double_t total_tr_err = 0; // Total for errors
    Double_t total_sh_err = 0;
    Int_t total_tr = 0; // Total for number of bins considered
    Int_t total_sh = 0; // Total for number of bins considered
  
    for (int i = 0; i < n_pid; i++) {
      TString input_file = filefolder + TString::Format("asym_pid_%.2f.root", PID_step * i);
      
      TFile *f = TFile::Open(input_file, "READ");
  
      TH2D* h_asym_tr = (TH2D*)f->Get(TString::Format("asymmetry_track_%.2f", PID_step * i));
      TH2D* h_asym_sh = (TH2D*)f->Get(TString::Format("asymmetry_shower_%.2f", PID_step * i));
  
      Double_t tr_error = h_asym_tr->GetBinError(bin);
      Double_t sh_error = h_asym_sh->GetBinError(bin);
      if (PID_step * i < Q_CUTOFF ) { // If below 0.6 take shower, else track.
        asym_sh_err.push_back(sh_error);
        total_sh_err = total_sh_err + sh_error;
        total_tr++;
      } else {
        asym_tr_err.push_back(tr_error);
        total_tr_err = total_tr_err + tr_error;
        total_sh++;
      }
    }
    average_tr_errors.push_back(total_tr_err / total_tr);
    average_sh_errors.push_back(total_sh_err / total_sh);
  }
  if (print) {
    for (unsigned int i = 0; i < N_PID_CLASSES.size(); i++) {
      cout << average_tr_errors[i] << endl;
    }
    for (unsigned int i = 0; i < N_PID_CLASSES.size(); i++) {
      cout << average_sh_errors[i] << endl;
    }
  }
  if (normalize) {
    Double_t tr_start = average_tr_errors[0];
    Double_t sh_start = average_sh_errors[0];
    for (unsigned int i = 0; i < N_PID_CLASSES.size(); i++) {
      average_tr_errors[i] = average_tr_errors[i] / tr_start;
      average_sh_errors[i] = average_sh_errors[i] / sh_start;
    }
  }


  // Plotting
  TCanvas *c1 = new TCanvas("c1", "c1", 200, 10, 700, 500);
  c1->SetGrid();

  const int n = N_PID_CLASSES_D.size();
  // track
  TGraph *gtr = new TGraph(n);
  for (int i=0; i < n; ++i) {
    gtr->SetPoint(i, N_PID_CLASSES_D.at(i), average_tr_errors.at(i));
  }
  gtr->Draw("AL");
  gtr->SetTitle(TString::Format("Error on bin %i as function of number of PID bins", bin));
  gtr->GetXaxis()->SetTitle("Number of bins used");
  gtr->GetYaxis()->SetTitle("Asymmetry error");
  gtr->SetLineColor(kYellow+2);
  gtr->SetMarkerColor(kYellow+2);
  gtr->SetMarkerStyle(4);
  gtr->SetMarkerSize(1);
  //gtr->GetYaxis()->SetRangeUser(0, 5e-3);
  gtr->GetYaxis()->SetRangeUser(0, 4);
  c1->Update();

  // shower
  TGraph *gsh = new TGraph(n);
  for (int i=0; i < n; ++i) {
    gsh->SetPoint(i, N_PID_CLASSES_D.at(i), average_sh_errors.at(i));
  }
  gsh->Draw("L");
  gsh->SetLineColor(kBlue+2);
  gsh->SetMarkerColor(kBlue+2);
  gsh->SetMarkerStyle(4);
  c1->Update();

  auto legend = new TLegend(0.1,0.7,0.48,0.9);
  legend->SetHeader("Channel","C"); // option "C" allows to center the header
  legend->AddEntry(gtr, "Track", "p");
  legend->AddEntry(gsh, "Shower", "p");
  legend->Draw();

}
