#include "TSystem.h"
#include "TROOT.h"
#include "TMath.h"
#include "TH2.h"
#include "TH3.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TStyle.h"

#include "DetResponse.h"
#include "FitFunction.h"
#include "NMHUtils.h"
#include "SummaryParser.h"
#include "SummaryEvent.h"
#include "HelperFunctions.C"

#include <iostream>
using namespace std;

void plot_relative_difference_asym_error() {
    TString filefolder = "./default_detres/";
    TString file_correlated_error = filefolder + "asym_pid.root";
    TString file_simple_error     = filefolder + "simplified_error_calc/asym_pid.root";
    TString output   = filefolder + "plot_relative_difference_asym_error.root";
    
    TFile *f_c_err = TFile::Open(file_correlated_error, "READ"); // Correlated error
    TFile *f_s_err  = TFile::Open(file_simple_error, "READ"); // Simplified error

    TH2D* h_asym_c_err_tr = (TH2D*)f_c_err->Get("asymmetry_track");
    TH2D* h_asym_s_err_tr = (TH2D*)f_s_err->Get("asymmetry_track");
    TH2D* h_asym_c_err_sh = (TH2D*)f_c_err->Get("asymmetry_shower");
    TH2D* h_asym_s_err_sh = (TH2D*)f_s_err->Get("asymmetry_shower");

//    Int_t nbinsX = h_t_all->GetXaxis()->GetNbins() + 1;
//    Int_t nbinsY = h_t_all->GetYaxis()->GetNbins() + 1;
//    Int_t nbins = nbinsX * nbinsY;

    TH2D *h_err_asym_c_err_tr  = GetErrorHisto2D(h_asym_c_err_tr);
    TH2D *h_err_asym_s_err_tr  = GetErrorHisto2D(h_asym_s_err_tr);
    TH2D *h_err_asym_c_err_sh  = GetErrorHisto2D(h_asym_c_err_sh);
    TH2D *h_err_asym_s_err_sh  = GetErrorHisto2D(h_asym_s_err_sh);
    std::vector<TH1D*> h1_proj_c_err_tr = VectorProjectionY2D(h_err_asym_c_err_tr); // Histogram containing errors in the bins for correlated error calculation for tracks: h_err_asym_c_err_tr
    std::vector<TH1D*> h1_proj_s_err_tr = VectorProjectionY2D(h_err_asym_s_err_tr);
    std::vector<TH1D*> h1_proj_c_err_sh = VectorProjectionY2D(h_err_asym_c_err_sh);
    std::vector<TH1D*> h1_proj_s_err_sh = VectorProjectionY2D(h_err_asym_s_err_sh);
    
    gStyle->SetOptStat(0);
//    gStyle->SetTitleSize(0.3);
    TCanvas *c1 = new TCanvas("c1", "c1", 1200, 500); // has to be outside of the if, so that the writing below works.
    c1->Divide(2,1);
    Double_t alpha = 0.5;
    // tracks
    for(unsigned int i = 1; i < h1_proj_c_err_tr.size(); i++) {
      c1->cd(1);
      h1_proj_s_err_tr[i]->Divide(h1_proj_c_err_tr[i]);
      h1_proj_s_err_tr[i]->SetLineColorAlpha(kBlue+2, alpha);
      if (i == 0) { 
        h1_proj_s_err_tr[i]->Draw(""); 
      }
      h1_proj_s_err_tr[i]->Draw("same");
      h1_proj_s_err_tr[i]->GetYaxis()->SetRangeUser(0, 10);
      h1_proj_s_err_tr[i]->GetXaxis()->SetRangeUser(-1, 0);
      h1_proj_s_err_tr[i]->SetTitle("#splitline{Ratio of the error in asymmetry calculation per bin}{track [simple error / correlated error]}");
      h1_proj_s_err_tr[i]->SetTitleSize(12);
      h1_proj_s_err_tr[i]->GetXaxis()->SetTitle("ct");
      h1_proj_s_err_tr[i]->GetYaxis()->SetTitle("Ratio of Simple error / Correlated error");
    }
    // showers
    for(unsigned int i = 1; i < h1_proj_c_err_sh.size(); i++) {
      c1->cd(2);
      h1_proj_s_err_sh[i]->Divide(h1_proj_c_err_sh[i]);
      h1_proj_s_err_sh[i]->SetLineColorAlpha(kBlue+2, alpha);
      if (i == 0) { 
        h1_proj_s_err_sh[i]->Draw(""); 
      }
      h1_proj_s_err_sh[i]->Draw("same");
      h1_proj_s_err_sh[i]->GetYaxis()->SetRangeUser(0, 10);
      h1_proj_s_err_sh[i]->GetXaxis()->SetRangeUser(-1, 0);
      h1_proj_s_err_sh[i]->SetTitle("#splitline{Ratio of the error in asymmetry calculation per bin}{shower [simple error / correlated error]}");
      h1_proj_s_err_sh[i]->SetTitleSize(12);
      h1_proj_s_err_sh[i]->GetXaxis()->SetTitle("ct");
      h1_proj_s_err_sh[i]->GetYaxis()->SetTitle("Ratio of Simple error / Correlated error");
    }
    c1->SaveAs(filefolder + "simplified_error_calc/relative_change_error_asym_per_bin.pdf");
}
