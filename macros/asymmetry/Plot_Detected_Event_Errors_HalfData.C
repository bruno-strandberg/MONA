#include "TSystem.h"
#include "TROOT.h"
#include "TMath.h"
#include "TH2.h"
#include "TH3.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TStyle.h"

#include "DetResponse.h"
#include "NMHUtils.h"
#include "SummaryParser.h"
#include "SummaryEvent.h"
#include "HelperFunctions.C"

#include <iostream>
using namespace std;

void Plot_Detected_Event_Errors_HalfData() {
  TString filefolder = "./default_detres/";
  TString file_all = filefolder + "default_expectated_evts_NO.root";
  TString file_gt  = filefolder + "/use_half_data_gt300/default_expectated_evts_NO.root";
  TString file_lt  = filefolder + "/use_half_data_lt300/default_expectated_evts_NO.root";
  TString output   = filefolder + "plot_expectation_errors.root";
  
  TFile *f_all = TFile::Open(file_all, "READ");
  TFile *f_gt  = TFile::Open(file_gt, "READ");
  TFile *f_lt  = TFile::Open(file_lt, "READ");

  std::tuple<TH2D*, TH2D*, TH2D*> h_tuple_all = ReadDetectorResponseFile(file_all);
  std::tuple<TH2D*, TH2D*, TH2D*> h_tuple_gt  = ReadDetectorResponseFile(file_gt);
  std::tuple<TH2D*, TH2D*, TH2D*> h_tuple_lt  = ReadDetectorResponseFile(file_lt);

  TH2D *h_t_all = std::get<0>(h_tuple_all);
  TH2D *h_t_gt  = std::get<0>(h_tuple_gt);
  TH2D *h_t_lt  = std::get<0>(h_tuple_lt);
  TH2D *h_s_all = std::get<1>(h_tuple_all);
  TH2D *h_s_gt  = std::get<1>(h_tuple_gt);
  TH2D *h_s_lt  = std::get<1>(h_tuple_lt);

  Int_t nbinsX = h_t_all->GetXaxis()->GetNbins() + 1;
  Int_t nbinsY = h_t_all->GetYaxis()->GetNbins() + 1;
  Int_t nbins = nbinsX * nbinsY;

  TH2D *h2_t_err_all = GetErrorHisto2D(h_t_all);
  TH2D *h2_t_err_gt  = GetErrorHisto2D(h_t_gt);
  TH2D *h2_t_err_lt  = GetErrorHisto2D(h_t_lt);
  std::vector<TH1D*> h1_proj_t_err = VectorProjectionY2D(h2_t_err_all);
  std::vector<TH1D*> h1_proj_t_gt  = VectorProjectionY2D(h2_t_err_gt );
  std::vector<TH1D*> h1_proj_t_lt  = VectorProjectionY2D(h2_t_err_lt );

  TH2D *h2_s_err_all = GetErrorHisto2D(h_s_all);
  TH2D *h2_s_err_gt  = GetErrorHisto2D(h_s_gt);
  TH2D *h2_s_err_lt  = GetErrorHisto2D(h_s_lt);
  std::vector<TH1D*> h1_proj_s_err = VectorProjectionY2D(h2_s_err_all);
  std::vector<TH1D*> h1_proj_s_gt  = VectorProjectionY2D(h2_s_err_gt );
  std::vector<TH1D*> h1_proj_s_lt  = VectorProjectionY2D(h2_s_err_lt );
  
  gStyle->SetOptStat(0);
//  gStyle->SetTitleSize(0.3);
  TCanvas *c1 = new TCanvas("c1", "c1", 1200, 1000); // has to be outside of the if, so that the writing below works.
  c1->Divide(2,2);
  Double_t alpha = 0.2;
  // tracks
  for(unsigned int i = 1; i < h1_proj_t_err.size(); i++) {
    c1->cd(1);
    h1_proj_t_gt[i]->Divide(h1_proj_t_err[i]);
    h1_proj_t_gt[i]->SetLineColorAlpha(kBlue+2, alpha);
    if (i == 0) { 
      h1_proj_t_gt[i]->Draw(""); 
    }
    h1_proj_t_gt[i]->Draw("same");
    h1_proj_t_gt[i]->GetYaxis()->SetRangeUser(0, 3);
    h1_proj_t_gt[i]->GetXaxis()->SetRangeUser(-1, 0);
    h1_proj_t_gt[i]->SetTitle("Ratio of the errors on a full data set and half data set [gt300 data / all data, track]");
    h1_proj_t_gt[i]->SetTitleSize(12);
  }
  for(unsigned int i = 1; i < h1_proj_t_err.size(); i++) {
    c1->cd(2);
    h1_proj_t_lt[i]->Divide(h1_proj_t_err[i]);
    h1_proj_t_lt[i]->SetLineColorAlpha(kYellow+2, alpha);
    if (i == 0) { 
      h1_proj_t_gt[i]->Draw(""); 
    }
    h1_proj_t_gt[i]->Draw("same");
    h1_proj_t_gt[i]->GetYaxis()->SetRangeUser(0, 3);
    h1_proj_t_gt[i]->GetXaxis()->SetRangeUser(-1, 0);
    h1_proj_t_gt[i]->SetTitle("Ratio of the errors on a full data set and half data set [lt300 data / all data, track]");
  }
  // showers
  for(unsigned int i = 1; i < h1_proj_t_err.size(); i++) {
    c1->cd(3);
    h1_proj_s_gt[i]->Divide(h1_proj_s_err[i]);
    h1_proj_s_gt[i]->SetLineColorAlpha(kBlue+2, alpha);
    if (i == 0) { 
      h1_proj_s_gt[i]->Draw(""); 
    }
    h1_proj_s_gt[i]->Draw("same");
    h1_proj_s_gt[i]->GetYaxis()->SetRangeUser(0, 3);
    h1_proj_s_gt[i]->GetXaxis()->SetRangeUser(-1, 0);
    h1_proj_s_gt[i]->SetTitle("Ratio of the errors on a full data set and half data set [gt300 data / all data, shower]");
  }
  for(unsigned int i = 1; i < h1_proj_t_err.size(); i++) {
    c1->cd(4);
    h1_proj_s_lt[i]->Divide(h1_proj_s_err[i]);
    h1_proj_s_lt[i]->SetLineColorAlpha(kYellow+2, alpha);
    if (i == 0) { 
      h1_proj_s_gt[i]->Draw(""); 
    }
    h1_proj_s_gt[i]->Draw("same");
    h1_proj_s_gt[i]->GetYaxis()->SetRangeUser(0, 3);
    h1_proj_s_gt[i]->GetXaxis()->SetRangeUser(-1, 0);
    h1_proj_s_gt[i]->SetTitle("Ratio of the errors on a full data set and half data set [lt300 data / all data, shower]");
  }

  c1->SaveAs(filefolder + "plot_expectation_error.pdf");
}
