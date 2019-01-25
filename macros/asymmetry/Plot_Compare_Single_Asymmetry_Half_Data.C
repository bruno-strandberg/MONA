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
#include "HelperFunctions.C"
#include "NMHUtils.h"
#include "SummaryParser.h"
#include "SummaryEvent.h"

#include <iostream>
using namespace std;

void Plot_Compare_Single_Asymmetry_Half_Data() {

  TString filefolder = "./default_detres/";
  TString s_alldata  = filefolder + "use_half_data_all/simple_error/asymmetry_default.root";
  TString s_data_gt  = filefolder + "use_half_data_gt300/simple_error/asymmetry_default.root";
  TString s_data_lt  = filefolder + "use_half_data_lt300/simple_error/asymmetry_default.root";
  
  TFile *f_all = TFile::Open(s_alldata, "READ");
  TFile *f_gt  = TFile::Open(s_data_gt, "READ");
  TFile *f_lt  = TFile::Open(s_data_lt, "READ");

  gStyle->SetPalette(kBird);
//gStyle->SetOptStat(0);

  TH2D *h_t_all = (TH2D*)f_all->Get("asymmetry_track");
  TH2D *h_t_gt  = (TH2D*)f_gt ->Get("asymmetry_track");
  TH2D *h_t_lt  = (TH2D*)f_lt ->Get("asymmetry_track");
  TH2D *h_s_all = (TH2D*)f_all->Get("asymmetry_shower");
  TH2D *h_s_gt  = (TH2D*)f_gt ->Get("asymmetry_shower");
  TH2D *h_s_lt  = (TH2D*)f_lt ->Get("asymmetry_shower");

  TH1D *h_t_gt_out = new TH1D("htgt", "#splitline{Relative error on track asymmerty}{compared to using half data}", 30, 0, 3);
  TH1D *h_t_lt_out = new TH1D("htlt", "#splitline{Relative error on track asymmerty}{compared to using half data}", 30, 0, 3);
  TH1D *h_s_gt_out = new TH1D("hsgt", "#splitline{Relative error on shower asymmerty}{compared to using half data}", 30, 0, 3);
  TH1D *h_s_lt_out = new TH1D("hslt", "#splitline{Relative error on shower asymmerty}{compared to using half data}", 30, 0, 3);

  for (Int_t xb = 1; xb <= h_s_all->GetXaxis()->GetNbins(); xb++) {
    for (Int_t yb = 1; yb <= h_s_all->GetYaxis()->GetNbins(); yb++) {
      Double_t error_all = h_t_all->GetBinError(xb,yb);
      Double_t error_gt  = h_t_gt ->GetBinError(xb,yb);
      Double_t error_lt  = h_t_lt ->GetBinError(xb,yb);

      if (error_all != 0) {
        error_gt = error_gt / error_all;
        error_lt = error_lt / error_all;
      } else {
        error_gt = 0;
        error_lt = 0;
      }

      //std::cout << "Error in bin " << xb << " " << yb << " is " << error_all << std::endl;
      //std::cout << "Rel. error in bin " << xb << " " << yb << " compared to gt is " << error_gt << std::endl;
      //std::cout << "Rel. error in bin " << xb << " " << yb << " compared to lt is " << error_lt << std::endl;
      
      if (error_gt != 0) h_t_gt_out->Fill(error_gt);
      if (error_lt != 0) h_t_lt_out->Fill(error_lt);
    }
  }
  for (Int_t xb = 1; xb <= h_s_all->GetXaxis()->GetNbins(); xb++) {
    for (Int_t yb = 1; yb <= h_s_all->GetYaxis()->GetNbins(); yb++) {
      Double_t error_all = h_s_all->GetBinError(xb,yb);
      Double_t error_gt  = h_s_gt ->GetBinError(xb,yb);
      Double_t error_lt  = h_s_lt ->GetBinError(xb,yb);

      if (error_all != 0) {
        error_gt = error_gt / error_all;
        error_lt = error_lt / error_all;
      } else {
        error_gt = 0;
        error_lt = 0;
      }

      //std::cout << "Error in bin " << xb << " " << yb << " is " << error_all << std::endl;
      //std::cout << "Rel. error in bin " << xb << " " << yb << " compared to gt is " << error_gt << std::endl;
      //std::cout << "Rel. error in bin " << xb << " " << yb << " compared to lt is " << error_lt << std::endl;
      
      if (error_gt != 0) h_s_gt_out->Fill(error_gt);
      if (error_lt != 0) h_s_lt_out->Fill(error_lt);
    }
  }

  TCanvas* c1 = new TCanvas("c1", "c1", 1200, 600);
  c1->Divide(2,1);
  c1->cd(1);
  h_t_gt_out->Draw();
  h_t_gt_out->GetXaxis()->SetTitle("Relative error");
  Double_t max_t = 1.1 * std::max(h_t_gt_out->GetMaximum(), h_t_lt_out->GetMaximum());
  h_t_gt_out->GetYaxis()->SetRangeUser(0, max_t);
  h_t_lt_out->SetLineStyle(2);
  h_t_lt_out->Draw("same");
  TLegend* leg_t = new TLegend(0.1, 0.7, 0.4, 0.8);
  leg_t->AddEntry(h_t_gt_out, "Data runID > 300", "l");
  leg_t->AddEntry(h_t_lt_out, "Data runID #leq 300", "l");
  leg_t->Draw();

  c1->cd(2);
  h_s_gt_out->Draw();
  h_s_gt_out->GetXaxis()->SetTitle("Relative error");
  Double_t max_s = 1.1 * std::max(h_s_gt_out->GetMaximum(), h_s_lt_out->GetMaximum());
  h_s_gt_out->GetYaxis()->SetRangeUser(0, max_s);
  h_s_lt_out->SetLineStyle(2);
  h_s_lt_out->Draw("same");
  TLegend* leg_s = new TLegend(0.1, 0.7, 0.4, 0.8);
  leg_s->AddEntry(h_s_gt_out, "Data runID > 300", "l");
  leg_s->AddEntry(h_s_lt_out, "Data runID #leq 300", "l");
  leg_s->Draw();
 
}
