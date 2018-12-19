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
#include "FitFunction.h"
#include "NMHUtils.h"
#include "SummaryParser.h"
#include "SummaryEvent.h"

#include <iostream>
using namespace std;

void error_plot(bool mass_ordering=true) {

  string order_string;
  if (mass_ordering) {
    order_string = "NO";
  }
  else {
    order_string = "IO";
  }

  TString input = Form("./default_detres/timing_%s.root", order_string.c_str());
  
  TFile *f_IO = TFile::Open(input, "READ");

  gStyle->SetPalette(kLightTemperature);
  gStyle->SetOptStat(0);

  TH3D *h1 = (TH3D*)f_IO->Get("detected_tracks");
  TCanvas *c1 = new TCanvas("c1", "c1", 1200, 1000); 
  c1->Divide(2,2);
  TH2D *h2 = (TH2D*)h1->Project3D("yx");
  c1->cd(1);
  h2->Draw("colz");

  h2->SetTitle(Form("Track events [%s]", order_string.c_str()));
  h2->SetXTitle("E_{reco} [GeV]");
  h2->SetYTitle("cos(#theta_{reco})");

  float size = 0.035;
  h2->GetXaxis()->SetTitleSize(size);
  h2->GetXaxis()->SetTitleOffset(1.2);
  h2->GetXaxis()->SetTitleFont(62);
  h2->GetXaxis()->SetLabelFont(62);
  h2->GetXaxis()->SetLabelSize(size);
  h2->GetYaxis()->SetTitleSize(size);
  h2->GetYaxis()->SetTitleOffset(1.2);
  h2->GetYaxis()->SetTitleFont(62);
  h2->GetYaxis()->SetLabelFont(62);
  h2->GetYaxis()->SetLabelSize(size);

  h2->GetXaxis()->SetRangeUser(3,100);
  h2->GetYaxis()->SetRangeUser(-1,0);
  c1->cd(1)->SetLogx();

  TH2D *h3 = (TH2D*)h2->Clone();
  h3->SetNameTitle("errorhist_track", "errorhist_track");
  h3->Reset();
  for (Int_t xb = 1; xb <= h3->GetXaxis()->GetNbins(); xb++) {
    for (Int_t yb = 1; yb <= h3->GetYaxis()->GetNbins(); yb++) {
      Double_t h2_err = h2->GetBinError(xb, yb);
      Double_t h2_cont=h2->GetBinContent(xb, yb);
      Double_t h2_err_relative;
      if (h2_cont != 0){
        h2_err_relative = h2_err / h2->GetBinContent(xb, yb);
      }
      else {
        h2_err_relative = 0;
      }
      h3->SetBinContent(xb, yb, h2_err_relative);
    }
  }
  c1->cd(2);
  h3->Draw("colz");

  h3->SetTitle(Form("Error track events [%s]", order_string.c_str()));
  h3->SetXTitle("E_{reco} [GeV]");
  h3->SetYTitle("cos(#theta_{reco})");

  h3->GetXaxis()->SetTitleSize(size);
  h3->GetXaxis()->SetTitleOffset(1.2);
  h3->GetXaxis()->SetTitleFont(62);
  h3->GetXaxis()->SetLabelFont(62);
  h3->GetXaxis()->SetLabelSize(size);
  h3->GetYaxis()->SetTitleSize(size);
  h3->GetYaxis()->SetTitleOffset(1.2);
  h3->GetYaxis()->SetTitleFont(62);
  h3->GetYaxis()->SetLabelFont(62);
  h3->GetYaxis()->SetLabelSize(size);

  h3->GetXaxis()->SetRangeUser(3,100);
  h3->GetYaxis()->SetRangeUser(-1,0);
  h3->GetZaxis()->SetRangeUser(0,0.2); 
  c1->cd(2)->SetLogx();
//  h1->Draw("colz");

  TH3D *h4 = (TH3D*)f_IO->Get("detected_showers");
  TH2D *h5 = (TH2D*)h4->Project3D("yx");
  c1->cd(3);
  h5->Draw("colz");

  h5->SetTitle(Form("Shower events [%s]", order_string.c_str()));
  h5->SetXTitle("E_{reco} [GeV]");
  h5->SetYTitle("cos(#theta_{reco})");

  h5->GetXaxis()->SetTitleSize(size);
  h5->GetXaxis()->SetTitleOffset(1.2);
  h5->GetXaxis()->SetTitleFont(62);
  h5->GetXaxis()->SetLabelFont(62);
  h5->GetXaxis()->SetLabelSize(size);
  h5->GetYaxis()->SetTitleSize(size);
  h5->GetYaxis()->SetTitleOffset(1.2);
  h5->GetYaxis()->SetTitleFont(62);
  h5->GetYaxis()->SetLabelFont(62);
  h5->GetYaxis()->SetLabelSize(size);

  h5->GetXaxis()->SetRangeUser(3,100);
  h5->GetYaxis()->SetRangeUser(-1,0);
  c1->cd(3)->SetLogx();

  TH2D *h6 = (TH2D*)h5->Clone();
  h6->SetNameTitle("errorhist_shower", "errorhist_shower");
  h6->Reset();
  for (Int_t xb = 1; xb <= h5->GetXaxis()->GetNbins(); xb++) {
    for (Int_t yb = 1; yb <= h5->GetYaxis()->GetNbins(); yb++) {
      Double_t h5_err  = h5->GetBinError(xb, yb);
      Double_t h5_cont = h5->GetBinContent(xb, yb);
      Double_t h5_err_relative;
      if (h5_cont != 0){
        h5_err_relative = h5_err / h5->GetBinContent(xb, yb);
      }
      else {
        h5_err_relative = 0;
      }
      h6->SetBinContent(xb, yb, h5_err_relative);
    }
  }
  c1->cd(4);
  h6->Draw("colz");

  h6->SetTitle(Form("Error track events [%s]", order_string.c_str()));
  h6->SetXTitle("E_{reco} [GeV]");
  h6->SetYTitle("cos(#theta_{reco})");

  h6->GetXaxis()->SetTitleSize(size);
  h6->GetXaxis()->SetTitleOffset(1.2);
  h6->GetXaxis()->SetTitleFont(62);
  h6->GetXaxis()->SetLabelFont(62);
  h6->GetXaxis()->SetLabelSize(size);
  h6->GetYaxis()->SetTitleSize(size);
  h6->GetYaxis()->SetTitleOffset(1.2);
  h6->GetYaxis()->SetTitleFont(62);
  h6->GetYaxis()->SetLabelFont(62);
  h6->GetYaxis()->SetLabelSize(size);

  h6->GetXaxis()->SetRangeUser(3,100);
  h6->GetYaxis()->SetRangeUser(-1,0);
  h6->GetZaxis()->SetRangeUser(0,0.2); 
  c1->cd(4)->SetLogx();

  TFile fout(Form("./default_detres/error_%s.root", order_string.c_str()), "RECREATE");
  h2->Write();
  h3->Write();
  h5->Write();
  h6->Write();
  fout.Close();
  c1->SaveAs(Form("./default_detres/error_plots_%s.pdf", order_string.c_str()));
}
