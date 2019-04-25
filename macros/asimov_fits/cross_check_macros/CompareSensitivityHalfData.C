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
#include "HelperFunctions.h"
#include "NMHUtils.h"
#include "SummaryParser.h"
#include "SummaryEvent.h"

#include <iostream>
using namespace std;


void CompareSensitivityHalfData() {

  bool print = false;
  bool plot = true;
  
  TFile *full = TFile::Open("output/root/AsimovFitIOTh23Range.root", "READ");
  TFile *gt   = TFile::Open("output/root/AsimovFitIOTh23Range_half_data_gt300.root", "READ");

  gStyle->SetPalette(kLightTemperature);
  gStyle->SetOptStat(0);

  TH2D *h_sens_full = (TH2D*)full->Get("sens_track_41");
  TH2D *h_sens_gt   = (TH2D*)gt->Get("sens_track_41");

  TH2D* h_sens_div = (TH2D*)h_sens_gt->Clone();
  h_sens_div->Divide(h_sens_full);

  for (Int_t xb = 1; xb <= h_sens_full->GetXaxis()->GetNbins(); xb++) {
    for (Int_t yb = 1; yb <= h_sens_full->GetYaxis()->GetNbins(); yb++) {

      Double_t S_full = h_sens_full->GetBinContent(xb, yb);
      Double_t S_gt   = h_sens_gt->GetBinContent(xb, yb);

      if ((S_gt / S_full > 1.03) or (S_gt / S_full < 0.97)) {
        if (print) {
          cout << "S track full " << S_full << endl;
          cout << "S track gt   " << S_gt << endl;
          cout << "Ratio track  " << S_gt / S_full << endl;
        }
      }
    }
  }

  TCanvas *c1 = new TCanvas("c1", "c1", 1200, 1000); 
  c1->Divide(2,2);
  c1->cd(1);
  h_sens_full->Draw("colz");
  h_sens_full->GetZaxis()->SetRangeUser(0,0.35); 

  c1->cd(2);
  h_sens_gt->Draw("colz");
  h_sens_gt->GetZaxis()->SetRangeUser(0,0.35); 

  c1->cd(3);
  h_sens_div->Draw("colz");
  h_sens_div->GetZaxis()->SetRangeUser(0,2); 
  

//  h2->SetTitle(Form("Track events [%s]", order_string.c_str()));
//  h2->SetXTitle("E_{reco} [GeV]");
//  h2->SetYTitle("cos(#theta_{reco})");
//
//  MakePlotNice(0.035, h2);
//
//  h2->GetXaxis()->SetRangeUser(3,100);
//  h2->GetYaxis()->SetRangeUser(-1,0);
//  c1->cd(1)->SetLogx();
//
//  TH2D* h3 = GetRelativeErrorHistogram(h2);
//  c1->cd(2);
//  h3->Draw("colz");
//
//  h3->SetTitle(Form("Error track events [%s]", order_string.c_str()));
//  h3->SetXTitle("E_{reco} [GeV]");
//  h3->SetYTitle("cos(#theta_{reco})");
//
//  MakePlotNice(0.035, h3);
//
//  h3->GetXaxis()->SetRangeUser(3,100);
//  h3->GetYaxis()->SetRangeUser(-1,0);
//  h3->GetZaxis()->SetRangeUser(0,0.2); 
//  c1->cd(2)->SetLogx();
//
//  TH3D *h4 = (TH3D*)f_IO->Get("detected_showers");
//  TH2D *h5 = (TH2D*)h4->Project3D("yx");
//  c1->cd(3);
//  h5->Draw("colz");
//
//  h5->SetTitle(Form("Shower events [%s]", order_string.c_str()));
//  h5->SetXTitle("E_{reco} [GeV]");
//  h5->SetYTitle("cos(#theta_{reco})");
//
//  MakePlotNice(0.035, h5);
//
//  h5->GetXaxis()->SetRangeUser(3,100);
//  h5->GetYaxis()->SetRangeUser(-1,0);
//  c1->cd(3)->SetLogx();
//
//  TH2D* h6 = GetRelativeErrorHistogram(h5);
//  c1->cd(4);
//  h6->Draw("colz");
//
//  h6->SetTitle(Form("Error track events [%s]", order_string.c_str()));
//  h6->SetXTitle("E_{reco} [GeV]");
//  h6->SetYTitle("cos(#theta_{reco})");
//
//  MakePlotNice(0.035, h6);
//
//  h6->GetXaxis()->SetRangeUser(3,100);
//  h6->GetYaxis()->SetRangeUser(-1,0);
//  h6->GetZaxis()->SetRangeUser(0,0.2); 
//  c1->cd(4)->SetLogx();
//
//  TFile fout(Form("./default_detres/plot_event_errors_%s.root", order_string.c_str()), "RECREATE");
//  h2->Write();
//  h3->Write();
//  h5->Write();
//  h6->Write();
//  fout.Close();
//  c1->SaveAs(Form("./default_detres/plot_event_errors_plots_%s.pdf", order_string.c_str()));
}
