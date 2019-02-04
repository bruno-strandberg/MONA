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

TH2D* GetRelativeErrorHistogram(TH2D* h1) { 
  TH2D* h2 = (TH2D*)h1->Clone();
  TString name = h1->GetName();
  h2->SetNameTitle("errorhist_"+name, "errorhist_"+name);
  h2->Reset();
  for (Int_t xb = 1; xb <= h1->GetXaxis()->GetNbins(); xb++) {
    for (Int_t yb = 1; yb <= h1->GetYaxis()->GetNbins(); yb++) {
      Double_t h1_err = h1->GetBinError(xb, yb);
      Double_t h1_cont= h1->GetBinContent(xb, yb);
      Double_t h1_err_relative;
      if (h1_cont != 0){
        h1_err_relative = h1_err / h1_cont;
      }
      else {
        h1_err_relative = 0;
      }
      h2->SetBinContent(xb, yb, h1_err_relative);
    }
  }
  return h2;
}

void PlotDetectedEventErrors(bool mass_ordering=true) {

  string order_string;
  if (mass_ordering) {
    order_string = "NO";
  }
  else {
    order_string = "IO";
  }

  TString input = Form("./default_detres/default_expected_evts_%s.root", order_string.c_str());
  
  TFile *f_IO = TFile::Open(input, "READ");

  gStyle->SetPalette(kBird);
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

  MakePlotNice(0.035, h2);

  h2->GetXaxis()->SetRangeUser(3,100);
  h2->GetYaxis()->SetRangeUser(-1,0);
  c1->cd(1)->SetLogx();

  TH2D* h3 = GetRelativeErrorHistogram(h2);
  c1->cd(2);
  h3->Draw("colz");

  h3->SetTitle(Form("Error track events [%s]", order_string.c_str()));
  h3->SetXTitle("E_{reco} [GeV]");
  h3->SetYTitle("cos(#theta_{reco})");

  MakePlotNice(0.035, h3);

  h3->GetXaxis()->SetRangeUser(3,100);
  h3->GetYaxis()->SetRangeUser(-1,0);
  h3->GetZaxis()->SetRangeUser(0,0.2); 
  c1->cd(2)->SetLogx();

  TH3D *h4 = (TH3D*)f_IO->Get("detected_showers");
  TH2D *h5 = (TH2D*)h4->Project3D("yx");
  c1->cd(3);
  h5->Draw("colz");

  h5->SetTitle(Form("Shower events [%s]", order_string.c_str()));
  h5->SetXTitle("E_{reco} [GeV]");
  h5->SetYTitle("cos(#theta_{reco})");

  MakePlotNice(0.035, h5);

  h5->GetXaxis()->SetRangeUser(3,100);
  h5->GetYaxis()->SetRangeUser(-1,0);
  c1->cd(3)->SetLogx();

  TH2D* h6 = GetRelativeErrorHistogram(h5);
  c1->cd(4);
  h6->Draw("colz");

  h6->SetTitle(Form("Error track events [%s]", order_string.c_str()));
  h6->SetXTitle("E_{reco} [GeV]");
  h6->SetYTitle("cos(#theta_{reco})");

  MakePlotNice(0.035, h6);

  h6->GetXaxis()->SetRangeUser(3,100);
  h6->GetYaxis()->SetRangeUser(-1,0);
  h6->GetZaxis()->SetRangeUser(0,0.2); 
  c1->cd(4)->SetLogx();

  TFile fout(Form("./default_detres/plot_event_errors_%s.root", order_string.c_str()), "RECREATE");
  h2->Write();
  h3->Write();
  h5->Write();
  h6->Write();
  fout.Close();
  c1->SaveAs(Form("./default_detres/plot_event_errors_plots_%s.pdf", order_string.c_str()));
}
