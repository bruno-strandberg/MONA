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

void asymmetry_plots_10bins() {

  for (int i = 0; i < 10; i++){
    TString input = Form("./pid_detres/asym_pid_0.%i.root", i);
    
    TFile *f = TFile::Open(input, "READ");

    gStyle->SetPalette(kLightTemperature);
    gStyle->SetOptStat(0);

    TH2D *h1 = (TH2D*)f->Get(Form("asymmetry_track_%i", i));
    TCanvas *c1 = new TCanvas(Form("c%i", i), Form("c%i", i), 630, 500); 
    h1->Draw("colz");

    h1->SetTitle(Form("Asymmetry tracks q_0.%i", i));
    h1->SetXTitle("E_{reco} [GeV]");
    h1->SetYTitle("cos(#theta_{reco})");

    float size = 0.035;
    h1->GetXaxis()->SetTitleSize(size);
    h1->GetXaxis()->SetTitleOffset(1.2);
    h1->GetXaxis()->SetTitleFont(62);
    h1->GetXaxis()->SetLabelFont(62);
    h1->GetXaxis()->SetLabelSize(size);
    h1->GetYaxis()->SetTitleSize(size);
    h1->GetYaxis()->SetTitleOffset(1.2);
    h1->GetYaxis()->SetTitleFont(62);
    h1->GetYaxis()->SetLabelFont(62);
    h1->GetYaxis()->SetLabelSize(size);

    h1->GetXaxis()->SetRangeUser(3,100);
    h1->GetYaxis()->SetRangeUser(-1,0);
    h1->GetZaxis()->SetRangeUser(-1,1); 
    h1->GetZaxis()->SetRangeUser(-1,1);
    c1->SetLogx();
    h1->Draw("colz");
    c1->SaveAs(Form("./pid_detres/asymmetry_tracks_0.%i.pdf", i));
    c1->SaveAs(Form("./pid_detres/png/asymmetry_tracks_0.%i.png", i));

    TH2D *h2 = (TH2D*)f->Get(Form("asymmetry_shower_%i", i));
    TCanvas *c2 = new TCanvas(Form("c%i", i+10), Form("c%i", i+10), 630, 500); 
    h2->Draw("colz");

    h2->SetTitle(Form("Asymmetry showers q_0.%i", i));
    h2->SetXTitle("E_{reco} [GeV]");
    h2->SetYTitle("cos(#theta_{reco})");

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
    h2->GetZaxis()->SetRangeUser(-1,1); 
    h2->GetZaxis()->SetRangeUser(-1,1);
    c2->SetLogx();
    h2->Draw("colz");
    c2->SaveAs(Form("./pid_detres/asymmetry_showers_0.%i.pdf", i));
    c2->SaveAs(Form("./pid_detres/png/asymmetry_showers_0.%i.png", i));
  }
}
