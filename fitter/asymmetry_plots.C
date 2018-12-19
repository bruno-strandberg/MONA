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

void asymmetry_plots() {

  TString input   = "./default_detres/asym_pid.root";
  
  TFile *f_IO = TFile::Open(input, "READ");

  gStyle->SetPalette(kLightTemperature);
  gStyle->SetOptStat(0);

  TH2D *h1 = (TH2D*)f_IO->Get("detected_showers_NO");
  TCanvas *c1 = new TCanvas("c1", "c1", 600, 500); 
  h1->Draw("colz");

  h1->SetTitle("Shower events [NO]");
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
//  h1->GetZaxis()->SetRangeUser(-1,1); 
  c1->SetLogx();
  h1->Draw("colz");
  c1->SaveAs("./default_detres/showers_NO.pdf");


  TH2D *h2 = (TH2D*)f_IO->Get("detected_showers_IO");
  TCanvas *c2 = new TCanvas("c2", "c2", 600, 500); 
  h2->Draw("colz");

  h2->SetTitle("Shower events [IO]");
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
//  h2->GetZaxis()->SetRangeUser(-1,1); 
  c2->SetLogx();
  h2->Draw("colz");
  c2->SaveAs("./default_detres/showers_IO.pdf");


  TH2D *h3 = (TH2D*)f_IO->Get("asymmetry_shower");
  TCanvas *c3 = new TCanvas("c3", "c3", 630, 500); 
  h3->Draw("colz");

  h3->SetTitle("Asymmetry showers");
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
  h3->GetZaxis()->SetRangeUser(-1,1); 
  c3->SetLogx();
  h3->Draw("colz");
  c3->SaveAs("./default_detres/asymmetry_showers.pdf");


  TH2D *h4 = (TH2D*)f_IO->Get("detected_tracks_NO");
  TCanvas *c4 = new TCanvas("c4", "c4", 600, 500); 
  h4->Draw("colz");

  h4->SetTitle("Track events [NO]");
  h4->SetXTitle("E_{reco} [GeV]");
  h4->SetYTitle("cos(#theta_{reco})");

  h4->GetXaxis()->SetTitleSize(size);
  h4->GetXaxis()->SetTitleOffset(1.2);
  h4->GetXaxis()->SetTitleFont(62);
  h4->GetXaxis()->SetLabelFont(62);
  h4->GetXaxis()->SetLabelSize(size);
  h4->GetYaxis()->SetTitleSize(size);
  h4->GetYaxis()->SetTitleOffset(1.2);
  h4->GetYaxis()->SetTitleFont(62);
  h4->GetYaxis()->SetLabelFont(62);
  h4->GetYaxis()->SetLabelSize(size);

  h4->GetXaxis()->SetRangeUser(3,100);
  h4->GetYaxis()->SetRangeUser(-1,0);
//  h4->GetZaxis()->SetRangeUser(-1,1); 
  c4->SetLogx();
  h4->Draw("colz");
  c4->SaveAs("./default_detres/tracks_NO.pdf");


  TH2D *h5 = (TH2D*)f_IO->Get("detected_tracks_IO");
  TCanvas *c5 = new TCanvas("c5", "c5", 600, 500); 
  h5->Draw("colz");

  h5->SetTitle("Track events [IO]");
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
//  h5->GetZaxis()->SetRangeUser(-1,1); 
  c5->SetLogx();
  h5->Draw("colz");
  c5->SaveAs("./default_detres/tracks_IO.pdf");


  TH2D *h6 = (TH2D*)f_IO->Get("asymmetry_track");
  TCanvas *c6 = new TCanvas("c6", "c6", 630, 500); 
  h6->Draw("colz");

  h6->SetTitle("Asymmetry tracks");
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
  h6->GetZaxis()->SetRangeUser(-1,1); 
  c6->SetLogx();
  h6->Draw("colz");
  c6->SaveAs("./default_detres/asymmetry_tracks.pdf");

}
