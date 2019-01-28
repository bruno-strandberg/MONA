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

  /* Script to plot the histogram for NO, IO and the asymmetry
   * for tracks and showers.
   */


void Plot_Asymmetry_Default() {

  TString input   = "./default_detres/asymmetry_default.root";
  
  TFile *f_IO = TFile::Open(input, "READ");

  gStyle->SetPalette(kBird);
  gStyle->SetOptStat(0);

  TH2D *h1 = (TH2D*)f_IO->Get("detected_showers_NO");
  TCanvas *c1 = new TCanvas("c1", "c1", 600, 500); 
  h1->Draw("colz");

  SetPlotTitle(h1, "Shower events [NO]");
  SetLabelSize(h1, 0.035);

  h1->GetXaxis()->SetRangeUser(3,100);
  h1->GetYaxis()->SetRangeUser(-1,0);
  c1->SetLogx();
  h1->Draw("colz");
  c1->SaveAs("./default_detres/showers_NO.pdf");


  TH2D *h2 = (TH2D*)f_IO->Get("detected_showers_IO");
  TCanvas *c2 = new TCanvas("c2", "c2", 600, 500); 
  h2->Draw("colz");

  SetPlotTitle(h2, "Shower events [IO]");
  SetLabelSize(h2, 0.035);

  h2->GetXaxis()->SetRangeUser(3,100);
  h2->GetYaxis()->SetRangeUser(-1,0);
  c2->SetLogx();
  h2->Draw("colz");
  c2->SaveAs("./default_detres/showers_IO.pdf");


  TH2D *h3 = (TH2D*)f_IO->Get("asymmetry_shower");
  TCanvas *c3 = new TCanvas("c3", "c3", 630, 500); 
  h3->Draw("colz");

  SetPlotTitle(h3, "Asymmetry showers");
  SetLabelSize(h3, 0.035);

  h3->GetXaxis()->SetRangeUser(3,100);
  h3->GetYaxis()->SetRangeUser(-1,0);
  h3->GetZaxis()->SetRangeUser(-1,1); 
  c3->SetLogx();
  h3->Draw("colz");
  c3->SaveAs("./default_detres/asymmetry_showers.pdf");


  TH2D *h4 = (TH2D*)f_IO->Get("detected_tracks_NO");
  TCanvas *c4 = new TCanvas("c4", "c4", 600, 500); 
  h4->Draw("colz");

  SetPlotTitle(h4, "Track events [NO]");
  SetLabelSize(h4, 0.035);

  h4->GetXaxis()->SetRangeUser(3,100);
  h4->GetYaxis()->SetRangeUser(-1,0);
  c4->SetLogx();
  h4->Draw("colz");
  c4->SaveAs("./default_detres/tracks_NO.pdf");


  TH2D *h5 = (TH2D*)f_IO->Get("detected_tracks_IO");
  TCanvas *c5 = new TCanvas("c5", "c5", 600, 500); 
  h5->Draw("colz");

  SetPlotTitle(h5, "Track events [IO]");
  SetLabelSize(h5, 0.035);

  h5->GetXaxis()->SetRangeUser(3,100);
  h5->GetYaxis()->SetRangeUser(-1,0);
  c5->SetLogx();
  h5->Draw("colz");
  c5->SaveAs("./default_detres/tracks_IO.pdf");


  TH2D *h6 = (TH2D*)f_IO->Get("asymmetry_track");
  TCanvas *c6 = new TCanvas("c6", "c6", 630, 500); 
  h6->Draw("colz");

  SetPlotTitle(h6, "Asymmetry tracks");
  SetLabelSize(h6, 0.035);

  h6->GetXaxis()->SetRangeUser(3,100);
  h6->GetYaxis()->SetRangeUser(-1,0);
  h6->GetZaxis()->SetRangeUser(-1,1); 
  c6->SetLogx();
  h6->Draw("colz");
  c6->SaveAs("./default_detres/asymmetry_tracks.pdf");

}
