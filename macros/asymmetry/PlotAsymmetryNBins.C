#include "TSystem.h"
#include "TROOT.h"
#include "TDatime.h"
#include "TH2.h"
#include "TH3.h"
#include "TFile.h"
#include "TStopwatch.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TString.h"

#include "DetResponse.h"
#include "HelperFunctions.C"
#include "NMHUtils.h"
#include "SummaryParser.h"
#include "SummaryEvent.h"

#include <iostream>
using namespace std;

void Plot_Asymmetry_NBins() {

  const int N_PID_CLASSES = 10;
  Double_t PID_step = 1 / float(N_PID_CLASSES);
  for (int i = 0; i < N_PID_CLASSES; i++){
    TString input = TString::Format("./pid_detres/pid_binning_%i/asymmetry_split_%.2f.root", N_PID_CLASSES, i * PID_step);

    TFile *f = TFile::Open(input, "READ");

    gStyle->SetPalette(kBird);
    gStyle->SetOptStat(0);

    TH2D *h1 = (TH2D*)f->Get(TString::Format("asymmetry_track_%.2f", i * PID_step));
    TCanvas *c1 = new TCanvas(TString::Format("c%i", i), TString::Format("c%i", i), 630, 500); 
    h1->Draw("colz");

    SetPlotTitle(h1, TString::Format("Asymmetry tracks q_%.2f", i * PID_step));
    SetLabelSize(h1, 0.035);

    h1->GetXaxis()->SetRangeUser(3,100);
    h1->GetYaxis()->SetRangeUser(-1,0);
    h1->GetZaxis()->SetRangeUser(-1,1); 
    h1->GetZaxis()->SetRangeUser(-1,1);
    c1->SetLogx();
    h1->Draw("colz");
    c1->SaveAs(TString::Format("./pid_detres/asymmetry_tracks_%.2f.pdf", i * PID_step));
    c1->SaveAs(TString::Format("./pid_detres/png/asymmetry_tracks_%.2f.png", i * PID_step));

    TH2D *h2 = (TH2D*)f->Get(TString::Format("asymmetry_shower_%.2f", i * PID_step));
    TCanvas *c2 = new TCanvas(TString::Format("c%i", i+10), TString::Format("c%i", i+10), 630, 500); 
    h2->Draw("colz");

    SetPlotTitle(h2, TString::Format("Asymmetry showers q_%.2f", i * PID_step));
    SetLabelSize(h2, 0.035);

    h2->GetXaxis()->SetRangeUser(3,100);
    h2->GetYaxis()->SetRangeUser(-1,0);
    h2->GetZaxis()->SetRangeUser(-1,1); 
    h2->GetZaxis()->SetRangeUser(-1,1);
    c2->SetLogx();
    h2->Draw("colz");
    c2->SaveAs(TString::Format("./pid_detres/asymmetry_showers_%.2f.pdf", i * PID_step));
    c2->SaveAs(TString::Format("./pid_detres/png/asymmetry_showers_%.2f.png", i * PID_step));
  }
}
