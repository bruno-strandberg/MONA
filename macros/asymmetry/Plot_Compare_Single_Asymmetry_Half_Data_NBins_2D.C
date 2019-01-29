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

void Plot_Compare_Single_Asymmetry_Half_Data_NBins_2D() {

  const int N_PID_CLASSES = 10;
  Double_t PID_step = 1 / float(N_PID_CLASSES);

  TString filefolder = TString::Format("./pid_detres/pid_binning_%i/", N_PID_CLASSES);

  TCanvas* ct = new TCanvas("ctrack", "ctrack", 2000, 800);
  ct->Divide(5,2);
  TCanvas* cs = new TCanvas("cshower", "cshower", 2000, 800);
  cs->Divide(5,2);

  for (int i = 0; i < N_PID_CLASSES; i++) {
    TString s_alldata  = filefolder + TString::Format("use_half_data_all/full_error/asymmetry_split_%.2f.root", PID_step * i);
    TString s_data_gt  = filefolder + TString::Format("use_half_data_gt300/full_error/asymmetry_split_%.2f.root", PID_step * i);
    TString s_data_lt  = filefolder + TString::Format("use_half_data_lt300/full_error/asymmetry_split_%.2f.root", PID_step * i);
    
    TFile *f_all = TFile::Open(s_alldata, "READ");
    TFile *f_gt  = TFile::Open(s_data_gt, "READ");
    TFile *f_lt  = TFile::Open(s_data_lt, "READ");

    gStyle->SetPalette(kBird);
//  gStyle->SetOptStat(0);

    TH2D *h_t_all = (TH2D*)f_all->Get( TString::Format("asymmetry_track_%.2f", PID_step * i) );
    TH2D *h_t_gt  = (TH2D*)f_gt ->Get( TString::Format("asymmetry_track_%.2f", PID_step * i) );
    TH2D *h_t_lt  = (TH2D*)f_lt ->Get( TString::Format("asymmetry_track_%.2f", PID_step * i) );
    TH2D *h_s_all = (TH2D*)f_all->Get( TString::Format("asymmetry_shower_%.2f", PID_step * i) );
    TH2D *h_s_gt  = (TH2D*)f_gt ->Get( TString::Format("asymmetry_shower_%.2f", PID_step * i) );
    TH2D *h_s_lt  = (TH2D*)f_lt ->Get( TString::Format("asymmetry_shower_%.2f", PID_step * i) );

    TH2D *h_t_nev_all = (TH2D*)f_all->Get( "detected_tracks_NO" );
    TH2D *h_t_nev_gt  = (TH2D*)f_gt ->Get( "detected_tracks_NO" );
    TH2D *h_t_nev_lt  = (TH2D*)f_lt ->Get( "detected_tracks_NO" );
    TH2D *h_s_nev_all = (TH2D*)f_all->Get( "detected_showers_NO" );
    TH2D *h_s_nev_gt  = (TH2D*)f_gt ->Get( "detected_showers_NO" );
    TH2D *h_s_nev_lt  = (TH2D*)f_lt ->Get( "detected_showers_NO" );

    TH2D *h_t_gt_out = std::get<0>(ErrorPlotWithHalfData2D(h_t_all, h_t_gt, h_t_lt, h_t_nev_all, h_t_nev_gt, h_t_nev_lt));
    h_t_gt_out->SetName(TString::Format("h_t_gt_%i", i));
    TH2D *h_t_lt_out = std::get<1>(ErrorPlotWithHalfData2D(h_t_all, h_t_gt, h_t_lt, h_t_nev_all, h_t_nev_gt, h_t_nev_lt));
    h_t_lt_out->SetName(TString::Format("h_t_lt_%i", i));

    TH2D *h_s_gt_out = std::get<0>(ErrorPlotWithHalfData2D(h_s_all, h_s_gt, h_s_lt, h_s_nev_all, h_s_nev_gt, h_s_nev_lt));
    h_s_gt_out->SetName(TString::Format("h_s_gt_%i", i));
    TH2D *h_s_lt_out = std::get<1>(ErrorPlotWithHalfData2D(h_s_all, h_s_gt, h_s_lt, h_s_nev_all, h_s_nev_gt, h_s_nev_lt));
    h_s_lt_out->SetName(TString::Format("h_s_lt_%i", i));

    ct->cd(i+1);
    h_t_gt_out->Draw("colz");
    h_t_gt_out->GetXaxis()->SetTitle("Relative error");
    h_t_gt_out->GetYaxis()->SetTitle("Number of events in NO");
    h_t_gt_out->SetTitle(TString::Format("#splitline{Relative error on track asymmetry Q=0.%i}{compared to using half data}", i));
//    Double_t max_t = 1.1 * std::max(h_t_gt_out->GetMaximum(), h_t_lt_out->GetMaximum());
//    h_t_gt_out->GetYaxis()->SetRangeUser(0, max_t);
//    h_t_lt_out->SetLineStyle(2);
//    h_t_lt_out->Draw("same");

//    TLegend* leg_t = new TLegend(0.1, 0.75, 0.4, 0.85);
//    leg_t->AddEntry(h_t_gt_out, "Data runID > 300", "l");
//    leg_t->AddEntry(h_t_lt_out, "Data runID #leq 300", "l");
//    leg_t->Draw();

    cs->cd(i+1);
    h_s_gt_out->Draw("colz");
    h_s_gt_out->GetXaxis()->SetTitle("Relative error");
    h_s_gt_out->GetYaxis()->SetTitle("Number of events in NO");
    h_s_gt_out->SetTitle(TString::Format("#splitline{Relative error on shower asymmetry Q=0.%i}{compared to using half data}", i));
//    Double_t max_s = 1.1 * std::max(h_s_gt_out->GetMaximum(), h_s_lt_out->GetMaximum());
//    h_s_gt_out->GetYaxis()->SetRangeUser(0, max_s);
//    h_s_lt_out->SetLineStyle(2);
//    h_s_lt_out->Draw("same");

//    TLegend* leg_s = new TLegend(0.1, 0.75, 0.4, 0.85);
//    leg_s->AddEntry(h_s_gt_out, "Data runID > 300", "l");
//    leg_s->AddEntry(h_s_lt_out, "Data runID #leq 300", "l");
//    leg_s->Draw();
//    c1->SaveAs(filefolder + TString::Format("use_half_data_all/full_error/plot_error_half_data_%i.pdf", i));
 
  }
}
