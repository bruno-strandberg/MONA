#include "TSystem.h"
#include "TROOT.h"
#include "TDatime.h"
#include "TH2.h"
#include "TH3.h"
#include "TFile.h"
#include "TLine.h"
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

  /* Script to plot the relative error on events for NO, IO for tracks and showers.
   */


void PlotRelativeErrorNEvents3Bins() {

  const int N_PID_CLASSES = 3;
  const Double_t PID_STEP = 1 / Double_t(N_PID_CLASSES);
  const Double_t PID_CUT = 0.6;
  const Double_t PID_EDGE = PID_CUT * N_PID_CLASSES;

  std::map<Int_t, Double_t> pid_map;
  pid_map.insert(std::make_pair(0, 0.0)); // shower
  pid_map.insert(std::make_pair(1, 0.2)); // middle group: shower
  pid_map.insert(std::make_pair(2, 0.6)); // track
  pid_map.insert(std::make_pair(3, 1.0)); // upper limit 

  TString filefolder = Form("./pid_detres/pid_binning_%i/qcut_0.2/", N_PID_CLASSES);

  std::vector<TH2D*> h_t_err_no;
  std::vector<TH2D*> h_s_err_no;
  std::vector<TH2D*> h_t_err_io;
  std::vector<TH2D*> h_s_err_io;

  for (Int_t i = 0; i < N_PID_CLASSES; i++) {
    TString s_file_NO = filefolder + Form("split_expected_evts_NO_%.2f.root", pid_map[i]);
  
    TFile *f_NO = TFile::Open(s_file_NO, "READ");

    h_t_err_no.push_back( (TH2D*)f_NO->Get( Form("track_response_%.2f_expct_err_yx", pid_map[i])  ));
    h_s_err_no.push_back( (TH2D*)f_NO->Get( Form("shower_response_%.2f_expct_err_yx", pid_map[i]) ));
  }

  gStyle->SetPalette(kLightTemperature);
  gStyle->SetOptStat(0);

  TCanvas *c1 = new TCanvas("c1", "c1", 1800,  600); 
  c1->Divide(3, 1);

  std::vector<std::pair<Double_t, Double_t>> lines;
  lines.push_back(std::make_pair(3  , 90));
  lines.push_back(std::make_pair(1.5, 50));
  lines.push_back(std::make_pair(3, 90));


  for(Int_t i = 0; i < N_PID_CLASSES; i++) {
    c1->cd(i+1);
    if (pid_map[i] < PID_CUT) {
      h_s_err_no[i]->Draw("colz");
      SetPlotTitle(h_s_err_no[i], "Shower events relative error [NO]");
      SetLabelSizes(h_s_err_no[i], 0.035);
      h_s_err_no[i]->GetXaxis()->SetRangeUser(1,100);
      h_s_err_no[i]->GetYaxis()->SetRangeUser(-1,0);
      h_s_err_no[i]->GetZaxis()->SetRangeUser(1e-2, 1);
      TLine *line1 = new TLine(lines[i].first, -1, lines[i].first, 0);
      line1->SetLineColor(kRed);
      line1->Draw();
      TLine *line2 = new TLine(lines[i].second, -1, lines[i].second, 0);
      line2->SetLineColor(kRed);
      line2->Draw();
    }
    else {
      h_t_err_no[i]->Draw("colz");
      SetPlotTitle(h_t_err_no[i], "Track events relative error [NO]");
      SetLabelSizes(h_t_err_no[i], 0.035);
      h_t_err_no[i]->GetXaxis()->SetRangeUser(1,100);
      h_t_err_no[i]->GetYaxis()->SetRangeUser(-1,0);
      h_t_err_no[i]->GetZaxis()->SetRangeUser(1e-2, 1);
      TLine *line1 = new TLine(lines[i].first, -1, lines[i].first, 0);
      line1->SetLineColor(kRed);
      line1->Draw();
      TLine *line2 = new TLine(lines[i].second, -1, lines[i].second, 0);
      line2->SetLineColor(kRed);
      line2->Draw();
    }
    c1->cd(i+1)->SetRightMargin(0.15);
    c1->cd(i+1)->SetLogx();
    c1->cd(i+1)->SetLogz();

  }

  c1->SaveAs(filefolder + "relative_errors_no.pdf");

}
