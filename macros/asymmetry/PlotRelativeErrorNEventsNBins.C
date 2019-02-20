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

  /* Script to plot the relative error on events for NO, IO for tracks and showers.
   */


void PlotRelativeErrorNEventsNBins() {

  const int N_PID_CLASSES = 10;
  const Double_t PID_STEP = 1 / float(N_PID_CLASSES);
  const Double_t PID_EDGE = 0.6 * N_PID_CLASSES;
  TString filefolder = Form("./pid_detres/pid_binning_%i/track_hybrid_energy/", N_PID_CLASSES);

  std::vector<TH2D*> h_t_err_no;
  std::vector<TH2D*> h_s_err_no;
  std::vector<TH2D*> h_t_err_io;
  std::vector<TH2D*> h_s_err_io;

  for (Int_t i = 0; i < N_PID_CLASSES; i++) {
    TString s_file_NO = filefolder + Form("split_expected_evts_NO_%.2f.root", i * PID_STEP);
    TString s_file_IO = filefolder + Form("split_expected_evts_IO_%.2f.root", i * PID_STEP);
  
    TFile *f_NO = TFile::Open(s_file_NO, "READ");
    TFile *f_IO = TFile::Open(s_file_IO, "READ");

    h_t_err_no.push_back( (TH2D*)f_NO->Get( Form("track_response_%.2f_expct_err_yx", i*PID_STEP)  ));
    h_s_err_no.push_back( (TH2D*)f_NO->Get( Form("shower_response_%.2f_expct_err_yx", i*PID_STEP) ));
    h_t_err_io.push_back( (TH2D*)f_IO->Get( Form("track_response_%.2f_expct_err_yx", i*PID_STEP)  ));
    h_s_err_io.push_back( (TH2D*)f_IO->Get( Form("shower_response_%.2f_expct_err_yx", i*PID_STEP) ));
  }

  gStyle->SetPalette(kLightTemperature);
  gStyle->SetOptStat(0);

  TCanvas *c1 = new TCanvas("c1", "c1", 2000,  800); 
  c1->Divide(5, 2);
  TCanvas *c2 = new TCanvas("c2", "c2", 2000,  800); 
  c2->Divide(5, 2);

  for(Int_t i = 0; i < N_PID_CLASSES; i++) {
    c1->cd(i+1);
    if (i < PID_EDGE) {
      h_s_err_no[i]->Draw("colz");
      SetPlotTitle(h_s_err_no[i], "Shower events relative error [NO]");
      SetLabelSizes(h_s_err_no[i], 0.035);
      h_s_err_no[i]->GetXaxis()->SetRangeUser(1,100);
      h_s_err_no[i]->GetYaxis()->SetRangeUser(-1,0);
      h_s_err_no[i]->GetZaxis()->SetRangeUser(1e-2, 1);
    }
    else {
      h_t_err_no[i]->Draw("colz");
      SetPlotTitle(h_t_err_no[i], "Track events relative error [NO]");
      SetLabelSizes(h_t_err_no[i], 0.035);
      h_t_err_no[i]->GetXaxis()->SetRangeUser(1,100);
      h_t_err_no[i]->GetYaxis()->SetRangeUser(-1,0);
      h_t_err_no[i]->GetZaxis()->SetRangeUser(1e-2, 1);
    }
    c1->cd(i+1)->SetRightMargin(0.15);
    c1->cd(i+1)->SetLogx();
    c1->cd(i+1)->SetLogz();


    c2->cd(i+1);
    if (i < PID_EDGE) {
      h_s_err_io[i]->Draw("colz");
      SetPlotTitle(h_s_err_io[i], "Shower events relative error [IO]");
      SetLabelSizes(h_s_err_io[i], 0.035);
      h_s_err_io[i]->GetXaxis()->SetRangeUser(1,100);
      h_s_err_io[i]->GetYaxis()->SetRangeUser(-1,0);
      h_s_err_io[i]->GetZaxis()->SetRangeUser(1e-2, 1);
    }
    else {
      h_t_err_io[i]->Draw("colz");
      SetPlotTitle(h_t_err_io[i], "Track events relative error [IO]");
      SetLabelSizes(h_t_err_io[i], 0.035);
      h_t_err_io[i]->GetXaxis()->SetRangeUser(1,100);
      h_t_err_io[i]->GetYaxis()->SetRangeUser(-1,0);
      h_t_err_io[i]->GetZaxis()->SetRangeUser(1e-2, 1);
    }
    c2->cd(i+1)->SetRightMargin(0.15);
    c2->cd(i+1)->SetLogx();
    c2->cd(i+1)->SetLogz();
  }

  c1->SaveAs(filefolder + "relative_errors_no.pdf");
  c2->SaveAs(filefolder + "relative_errors_io.pdf");

}
