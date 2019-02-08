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


void PlotRelativeErrorNEventsDefault() {

  TString filefolder = "./default_detres/";
  TString s_file_NO = filefolder + "default_expected_evts_NO.root";
  TString s_file_IO = filefolder + "default_expected_evts_IO.root";
  
  TFile *f_NO = TFile::Open(s_file_NO, "READ");
  TFile *f_IO = TFile::Open(s_file_IO, "READ");

  gStyle->SetPalette(kLightTemperature);
  gStyle->SetOptStat(0);

  TH2D *h_t_err_no = (TH2D*)f_NO->Get("track_response_expct_err_yx");
  TH2D *h_s_err_no = (TH2D*)f_NO->Get("shower_response_expct_err_yx");
  TH2D *h_t_err_io = (TH2D*)f_IO->Get("track_response_expct_err_yx");
  TH2D *h_s_err_io = (TH2D*)f_IO->Get("shower_response_expct_err_yx");

  TCanvas *c1 = new TCanvas("c1", "c1", 1200, 1000); 
  c1->Divide(2, 2);

  c1->cd(1);
  h_t_err_no->Draw("colz");
  SetPlotTitle(h_t_err_no, "Track events relative error [NO]");
  SetLabelSizes(h_t_err_no, 0.035);
  h_t_err_no->GetXaxis()->SetRangeUser(1,100);
  h_t_err_no->GetYaxis()->SetRangeUser(-1,0);
  h_t_err_no->GetZaxis()->SetRangeUser(1e-2 ,1);
  c1->cd(1)->SetRightMargin(0.15);
  c1->cd(1)->SetLogx();
  c1->cd(1)->SetLogz();

  c1->cd(2);
  h_s_err_no->Draw("colz");
  SetPlotTitle(h_s_err_no, "Shower events relative error [NO]");
  SetLabelSizes(h_s_err_no, 0.035);
  h_s_err_no->GetXaxis()->SetRangeUser(1,100);
  h_s_err_no->GetYaxis()->SetRangeUser(-1,0);
  h_s_err_no->GetZaxis()->SetRangeUser(1e-2, 1);
  c1->cd(2)->SetRightMargin(0.15);
  c1->cd(2)->SetLogx();
  c1->cd(2)->SetLogz();

  c1->cd(3);
  h_t_err_io->Draw("colz");
  SetPlotTitle(h_t_err_io, "Track events relative error [IO]");
  SetLabelSizes(h_t_err_io, 0.035);
  h_t_err_io->GetXaxis()->SetRangeUser(1,100);
  h_t_err_io->GetYaxis()->SetRangeUser(-1,0);
  h_t_err_io->GetZaxis()->SetRangeUser(1e-2, 1);
  c1->cd(3)->SetRightMargin(0.15);
  c1->cd(3)->SetLogx();
  c1->cd(3)->SetLogz();

  c1->cd(4);
  h_s_err_io->Draw("colz");
  SetPlotTitle(h_s_err_io, "Shower events relative error [IO]");
  SetLabelSizes(h_s_err_io, 0.035);
  h_s_err_io->GetXaxis()->SetRangeUser(1,100);
  h_s_err_io->GetYaxis()->SetRangeUser(-1,0);
  h_s_err_io->GetZaxis()->SetRangeUser(1e-2, 1);
  c1->cd(4)->SetRightMargin(0.15);
  c1->cd(4)->SetLogx();
  c1->cd(4)->SetLogz();

  c1->SaveAs(filefolder + "relative_errors.pdf");
}
