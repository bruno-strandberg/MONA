#include <iomanip>

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

void histo_MC_total_entries() {

  const int N_PID_CLASSES = 15;
  Double_t PID_step = 1 / float(N_PID_CLASSES);

  TString input_def_t_NO = "./default_detres/track_response_timing_NO.root";
  TString input_def_s_NO = "./default_detres/shower_response_timing_NO.root";
  TString input_def_t_IO = "./default_detres/track_response_timing_IO.root";
  TString input_def_s_IO = "./default_detres/shower_response_timing_IO.root";
  
  TFile *f_def_t_NO = TFile::Open(input_def_t_NO, "READ");
  TFile *f_def_s_NO = TFile::Open(input_def_s_NO, "READ");
  TFile *f_def_t_IO = TFile::Open(input_def_t_IO, "READ");
  TFile *f_def_s_IO = TFile::Open(input_def_s_IO, "READ");

  TH3D *h_t_NO_def = (TH3D*)f_def_t_NO->Get("hresp");
  TH3D *h_s_NO_def = (TH3D*)f_def_s_NO->Get("hresp");
  TH3D *h_t_IO_def = (TH3D*)f_def_t_IO->Get("hresp");
  TH3D *h_s_IO_def = (TH3D*)f_def_s_IO->Get("hresp");

  Double_t n_t_NO_def = h_t_NO_def->Integral();
  Double_t n_s_NO_def = h_s_NO_def->Integral();
  Double_t n_t_IO_def = h_t_IO_def->Integral();
  Double_t n_s_IO_def = h_s_IO_def->Integral();
  
  cout.precision(8);
  cout << "N track NO" << " | " << "N track IO" << " | " << "N shower NO" << " | " << "N shower IO" << " | " << endl;
  cout << n_t_NO_def << " | " << n_t_IO_def << " | " << n_s_NO_def << " | " << n_s_IO_def << " | " << endl;
  

  for (int i = 0; i < N_PID_CLASSES; i++){
    TString input_t_NO = Form("./pid_detres/track_response_NO_%.2f.root", i * N_PID_CLASSES);
    TString input_t_IO = Form("./pid_detres/track_response_IO_%.2f.root", i * N_PID_CLASSES);
    TString input_s_NO = Form("./pid_detres/shower_response_NO_%.2f.root", i * N_PID_CLASSES);
    TString input_s_IO = Form("./pid_detres/shower_response_IO_%.2f.root", i * N_PID_CLASSES);
    
    TFile *f_t_NO = TFile::Open(input_t_NO, "READ");
    TFile *f_t_IO = TFile::Open(input_t_IO, "READ");
    TFile *f_s_NO = TFile::Open(input_s_NO, "READ");
    TFile *f_s_IO = TFile::Open(input_s_IO, "READ");

    TH3D *h_t_NO = (TH3D*)f_t_NO->Get("hresp");
    TH3D *h_t_IO = (TH3D*)f_t_IO->Get("hresp");
    TH3D *h_s_NO = (TH3D*)f_s_NO->Get("hresp");
    TH3D *h_s_IO = (TH3D*)f_s_IO->Get("hresp");

    Double_t n_t_NO = h_t_NO->Integral();
    Double_t n_t_IO = h_t_IO->Integral();
    Double_t n_s_NO = h_s_NO->Integral();
    Double_t n_s_IO = h_s_IO->Integral();
    
    if (i == 0) {
        cout << cout.width(25) << "N track NO" << " | " << "N track IO" << " | " << "N shower NO" << " | " << "N shower IO" << " | " << endl;
    }
    cout.precision(7);
    cout << "Values for " << Form("%.2f", i * PID_step) << " < q < " << Form("%.2f", (i+1) * PID_step) << " | " 
         << n_t_NO << " | " << n_t_IO << " | " << n_s_NO << " | " << n_s_IO << " | " << endl;
         //<< cout.width(12) << n_t_NO << " | " << cout.width(12) << n_t_IO << " | " << cout.width(12) << n_s_NO << " | " << cout.width(12) << n_s_IO << " | " << endl;

  }
}
