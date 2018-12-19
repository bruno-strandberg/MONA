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


// THIS SCRIPT IS DEPRECATED, I AM NOT SURE WHERE I GOT THE NUMBERS FROM FOR EXPECTED ENTRIES. DO NOT USE THIS UNTIL FIXED.

void histo_detected_total_entries() {

  const int N_PID_CLASSES = 15;
  Double_t PID_step = 1 / float(N_PID_CLASSES);

  TString input_def = "./default_detres/asym_pid.root";
  
  TFile *f_def = TFile::Open(input_def, "READ");

  TH2D *h_t_NO_def = (TH2D*)f_def->Get("detected_tracks_NO");
  TH2D *h_s_NO_def = (TH2D*)f_def->Get("detected_showers_NO");
  TH2D *h_t_IO_def = (TH2D*)f_def->Get("detected_tracks_IO");
  TH2D *h_s_IO_def = (TH2D*)f_def->Get("detected_showers_IO");

  Double_t n_t_NO_def = h_t_NO_def->Integral();
  Double_t n_s_NO_def = h_s_NO_def->Integral();
  Double_t n_t_IO_def = h_t_IO_def->Integral();
  Double_t n_s_IO_def = h_s_IO_def->Integral();
  
  cout.precision(8);
  cout << "N track NO" << " | " << "N track IO" << " | " << "N shower NO" << " | " << "N shower IO" << " | " << endl;
  cout << n_t_NO_def << " | " << n_t_IO_def << " | " << n_s_NO_def << " | " << n_s_IO_def << " | " << endl;
  

  for (int i = 0; i < 10; i++){
    TString input = Form("./pid_detres/asym_pid_%.2f.root", i * PID_step);
    
    TFile *f = TFile::Open(input, "READ");

    TH2D *h_t_NO = (TH2D*)f->Get("detected_tracks_NO");
    TH2D *h_t_IO = (TH2D*)f->Get("detected_tracks_IO");
    TH2D *h_s_NO = (TH2D*)f->Get("detected_showers_NO");
    TH2D *h_s_IO = (TH2D*)f->Get("detected_showers_IO");

    Double_t n_t_NO = h_t_NO->Integral();
    Double_t n_t_IO = h_t_IO->Integral();
    Double_t n_s_NO = h_s_NO->Integral();
    Double_t n_s_IO = h_s_IO->Integral();
    
    if (i == 0) {
        cout << cout.width(25) << "N track NO" << " | " << "N track IO" << " | " << "N shower NO" << " | " << "N shower IO" << " | " << endl;
    }
    cout.precision(7);
    cout << "Values for " << Form("%.2f", i * 0.1) << " < q < " << Form("%.2f", (i+1) * 0.1) << " | " 
         << n_t_NO << " | " << n_t_IO << " | " << n_s_NO << " | " << n_s_IO << " | " << endl;
         //<< cout.width(12) << n_t_NO << " | " << cout.width(12) << n_t_IO << " | " << cout.width(12) << n_s_NO << " | " << cout.width(12) << n_s_IO << " | " << endl;

  }
}
