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
#include "NMHUtils.h"
#include "SummaryParser.h"
#include "SummaryEvent.h"

#include <iostream>
using namespace std;

void Print_Detected_Events_NBins() {

  const int N_PID_CLASSES = 10;
  Double_t PID_step = 1 / float(N_PID_CLASSES);

  TString input_def = "./default_detres/asymmetry_default.root";
  
  TFile *f_def = TFile::Open(input_def, "READ");

  TH3D *h_t_NO_def = (TH3D*)f_def->Get("detected_tracks_NO");
  TH3D *h_s_NO_def = (TH3D*)f_def->Get("detected_showers_NO");
  TH3D *h_t_IO_def = (TH3D*)f_def->Get("detected_tracks_IO");
  TH3D *h_s_IO_def = (TH3D*)f_def->Get("detected_showers_IO");

  Double_t n_t_NO_def = h_t_NO_def->Integral();
  Double_t n_s_NO_def = h_s_NO_def->Integral();
  Double_t n_t_IO_def = h_t_IO_def->Integral();
  Double_t n_s_IO_def = h_s_IO_def->Integral();
  
  cout.precision(8);
  cout << "N track NO" << " | " << "N track IO" << " | " << "N shower NO" << " | " << "N shower IO" << " | " << endl;
  cout << n_t_NO_def << " | " << n_t_IO_def << " | " << n_s_NO_def << " | " << n_s_IO_def << " | " << endl;
  

  for (int i = 0; i < N_PID_CLASSES; i++){
    TString input = Form("./pid_detres/pid_binning_%i/asymmetry_split_%.2f.root", N_PID_CLASSES, PID_step * i);
    
    TFile *f = TFile::Open(input, "READ");

    TH3D *h_t_NO = (TH3D*)f->Get("detected_tracks_NO");
    TH3D *h_t_IO = (TH3D*)f->Get("detected_tracks_IO");
    TH3D *h_s_NO = (TH3D*)f->Get("detected_showers_NO");
    TH3D *h_s_IO = (TH3D*)f->Get("detected_showers_IO");

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
  }
}
