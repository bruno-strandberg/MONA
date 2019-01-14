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

void Print_ResponseFunction_Events() {

  const int N_PID_CLASSES = 10;
  Double_t PID_step = 1 / float(N_PID_CLASSES);

  TString input_def_t = "./default_detres/track_response.root";
  TString input_def_s = "./default_detres/shower_response.root";
  
  TFile *f_def_t = TFile::Open(input_def_t, "READ");
  TFile *f_def_s = TFile::Open(input_def_s, "READ");

  TH3D *h_t_def = (TH3D*)f_def_t->Get("hresp");
  TH3D *h_s_def = (TH3D*)f_def_s->Get("hresp");

  Double_t n_t_def = h_t_def->Integral();
  Double_t n_s_def = h_s_def->Integral();
  
  cout.precision(8);
  cout << "N track" << " | " << "N shower" << " | " << endl;
  cout << n_t_def << " | " << n_s_def << " | " << endl;
  

  for (int i = 0; i < N_PID_CLASSES; i++){
    TString input_t = Form("./pid_detres/pid_binning_%i/track_response_%.2f.root",  N_PID_CLASSES, i * PID_step);
    TString input_s = Form("./pid_detres/pid_binning_%i/shower_response_%.2f.root", N_PID_CLASSES, i * PID_step);
    
    TFile *f_t = TFile::Open(input_t, "READ");
    TFile *f_s = TFile::Open(input_s, "READ");

    TH3D *h_t = (TH3D*)f_t->Get("hresp");
    TH3D *h_s = (TH3D*)f_s->Get("hresp");

    Double_t n_t = h_t->Integral();
    Double_t n_s = h_s->Integral();
    
    if (i == 0) {
        cout << cout.width(25) << "N track" << " | " << "N shower" << " | " << endl;
    }
    cout.precision(7);
    cout << "Values for " << Form("%.2f", i * PID_step) << " < q < " << Form("%.2f", (i+1) * PID_step) << " | " 
         << n_t << " | " << n_s << " | " << endl;
  }
}
