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

 /* Script to print the relative errors of an asymmetry histogram for one file and one histogram
  */

void Plot_Single_Asymmetry_Errors(TString filepath="./pid_detres/pid_binning_10/asymmetry_split_0.90.root",
                           TString name="asymmetry_track_0.90") {

    TFile *f = TFile::Open(filepath, "READ");
    TH2D *h_asym = (TH2D*)f->Get(name);
    h_asym->SetDirectory(0);

    TH2D *h_asym_err = (TH2D*)h_asym->Clone("h_a_rel_err");
    h_asym_err->SetNameTitle("h_a_rel_err", "h_a_rel_err");
    h_asym_err->Reset();
    h_asym_err->SetDirectory(0);

    for (Int_t xb = 1; xb <= h_asym->GetXaxis()->GetNbins(); xb++) {
      for (Int_t yb = 1; yb <= h_asym->GetYaxis()->GetNbins(); yb++) {
	
        Double_t N_h_asym     = h_asym->GetBinContent(xb, yb);
        Double_t N_h_asym_err = h_asym->GetBinError(xb, yb);
        Double_t relative_err = 0;
        if (N_h_asym != 0) {
          relative_err = N_h_asym_err / N_h_asym;
        } else {
          relative_err = 0;
        }
        h_asym_err->SetBinContent(xb, yb, std::abs(rel_err));
      }
    }
    gStyle->SetOptStat(0);
    TCanvas *c1 = new TCanvas("c1", "c1", 800, 600); 
    gStyle->SetPalette(kBird);
    h_asym_err->Draw("colz");
    h_asym_err->GetZaxis()->UnZoom();
    c1->SetLogx();
    c1->SetLogz();
}
