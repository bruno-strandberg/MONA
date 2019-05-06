#include "TCanvas.h"
#include "TF1.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TTree.h"
#include "TFile.h"
#include "TRandom.h"
#include "TStyle.h"


void FitOverEstimationNBins(Int_t nbins=2) {

  TFile file("data_chi2.root", "READ");

  TTree* tree_in = (TTree*)file.Get( Form("data_tree_%i_no", nbins ) );

  TCanvas* c1 = new TCanvas("c1","c1");
  Int_t n = tree_in->Draw("percentage:fit_chi2", "", "goff"); 
  TGraph *g = new TGraph(n, tree_in->GetV1(), tree_in->GetV2()); 
  c1->SetLogx();
  g->SetMinimum(0);
 // g->SetMaximum(120);
  g->Draw("ap"); 

  TF1 *f_tr = new TF1("f_tr", "[0] + ([1] / x)", 1e-3, 10);
  f_tr->SetParameters(1, 1); // Initial values of the params.
  g->Fit(f_tr);
  g->SetTitle("Fit: p0 + (p1 / x)");
  g->GetXaxis()->SetTitle("Fraction of total #MC events");
  g->GetYaxis()->SetTitle("#Delta #chi^{2}");
  f_tr = g->GetFunction("f_tr");
  g->SetMarkerStyle(6);
  gStyle->SetOptFit(kTRUE);

}
