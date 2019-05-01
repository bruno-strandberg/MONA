#include "TCanvas.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TTree.h"
#include "TFile.h"
#include "TRandom.h"
#include "TStyle.h"


TTree* CalculateChi2Infinite(TTree* input_ttree, TString tree_nametitle);

void CreateRootFilesForPlot() {

  std::vector<Int_t> pid_cats = {2,3,4,5,10};

  TFile* f = TFile::Open("./data_chi2_extraplotation.root", "RECREATE");
  TFile* f_read = TFile::Open("./data_chi2_th23.root", "READ");
  for (auto pid: pid_cats) {

    TTree* tree_in = (TTree*)f_read->Get( Form("data_tree_%i_no", pid) );

    TTree* t_fit_chi2 = CalculateChi2Infinite(tree_in, Form("%ibins_no", pid));
    delete tree_in;

    f->cd();
    t_fit_chi2->Write();
  }
  for (auto pid: pid_cats) {

    TTree* tree_in = (TTree*)f_read->Get( Form("data_tree_%i_io", pid) );

    TTree* t_fit_chi2 = CalculateChi2Infinite(tree_in, Form("%ibins_io", pid));
    delete tree_in;

    f->cd();
    t_fit_chi2->Write();
  }
  f->Write();
  f->Close();
}


TTree* CalculateChi2Infinite(TTree* input_ttree, TString tree_nametitle) {
  
  Int_t th23_min = input_ttree->GetMinimum("th23");
  Int_t th23_max = input_ttree->GetMaximum("th23");

  Int_t th23;
  Double_t fit_chi2;

  TTree* t_fit_chi2 = new TTree("fit_tree_" + tree_nametitle, "Chi2_inf values from fit " + tree_nametitle);
  TBranch* b_th23     = t_fit_chi2->Branch("th23", &th23, "th23/I");
  TBranch* b_fit_chi2 = t_fit_chi2->Branch("fit_chi2", &fit_chi2, "fit_chi2/D");

  for (Int_t i = th23_min; i <= th23_max; i++) {
    Int_t n = input_ttree->Draw("percentage:fit_chi2", Form("th23==%i", i), "goff"); 
    TGraph *g = new TGraph(n, input_ttree->GetV1(), input_ttree->GetV2()); 

    TF1 *f_tr = new TF1("f_tr", "[0] + ([1] / x)", 1e-3, 10);
    f_tr->SetParameters(1, 1); // Initial values of the params.
    TFitResultPtr result = g->Fit(f_tr, "S"); // Q = quiet, S = result in pointer
    fit_chi2 = result->Value(0);
    th23 = i;

    gStyle->SetOptFit(kTRUE);
    g->Write(Form("scatter_chi2_%i_", i) + tree_nametitle);

    t_fit_chi2->Fill();
  }

  return t_fit_chi2;
}
