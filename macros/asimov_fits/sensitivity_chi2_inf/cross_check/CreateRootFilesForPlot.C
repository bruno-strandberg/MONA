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


TTree* CalculateChi2Infinite(TTree* input_ttree, TString tree_nametitle, Bool_t write_scatter=kFALSE);
TTree* GetChi2(TTree* input_ttree, TString tree_nametitle);
TTree* GetChi2QRand(TTree* input_ttree, TString tree_nametitle);

void CreateRootFilesForPlot() {

  std::vector<Int_t> pid_cats = {2,3,4,5,10};

  TFile* f_read = TFile::Open("./data_chi2_th23.root", "READ");
  TFile* f = TFile::Open("./data_chi2_th23_extraplotation.root", "RECREATE");
  for (auto pid: pid_cats) {

    TTree* tree_in = (TTree*)f_read->Get( Form("data_tree_%i_no", pid) );

    TTree* t_inf_chi2 = CalculateChi2Infinite(tree_in, Form("%ibins_no", pid));
    TTree* t_finite_chi2 = GetChi2QRand(tree_in, Form("%ibins_no", pid));

    f->cd();
    t_inf_chi2->Write();
    t_finite_chi2->Write();
  }
  for (auto pid: pid_cats) {

    TTree* tree_in = (TTree*)f_read->Get( Form("data_tree_%i_io", pid) );

    TTree* t_inf_chi2 = CalculateChi2Infinite(tree_in, Form("%ibins_io", pid));
    TTree* t_finite_chi2 = GetChi2QRand(tree_in, Form("%ibins_io", pid));

    f->cd();
    t_inf_chi2->Write();
    t_finite_chi2->Write();
  }
  f->Close();
}


TTree* CalculateChi2Infinite(TTree* input_ttree, TString tree_nametitle, Bool_t write_scatter) {
  
  Int_t th23_min = input_ttree->GetMinimum("th23");
  Int_t th23_max = input_ttree->GetMaximum("th23");

  Int_t th23;
  Double_t inf_chi2;
  Double_t inf_err;
  Double_t sqrt_inf_chi2;
  Double_t sqrt_inf_err;
  Double_t kParameter;

  TTree* t_inf_chi2 = new TTree("inf_" + tree_nametitle, "Chi2_inf values from fit " + tree_nametitle);
  TBranch* b_th23          = t_inf_chi2->Branch("th23", &th23, "th23/I");
  TBranch* b_inf_chi2      = t_inf_chi2->Branch("inf_chi2", &inf_chi2, "inf_chi2/D");
  TBranch* b_inf_err       = t_inf_chi2->Branch("inf_err", &inf_err, "inf_err/D");
  TBranch* b_sqrt_inf_chi2 = t_inf_chi2->Branch("sqrt_inf_chi2", &sqrt_inf_chi2, "sqrt_inf_chi2/D");
  TBranch* b_sqrt_inf_err  = t_inf_chi2->Branch("sqrt_inf_err", &sqrt_inf_err, "sqrt_inf_err/D");
  TBranch* b_kParameter    = t_inf_chi2->Branch("kParameter", &kParameter, "kParameter/D");

  for (Int_t i = th23_min; i <= th23_max; i++) {
    Int_t n = input_ttree->Draw("percentage:fit_chi2", Form("th23==%i", i), "goff"); 
    TGraph *g = new TGraph(n, input_ttree->GetV1(), input_ttree->GetV2()); 

    TF1 *f_tr = new TF1("f_tr", "[0] + ([1] / x)", 1e-3, 10);
    f_tr->SetParameters(1, 1); // Initial values of the params.
    TFitResultPtr result = g->Fit(f_tr, "QS"); // Q = quiet, S = result in pointer

    Int_t NFreeParams = result->NFreeParameters();
    std::vector<double> pars = result->Parameters();
    std::vector<double> errs = result->Errors();

    th23 = i;
    inf_chi2 = pars[0];
    sqrt_inf_chi2 = TMath::Sqrt(inf_chi2);
    inf_err = errs[0];
    sqrt_inf_err = inf_err / (2. * sqrt_inf_chi2); // Standard error propagation.
    kParameter = pars[1];

    gStyle->SetOptFit(kTRUE);
    if (write_scatter) g->Write(Form("scatter_chi2_%i_", i) + tree_nametitle);

    t_inf_chi2->Fill();
  }

  return t_inf_chi2;
}

TTree* GetChi2(TTree* input_ttree, TString tree_nametitle) {
  
  // CloneTree does not work as expected...
  TTree* t_fit_chi2 = input_ttree->CopyTree("percentage==1");
  t_fit_chi2->SetName("chi2_" + tree_nametitle);
  t_fit_chi2->SetTitle("Chi2 values " + tree_nametitle);

  return t_fit_chi2;
}


/* This function is necessary because for percentage==1 you expect a unique value
 * for the sensitivity. However, when one randomizes the Q value, for 3+ bins the
 * sensitivity is non-constant anymore, depending on the randomization. This is 
 * because signal events might end up in signal regions, and you have no control
 * over this. For 2 bins, they always end up in the same bin, nothing happens, and
 * the sensitivity is constant, whether you randomize or not.
 * For 3+ bins, you can move events to specific bins (accidentally) and for 100%
 * of event usage, the value of S increases/decreases.
 */
TTree* GetChi2QRand(TTree* input_ttree, TString tree_nametitle) {
  
  // Trees will be added this this list and merged.
  TList *treeList = new TList;

  TTree* tFitChi2 = input_ttree->CopyTree("percentage==1");

  Int_t th23Min = tFitChi2->GetMinimum("th23");
  Int_t th23Max = tFitChi2->GetMaximum("th23");

  for (Int_t th23 = th23Min; th23 <= th23Max; th23++) {
    TTree* tFitChi2Th23Cut = tFitChi2->CopyTree(Form("th23==%i", th23));
    Double_t sMin = tFitChi2Th23Cut->GetMinimum("fit_chi2");

    TTree* tFitChi2Min = tFitChi2Th23Cut->CopyTree(Form("fit_chi2==%f", sMin));

    treeList->Add(tFitChi2Min);
  }

  TTree* mergedTree = TTree::MergeTrees(treeList);
  mergedTree->SetName("chi2_" + tree_nametitle);
  mergedTree->SetTitle("Chi2 values " + tree_nametitle);

  return mergedTree;
}
