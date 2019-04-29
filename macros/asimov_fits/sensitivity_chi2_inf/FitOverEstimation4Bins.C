#include "TCanvas.h"
#include "TF1.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TTree.h"
#include "TFile.h"
#include "TRandom.h"


void FitOverEstimation4Bins() {

  TString MONADIR = (TString)getenv("MONADIR") + "/macros/asimov_fits/";
  // Create TFile with the csv data and a slice of TTree with the data needed for the fit
  TFile file_out("FitOverEstimation.root", "RECREATE");
  TTree tree_out("data_tree", "Data for unbinned fit");

  tree_out.ReadFile(MONADIR + "/output/csv/SensChi2Inf/AsimovFit4BinsNO_PercentageOfMC/AsimovFit4BinsNO_PercentageOfMC_0-99.csv", 
                    "percentage/D:Ebins/I:ctBins/I:sens_0/D:sens_1/D:sens_2/D:sens_3/D:fit_chi2/D", ',');

  tree_out.Write();
  file_out.Close();

//  TFile file_in("UnbinnedFit.root", "READ");
//  TTree* tree_in = (TTree*)file_in.Get("data_tree");
//  tree_in->SetDirectory(0); // Decouple the tree from the file and keep it in memory
//  file_in.Close(); // Close the file so that RooFit can make plots on screen. (bug?)

  TTree* tree_in = new TTree("data_tree", "Data for unbinned fit");
  tree_in->ReadFile(MONADIR + "/output/csv/SensChi2Inf/AsimovFit4BinsNO_PercentageOfMC/AsimovFit4BinsNO_PercentageOfMC_0-99.csv", 
                    "percentage/D:Ebins/I:ctBins/I:sens_0/D:sens_1/D:sens_2/D:sens_3/D:fit_chi2/D", ',');

//  Double_t inverted_percentage;
//  Double_t percentage;
//  tree_in->SetBranchAddress("percentage", &percentage);
//  TBranch* new_branch = tree_in->Branch("inverted_percentage", &inverted_percentage, "inverted_percentage/D");
//  Int_t nentries = tree_in->GetEntries(); // read the number of entries in tree
//  for (Int_t i = 0; i < nentries; i++) {
//    tree_in->GetEntry(i);
//    inverted_percentage = 1 / percentage;
//    new_branch->Fill();
//  }

//  TCanvas* c1 = new TCanvas("c1","c1");
//  Int_t n = tree_in->Draw("inverted_percentage:fit_chi2", "", "goff"); 
//  TGraph *g = new TGraph(n, tree_in->GetV1(), tree_in->GetV2()); 
//  g->SetMinimum(0);
// // g->SetMaximum(120);
//  g->Draw("ap"); 
// 
//  TF1 *f_tr = new TF1("f_tr", "[0] + (x / [1])^[2]", 1e-3, 10);
//  f_tr->SetParameters(1, 1, 1); // Initial values of the params.
//  g->Fit(f_tr);
//  //g->FitPanel();
//  f_tr = g->GetFunction("f_tr");

  TCanvas* c1 = new TCanvas("c1","c1");
  Int_t n = tree_in->Draw("percentage:fit_chi2", "", "goff"); 
  TGraph *g = new TGraph(n, tree_in->GetV1(), tree_in->GetV2()); 
  c1->SetLogx();
  g->SetMinimum(0);
 // g->SetMaximum(120);
  g->Draw("ap"); 

  TF1 *f_tr = new TF1("f_tr", "[0] + ([1] / x)^[2]", 1e-3, 10);
  f_tr->SetParameters(1, 1, 1); // Initial values of the params.
  g->Fit(f_tr);
  //g->FitPanel();
  f_tr = g->GetFunction("f_tr");


}
