#include "TCanvas.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TTree.h"
#include "TFile.h"
#include "TRandom.h"

TTree* CalculateChi2Infinite(TTree* input_ttree, TString tree_nametitle);

void FitOverEstimation2BinsTh23() {

  TString MONADIR = (TString)getenv("MONADIR") + "/macros/asimov_fits/";

  TTree* tree_in_no = new TTree("data_tree_no", "Data for unbinned fit");
  tree_in_no->ReadFile(MONADIR + "/output/csv/SensChi2Inf/AsimovFitNOTh23Range_PercentageOfMC.csv", 
                    "percentage/D:Ebins/I:ctBins/I:th23/D:sinSqTh23/D:sens_0/D:sens_err_0/D:sens_1/D:sens_err_1/D:fit_chi2/D", ',');
  TTree* tree_in_io = new TTree("data_tree_io", "Data for unbinned fit");
  tree_in_io->ReadFile(MONADIR + "/output/csv/SensChi2Inf/AsimovFitIOTh23Range_PercentageOfMC.csv", 
                    "percentage/D:Ebins/I:ctBins/I:th23/D:sinSqTh23/D:sens_0/D:sens_err_0/D:sens_1/D:sens_err_1/D:fit_chi2/D", ',');

  TTree* t_fit_chi2_2bins_no = CalculateChi2Infinite(tree_in_no, "2BinsNO");
  TTree* t_fit_chi2_2bins_io = CalculateChi2Infinite(tree_in_io, "2BinsIO");

//  Int_t th23 = 0;
//  Double_t fit_chi2 = 0;
//
//  TTree* t_fit_chi2 = new TTree("fit_tree", "Chi2_inf values from fit");
//  TBranch* b_th23     = t_fit_chi2->Branch("th23", &th23, "th23/I");
//  TBranch* b_fit_chi2 = t_fit_chi2->Branch("fit_chi2", &fit_chi2, "fit_chi2/D");
//
//  for (Int_t i = 40; i <= 50; i++) {
//    Int_t n = tree_in->Draw("percentage:fit_chi2", Form("th23==%i", i), "goff"); 
//    TGraph *g = new TGraph(n, tree_in->GetV1(), tree_in->GetV2()); 
//
//    TF1 *f_tr = new TF1("f_tr", "[0] + ([1] / x)^[2]", 1e-3, 10);
//    f_tr->SetParameters(1, 1, 1); // Initial values of the params.
//    TFitResultPtr result = g->Fit(f_tr, "QS");
//    fit_chi2 = result->Value(0);
//    th23 = i;
//
//    t_fit_chi2->Fill();
//  }

  TCanvas* c1 = new TCanvas("c1", "c1");

  Int_t n = t_fit_chi2_2bins_no->Draw("th23:fit_chi2", "");
  TGraph *g = new TGraph(n, t_fit_chi2_2bins_no->GetV1(), t_fit_chi2_2bins_no->GetV2()); 

  Int_t m = t_fit_chi2_2bins_io->Draw("th23:fit_chi2", "");
  TGraph *h = new TGraph(n, t_fit_chi2_2bins_io->GetV1(), t_fit_chi2_2bins_io->GetV2()); 

  g->SetLineColor(kBlue+1);
  g->SetMarkerColor(kBlue+1);
  g->SetMarkerStyle(21);
  h->SetLineColor(kRed+2);
  h->SetMarkerColor(kRed+2);
  h->SetMarkerStyle(21);

  g->SetMinimum(0);
  g->SetMaximum(20);
  g->Draw("apl"); 

  h->Draw("pl");
  
  g->SetTitle("#Delta #chi^{2} at infinite statistics");
  g->GetXaxis()->SetTitle("#theta_{23}");
  g->GetYaxis()->SetTitle("#Delta #chi^{2}");


  TLegend* leg = new TLegend(0.7, 0.75, 0.9, 0.9);
  leg->AddEntry(g, "#Delta #chi^{2} IO", "lp");
  leg->AddEntry(h, "#Delta #chi^{2} NO", "lp");
  leg->Draw();
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

    TF1 *f_tr = new TF1("f_tr", "[0] + ([1] / x)^[2]", 1e-3, 10);
    f_tr->SetParameters(1, 1, 1); // Initial values of the params.
    TFitResultPtr result = g->Fit(f_tr, "QS");
    fit_chi2 = result->Value(0);
    th23 = i;

    t_fit_chi2->Fill();
  }

  return t_fit_chi2;
}
