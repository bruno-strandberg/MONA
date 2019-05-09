#include "TCanvas.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TTree.h"
#include "TFile.h"
#include "TRandom.h"

TTree* CalculateChi2Infinite(TTree* input_ttree, TString tree_nametitle);
TString CsvFilePath(Int_t n, string ordering);
TString CsvColumns(Int_t n);

void FitOverEstimationNBinsTh23(Int_t nbins=2) {

  TFile file("data_chi2_th23.root", "READ");

  TTree* tree_in_no = (TTree*)file.Get( Form("data_tree_%i_no", nbins ) );
  TTree* tree_in_io = (TTree*)file.Get( Form("data_tree_%i_io", nbins ) );

  TTree* t_fit_chi2_no = CalculateChi2Infinite(tree_in_no, Form("%iBinsNO", nbins));
  TTree* t_fit_chi2_io = CalculateChi2Infinite(tree_in_io, Form("%iBinsIO", nbins));

  TCanvas* c1 = new TCanvas("c1", "c1");

  // Draw the regular chi2, no/io
  Int_t n1 = tree_in_no->Draw("th23:sqrt_fit_chi2", "percentage==1", "goff"); 
  TGraph *g_no = new TGraph(n1, tree_in_no->GetV1(), tree_in_no->GetV2()); 
  g_no->Sort(&TGraph::CompareX);

  Int_t n2 = tree_in_io->Draw("th23:sqrt_fit_chi2", "percentage==1", "goff"); 
  TGraph *g_io = new TGraph(n2, tree_in_io->GetV1(), tree_in_io->GetV2()); 
  g_io->Sort(&TGraph::CompareX);


  // Draw the fitted chi2_inf for N bins, no/io
  Int_t n = t_fit_chi2_no->Draw("th23:sqrt_fit_chi2:sqrt_chi2_err", "", "goff");
  TGraphErrors *g = new TGraphErrors(n, t_fit_chi2_no->GetV1(), t_fit_chi2_no->GetV2(), 0, t_fit_chi2_no->GetV3()); 

  Int_t m = t_fit_chi2_io->Draw("th23:sqrt_fit_chi2:sqrt_chi2_err", "", "goff");
  TGraphErrors *h = new TGraphErrors(n, t_fit_chi2_io->GetV1(), t_fit_chi2_io->GetV2(), 0, t_fit_chi2_io->GetV3()); 

  g->SetLineColor(kBlue+1);
  g->SetMarkerColor(kBlue+1);
//  g->SetMarkerStyle(21);
  h->SetLineColor(kRed+2);
  h->SetMarkerColor(kRed+2);
//  h->SetMarkerStyle(21);

  g->SetMinimum(0);
  g->SetMaximum(6);
  g->Draw("apl"); 
  h->Draw("pl");

  g_no->SetLineStyle(9);
  g_no->SetLineColor(kBlue+1);
  g_io->SetLineStyle(9);
  g_io->SetLineColor(kRed+2);
  g_no->Draw("pl");
  g_io->Draw("pl");
  
  g->SetTitle((TString)"#LT #Delta #chi^{2} #GT at infinite statistics" + Form(" %i PID categories", nbins));
  g->GetXaxis()->SetTitle("#theta_{23}");
  g->GetYaxis()->SetTitle("#sqrt{ #Delta #chi^{2} }");


  TLegend* leg = new TLegend(0.1, 0.75, 0.3, 0.9);
  leg->AddEntry(g, "#LT #Delta #chi^{2} #GT IO at #infty", "lp");
  leg->AddEntry(h, "#LT #Delta #chi^{2} #GT NO at #infty", "lp");
  leg->AddEntry(g_no, "#Delta #chi^{2} IO", "lpe");
  leg->AddEntry(g_io, "#Delta #chi^{2} NO", "lpe");
  leg->Draw();
}


TTree* CalculateChi2Infinite(TTree* input_ttree, TString tree_nametitle) {
  
  Int_t th23_min = input_ttree->GetMinimum("th23");
  Int_t th23_max = input_ttree->GetMaximum("th23");

  Int_t th23;
  Double_t fit_chi2;
  Double_t fit_err;
  Double_t sqrt_fit_chi2;
  Double_t sqrt_chi2_err;

  TTree* t_fit_chi2 = new TTree("fit_tree_" + tree_nametitle, "Chi2_inf values from fit " + tree_nametitle);
  TBranch* b_th23          = t_fit_chi2->Branch("th23", &th23, "th23/I");
  TBranch* b_fit_chi2      = t_fit_chi2->Branch("fit_chi2", &fit_chi2, "fit_chi2/D");
  TBranch* b_fit_err       = t_fit_chi2->Branch("fit_err", &fit_err, "fit_err/D");
  TBranch* b_sqrt_fit_chi2 = t_fit_chi2->Branch("sqrt_fit_chi2", &sqrt_fit_chi2, "sqrt_fit_chi2/D");
  TBranch* b_sqrt_chi2_err = t_fit_chi2->Branch("sqrt_chi2_err", &sqrt_chi2_err, "sqrt_chi2_err/D");

  for (Int_t i = th23_min; i <= th23_max; i++) {
    Int_t n = input_ttree->Draw("percentage:fit_chi2", Form("th23==%i", i), "goff"); 
    TGraph *g = new TGraph(n, input_ttree->GetV1(), input_ttree->GetV2()); 

    TF1 *f_tr = new TF1("f_tr", "[0] + ([1] / x)^[2]", 1e-3, 10);
    f_tr->SetParameters(1, 1, 1); // Initial values of the params.
    TFitResultPtr result = g->Fit(f_tr, "S"); // Q = quiet, S = result in pointer

    Int_t NFreeParams = result->NFreeParameters();
    std::vector<double> pars = result->Parameters();
    std::vector<double> errs = result->Errors();

    th23 = i;
    fit_chi2 = pars[0];
    sqrt_fit_chi2 = TMath::Sqrt(fit_chi2);
    fit_err = errs[0];
    sqrt_chi2_err = fit_err / (2. * sqrt_fit_chi2); // Standard error propagation.

    t_fit_chi2->Fill();
  }

  return t_fit_chi2;
}
