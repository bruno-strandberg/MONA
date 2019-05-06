#include "TCanvas.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TTree.h"
#include "TFile.h"
#include "TRandom.h"

void FitOverEstimationCompareGraphsIO() {

  std::vector<Int_t> pid_cats = {2,3,4,5,10};

  TCanvas* c1 = new TCanvas("c1", "c1");
  TLegend* leg = new TLegend(0.1, 0.7, 0.5, 0.9);


  TFile *f = TFile::Open("data_chi2_th23_extraplotation.root", "READ");
     
  TTree* t_2no  = (TTree*)f->Get("fit_tree_2bins_no");
  TTree* t_3no  = (TTree*)f->Get("fit_tree_3bins_no");
  TTree* t_4no  = (TTree*)f->Get("fit_tree_4bins_no");
  TTree* t_5no  = (TTree*)f->Get("fit_tree_5bins_no");
  TTree* t_10no = (TTree*)f->Get("fit_tree_10bins_no");

  // Draw the fitted chi2_inf
  Int_t n2 = t_2no->Draw("th23:sqrt_fit_chi2", "");
  TGraph *g2 = new TGraph(n2, t_2no->GetV1(), t_2no->GetV2()); 
  Int_t n3 = t_3no->Draw("th23:sqrt_fit_chi2", "");
  TGraph *g3 = new TGraph(n3, t_3no->GetV1(), t_3no->GetV2()); 
  Int_t n4 = t_4no->Draw("th23:sqrt_fit_chi2", "");
  TGraph *g4 = new TGraph(n4, t_4no->GetV1(), t_4no->GetV2()); 
  Int_t n5 = t_5no->Draw("th23:sqrt_fit_chi2", "");
  TGraph *g5 = new TGraph(n5, t_5no->GetV1(), t_5no->GetV2()); 
  Int_t n10 = t_10no->Draw("th23:sqrt_fit_chi2", "");
  TGraph *g10 = new TGraph(n10, t_10no->GetV1(), t_10no->GetV2()); 

  g2->SetLineColor(kBlue+1);
  g2->SetMarkerColor(kBlue+1);
  g2->SetMarkerStyle(21);

  g3->SetLineColor(kGreen+3);
  g3->SetMarkerColor(kGreen+3);
  g3->SetMarkerStyle(21);

  g4->SetLineColor(kRed+2);
  g4->SetMarkerColor(kRed+2);
  g4->SetMarkerStyle(21);

  g5->SetLineColor(kAzure+1);
  g5->SetMarkerColor(kAzure+1);
  g5->SetMarkerStyle(21);

  g10->SetLineColor(kBlack);
  g10->SetMarkerColor(kBlack);
  g10->SetMarkerStyle(21);

  g2->SetMinimum(0);
  g2->SetMaximum(7);
  g2->Draw("apl"); 
  g3->Draw("pl"); 
  g4->Draw("pl"); 
  g5->Draw("pl"); 
  g10->Draw("pl"); 

    
  g2->SetTitle("#LT #Delta #chi^{2} #GT at infinite statistics ");
  g2->GetXaxis()->SetTitle("#theta_{23}");
  g2->GetYaxis()->SetTitle("#sqrt{ #Delta #chi^{2} }");


  leg->AddEntry(g2, "#LT #Delta #chi^{2} #GT IO at #infty statistics 2 PID categories", "lp");
  leg->AddEntry(g3, "#LT #Delta #chi^{2} #GT IO at #infty statistics 3 PID categories", "lp");
  leg->AddEntry(g4, "#LT #Delta #chi^{2} #GT IO at #infty statistics 4 PID categories", "lp");
  leg->AddEntry(g5, "#LT #Delta #chi^{2} #GT IO at #infty statistics 5 PID categories", "lp");
  leg->AddEntry(g10, "#LT #Delta #chi^{2} #GT IO at #infty statistics 10 PID categories", "lp");
  leg->Draw();
}

