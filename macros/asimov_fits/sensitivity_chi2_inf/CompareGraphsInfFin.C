#include "TCanvas.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TTree.h"
#include "TFile.h"
#include "TRandom.h"

TCanvas* ComparisonGraph(Int_t nPidCategories, TString ordering, Double_t yMin, Double_t yMax);

void CompareGraphsInfFin() {

  std::vector<Int_t> pid_cats = {2,3,4,5,10};

  TFile* fout = new TFile("CompareGraphsInfFin.root", "UPDATE");

  TCanvas* c2 = ComparisonGraph(2, "io", 2.5, 5);
  fout->cd();
  c2->Write("comaprison_2bins_io");

  TCanvas* c3 = ComparisonGraph(3, "io", 2.5, 5);
  fout->cd();
  c3->Write("comaprison_3bins_io");

  TCanvas* c4 = ComparisonGraph(4, "io", 2.5, 5);
  fout->cd();
  c4->Write("comaprison_4bins_io");

  TCanvas* c5 = ComparisonGraph(5, "io", 2.5, 5);
  fout->cd();
  c5->Write("comaprison_5bins_io");

  for (auto c: {c2,c3,c4,c5}) delete c;

  TCanvas* c6 = ComparisonGraph(2, "no", 2.5, 5.5);
  fout->cd();
  c6->Write("comaprison_2bins_no");

  TCanvas* c7 = ComparisonGraph(3, "no", 2, 6);
  fout->cd();
  c7->Write("comaprison_3bins_no");

  TCanvas* c8 = ComparisonGraph(4, "no", 2.5, 6);
  fout->cd();
  c8->Write("comaprison_4bins_no");

  TCanvas* c9 = ComparisonGraph(5, "no", 2.5, 6);
  fout->cd();
  c9->Write("comaprison_5bins_no");

  fout->Close();

}


TCanvas* ComparisonGraph(Int_t nPidCategories, TString ordering, Double_t yMin, Double_t yMax) {
  TCanvas* c = new TCanvas(Form("c%i", nPidCategories), Form("c%i", nPidCategories));
  TLegend* leg = new TLegend(0.1, 0.7, 0.5, 0.9);

  TFile *fNormalQ = TFile::Open("FitOverEstimationCompareGraphs.root", "READ");
  TFile *fRandomQ = TFile::Open("./cross_check/FitOverEstimationCompareGraphs.root", "READ");

  TString getString = Form("comaprison_%ibins_", nPidCategories) + ordering;

  TGraphErrors* gBinNormQFinIO = (TGraphErrors*)fNormalQ->Get(getString + "_fin"); // N Bins, Normal Q, Finite statistics, IO
  gBinNormQFinIO->SetLineColor(kBlue+1);

  TGraphErrors* gBinNormQInfIO = (TGraphErrors*)fNormalQ->Get(getString + "_inf");
  gBinNormQInfIO->SetLineColor(kBlue+1);
  gBinNormQInfIO->SetLineStyle(7);

  TGraphErrors* gBinRandQFinIO = (TGraphErrors*)fRandomQ->Get(getString + "_fin");
  gBinRandQFinIO->SetLineColor(kGreen-3);

  TGraphErrors* gBinRandQInfIO = (TGraphErrors*)fRandomQ->Get(getString + "_inf");
  gBinRandQInfIO->SetLineColor(kGreen-3);
  gBinRandQInfIO->SetLineStyle(7);
    

  gBinNormQFinIO->SetMinimum(yMin);
  gBinNormQFinIO->SetMaximum(yMax);

  gBinNormQFinIO->Draw();
  gBinNormQInfIO->Draw("same");
  gBinRandQFinIO->Draw("same");
  gBinRandQInfIO->Draw("same");

  gBinNormQFinIO->SetTitle(Form("Sensitivity comparison for %i PID categories", nPidCategories));
  gBinNormQFinIO->GetXaxis()->SetTitle("#theta_{23}");
  gBinNormQFinIO->GetYaxis()->SetTitle("#sqrt{ #Delta #chi^{2} }");

  ordering.ToUpper();
  leg->AddEntry(gBinNormQFinIO, "#Delta #chi^{2}, normal Q, " + ordering, "lp");
  leg->AddEntry(gBinNormQInfIO, "#LT #Delta #chi^{2} #GT at #infty, normal Q, " + ordering, "lpe");
  leg->AddEntry(gBinRandQFinIO, "#Delta #chi^{2}, random Q, " + ordering, "lp");
  leg->AddEntry(gBinRandQInfIO, "#LT #Delta #chi^{2} #GT at #infty, random Q, " + ordering, "lpe");
  leg->Draw();

  return c;
};
