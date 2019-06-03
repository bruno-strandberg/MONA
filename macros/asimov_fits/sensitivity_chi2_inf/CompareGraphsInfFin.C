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
  c2->SaveAs("CompareInfFin2IO.pdf");

  TCanvas* c3 = ComparisonGraph(3, "io", 2.5, 5);
  fout->cd();
  c3->Write("comaprison_3bins_io");
  c3->SaveAs("CompareInfFin3IO.pdf");

  TCanvas* c4 = ComparisonGraph(4, "io", 2.5, 5);
  fout->cd();
  c4->Write("comaprison_4bins_io");
  c4->SaveAs("CompareInfFin4IO.pdf");

  TCanvas* c5 = ComparisonGraph(5, "io", 2.5, 5);
  fout->cd();
  c5->Write("comaprison_5bins_io");
  c5->SaveAs("CompareInfFin5IO.pdf");

  for (auto c: {c2,c3,c4,c5}) delete c;

  TCanvas* c6 = ComparisonGraph(2, "no", 2.5, 6);
  fout->cd();
  c6->Write("comaprison_2bins_no");
  c6->SaveAs("CompareInfFin2NO.pdf");

  TCanvas* c7 = ComparisonGraph(3, "no", 2.5, 6);
  fout->cd();
  c7->Write("comaprison_3bins_no");
  c7->SaveAs("CompareInfFin3NO.pdf");

  TCanvas* c8 = ComparisonGraph(4, "no", 2.5, 6);
  fout->cd();
  c8->Write("comaprison_4bins_no");
  c8->SaveAs("CompareInfFin4NO.pdf");

  TCanvas* c9 = ComparisonGraph(5, "no", 2.5, 6);
  fout->cd();
  c9->Write("comaprison_5bins_no");
  c9->SaveAs("CompareInfFin5NO.pdf");

  fout->Close();


}


TCanvas* ComparisonGraph(Int_t nPidCategories, TString ordering, Double_t yMin, Double_t yMax) {
  TCanvas* c = new TCanvas(Form("c%i", nPidCategories), Form("c%i", nPidCategories));
  TLegend* leg = new TLegend(0.1, 0.7, 0.5, 0.9);

  TFile *fNormalQ = TFile::Open("FitOverEstimationCompareGraphs.root", "READ");
  TFile *fRandomQ = TFile::Open("./cross_check/FitOverEstimationCompareGraphs.root", "READ");

  TString getString = Form("comaprison_%ibins_", nPidCategories) + ordering;

  TGraphErrors* gBinNormQFinIO = (TGraphErrors*)fNormalQ->Get(getString + "_fin"); // N Bins, Normal Q, Finite statistics, IO
  gBinNormQFinIO->SetLineColor(kRed+2);
  gBinNormQFinIO->SetLineStyle(7);

  TGraphErrors* gBinNormQInfIO = (TGraphErrors*)fNormalQ->Get(getString + "_inf");
  gBinNormQInfIO->SetLineColor(kRed+2);
  gBinNormQInfIO->SetLineWidth(2);

  TGraphErrors* gBinRandQFinIO = (TGraphErrors*)fRandomQ->Get(getString + "_fin");
  gBinRandQFinIO->SetLineColor(kBlue+1);
  gBinRandQFinIO->SetLineStyle(7);

  TGraphErrors* gBinRandQInfIO = (TGraphErrors*)fRandomQ->Get(getString + "_inf");
  gBinRandQInfIO->SetLineColor(kBlue+1);
  gBinRandQInfIO->SetLineWidth(2);
    

  gBinNormQFinIO->SetMinimum(yMin);
  gBinNormQFinIO->SetMaximum(yMax);

  gBinNormQFinIO->Draw();
  gBinNormQInfIO->Draw("same");
  gBinRandQFinIO->Draw("same");
  gBinRandQInfIO->Draw("same");

  ordering.ToUpper();
  gBinNormQFinIO->SetTitle(Form("Sensitivity comparison: %i PID categories " + ordering, nPidCategories));
  gBinNormQFinIO->GetXaxis()->SetTitle("#theta_{23}");
  gBinNormQFinIO->GetYaxis()->SetTitle("#sqrt{ #Delta #chi^{2} }");

  leg->AddEntry(gBinNormQFinIO, "#Delta #chi^{2}, real track score ", "lp");
  leg->AddEntry(gBinNormQInfIO, "#Delta #chi^{2} #infty, real track score" , "lpe");
  leg->AddEntry(gBinRandQFinIO, "#Delta #chi^{2}, random track score", "lp");
  leg->AddEntry(gBinRandQInfIO, "#Delta #chi^{2} #infty, random track score", "lpe");
  leg->Draw();

  return c;
};
