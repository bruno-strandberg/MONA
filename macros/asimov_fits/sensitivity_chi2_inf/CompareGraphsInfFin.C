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

  TCanvas* c21 = ComparisonGraph(21, "io", 2.5, 5);
  fout->cd();
  c21->Write("comaprison_21bins_io");
  c21->SaveAs("CompareInfFin21IO.pdf");

  TCanvas* c22 = ComparisonGraph(22, "io", 2.5, 5);
  fout->cd();
  c22->Write("comaprison_22bins_io");
  c22->SaveAs("CompareInfFin22IO.pdf");

  TCanvas* c31 = ComparisonGraph(31, "io", 2.5, 5);
  fout->cd();
  c31->Write("comaprison_31bins_io");
  c31->SaveAs("CompareInfFin31IO.pdf");

  TCanvas* c32 = ComparisonGraph(32, "io", 2.5, 5);
  fout->cd();
  c32->Write("comaprison_32bins_io");
  c32->SaveAs("CompareInfFin32IO.pdf");

  TCanvas* c41 = ComparisonGraph(41, "io", 2.5, 5);
  fout->cd();
  c41->Write("comaprison_41bins_io");
  c41->SaveAs("CompareInfFin41IO.pdf");

  for (auto c: {c2,c3,c4,c5,c21,c22,c31,c32,c41}) delete c;

  TCanvas* k2 = ComparisonGraph(2, "no", 2.5, 6);
  fout->cd();
  k2->Write("comaprison_2bins_no");
  k2->SaveAs("CompareInfFin2NO.pdf");

  TCanvas* k3 = ComparisonGraph(3, "no", 2.5, 6);
  fout->cd();
  k3->Write("comaprison_3bins_no");
  k3->SaveAs("CompareInfFin3NO.pdf");

  TCanvas* k4 = ComparisonGraph(4, "no", 2.5, 6);
  fout->cd();
  k4->Write("comaprison_4bins_no");
  k4->SaveAs("CompareInfFin4NO.pdf");

  TCanvas* k5 = ComparisonGraph(5, "no", 2.5, 6);
  fout->cd();
  k5->Write("comaprison_5bins_no");
  k5->SaveAs("CompareInfFin5NO.pdf");

  TCanvas* k21 = ComparisonGraph(21, "no", 2.5, 6);
  fout->cd();
  k21->Write("comaprison_21bins_no");
  k21->SaveAs("CompareInfFin21NO.pdf");

  TCanvas* k22 = ComparisonGraph(22, "no", 2.5, 6);
  fout->cd();
  k22->Write("comaprison_22bins_no");
  k22->SaveAs("CompareInfFin22NO.pdf");

  TCanvas* k31 = ComparisonGraph(31, "no", 2.5, 6);
  fout->cd();
  k31->Write("comaprison_31bins_no");
  k31->SaveAs("CompareInfFin31NO.pdf");

  TCanvas* k32 = ComparisonGraph(32, "no", 2.5, 6);
  fout->cd();
  k32->Write("comaprison_32bins_no");
  k32->SaveAs("CompareInfFin32NO.pdf");

  TCanvas* k41 = ComparisonGraph(41, "no", 2.5, 6);
  fout->cd();
  k41->Write("comaprison_41bins_no");
  k41->SaveAs("CompareInfFin41NO.pdf");


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
  gBinNormQFinIO->SetLineStyle(1);
  gBinNormQFinIO->SetLineWidth(2);

  TGraphErrors* gBinNormQInfIO = (TGraphErrors*)fNormalQ->Get(getString + "_inf");
  gBinNormQInfIO->SetLineColor(kRed+2);
  gBinNormQInfIO->SetLineWidth(2);

  TGraphErrors* gBinRandQFinIO = (TGraphErrors*)fRandomQ->Get(getString + "_fin");
  gBinRandQFinIO->SetLineColor(kBlue+1);
  gBinRandQFinIO->SetLineStyle(1);
  gBinRandQFinIO->SetLineWidth(2);

  TGraphErrors* gBinRandQInfIO = (TGraphErrors*)fRandomQ->Get(getString + "_inf");
  gBinRandQInfIO->SetLineColor(kBlue+1);
  gBinRandQInfIO->SetLineWidth(2);
    

  gBinNormQFinIO->SetMinimum(yMin);
  gBinNormQFinIO->SetMaximum(yMax);

  gBinNormQFinIO->Draw();
  //gBinNormQInfIO->Draw("same");
  gBinRandQFinIO->Draw("same");
  //gBinRandQInfIO->Draw("same");

  ordering.ToUpper();
  gBinNormQFinIO->SetTitle(Form("Sensitivity comparison: %i PID categories " + ordering, nPidCategories));
  gBinNormQFinIO->GetXaxis()->SetTitle("#theta_{23}");
  gBinNormQFinIO->GetYaxis()->SetTitle("#sqrt{ #Delta #chi^{2} }");

  leg->AddEntry(gBinNormQFinIO, "#Delta #chi^{2}, real track score ", "lp");
  //leg->AddEntry(gBinNormQInfIO, "#Delta #chi^{2} #infty, real track score" , "lpe");
  leg->AddEntry(gBinRandQFinIO, "#Delta #chi^{2}, random track score", "lp");
  //leg->AddEntry(gBinRandQInfIO, "#Delta #chi^{2} #infty, random track score", "lpe");
  leg->Draw();

  return c;
};
