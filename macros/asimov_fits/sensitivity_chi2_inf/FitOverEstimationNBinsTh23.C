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
TString CsvFilePath(Int_t n, string ordering);
TString CsvColumns(Int_t n);

void FitOverEstimation2BinsTh23(Int_t nbins=2) {

  TString MONADIR = (TString)getenv("MONADIR") + "/macros/asimov_fits/";
  TString file_no = MONADIR + CsvFilePath(nbins, "NO");
  TString file_io = MONADIR + CsvFilePath(nbins, "IO");
  TString columns_no = CsvColumns(nbins);
  TString columns_io = CsvColumns(nbins);

  TTree* tree_in_no = new TTree("data_tree_no", "Data for unbinned fit");
  tree_in_no->ReadFile(file_no, columns_no, ',');
  TTree* tree_in_io = new TTree("data_tree_io", "Data for unbinned fit");
  tree_in_io->ReadFile(file_io, columns_io, ',');

  TTree* t_fit_chi2_2bins_no = CalculateChi2Infinite(tree_in_no, "2BinsNO");
  TTree* t_fit_chi2_2bins_io = CalculateChi2Infinite(tree_in_io, "2BinsIO");

  TCanvas* c1 = new TCanvas("c1", "c1");

  // Draw the regular chi2, no/io
  Int_t n1 = tree_in_no->Draw("th23:fit_chi2", "percentage==1", "goff"); 
  TGraph *g_no = new TGraph(n1, tree_in_no->GetV1(), tree_in_no->GetV2()); 
  g_no->Sort(&TGraph::CompareX);

  Int_t n2 = tree_in_io->Draw("th23:fit_chi2", "percentage==1", "goff"); 
  TGraph *g_io = new TGraph(n2, tree_in_io->GetV1(), tree_in_io->GetV2()); 
  g_io->Sort(&TGraph::CompareX);


  // Draw the fitted chi2_inf for 2 bins, no/io
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
  g->SetMaximum(30);
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
  g->GetYaxis()->SetTitle("#Delta #chi^{2}");


  TLegend* leg = new TLegend(0.1, 0.75, 0.3, 0.9);
  leg->AddEntry(g, "#LT #Delta #chi^{2} #GT IO at #infty", "lp");
  leg->AddEntry(h, "#LT #Delta #chi^{2} #GT NO at #infty", "lp");
  leg->AddEntry(g_no, "#Delta #chi^{2} NO", "lp");
  leg->AddEntry(g_io, "#Delta #chi^{2} IO", "lp");
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
    TFitResultPtr result = g->Fit(f_tr, "S"); // Q = quiet, S = result in pointer
    fit_chi2 = result->Value(0);
    th23 = i;

    t_fit_chi2->Fill();
  }

  return t_fit_chi2;
}


TString CsvFilePath(Int_t n, string ordering) {
  TString filepath;

  switch (n) {
    case 2:
      filepath = Form("/output/csv/SensChi2Inf/AsimovFit%sTh23Range_PercentageOfMC/AsimovFit%sTh23Range_PercentageOfMC_0-99.csv", ordering.c_str(), ordering.c_str());
      break;
    case 3:
      filepath = Form("/output/csv/SensChi2Inf/AsimovFit%iBins%sTh23Range_PercentageOfMC/AsimovFit%iBins%sTh23Range_PercentageOfMC_0-99.csv", n, ordering.c_str(), n, ordering.c_str());
      break;
    case 4:
      filepath = Form("/output/csv/SensChi2Inf/AsimovFit%iBins%sTh23Range_PercentageOfMC/AsimovFit%iBins%sTh23Range_PercentageOfMC_0-99.csv", n, ordering.c_str(), n, ordering.c_str());
      break;
    case 5:
      filepath = Form("/output/csv/SensChi2Inf/AsimovFit%iBins%sTh23Range_PercentageOfMC/AsimovFit%iBins%sTh23Range_PercentageOfMC_0-99.csv", n, ordering.c_str(), n, ordering.c_str());
      break;
    case 10:
      filepath = Form("/output/csv/SensChi2Inf/AsimovFit%iBins%sTh23Range_PercentageOfMC/AsimovFit%iBins%sTh23Range_PercentageOfMC_0-99.csv", n, ordering.c_str(), n, ordering.c_str());
      break;
  }

  return filepath;
}


TString CsvColumns(Int_t n) {
  TString columns;

  switch (n) {
    case 2:
      columns = "percentage/D:Ebins/I:ctBins/I:th23/D:sinSqTh23/D:sens_0/D:sens_err_0/D:sens_1/D:sens_err_1/D:fit_chi2/D";
      break;
    case 3:
      columns = "percentage/D:Ebins/I:ctBins/I:th23/D:sinSqTh23/D:sens_0/D:sens_err_0/D:sens_1/D:sens_err_1/D:sens_2/D:sens_err_2/D:fit_chi2/D";
      break;
    case 4:
      columns = "percentage/D:Ebins/I:ctBins/I:th23/D:sinSqTh23/D:sens_0/D:sens_err_0/D:sens_1/D:sens_err_1/D:sens_2/D:sens_err_2/D:sens_3/D:sens_err_3/D:fit_chi2/D";
      break;
    case 5:
      columns = "percentage/D:Ebins/I:ctBins/I:th23/D:sinSqTh23/D:sens_0/D:sens_err_0/D:sens_1/D:sens_err_1/D:sens_2/D:sens_err_2/D:sens_3/D:sens_err_3/D:sens_4/D:sens_err_4/D:fit_chi2/D";
      break;
    case 10:
      columns = "percentage/D:Ebins/I:ctBins/I:th23/D:sinSqTh23/D:sens_0/D:sens_err_0/D:sens_1/D:sens_err_1/D:sens_2/D:sens_err_2/D:sens_3/D:sens_err_3/D:sens_4/D:sens_err_4/D:sens_5/D:sens_err_5/D:sens_6/D:sens_err_6/D:sens_7/D:sens_err_7/D:sens_8/D:sens_err_8/D:sens_9/D:sens_err_9/D:fit_chi2/D";
      break;
  }

  return columns;
}
