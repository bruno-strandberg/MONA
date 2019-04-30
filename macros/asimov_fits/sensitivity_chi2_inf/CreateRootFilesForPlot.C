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

void CreateRootFilesForPlot() {

  std::vector<Int_t> pid_cats = {2,3,4,5,10};

  TFile *f = TFile::Open("./Chi2Extraplotation.root", "RECREATE");
  for (auto pid: pid_cats) {

    TString MONADIR = (TString)getenv("MONADIR") + "/macros/asimov_fits/";
    TString file = MONADIR + CsvFilePath(pid, "NO");
    TString columns = CsvColumns(pid);

    TTree* tree_in = new TTree("data_tree", "data for fit");
    tree_in->ReadFile(file, columns, ',');

    TTree* t_fit_chi2 = CalculateChi2Infinite(tree_in, Form("%ibins_no", pid));
    delete tree_in;

    t_fit_chi2->Write();
  }
  for (auto pid: pid_cats) {

    TString MONADIR = (TString)getenv("MONADIR") + "/macros/asimov_fits/";
    TString file = MONADIR + CsvFilePath(pid, "IO");
    TString columns = CsvColumns(pid);

    TTree* tree_in = new TTree("data_tree", "data for fit");
    tree_in->ReadFile(file, columns, ',');

    TTree* t_fit_chi2 = CalculateChi2Infinite(tree_in, Form("%ibins_io", pid));
    delete tree_in;

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
