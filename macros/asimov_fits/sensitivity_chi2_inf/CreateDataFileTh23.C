#include "TCanvas.h"
#include "TF1.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TTree.h"
#include "TFile.h"
#include "TRandom.h"

void WriteSqrtChi2(TTree* tree);
TString CsvFilePath(Int_t n, string ordering);
TString CsvColumns(Int_t n);

void CreateDataFileTh23() {

  TString MONADIR = (TString)getenv("MONADIR") + "/macros/asimov_fits/";
  // Create TFile with the csv data and a slice of TTree with the data needed for the fit
  TFile file_out("data_chi2_th23.root", "RECREATE");

  std::vector<Int_t> pid_cats = {2,3,4,5,10};

  for (auto order: {"no", "io"}) {
    for (auto pid: pid_cats) {
      TTree* tree_out = new TTree( Form("data_tree_%i_%s", pid, order), "Data for chi2_inf fit");

      tree_out->ReadFile(MONADIR + CsvFilePath(pid, order), CsvColumns(pid), ',');

      WriteSqrtChi2(tree_out);

      tree_out->Write();
    }
  }

  file_out.Close();
}

void WriteSqrtChi2(TTree* tree) {
  // Branch to read from
  Double_t fit_chi2;
  tree->SetBranchAddress("fit_chi2", &fit_chi2);

  // Branch to write to
  Double_t sqrt_fit_chi2;
  TBranch* sqrt_branch = tree->Branch("sqrt_fit_chi2", &sqrt_fit_chi2, "sqrt_fit_chi2/D");

  Int_t nentries = tree->GetEntries(); // read the number of entries in tree
  for (Int_t i = 0; i < nentries; i++) {
    tree->GetEntry(i);
    sqrt_fit_chi2 = TMath::Sqrt(fit_chi2);
    sqrt_branch->Fill();
  }
}

TString CsvFilePath(Int_t n, string ordering) {
  TString filepath;
  for (auto & c: ordering) c = std::toupper(c); // Make the string upper case

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
