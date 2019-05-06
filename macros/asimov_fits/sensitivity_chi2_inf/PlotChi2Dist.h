
#include "TBranch.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1D.h"
#include "TTree.h"
#include "TStyle.h"

#include <iostream>

#include "EventFilter.h"
#include "NMHUtils.h"

using namespace std;

std::vector<TH1D*> PlotChi2Distribution(TTree* tree, Int_t use_th23=-1, Int_t n_slices_in_graph=10) {

  const Int_t N_SLICES = n_slices_in_graph;

  TTree* tree_in = tree;
  TBranch* b_percentage    = tree_in->GetBranch("percentage");
  TBranch* b_th23          = tree_in->GetBranch("th23");
  TBranch* b_fit_chi2      = tree_in->GetBranch("fit_chi2");
  TBranch* b_sqrt_fit_chi2 = tree_in->GetBranch("sqrt_fit_chi2");

  Double_t percentage;
  Double_t th23;
  Double_t fit_chi2;
  Double_t sqrt_fit_chi2;

  b_percentage->SetAddress(&percentage);
  b_th23->SetAddress(&th23);
  b_fit_chi2->SetAddress(&fit_chi2);
  b_sqrt_fit_chi2->SetAddress(&sqrt_fit_chi2);

  // Percentages
  std::vector<Double_t> PERCENTAGES = NMHUtils::GetLogBins(N_SLICES, 1, 100);
  std::map<Double_t, Double_t> min_map;
  std::map<Double_t, Double_t> max_map;

  // Make the maps for min and max per histogram
  for (Int_t i = 0; i <= N_SLICES; i++) {
    Double_t _min = 1E30;
    Double_t _max = 0;
    for (Int_t j = 0; j < tree_in->GetEntries(); j++) {
      tree_in->GetEntry(j);
      if (th23 == use_th23) continue;

      // For some reason percentage==PERCENTAGES[i] does not work.
      if (std::abs((percentage - PERCENTAGES[i]/100.) / percentage) < 1e-3) {
        _min = std::min(_min, fit_chi2);
        _max = std::max(_max, fit_chi2);
      }
    }

    // Debugging
    //cout << "finished looping through file, inserting at " << frac << "  min " << _min << " max " << _max << endl;
    min_map[PERCENTAGES[i]] = _min;
    max_map[PERCENTAGES[i]] = _max;
  }

  // Create histogram vector
  std::vector<TH1D*> th1d_vector;
  for (Int_t i = 0; i <= N_SLICES; i++) {
    Double_t perc = PERCENTAGES[i];
    TH1D* g = new TH1D( Form("h_chi2_dist_%i", i), Form("Chi2 dist at %.2f percent MC in RM", perc),
                        30, 0.9 * min_map[perc], 1.1 * max_map[perc]);
    th1d_vector.push_back(g);
  }

  // Fill histograms
  for (Int_t i = 0; i <= N_SLICES; i++) {
    for (Int_t j = 0; j < tree_in->GetEntries(); j++) {
      tree_in->GetEntry(j);
      if (th23 == use_th23) continue;

      if (std::abs((percentage - PERCENTAGES[i]/100.) / percentage) < 1e-3) {
        th1d_vector[i]->Fill(fit_chi2);
      }
    }
  }

  // Plot histograms into one canvas
  TCanvas* c1 = new TCanvas("c1", "c1", 1200, 800);

  c1->Divide(4,3);

  for (Int_t i = 0; i <= N_SLICES; i++) {
    c1->cd(i);
    th1d_vector[i]->Draw("same"); // Same needed for Divide

    TF1 *f = new TF1("f", "gaus", 1, 100);
    f->SetParameters(1, 1); // Initial values of the params.
    th1d_vector[i]->Fit(f);
    gStyle->SetOptFit(kTRUE);
  }

  return th1d_vector;

}
