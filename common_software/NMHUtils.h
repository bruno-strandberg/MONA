#ifndef NMHUtils_h
#define NMHUtils_h

#include "TH2.h"
#include <tuple>
#include<vector>

using namespace std;

/**
 * A namespace that collects miscellaneous useful functions. 
 */
namespace NMHUtils {

  vector<Double_t> GetLogBins(Int_t nbins, Double_t low, Double_t high);
  vector<Double_t> GetBins(Int_t nbins, Double_t low, Double_t high);
  vector<TString>  ReadLines(TString input_file);
  Bool_t           FileExists(TString filename, Double_t size = 0);
  TString          Getcwd();
  Bool_t           BinsMatch(TH1 *h1, TH1 *h2);

  std::tuple<TH2D*, Double_t, Double_t> 
    Asymmetry(TH2D *h1, TH2D* h2, TString nametitle, 
	      Double_t xlow = -1e10, Double_t xhigh = 1e10,
	      Double_t ylow = -1e10, Double_t yhigh = 1e10);
  std::tuple<Double_t, Double_t> 
    SquaredSumErrorProp(std::vector<double> values, std::vector<double> errors);
};

#endif
