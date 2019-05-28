#ifndef NMHUtils_h
#define NMHUtils_h

#include "TH1.h"
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
  Bool_t           HeaderParameterMatch(TString file1, TString file2, TString parameter_name);
  Bool_t           DatatagMatch(TString file1, TString file2);

  std::tuple<TH1*, Double_t, Double_t> 
    Asymmetry(TH1 *h1, TH1* h2, TString nametitle, 
          Double_t xlow = -1e10, Double_t xhigh = 1e10,
          Double_t ylow = -1e10, Double_t yhigh = 1e10,
          Double_t zlow = -1e10, Double_t zhigh = 1e10,
          Bool_t simple_error=false);
  std::tuple<Double_t, Double_t> 
    SquaredSumErrorProp(std::vector<std::pair<Double_t, Double_t>> value_and_error);
};

#endif
