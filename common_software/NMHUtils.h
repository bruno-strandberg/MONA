#include "TAxis.h"
#include "TMath.h"

/**
 * A namespace that collects miscellaneous useful functions. 
 */
namespace NMHUtils {

  /** Function to convert equidistant bins on a log scale to linear scale bins.
   *
   *  For example, say one wants 40 bins of equal width on a log scale from 0.1 to 100. Then call
   *  this function as vector<Double_t> low_edges = GetLogBins(40, 0.1, 100). One can then create a
   *  histogram TH1D h( "h","h", 40, &low_edges[0] ) that will have 40 bins on the x axis that are
   *  of equal size on a logarithmic scale and of different size on a linear scale.
   *
   *  \param nbins   Number of bins
   *  \param low     Lowest value on the axis
   *  \param high    Highest value on the axis
   *  \return        Vector with bin low edges on linear scale
   */
  vector<Double_t> GetLogBins(Int_t nbins, Double_t low, Double_t high) {

    // greate an axis with equal bins on log scale
    TAxis log_axis(nbins, TMath::Log10(low),  TMath::Log10(high) );

    // convert the bin low edges back to linear scale
    vector<Double_t> low_edges;
    for (Int_t bin = 1; bin <= log_axis.GetNbins() + 1; bin++) {
      low_edges.push_back( TMath::Power( 10, log_axis.GetBinLowEdge(bin) ) );
    }

    //return the vector with bin edges
    return low_edges;
    
  }

}
