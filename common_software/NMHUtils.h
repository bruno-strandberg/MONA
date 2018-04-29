#include "TAxis.h"
#include "TMath.h"
#include <fstream>

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

  //****************************************************************************

  /** A function that reads the lines of a plain-text file to a vector.
   *
   *  Empty lines and lines starting with # are ignored.
   *
   * \param input_file   input plain-text file
   * \return             a vector of TString's, each element corresponds to one line
   *
   */
  vector<TString> ReadLines(TString input_file) {
    
    vector<TString> FileLines;

    string line;
    ifstream inputfile((const Char_t*)input_file);

    if ( inputfile.is_open() ) {

      while ( getline (inputfile, line) ) {

	//ignore hashtagged and empty lines
	if ( line[0] == '#' || line == "") continue;

	FileLines.push_back( (TString)line );

      }

      inputfile.close();

    }
    else {
      cout << "ERROR: NMHUtils::ReadLines() did not find file " << input_file << endl;
    }

    return FileLines;

  }
  
}
