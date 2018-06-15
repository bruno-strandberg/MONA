#include "TAxis.h"
#include "TMath.h"
#include <fstream>
#include <sys/stat.h>
#include <stdexcept>
#include "TH2.h"
#include <tuple>

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
  
  //****************************************************************************

  /**
   *  Function to check whether file exists.
   * \param filename Path to the file
   * \param size     File size in MB; If specified, check that file is larger than size (default: 0)
   * \return False if file does not exist or is smaller than size; True otherwise.
   */
  Bool_t FileExists(TString filename, Double_t size = 0) {
    
    struct stat buf;
    if ( (stat(filename, &buf) != 0) ) {
      return false;
    }
    else {
      return ( (Double_t)buf.st_size*1e-6 > size );
    }

  }

  //****************************************************************************

  /**
   *  Function to calculate bin-by-bin 'asymmetry' between two histograms.
   *
   *  Asymmetry is defined as \f$ A(i) = (N_{h1}^{bin i} - N_{h2}^{bin i})/\sqrt{N_{h1}^{bin i}} \f$.
   *
   * \param h1           First histogram
   * \param h2           Second histogram
   * \param nametitle    String used as name and title for the created asymmetry histogram
   * \param xlow         X-bins with centers below xlow are excluded
   * \param xlow         X-bins with centers above xhigh are excluded
   * \param ylow         Y-bins with centers below ylow are excluded
   * \param ylow         Y-bins with centers above yhigh are excluded
   * \param ReverseSign  The sign of the numerator is reversed, thus becoming 
   *                     \f$ (N_{h2}^{bin i} - N_{h1}^{bin i}) \f$.
   * \param BothDenoms   If false, the asymmetry in bin i is zero if \f$ N_{h1}^{bin i} = 0\f$.
   *                     If true, the asymmetry in bin i is \f$ -\sqrt{N_{h2}^{bin i}} \f$ if 
   *                     \f$ N_{h2}^{bin i} > 0\f$ and \f$ N_{h1}^{bin i} = 0\f$. If both are
   *                     0 the asymmetry is 0.
   * \return             A std::tuple with elements: 
   *                     0) a pointer to a histgram with bin-by-bin asymmetries;
   *                     1) the quantity \f$ \sqrt{\sum_{i} A_i^2} \f$ (combined asymmetry); 
   *                     3) the quantity \f$ ( \sum_{bin i} (N_{h1}^i - N_{h2}^i)^2 ) \f$
   *                       (the actual chi2 between the two histograms)
   *                     4) the number of considered bins (degress of freedom)
   *                     
   */
  std::tuple<TH2D*, Double_t, Double_t, Double_t> 
    Asymmetry(TH2D *h1, TH2D* h2, TString nametitle, 
	      Double_t xlow = -1e10, Double_t xhigh = 1e10,
	      Double_t ylow = -1e10, Double_t yhigh = 1e10,
	      Bool_t ReverseSign = kFALSE, Bool_t BothDenoms  = kFALSE) {

    //------------------------------------------------------------
    // check that both histograms have the same binning
    //------------------------------------------------------------

    if ( ( h1->GetXaxis()->GetNbins() != h2->GetXaxis()->GetNbins() ) || 
	 ( h1->GetYaxis()->GetNbins() != h2->GetYaxis()->GetNbins() ) ) {
      throw std::invalid_argument( "ERROR! NMHUtils::Asymmetry() input histograms bin count mismatch.");
    }

    for (Int_t i = 1; i <= h1->GetXaxis()->GetNbins(); i++) {
      if ( h1->GetXaxis()->GetBinLowEdge(i) != h2->GetXaxis()->GetBinLowEdge(i) ) {
	throw std::invalid_argument( "ERROR! NMHUtils::Asymmetry() input histograms bin low edge mismatch on X axis" );
      }
    }

    for (Int_t i = 1; i <= h1->GetYaxis()->GetNbins(); i++) {
      if ( h1->GetYaxis()->GetBinLowEdge(i) != h2->GetYaxis()->GetBinLowEdge(i) ) {
	throw std::invalid_argument( "ERROR! NMHUtils::Asymmetry() input histograms bin low edge mismatch on Y axis" );
      }
    }

    //------------------------------------------------------------
    // calculate the asymmetry
    //------------------------------------------------------------

    TH2D *h_asym = (TH2D*)h1->Clone();
    h_asym->SetNameTitle(nametitle,nametitle);
    h_asym->Reset();
    h_asym->SetDirectory(0);

    Double_t asym  = 0.;
    Double_t chi2  = 0.;
    Double_t Nbins = 0.;

    for (Int_t xb = 1; xb <= h_asym->GetXaxis()->GetNbins(); xb++) {
      for (Int_t yb = 1; yb <= h_asym->GetYaxis()->GetNbins(); yb++) {
	
	Double_t N_h1 = h1->GetBinContent(xb, yb);
	Double_t N_h2 = h2->GetBinContent(xb, yb);	

	chi2  += (N_h1 - N_h2) * (N_h1 - N_h2);
	Nbins += 1;
	
	Double_t A    = 0;

	if      (   N_h1 > 0                 ) { A = (N_h1 - N_h2)/TMath::Sqrt(N_h1); }
	else if ( ( N_h2 > 0 ) && BothDenoms ) { A = (N_h1 - N_h2)/TMath::Sqrt(N_h2); }
	else                                   { A = 0.;                              }
	
	if (ReverseSign) A = -A;
	
	Double_t xc = h_asym->GetXaxis()->GetBinCenter(xb);
	Double_t yc = h_asym->GetYaxis()->GetBinCenter(yb);

	if ( (xc < xlow) || (xc > xhigh) || ( yc < ylow) || ( yc > yhigh ) ) {
	  A = 0.;
	}

	h_asym->SetBinContent(xb, yb, A);

	asym += A * A;

      }
    }

    asym = TMath::Sqrt( asym );

    return std::make_tuple(h_asym, asym, chi2, Nbins);

  }

}
