#include "NMHUtils.h"

#include "TAxis.h"
#include "TMath.h"
#include <fstream>
#include <sys/stat.h>
#include <stdexcept>
#include<iostream>

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
vector<Double_t> NMHUtils::GetLogBins(Int_t nbins, Double_t low, Double_t high) {

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

/** 
    Function to get bin edges that can be input to TH1/2/3.
    
    Useful when constructing TH3 histgrams with E-axis on a log scale (see `NMHUtils::GetLogBins`) and cos-theta and bjorken-y on linear scales.

    \param nbins   Number of bins
    \param low     Lowest value on the axis
    \param high    Highest value on the axis
    \return        Vector with bin low edges

*/
vector<Double_t> NMHUtils::GetBins(Int_t nbins, Double_t low, Double_t high) {

  TAxis axis( nbins, low, high );
  vector<Double_t> edges;
  
  for (Int_t bin = 1; bin <= axis.GetNbins() + 1; bin++) {
    edges.push_back( axis.GetBinLowEdge(bin) );
  }

  return edges;

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
vector<TString> NMHUtils::ReadLines(TString input_file) {
    
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
 * \return         False if file does not exist or is smaller than size; True otherwise.
 */
Bool_t NMHUtils::FileExists(TString filename, Double_t size) {
    
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
 * \param xhigh        X-bins with centers above xhigh are excluded
 * \param ylow         Y-bins with centers below ylow are excluded
 * \param yhigh        Y-bins with centers above yhigh are excluded
 * \return             A std::tuple with elements: 
 *                     0) a pointer to a histgram with bin-by-bin asymmetries;
 *                     1) the quantity \f$ \sqrt{\sum_{i} A_i^2} \f$ (combined asymmetry); 
 *                     2) the quantity \f$ \Delta A$ (error on combined asymmetry);
 *                     3) the number of considered bins (degress of freedom)
 *                     
 */
std::tuple<TH2D*, Double_t, Double_t, Double_t>
NMHUtils::Asymmetry(TH2D *h1, TH2D* h2, TH2D* hcorr, TString nametitle, 
		    Double_t xlow, Double_t xhigh,
		    Double_t ylow, Double_t yhigh) {
  
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

  Double_t asym     = 0.;
  Double_t asym_err = 0.;
  Double_t Nbins    = 0.;

  for (Int_t xb = 1; xb <= h_asym->GetXaxis()->GetNbins(); xb++) {
    for (Int_t yb = 1; yb <= h_asym->GetYaxis()->GetNbins(); yb++) {
	
      Double_t N_h1 = h1->GetBinContent(xb, yb);
      Double_t N_h2 = h2->GetBinContent(xb, yb);	
      Double_t N_h1_err = h1->GetBinError(xb, yb);
      Double_t N_h2_err = h2->GetBinError(xb, yb);	

      Double_t correlation_coef = hcorr->GetBinContent(xb, yb);
	
      Double_t A     = 0;
      Double_t A_err = 0;

      if ( N_h1 > 0 ) { 
        A = (N_h1 - N_h2)/TMath::Sqrt(N_h1); 
        A_err = std::pow(0.5*(N_h1 + N_h2) / std::pow(N_h1, 1.5), 2.) * std::pow(N_h1_err, 2.) +
                std::pow(-1 / std::sqrt(N_h1), 2.) * std::pow(N_h2_err, 2.) -
                2 * correlation_coef * N_h1_err * N_h2_err * (0.5*(N_h1 + N_h2) / std::pow(N_h1, 2.));
        //Double_t Bla_err = A_err;
        if ((A_err < 0) and (std::abs(A_err) < 1e-10)) { A_err = 0.; } // This happened on day 1, has disappeared since...
        A_err = std::sqrt(A_err);
        
        //if (A_err != A_err) { cout << Bla_err << " " << N_h1_err << " " << N_h2_err << " " << N_h1 << " " << N_h2 << endl; }
      }
      else { 
        A = 0.; 
        A_err = 0.;
      }
	
      Double_t xc = h_asym->GetXaxis()->GetBinCenter(xb);
      Double_t yc = h_asym->GetYaxis()->GetBinCenter(yb);

      if ( (xc < xlow) || (xc > xhigh) || ( yc < ylow) || ( yc > yhigh ) ) {
	A = 0.;
      }
      else {
	Nbins += 1;
      }

    //  if ((xb == 11) and (yb == 7)) { cout << A << " " << A_err << endl; }
      h_asym->SetBinContent(xb, yb, A);
      h_asym->SetBinError(xb, yb, A_err);

      asym += A * A;
      asym_err += A * A * A_err * A_err; // The denominator is a function of all bins, so it is done outside of the for-loop, below. 

    }
  }

  asym = TMath::Sqrt( asym );
  asym_err = 1. / asym * TMath::Sqrt( asym_err );

  return std::make_tuple(h_asym, asym, asym_err, Nbins);

}


//****************************************************************************

/**
 *  Function calculate the value and error propagated on it for N values that are 
 *  summed quadratically.
 *
 *  Sum is defined as \f$ S = \sqrt{ \Sigma_i x_i^2 \f$ for values \f${x_1, ..., x_N} }\f$
 *  with errors \f${\Delta x_1, ..., \Delta x_N } \f$.
 *
 * \param values       Vector of length N with the values to be summed
 * \param errors       Vector of length N with the errors to be propagated
 * \return             A std::tuple with elements:
 *                     0) The square-root of the sum of squares of the values
 *                     1) The propagated error of the sum.
 */
std::tuple<Double_t, Double_t> NMHUtils::SquaredSumErrorProp(std::vector<Double_t> values, std::vector<Double_t> errors) {
  Double_t size = values.size();
  Double_t total_value = 0;
  Double_t total_error = 0;
  for (int i = 0; i < size; i++) {
    total_value += values[i] * values[i];
    total_error += values[i] * values[i] * errors[i] * errors[i];
  }

  total_value = std::sqrt(total_value);
  total_error = 1 / total_value * std::sqrt(total_error);

  return std::make_tuple(total_value, total_error);
}
