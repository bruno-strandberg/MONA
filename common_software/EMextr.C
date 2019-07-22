#include "EMextr.h"
#include "NMHUtils.h"

#include "TFile.h"
#include "TGraph.h"

#include <stdexcept>
#include <iostream>

using namespace std;

/** Constructor.

    \param tb       Pointer to a `TH3D` histogram with the desired binning configuration.
    \param effmfile MONA format effective mass file
 */
EMextr::EMextr(TH3D* tb, TString effmfile) {  

  //-----------------------------------------------------------------------
  // init histograms
  //-----------------------------------------------------------------------

  for (auto f: fFlavMap) {
    for (auto i: fIntMap) {
      for (auto p: fPolMap) {
	
	TString hname = "meffextr_" + f.second + "_" + i.second + "_" + p.second;
	fhMeff[f.first][i.first][p.first] = (TH3D*)tb->Clone(hname);
	fhMeff[f.first][i.first][p.first]->SetDirectory(0);
	fhMeff[f.first][i.first][p.first]->SetNameTitle(hname, hname);
	fhMeff[f.first][i.first][p.first]->Reset();

      }
    }
  }

  //-----------------------------------------------------------------------
  // get template histogram from effective mass file & calculate the number of bins for `EffMass` constructor
  //-----------------------------------------------------------------------

  TFile ftemp(effmfile, "READ");

  if ( !ftemp.IsOpen() ) {
    throw std::invalid_argument("ERROR! EMextr::EMextr() cannot open file " + (string)effmfile);
  }

  TH3D* htemp = (TH3D*)ftemp.Get("hgen_elec_nc_nu");

  if ( htemp == NULL ) {
    throw std::invalid_argument("ERROR! EMextr::EMextr() cannot find 'hgen_elec_nc_nu' in input file");
  }

  Int_t ebins  = GetEMbins( tb->GetXaxis(), htemp->GetXaxis(), fEBinLims.first , fEBinLims.second );
  Int_t ctbins = GetEMbins( tb->GetYaxis(), htemp->GetYaxis(), fCtBinLims.first, fCtBinLims.second);
  Int_t bybins = GetEMbins( tb->GetZaxis(), htemp->GetZaxis(), fByBinLims.first, fByBinLims.second);

  //-----------------------------------------------------------------------
  // no interpolation is performed in bjorken-y, check that input binning matches the number used
  // for the initialisation of EffMass
  //-----------------------------------------------------------------------

  if ( bybins != tb->GetZaxis()->GetNbins() ) {
    throw std::invalid_argument("ERROR! EMextr::EMextr() the requested " + to_string( tb->GetZaxis()->GetNbins() ) + " bjorken-y bins is not supported, please choose a number between 1 and 4." );
  }

  ftemp.Close();

  //-----------------------------------------------------------------------
  // init the effective mass, fill the data in MC range and finally extrapolate
  //-----------------------------------------------------------------------

  fEM = new EffMass( effmfile, ebins, ctbins, bybins );

  // only interpolate and extrapolate if binnings mismatch
  Bool_t intexpolate = !NMHUtils::BinsMatch( tb, fEM->GetMeff3DH(0,0,0) );

  FillDataMCrange(intexpolate);
  Extrapolate(intexpolate);

}

//====================================================================================

/** Destructor. */
EMextr::~EMextr() {

  for (auto f: fFlavMap) {
    for (auto i: fIntMap) {
      for (auto p: fPolMap) {
	
	delete fhMeff[f.first][i.first][p.first];

      }
    }
  }

  if ( fEM ) delete fEM;

}

//====================================================================================

/** This function returns the number of bins the axis `R` has within the minimum to maximum range of axis `E`.

    The function expects R_min < E_min and R_max > E_max

    \param R  Pointer to `TAxis` of `EMextr::fhMeff`
    \param E  Pointer to `TAxis` of `EffMass::fhMeff`
    \bins_min If `R` has less than `bins_min` in the range of `E`, `bins_min` is returned
    \bins_max If `R` has more than `bins_max` in the range of `E`, `bins_max` is returned
*/
Int_t EMextr::GetEMbins(TAxis *R, TAxis *E, Int_t bins_min, Int_t bins_max) {

  if ( R->GetXmin() > E->GetXmin() || R->GetXmax() < E->GetXmax() ) {
    throw std::invalid_argument("ERROR! EMextr::GetEMbins() expecting the requested range to be larger than the existing range.");
  }

  // MC limits in the effective mass histogram axis, fetched from the effective mass file
  Double_t min  = E->GetXmin();
  Double_t max  = E->GetXmax();
  Int_t    bins = E->GetNbins();

  // nr of bins requested in the MC range
  Int_t bins_req = R->FindBin( max ) - R->FindBin( min );

  //----------------------------------------------------------------
  // if more or less bins is requested than the input limits, return the corresponding limit
  //----------------------------------------------------------------

  if ( bins_req > bins_max ) {
    return bins_max;
  }

  if ( bins_req < bins_min ) {
    return bins_min;
  }

  //----------------------------------------------------------------
  // return the number of bins; if bins cannot be divided to bins_req, find the first suitable divisor
  //----------------------------------------------------------------

  if ( bins % bins_req == 0 ) {
    return bins_req;
  }
  else {

    while ( bins_req > bins_min && bins % bins_req != 0 ) { bins_req--; }

    return bins_req;

  }
 
}

//====================================================================================

/** This method uses the member `fEM` to interpolate the effective mass values at the bin centers of the member histograms `fhMeff`.
    \param interpolate  true - use interpolation, false - use bin center values of `EffMass`
 */
void EMextr::FillDataMCrange(Bool_t interpolate) {

  for (auto f: fFlavMap) {
    for (auto i: fIntMap) {
      for (auto p: fPolMap) {
	
	TH3D*  hmeff = fhMeff[f.first][i.first][p.first];

	// loop over the bins and fill MC data
	for (Int_t xbin = 1; xbin <= hmeff->GetXaxis()->GetNbins(); xbin++) {
	  for (Int_t ybin = 1; ybin <= hmeff->GetYaxis()->GetNbins(); ybin++) {
	    for (Int_t zbin = 1; zbin <= hmeff->GetZaxis()->GetNbins(); zbin++) {

	      Double_t E  = hmeff->GetXaxis()->GetBinCenter( xbin );
	      Double_t ct = hmeff->GetYaxis()->GetBinCenter( ybin );
	      Double_t by = hmeff->GetZaxis()->GetBinCenter( zbin );
	      
	      Double_t bc = 0.0;

	      // try calculating the effective mass, catch the error if the point is outside the
	      // range and assign that bin content to 0 (this point is extrapolated later)
	      try { 
		bc = fEM->GetMeff(f.first, i.first, p.first, E, ct, by, interpolate);
	      }
	      catch (std::invalid_argument& ia) {
		bc = 0.0;
	      }

	      hmeff->SetBinContent( xbin, ybin, zbin, bc );

	    }
	  }
	}

      }
    }
  }
  
}

//====================================================================================

/** This function extrapolates the effective mass values to energies outside the Monte-Carlo range. 
    \param extrapolate  Perform extrapolation
 */
void EMextr::Extrapolate(Bool_t extrapolate) {

  if ( !extrapolate ) return;

  for (auto f: fFlavMap) {
    for (auto i: fIntMap) {
      for (auto p: fPolMap) {
	
	TH3D*  hmeff = fhMeff[f.first][i.first][p.first];

	for (Int_t ctbin = 1; ctbin <= hmeff->GetYaxis()->GetNbins(); ctbin++) {
	  for (Int_t bybin = 1; bybin <= hmeff->GetZaxis()->GetNbins(); bybin++) {

	    //------------------------------------------------------------------------------
	    // Get the 1D histogram in energy at given (cos-theta, bjorken-y) & put non-empty bin data to a TGraph
	    //------------------------------------------------------------------------------

	    TH1D* slice = hmeff->ProjectionX("slice", ctbin, ctbin, bybin, bybin);

	    TGraph graph;

	    for (Int_t ebin = 1; ebin <= slice->GetXaxis()->GetNbins(); ebin++) {
	      
	      Double_t E  = slice->GetXaxis()->GetBinCenter( ebin );
	      Double_t bc = slice->GetBinContent( ebin );
	      
	      if ( bc > 0.) graph.SetPoint( graph.GetN(), E, bc );

	    }

	    delete slice;

	    //------------------------------------------------------------------------------
	    // use the graph to perform simple extrapolation to ranges outside the MC range
	    // the extrapolation is based on 1D splines in energy only
	    //------------------------------------------------------------------------------

	    for (Int_t ebin = 1; ebin <= slice->GetXaxis()->GetNbins(); ebin++) {

	      Double_t E  = hmeff->GetXaxis()->GetBinCenter( ebin );
	      Double_t bc = hmeff->GetBinContent(ebin, ctbin, bybin);
	      
	      if ( bc == 0. ) bc = graph.Eval( E );  // extrapolate bin content if bin was empty
	      if ( bc < 0.  ) bc = 0;                // if extrapolated to below 0 set to 0

	      hmeff->SetBinContent( ebin, ctbin, bybin, bc );

	    }

	  } // end loop over bjorken-y
	} // end loop over cos-theta

      }
    }
  }

}

//====================================================================================

/** This function returns the effective mass value in Ton.
    \param flavor       Neutrino flavor (0 - elec, 1 - muon, 2 - tau)
    \param iscc         NC = 0, CC = 1
    \param isnb         nu = 0, nub = 1
    \param E_true       True neutrino energy
    \param Ct_true      True neutrino direction
    \param By_true      True bjorken-y
    \return effective mass in Ton
 */
Double_t EMextr::GetMeff(Int_t flavor, Bool_t iscc, Bool_t isnb, 
			 Double_t E_true, Double_t Ct_true, Double_t By_true) {

  if ( fFlavMap.find(flavor) == fFlavMap.end() || fIntMap.find(iscc) == fIntMap.end() ||
       fPolMap.find(isnb) == fPolMap.end() ) {
    throw std::invalid_argument("ERROR! EMextr::GetMeff() unknown neutrino type: " + to_string(flavor) + " " + to_string(iscc) + " " + to_string(isnb) );
  }

  TH3D* hmeff = fhMeff[flavor][iscc][isnb];

  if ( E_true < hmeff->GetXaxis()->GetXmin() || E_true > hmeff->GetXaxis()->GetXmax() ) {
    throw std::invalid_argument("ERROR! EMextr::GetMeff() energy outside the histogram range");
  }

  if ( Ct_true < hmeff->GetYaxis()->GetXmin() || Ct_true > hmeff->GetYaxis()->GetXmax() ) {
    throw std::invalid_argument("ERROR! EMextr::GetMeff() cos-theta outside the histogram range");
  }

  if ( By_true < hmeff->GetZaxis()->GetXmin() || By_true > hmeff->GetZaxis()->GetXmax() ) {
    throw std::invalid_argument("ERROR! EMextr::GetMeff() bjorken-y outside the histogram range");
  }

  Int_t ebin  = hmeff->GetXaxis()->FindBin(E_true);
  Int_t ctbin = hmeff->GetYaxis()->FindBin(Ct_true);
  Int_t bybin = hmeff->GetZaxis()->FindBin(By_true);

  return hmeff->GetBinContent( ebin, ctbin, bybin );

}

//====================================================================================

/** This function returns a pointer to the 3D histogram with effective mass data in Tons.
    \param flavor       Neutrino flavor (0 - elec, 1 - muon, 2 - tau)
    \param iscc         NC = 0, CC = 1
    \param isnb         nu = 0, nub = 1
    \return Pointer to a `TH3D` with effective mass data.
 */
TH3D* EMextr::GetMeff3DH(Int_t flavor, Bool_t iscc, Bool_t isnb) {

  if ( fFlavMap.find(flavor) == fFlavMap.end() || fIntMap.find(iscc) == fIntMap.end() ||
       fPolMap.find(isnb) == fPolMap.end() ) {
    throw std::invalid_argument("ERROR! EMextr::GetMeff3DH() unknown neutrino type: " + to_string(flavor) + " " + to_string(iscc) + " " + to_string(isnb) );
  }

  return fhMeff[flavor][iscc][isnb];

}

//====================================================================================

/** This function returns a pointer to TH1D with energy on the x-axis and eff mass on the y-axis.
    \param flavor       Neutrino flavor (0 - elec, 1 - muon, 2 - tau)
    \param iscc         NC = 0, CC = 1
    \param isnb         nu = 0, nub = 1
    \param ct           cos-theta value at which the slice is drawn
    \param by           bjorken-y value at which the slice is drawn
    \return `TH1D` with the effective mass curve for given input
*/
TH1D* EMextr::GetSlice(Int_t flavor, Bool_t iscc, Bool_t isnb, Double_t ct, Double_t by) {

  if ( fFlavMap.find(flavor) == fFlavMap.end() || fIntMap.find(iscc) == fIntMap.end() ||
       fPolMap.find(isnb) == fPolMap.end() ) {
    throw std::invalid_argument("ERROR! EMextr::GetSlice() unknown neutrino type: " + to_string(flavor) + " " + to_string(iscc) + " " + to_string(isnb) );
  }

  TH3D* hmeff = fhMeff[flavor][(UInt_t)iscc][(UInt_t)isnb];

  TString nametitle = "meff_ct=" + (TString)to_string(ct) + "_by=" + (TString)to_string(by);
  Int_t ybin = hmeff->GetYaxis()->FindBin( ct );
  Int_t zbin = hmeff->GetZaxis()->FindBin( by );

  TH1D* slice = hmeff->ProjectionX(nametitle, ybin, ybin, zbin, zbin);

  return slice;

}
