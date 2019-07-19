#include "NuXsec.h"
#include "TFile.h"
#include "TSystem.h"
#include <algorithm>
#include <iostream>
#include <stdexcept>

//*************************************************************************

/**
 * Constructor.

 * \param  bybins       Number of bjorken-y bins used in the analysis
 * \param  xsecfile     A root file from genie's `gspl2root` tool with neutrino xsec data, created by `MONA/apps/xsec_extractor`
 * \param  by_dist_file A roof file in specific format with energy vs bjorken-y distributions, created by `MONA/apps/bjorkeny_dists`
 */
NuXsec::NuXsec(UInt_t bybins, TString xsecfile, TString by_dist_file) {
  
  fXsecFile = xsecfile;
  fByFile   = by_dist_file;  

  //init the maps fGraphs and fByHists that holds the TGraphs with xsec data

  InitMaps(fXsecFile, fByFile, bybins);

  fByBinsSet = (bybins != 0);
  fGraphsSet = false;
  f_g_H1     = NULL;
  f_g_O16    = NULL;
  f_h_by     = NULL;

}

//*************************************************************************

/**
 * Destructor.
 */
NuXsec::~NuXsec() {

  f_g_H1  = NULL;
  f_g_O16 = NULL;
  f_h_by  = NULL;
  
  for (auto& n: fGraphs) {
    for (auto g: n.second) if (g) delete g;
  }

  for (auto& h: fByHists) {
    if (h.second) delete h.second;
  }

}

//*************************************************************************

/**
 * This function selects the interaction for which the fuction GetXsec(E) returns the xsec.
 *
 * \param  nu_flavor     Neutrino flavor
 * \param  is_cc         Interaction type (0 - nc, 1 - cc)
 * \param  is_nubar      0 - particle, 1 - antiparticle.
 *
 */
void NuXsec::SelectInteraction(Int_t nu_flavor, Bool_t is_cc, Bool_t is_nubar) {

  //--------------------------------------------------------------
  //check that the interaction type is supported (in the maps)
  //--------------------------------------------------------------  

  Bool_t supported = true;

  if ( fNu_flavs.find(nu_flavor) == fNu_flavs.end() ) {
    cout << "ERROR! NuXsec::SelectInteraction() neutrino flavor " << nu_flavor 
	 << " not supported" << endl;
    supported = false;
  }

  if ( fInt_types.find((Bool_t)is_cc) == fInt_types.end() ) {
    cout << "ERROR! NuXsec::SelectInteraction() interaction " << is_cc <<
      " not supported" << endl;
    supported = false;
  }

  if ( fP_types.find((Bool_t)is_nubar) == fP_types.end() ) {
    cout << "ERROR! NuXsec::SelectInteraction() particle type " << is_nubar <<
      " not supported" << endl;
    supported = false;
  }

  //--------------------------------------------------------------
  //If the interaction is supported, set the global graph pointers
  //to the xsec graphs of the selected interaction; set the bjorken-y distribution pointer
  //--------------------------------------------------------------  

  if (supported) {

    TString searchstr = CreateString(nu_flavor, is_cc, is_nubar);
    f_g_H1  = fGraphs[searchstr][0];
    f_g_O16 = fGraphs[searchstr][1];
    f_h_by  = fByHists[searchstr];

    fGraphsSet = true;
  }
  else {
    f_g_H1  = NULL;
    f_g_O16 = NULL;
    f_h_by  = NULL;
    fGraphsSet = false;
  }

}

//*************************************************************************

/**
 * This function returns the cross-section per nucleon in H20 in units m^2.
 *
 * The graphs store the cross-section in units 1e-38 cm^2. To convert this to m^2, the return
 * is multiplied by 1e-38 * 1e-4, the latter term coming from cm^2 = 0.0001 m^2.
 *
 * \param  E  Neutrino energy
 * \return    (2 x cs_proton + cs_oxygen)/18
 */
Double_t NuXsec::GetXsec(Double_t E) {

  if (!fGraphsSet || !f_g_H1 || !f_g_O16) {
    throw std::invalid_argument( "ERROR! NuXsec::GetXsec() graphs not set! Did you run SelectInteraction(...)?");
  }

  if ( E < fEmin || E > fEmax ) {
    throw std::invalid_argument( "ERROR! NuXsec::GetXsec() energy " + to_string(E) + " out of range, limited by the input file " + (string)fXsecFile );
  }

  Double_t cs_O16 = f_g_O16->Eval(E);
  Double_t cs_H1  = f_g_H1->Eval(E);

  return (2 * cs_H1 + cs_O16)/18. * 1e-42;

}

//*************************************************************************

/**
 * This function returns the cross-section per nucleon in H20 in units m^2.
 *
 * \param  nu_flavor     Neutrino flavor
 * \param  is_cc         Interaction type (0 - nc, 1 - cc)
 * \param  is_nubar      0 - particle, 1 - antiparticle.
 * \param                E  Neutrino energy
 * \return               (2 x cs_proton + cs_oxygen)/18
 */
Double_t NuXsec::GetXsec(Int_t nu_flavor, Bool_t is_cc, Bool_t is_nubar, Double_t E) {
  SelectInteraction(nu_flavor, is_cc, is_nubar);
  return GetXsec(E);
}

//*************************************************************************

/**
   This function specifies the fraction of events that lie in a bjorken-y bin corresponding to argument by at neutrino energy E.

   \param E   Neutrino energy
   \param by  Bjorken-y value
   \return    Fraction of events in a bjorken-y bin for given energy
 */
Double_t NuXsec::GetBYfrac(Double_t E, Double_t by) {

  if ( !fByBinsSet ) {
    throw std::invalid_argument( "ERROR! NuXsec::GetBYfrac() bjorken-y bins not specified. Please specify bjorken-y bins at construction, e.g. NuXsec(4)");
  }

  if ( !f_h_by ) {
    throw std::invalid_argument( "ERROR! NuXsec::GetBYfrac() histgram pointer not set! Did you run SelectInteraction(...)?");
  }

  // If 1 bjorken-y bin is used, return 1 without further checks. This step is necessary, as sometimes xsec data
  // is required outside the available MC data range that was used to create the bjorken-y distributions.
  if ( f_h_by->GetYaxis()->GetNbins() == 1 ) return 1.0;

  if ( E < f_h_by->GetXaxis()->GetXmin() || E >= f_h_by->GetXaxis()->GetXmax() ) {
    throw std::invalid_argument("ERROR! NuXsec::GetBYfrac() energy " + to_string(E) + " out of range, limited by the input file " + (string)fByFile );
  }

  if ( by < f_h_by->GetYaxis()->GetXmin() || by > f_h_by->GetYaxis()->GetXmax() ) {
    throw std::invalid_argument("ERROR! NuXsec::GetBYfrac() bjorken-y " + to_string(by) + " out of range." );
  }

  Int_t ebin  = f_h_by->GetXaxis()->FindBin(E);
  Int_t bybin = f_h_by->GetYaxis()->FindBin(by);

  // due to technical reasons ( TAxis range [0, 1) ) by = 1.0 fetches the over-flow bin value, which is 0
  // this does not make physical sense for bjorken-y. as by > 1.0 throws an error above, here the bybin is 
  // changed from overflow to the last bin
  if ( bybin > f_h_by->GetYaxis()->GetNbins() ) {
    bybin = bybin - 1;
  }

  return f_h_by->GetBinContent(ebin, bybin);

}

//*************************************************************************

/**
   This function specifies the fraction of events that lie in a bjorken-y bin corresponding to argument by at neutrino energy E.

   \param  nu_flavor     Neutrino flavor
   \param  is_cc         Interaction type (0 - nc, 1 - cc)
   \param  is_nubar      0 - particle, 1 - antiparticle.
   \param                E   Neutrino energy
   \param by             Bjorken-y value
   \return               Fraction of events in a bjorken-y bin for given energy
 */
Double_t NuXsec::GetBYfrac(Int_t nu_flavor, Bool_t is_cc, Bool_t is_nubar, Double_t E, Double_t by) {
  SelectInteraction(nu_flavor, is_cc, is_nubar);
  return GetBYfrac(E, by);
}

//*************************************************************************

/**
 * This function initialises the maps fGraphs and fByHists that contain TGraphs with xsec data and histograms with bjorken-y data.
 *
 * \param  xsecfile      File with xsec TGraphs.
 * \param  byfile        File with energy vs bjorken-y distribution histograms.
 * \param  bybins        Number of bjorken-y bins used in the analysis
 */
void NuXsec::InitMaps(TString xsecfile, TString byfile, UInt_t bybins) {

  TFile *fxsec = new TFile(xsecfile, "READ");
  TFile *fby   = new TFile(byfile,   "READ");
  if ( !fxsec->IsOpen() || !fby->IsOpen() ) {
    throw std::invalid_argument( "ERROR! NuXsec::InitMaps() could not find xsec or bjorken-y file!");
  }

  vector<double> maxima, minima; // vectors to define the valid energy range

  //loop over maps, use the CreateString() function to create a lookup string
  //for each combination of nu_flavor, int_type, p_type. Clone the xsec graphs/by hists to the maps.
  
  for (auto& nuflav: fNu_flavs) {
    for (auto& inttype: fInt_types) {
      for (auto& ptype: fP_types) {

	TString gname1 = "nu_" + nuflav.second + ptype.second + "_H1/tot" + inttype.second;
	TString gname2 = "nu_" + nuflav.second + ptype.second + "_O16/tot" + inttype.second;
	TString hname  = nuflav.second + ptype.second + inttype.second;      //name in file
	TString clonename = hname;                                           //name used in class
	if (inttype.first == 0) hname = "e" + ptype.second + inttype.second; //for NC only elec's
	
	vector<TGraph*> graphs;
	graphs.push_back( (TGraph*)( fxsec->Get(gname1)->Clone() ) );
	graphs.push_back( (TGraph*)( fxsec->Get(gname2)->Clone() ) );
	for (auto g: graphs) {
	  maxima.push_back( g->GetXaxis()->GetXmax() );
	  minima.push_back( g->GetXaxis()->GetXmin() );
	}

	TH2D *hist = (TH2D*)fby->Get(hname)->Clone();
	hist->SetDirectory(0);
	hist->SetNameTitle(clonename, clonename);

	// calculate the rebinning for bjorken-y histograms, such that the input histograms
	// are left with the desired number of bjorken-y bins. The constructor can be used without
	// specifying binning, in which case bybins = 0 and only 1 bjorken-y bin is used
	Int_t existing_bins = hist->GetYaxis()->GetNbins();
	Int_t rebinning     = existing_bins;
	Int_t modulus       = 0;
	if (bybins > 0) {
	  rebinning = existing_bins/(Int_t)bybins;
	  modulus   = existing_bins%(Int_t)bybins;
	}

	if ( modulus != 0 ) {
	  throw std::invalid_argument( "ERROR! NuXsec::InitMaps() the input file has " + to_string(existing_bins) + " bjorken-y bins, which cannot be re-binned to " + to_string(bybins) );
	}

	hist->Rebin2D(1, rebinning);

	TString namestr = CreateString( nuflav.first, (Bool_t)inttype.first, (Bool_t)ptype.first );
	
	fGraphs.insert( pair<TString, vector<TGraph*> > {namestr, graphs} );
	fByHists.insert( pair<TString, TH2D*> {namestr, hist} );
	
      }
    }
  }

  fxsec->Close();
  delete fxsec;

  fby->Close();
  delete fby;

  // define the maximum and minimum validity range
  fEmin = *std::max_element( minima.begin(), minima.end() );
  fEmax = *std::min_element( maxima.begin(), maxima.end() );

}

//*************************************************************************

/**
 * This function creates a string for searching/populating the map with xsec graphs.
 *
 * \param  nu_flavor     Neutrino flavor
 * \param  is_cc         Interaction type (0 - nc, 1 - cc)
 * \param  is_nubar      0 - particle, 1 - antiparticle.
 * \return               A search string for the fGraphs map.
 */
TString NuXsec::CreateString(Int_t nu_flavor, Bool_t is_cc, Bool_t is_nubar) {

  return fNu_flavs[nu_flavor] + fP_types[(Int_t)is_nubar] + fInt_types[(Int_t)is_cc];

}
