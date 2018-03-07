#include "CScalculator.h"
#include "TFile.h"
#include<iostream>

//*************************************************************************

/**
 * Constructor.
 * \param  xsecfile  A root file in specific format with neutrino xsec data.
 */
CScalculator::CScalculator(TString xsecfile) {

  //init the maps that are used throughout the code

  nu_flavs.insert( pair<Int_t, TString> {0, "e"}   );
  nu_flavs.insert( pair<Int_t, TString> {1, "mu"}  );
  nu_flavs.insert( pair<Int_t, TString> {2, "tau"} );

  int_types.insert( pair<Int_t, TString> {0, "_nc"} );
  int_types.insert( pair<Int_t, TString> {1, "_cc"} );

  p_types.insert( pair<Int_t, TString> {0, ""} );
  p_types.insert( pair<Int_t, TString> {1, "_bar"} );

  //init the map fGraphs that holds the TGraphs with xsec data

  if ( !InitGraphs(xsecfile) ) {
    cout << "ERROR! CSCalculator::CSCalculator(). Could not find xsec file!" << endl;
  }

  fGraphsSet = false;
  f_g_H1     = NULL;
  f_g_O16    = NULL;

}

//*************************************************************************

/**
 * Destructor.
 */
CScalculator::~CScalculator() {

  f_g_H1  = NULL;
  f_g_O16 = NULL;
  
  for (auto& n: fGraphs) {
    for (auto g: n.second) if (g) delete g;
  }

}

//*************************************************************************

/**
 * This function selects the interaction for which the fuction CalculateCS(E) returns the xsec.
 *
 * \param  nu_flavor     Neutrino flavor
 * \param  is_cc         Interaction type (0 - nc, 1 - cc)
 * \param  is_nubar      0 - particle, 1 - antiparticle.
 *
 */
void CScalculator::SelectInteraction(Int_t nu_flavor, Bool_t is_cc, Bool_t is_nubar) {

  //--------------------------------------------------------------
  //check that the interaction type is supported (in the maps)
  //--------------------------------------------------------------  

  Bool_t supported = true;

  if ( nu_flavs.find(nu_flavor) == nu_flavs.end() ) {
    cout << "ERROR! CScalculator::SelectInteraction() neutrino flavor " << nu_flavor 
	 << " not supported" << endl;
    supported = false;
  }

  if ( int_types.find((Bool_t)is_cc) == int_types.end() ) {
    cout << "ERROR! CScalculator::SelectInteraction() interaction " << is_cc <<
      " not supported" << endl;
    supported = false;
  }

  if ( p_types.find((Bool_t)is_nubar) == p_types.end() ) {
    cout << "ERROR! CScalculator::SelectInteraction() particle type " << is_nubar <<
      " not supported" << endl;
    supported = false;
  }

  //--------------------------------------------------------------
  //If the interaction is supported, set the global graph pointers
  //to the xsec graphs of the selected interaction
  //--------------------------------------------------------------  

  if (supported) {

    TString searchstr = CreateString(nu_flavor, is_cc, is_nubar);
    f_g_H1  = fGraphs[searchstr][0];
    f_g_O16 = fGraphs[searchstr][1];

    fGraphsSet = true;
  }
  else {
    f_g_H1  = NULL;
    f_g_O16 = NULL;
    fGraphsSet = false;
  }

}

//*************************************************************************

/**
 * This function returns the cross-section per nucleon in H20.
 *
 * \param  E  Neutrino energy
 * \return    (2 x cs_proton + cs_oxygen)/18
 */
Double_t CScalculator::CalculateCS(Double_t E) {

  if (!fGraphsSet || !f_g_H1 || !f_g_O16) {
    cout << "ERROR! CScalculator::CalculateCS() graphs not set, returning 0." << endl;
    return 0.;
  }

  Double_t cs_O16 = f_g_O16->Eval(E);
  Double_t cs_H1  = f_g_H1->Eval(E);

  return (2 * cs_H1 + cs_O16)/18.;

}

//*************************************************************************

/**
 * This function initialises the map fGraphs that contains TGraphs with xsec data.
 *
 * \param  xsecfile      File with xsec TGraphs.
 * \return               true if xsec file found, false otherwise.
 */
Bool_t CScalculator::InitGraphs(TString xsecfile) {

  TFile *f = new TFile(xsecfile,"READ");
  if ( !f->IsOpen() ) return false;

  //loop over maps, use the CreateString() function to create a lookup string
  //for each combination of nu_flavor, int_type, p_type. Clone the xsec graphs to the map.
  
  for (auto& nuflav: nu_flavs) {
    for (auto& inttype: int_types) {
      for (auto& ptype: p_types) {

	TString gname1 = "nu_" + nuflav.second + ptype.second + "_H1/tot" + inttype.second;
	TString gname2 = "nu_" + nuflav.second + ptype.second + "_O16/tot" + inttype.second;
	
	vector<TGraph*> graphs;
	graphs.push_back( (TGraph*)( f->Get(gname1)->Clone() ) );
	graphs.push_back( (TGraph*)( f->Get(gname2)->Clone() ) );
	
	TString namestr = CreateString( nuflav.first, (Bool_t)inttype.first, (Bool_t)ptype.first );
	
	fGraphs.insert( pair<TString, vector<TGraph*> > {namestr, graphs} );
	
      }
    }
  }

  f->Close();
  delete f;

  return true;

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
TString CScalculator::CreateString(Int_t nu_flavor, Bool_t is_cc, Bool_t is_nubar) {

  return nu_flavs[nu_flavor] + p_types[(Int_t)is_nubar] + int_types[(Int_t)is_cc];

}
