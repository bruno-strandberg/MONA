#ifndef AtmFlux_h
#define AtmFlux_h

#include<map>
#include<iostream>
#include "TGraph2D.h"
using namespace std;

/** Options for Atmospheric flux class.
 *
 *  h - honda; frj - frejus site; grn - gran sasso site; mtn - with mountain;
 *  min/max - solar minimum/maximum.
 */
enum AtmFluxOpt {h_frj_min, h_frj_max, h_frj_mtn_min, h_frj_mtn_max,
		 h_grn_min, h_grn_max, h_grn_mtn_min, h_grn_mtn_max};

/** A class that reads in atm neutrino flux data from files and provides a flux calculator. 
 *
 *  To use, init AtmFlux class and use GetHondaFlux function to get the atmospheric neutrino flux
 *  for a specified neutrino at a certain energy and angle. Different flux options for the constructor
 *  are listed in AtmFluxOpt. E.g. AtmFlux f(AtmFluxOpt::h_frj_min).
 *
 *  NB! There is a small approximation involved, see documentation for function ReadHondaFlux.
 * 
 */
class AtmFlux {

 public:
  AtmFlux(UInt_t opt = AtmFluxOpt::h_frj_min, Bool_t debug=false);
  ~AtmFlux();
  Double_t GetHondaFlux(UInt_t nu_flavor, Bool_t is_nubar, Double_t E, Double_t cosz);
  
 private:
  Bool_t   ReadHondaFlux(TString fname, Bool_t debug=false);
  
  map <UInt_t, TString> fFileMap;      //!< map linking options and flux file names
  static const int fNflavs = 2;        //!< number of flavors (e, mu, no tau)
  static const int fPtypes = 2;        //!< number of types (nu, nubar)
  TGraph2D* fGraphs[fNflavs][fPtypes]; //!< 2D graphs that hold the flux data
  
};

#endif
