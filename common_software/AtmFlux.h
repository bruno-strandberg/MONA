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
 *  To use, init AtmFlux class and use Flux_dE_dcosz function to get the differential atmospheric 
 *  neutrino flux for a specified neutrino at a certain energy and angle. Different flux options for
 *  the constructor are listed in AtmFluxOpt. E.g. AtmFlux f(AtmFluxOpt::h_frj_min).
 *
 *  NB! There is a small approximation involved, see documentation for function ReadHondaFlux.
 * 
 */
class AtmFlux {

 public:
  AtmFlux(UInt_t opt = AtmFluxOpt::h_frj_min, Bool_t debug=false);
  ~AtmFlux();
  Double_t Flux_dE_dcosz(UInt_t nu_flavor, Bool_t is_nubar, Double_t E, Double_t cosz);
  
 private:
  void ReadHondaFlux(TString fname, Bool_t debug=false);
  
  map <UInt_t, TString> fFileMap;       //!< map linking options and flux file names
  static const int fNflavs = 2;         //!< number of flavors (e, mu, no tau)
  static const int fPtypes = 2;         //!< number of types (nu, nubar)
  TGraph2D* fGraphs[fNflavs][fPtypes];  //!< 2D graphs that hold the flux data
  Double_t fEmin;                       //!< minimum energy range of the interpolation
  Double_t fEmax;                       //!< maximum energy range of the interpolation
  Double_t fCtmin;                      //!< minimum cos-theta range of the interpolation
  Double_t fCtmax;                      //!< maximum cos-theta range of the interpolation
  static constexpr double fCt_bw = 0.1; //!< bin width used in input honda tables
  static const Int_t fPower = 3;        //!< E power for more accurate interpolation
};

#endif
