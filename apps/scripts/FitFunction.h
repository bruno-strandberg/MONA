#ifndef FitFunction_h
#define FitFunction_h

#include "AtmFlux.h"
#include "NuXsec.h"
#include "DetResponse.h"

#include "PMNS_Fast.h"
#include "PremModel.h"

#include "TStopwatch.h" //debug/dev
#include "TH2.h"
#include "TH3.h"

/**
   This is a first attempt to implement a fit function class and will need severe re-desing/re-writing.
 */
class FitFunction {

 public:
  FitFunction(DetResponse *DR, Double_t op_time,
	      TString meff_elec_cc, TString meff_muon_cc, 
	      TString meff_tau_cc, TString meff_elec_nc);
  ~FitFunction();

  void ReadMeffHists(TString meff_elec_cc, TString meff_muon_cc, 
		     TString meff_tau_cc, TString meff_elec_nc);


  // this operator is implemented, such that `FitFunction` object can be used with TF1/TF2/TF3
  std::tuple<double, double> operator() (double *x, double *p);
  std::tuple<double, double> ResponseCount(double *x, double *p);
  std::vector<Double_t> Weights (double *x);

  Double_t TrueDetected (Int_t ebin_true, Int_t ctbin_true, Int_t bybin_true, 
			 Int_t flav, Int_t iscc, Int_t isnb, double *p);

  void     ProbCacher(Double_t th12, Double_t th13, Double_t th23, Double_t dcp, Double_t dm21, Double_t dm31);
  Double_t GetCachedProb(Int_t isnb, Int_t flav_in, Int_t flav_out, Int_t ebin, Int_t ctbin);
  Double_t GetCachedFlux(Int_t flav, Int_t isnb, Int_t ebin, Int_t ctbin);
  Double_t GetCachedXsec(Int_t flav, Int_t iscc, Int_t isnb, Int_t ebin);
  void     InitCacheHists(TH3D *h_template, AtmFlux *flux, NuXsec *xsec, Double_t op_time);

 private:

  TH3D               *fHB;             //!< pointer to the histogram from DetResponse for binning info
  AtmFlux            *fFlux;           //!< atm flux calculator
  NuXsec             *fXsec;           //!< xsec calculator
  OscProb::PMNS_Fast *fProb;           //!< oscillation probability calculator
  OscProb::PremModel *fPrem;           //!< earth model
  DetResponse        *fResponse;       //!< detector response
  TH3D               *fhMeff[3][2][2]; //!< effective mass histograms with structure [flav][is_cc][is_nub]

  Double_t fOpTime; //!< operation time

  const double fMp = 1.672621898e-27;                     //!< proton mass in kg
  const double fMn = 1.674927471e-27;                     //!< neutron mass in kg
  const double fMN = (fMp+fMn)/2;                         //!< nucleon-average mass
  const double fSec_per_y   = 365.2421897 * 24 * 60 * 60; //!< seconds in a tropical year
  const double fKg_per_Mton = 1e9;                        //!< kg per MTon (MTon = 1e6 Ton; Ton = 1e3 kg)

  TStopwatch *fPathCalc; //!< path calculation timer
  TStopwatch *fOscCalc;  //!< oscillation calculation timer
  TStopwatch *fAtmCalc;  //!< atmospheric flux calculation timer
  TStopwatch *fRestCalc; //!< calculation of ther steps timer
  Long64_t    fOscCalls; //!< counter to the number of oscillator calls

  // cached oscillation probabilities
  TH2D *fhOscProb[2][2][3]; //!< cached osc probs with structure [isnb][elec/muon][flavor]
  Double_t f_cache_th12;    //!< cached theta12 valeu
  Double_t f_cache_th13;    //!< cached theta13 value
  Double_t f_cache_th23;    //!< cached theta23 value
  Double_t f_cache_dcp ;    //!< cached delta-cp value
  Double_t f_cache_dm21;    //!< cached dm21 value
  Double_t f_cache_dm31;    //!< cached dm31 value

  // atm flux and xsec cache
  TH2D *fhAtmFluxCache[2][2]; //!< cached fluxes with structure [elec/muon][isnb]
  TH1D *fhXsecCache[3][2][2]; //!< cached xsec with structure [flav][iscc][isnb]

};

#endif
