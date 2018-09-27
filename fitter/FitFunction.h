#ifndef FitFunction_h
#define FitFunction_h

#include "AtmFlux.h"
#include "NuXsec.h"
#include "DetResponse.h"

#include "PMNS_Fast.h"
#include "PremModel.h"

#include "TStopwatch.h" //debug/dev

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
  double operator() (double *x, double *p);

  Double_t TrueDetected (Double_t e_true, Double_t ct_true, Double_t by_true, 
			 Int_t flav, Int_t iscc, Int_t isnb, double *p);

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

  // timers for development
  TStopwatch *fPathCalc;
  TStopwatch *fOscCalc;
  TStopwatch *fAtmCalc;
  TStopwatch *fRestCalc;
  Long64_t    fOscCalls;

};

#endif
