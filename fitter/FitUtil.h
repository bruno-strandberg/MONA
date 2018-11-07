#ifndef FitUtil_h
#define FitUtil_h

// NMH and OscProb
#include "DetResponse.h"
#include "AtmFlux.h"
#include "NuXsec.h"
#include "PMNS_Fast.h"
#include "PremModel.h"

// Root
#include "TH3.h"
#include "TStopwatch.h"

// RooFit
#include "RooRealVar.h"
#include "RooRealProxy.h"
#include "RooArgSet.h"
#include "RooArgList.h"

// standard cpp
#include <map>

class FitUtil {

 public:

  // constructors/destructors

  /** Default constructor. */
  FitUtil() {};

  FitUtil(Double_t op_time, TH3 *h_template,
	  Double_t emin, Double_t emax, Double_t ctmin, Double_t ctmax, Double_t bymin, Double_t bymax,
	  TString meffh_elec_cc, TString meffh_muon_cc, TString meffh_tau_cc, TString meffh_elec_nc);
  ~FitUtil();

  //------------------------------------------------------------------
  // public functions that are called in `FitPDF`
  //------------------------------------------------------------------
  Double_t PdfEvaluate (const std::map<TString, RooRealProxy*> &parmap, DetResponse *resp);
  std::pair<TH3D*, Double_t> PdfExpectation(const std::map<TString, RooRealProxy*> &parmap, DetResponse *resp, const char* rangeName);

  // setters/getters
  RooArgSet   GetSet()         { return fParSet;  }
  RooArgList  GetObs()         { return fObsList; }
  TH3D*       GetBinningHist() { return fHB;      }
  void        SetNOlims();
  void        SetIOlims();
  
  // DEV: this should become private; kept public at the moment for comparisons with ROOT
  std::pair<Double_t, Double_t> RecoEvts(DetResponse *resp,
					 Double_t E_reco, Double_t Ct_reco, Double_t By_reco,
					 Double_t SinsqTh12, Double_t SinsqTh13, Double_t SinsqTh23,
					 Double_t Dcp, Double_t Dm21, Double_t Dm31);
 private:

  //------------------------------------------------------------------
  // private functions
  //------------------------------------------------------------------

  std::pair<Double_t, Double_t> TrueEvts(Int_t ebin_true, Int_t ctbin_true, Int_t bybin_true,
					 UInt_t flav, UInt_t iscc, UInt_t isnb,
					 Double_t SinsqTh12, Double_t SinsqTh13, Double_t SinsqTh23, 
					 Double_t Dcp, Double_t Dm21, Double_t Dm31);

  void ProbCacher(Double_t SinsqTh12, Double_t SinsqTh13, Double_t SinsqTh23, 
		  Double_t Dcp, Double_t Dm21, Double_t Dm31);
  
  void InitFitVars(Double_t emin, Double_t emax, Double_t ctmin, Double_t ctmax, Double_t bymin, Double_t bymax);
  void InitCacheHists(TH3D *h_template);
  void FillFluxAndXsecCache(AtmFlux *flux, NuXsec *xsec, Double_t op_time);
  void ReadMeffHists(TH3D *h_template, TString meffh_elec_cc, TString meffh_muon_cc, 
		     TString meffh_tau_cc, TString meffh_elec_nc);
  std::tuple<Double_t, Double_t, Int_t, Int_t> GetRange(Double_t min, Double_t max, TAxis *axis);
  enum rangeret { MIN=0, MAX, MINBIN, MAXBIN };           //!< enum for function `GetRange` return
    
  //------------------------------------------------------------------
  // constants for detected neutrino count calculation
  //------------------------------------------------------------------

  const double fMp = 1.672621898e-27;                     //!< proton mass in kg
  const double fMn = 1.674927471e-27;                     //!< neutron mass in kg
  const double fMN = (fMp+fMn)/2;                         //!< nucleon-average mass
  const double fSec_per_y   = 365.2421897 * 24 * 60 * 60; //!< seconds in a tropical year
  const double fKg_per_MTon = 1e9;                        //!< kg per MTon (MTon = 1e6 Ton; Ton = 1e3 kg)

  static const UInt_t fFlavs = 3;                         //!< number of flavors
  static const UInt_t fInts  = 2;                         //!< number of interactions (CC/NC)
  static const UInt_t fPols  = 2;                         //!< number of polarisations (nu/nub)

  enum flavors {ELEC = 0, MUON, TAU};                     //!< enum for flavors

  //------------------------------------------------------------------
  // private members for detected neutrino count calculation
  //------------------------------------------------------------------

  Double_t            fOpTime;                      //!< operation time
  TH3D               *fHB;                          //!< a template histogram that defines the binning
  AtmFlux            *fFlux;                        //!< atm flux calculator
  NuXsec             *fXsec;                        //!< xsec calculator
  OscProb::PMNS_Fast *fProb;                        //!< oscillation probability calculator
  OscProb::PremModel *fPrem;                        //!< earth model

  //------------------------------------------------------------------
  // private members for caching
  //------------------------------------------------------------------

  TH2D *fhFluxCache[fFlavs][fPols];         //!< atm flux cache with struct. [flav][is_nub]
  TH2D *fhOscCache [fFlavs][fFlavs][fPols]; //!< osc prob cache with struct. [flav_in][flav_out][isnb]
  TH1D *fhXsecCache[fFlavs][fInts][fPols];  //!< xsec cache with struct. [flav][is_cc][isnb]
  TH3D *fhMeff     [fFlavs][fInts][fPols];  //!< effective mass hists with struct. [flav][is_cc][is_nub]

  Double_t f_cache_sinsqth12;    //!< cached theta12 value
  Double_t f_cache_sinsqth13;    //!< cached theta13 value
  Double_t f_cache_sinsqth23;    //!< cached theta23 value
  Double_t f_cache_dcp ;         //!< cached delta-cp value
  Double_t f_cache_dm21;         //!< cached dm21 value
  Double_t f_cache_dm31;         //!< cached dm31 value

  //------------------------------------------------------------------
  // private members for determining the integration range and oscillation calculation range
  //------------------------------------------------------------------
  Int_t fEbin_min;               //!< minimum energy bin number included
  Int_t fEbin_max;               //!< maximum energy bin number included
  Int_t fCtbin_min;              //!< minimum cos-theta bin number included
  Int_t fCtbin_max;              //!< minimum cos-theta bin number included
  Int_t fBybin_min;              //!< minimum bjorken-y bin number included
  Int_t fBybin_max;              //!< minimum bjorken-y bin number included
  
  //------------------------------------------------------------------
  // private members for monitoring performance
  //------------------------------------------------------------------
  Long64_t    fOscCalls;
  TStopwatch *fOscCalcTime;

  //------------------------------------------------------------------
  // variables that need to be recognized by the RooFit minimizer
  //------------------------------------------------------------------

  RooRealVar *fE_reco;    //!< reconstructed energy observable
  RooRealVar *fCt_reco;   //!< reconstructed cos-theta observable
  RooRealVar *fBy_reco;   //!< reconstructed bjorken-y observable
  RooRealVar *fSinsqTh12; //!< sin^2 theta-12 parameter
  RooRealVar *fSinsqTh13; //!< sin^2 theta-13 parameter
  RooRealVar *fSinsqTh23; //!< sin^2 theta-23 parameter
  RooRealVar *fDcp;      //!< delta-cp in \f$\pi$\f's, as given by the PDG group (e.g. 1.38) parameter
  RooRealVar *fDm21;     //!< small mass-splitting-square parameter in eV^2
  RooRealVar *fDm31;     //!< large mass-splitting-square parameter in eV^2

  RooArgSet   fParSet;   //!< set that includes all variables (observables+parameters)
  RooArgList  fObsList;  //!< list (ordered!) that includes only observables (Energy, cos-theta, bjorken-y)

};

#endif
