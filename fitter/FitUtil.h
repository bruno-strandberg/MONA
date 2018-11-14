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
#include "TMath.h"

// RooFit
#include "RooRealVar.h"
#include "RooRealProxy.h"
#include "RooArgSet.h"
#include "RooArgList.h"

// standard cpp
#include <map>

/** This class is used in conjunction with `FitPDF` to fit NMO data with `RooFit` 

    Main ideas: internal oscillation cache, such that i.e. with 10 histograms fitted simultaneously, the number of oscillator calls in minimized.
    Main ideas: expansion of the fit model is performed here, no modifications of the `FitPDF` class is required.
*/
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

  /** Get the `RooArgSet` with all parameters known to `RooFit`
      \return `RooArgSet` with known parameters.
   */
  RooArgSet   GetSet()         { return fParSet;  }
  /** Get the `RooArgList` with observables (Energy, cos-theta, bjorken-y)
      \return `RooArgList` with observables (Energy, cos-theta, bjorken-y)
   */
  RooArgList  GetObs()         { return fObsList; }
  /** Get the 3D histogram that stores the binning information
      \return `TH3D` with the binning used in the analysis.
   */
  TH3D*       GetBinningHist() { return fHB;      }
  void        SetNOlims();
  void        SetNOcentvals();
  void        SetIOlims();
  void        SetIOcentvals();
  
  // this function should be private, but is kept public for comparisons with ROOT
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
  // central values and limits for oscillation parameters, used throughout the class
  //------------------------------------------------------------------

  /** internal structure to store parameter central values and limits*/
  struct fpar {
    Double_t cv;   //!< central value
    Double_t min;  //!< minimum
    Double_t max;  //!< max

    /** constructor 
	\param _cv  central value
	\param _min minimum
	\param _max maximum
    */
    fpar(Double_t _cv, Double_t _min, Double_t _max) {

      if (_cv < _min) {
	throw std::invalid_argument( "ERROR! FitUtil::fpar::fpar() central value " + to_string(_cv) + " is smaller than the minimum " + to_string(_min) );
      }

      if (_cv > _max) {
	throw std::invalid_argument( "ERROR! FitUtil::fpar::fpar() central value " + to_string(_cv) + " is larger than the maximum " + to_string(_max) );
      }
      
      cv  = _cv;
      min = _min;
      max = _max;
    }
    
  };
  
  const fpar f_NO_sinsqth12 = fpar(0.297, 0.25, 0.354);    //!< central value and limits for \f$ sin^2\theta_{12} \f$, normal ordering
  const fpar f_IO_sinsqth12 = fpar(0.297, 0.25, 0.354);    //!< central value and limits for \f$ sin^2\theta_{12} \f$, inverse ordering

  const fpar f_NO_sinsqth13 = fpar(0.0215, 0.019, 0.024);  //!< central value and limits for \f$ sin^2\theta_{13} \f$, normal ordering
  const fpar f_IO_sinsqth13 = fpar(0.0216, 0.019, 0.0242); //!< central value and limits for \f$ sin^2\theta_{13} \f$, inverse ordering

  const fpar f_NO_sinsqth23 = fpar(0.425, 0.381, 0.615);   //!< central value and limits for \f$ sin^2\theta_{23} \f$, normal ordering
  const fpar f_IO_sinsqth23 = fpar(0.589, 0.384, 0.636);   //!< central value and limits for \f$ sin^2\theta_{23} \f$, inverse ordering

  const fpar f_NO_dcp = fpar(1.38, 0., TMath::Pi() );      //!< central value and limits for \f$ \delta_{CP} \f$, normal ordering
  const fpar f_IO_dcp = fpar(1.31, 0., TMath::Pi() );      //!< central value and limits for \f$ \delta_{CP} \f$, inverse ordering

  const fpar f_NO_dm21 = fpar(7.37e-5, 6.93e-5, 7.96e-5);  //!< central value and limits for \f$ \Delta m_{21}^2 \f$, normal ordering
  const fpar f_IO_dm21 = fpar(7.37e-5, 6.93e-5, 7.96e-5);  //!< central value and limits for \f$ \Delta m_{21}^2 \f$, inverse ordering

  const fpar f_NO_dm31 = fpar(2.56e-3, 2.45e-3, 2.69e-3);  //!< central value and limits for \f$ \Delta m_{31}^2 \f$, normal ordering
  /// central value and limits for \f$ \Delta m_{31}^2 \f$, inverse ordering
  const fpar f_IO_dm31 = fpar(-2.54e-3+f_IO_dm21.cv, -2.66e-3+f_IO_dm21.cv, -2.42e-3+f_IO_dm21.cv);
    
  //------------------------------------------------------------------
  // private members for detected neutrino count calculation
  //------------------------------------------------------------------

  Double_t            fOpTime;                      //!< operation time in years
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

  Double_t f_cache_sinsqth12;    //!< internally cached theta12 value
  Double_t f_cache_sinsqth13;    //!< internally cached theta13 value
  Double_t f_cache_sinsqth23;    //!< internally cached theta23 value
  Double_t f_cache_dcp ;         //!< internally cached delta-cp value
  Double_t f_cache_dm21;         //!< internally cached dm21 value
  Double_t f_cache_dm31;         //!< internally cached dm31 value

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
  Long64_t    fOscCalls;         //!< counts the number of oscillator calls
  TStopwatch *fOscCalcTime;      //!< counts the accumulated time for oscillation calculation

  //------------------------------------------------------------------
  // variables that need to be recognized by the RooFit minimizer
  //------------------------------------------------------------------

  RooRealVar *fE_reco;    //!< reconstructed energy observable in GeV
  RooRealVar *fCt_reco;   //!< reconstructed cos-theta observable
  RooRealVar *fBy_reco;   //!< reconstructed bjorken-y observable
  RooRealVar *fSinsqTh12; //!< \f$ sin^2\theta_{12} \f$ fit parameter
  RooRealVar *fSinsqTh13; //!< \f$ sin^2\theta_{13} \f$ fit parameter
  RooRealVar *fSinsqTh23; //!< \f$ sin^2\theta_{23} \f$ fit parameter
  RooRealVar *fDcp;       //!< \f$ \delta_{CP} \f$ fit parameter in \f$ \pi $\f's, as given by the PDG group (e.g. 1.38)
  RooRealVar *fDm21;      //!< \f$ \Delta m_{21}^2 \f$ fit parameter in eV^2
  RooRealVar *fDm31;      //!< \f$ \Delta m_{31}^2 \f$ fit parameter in eV^2

  RooArgSet   fParSet;   //!< set that includes all variables (observables+parameters)
  RooArgList  fObsList;  //!< list (ordered!) that includes only observables (Energy, cos-theta, bjorken-y)

};

#endif
