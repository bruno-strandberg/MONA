#ifndef FitUtil_h
#define FitUtil_h

// NMH and OscProb
#include "DetResponse.h"
#include "AtmFlux.h"
#include "NuXsec.h"
#include "EffMass.h"
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

/// type definition for passing arguments between RooFit wrapper class `FitPDF` and worker-functions in `FitUtil`.
typedef std::map<TString, RooRealProxy*> proxymap_t;

/** This class is used in conjunction with `FitPDF` to fit NMO data with `RooFit`.

    The class hosts the fit parameters (e.g. the 6 oscillation parameters) and uses elements of `common_software/` and `OscProb` to provide functions that predict the number of expected events in a true and reco (E, cos-theta, bjorken-y) bin. Additionally, it has functions that are to be called inside `FitPDF` class - together, `FitUtil` and `FitPDF` enable the use of `RooFit` for fitting NMO data. See `fitter_software/README.md` for more info.

    Each `FitPDF` class will require a pointer to a `FitUtil` instance. One `FitUtil` instance can (and should!) be shared between several `FitPDF` instances, if several `FitPDF` instances are required (usually yes - at least one for tracks and one for showers). `FitUtil` employs an internal cache for oscillation probabilities, which speeds up the fitting procedure when several histograms/data sets are fitted in parallel (e.g. when several PID bins are used).

    The `FitPDF` class that uses `FitUtil` is a modification of a class that was auto-generated with `RooClassFactory::makePdf`. The idea is that `FitPDF` acts merely as a wrapper class to get access to `RooFit` niceties (such as simultaneous fitting), the model calculations and parameters are defined here in `FitUtil` class.

    <B> How to extend? </B> Simply put, the `FitPDF::evaluate` method is called each time during the fitting (`pdf->fitTo`) when the fitter moves to a new \f$ (E_{\rm reco}, cos\theta_{\rm reco}, by_{\rm reco}) \f$ bin. `FitPDF::evaluate` calls `PdfEvaluate` of this class, which returns the expected event density in the reco bin. Let's say one now wishes to add another fit parameter that scales the atmospheric flux. To do that, one must: 

    -# create a new `RooRealVar*` member in this class, initialise it <B> and add it to `fParSet` </B>. Basically, one must do exactly what is done for `fDm31`. `FitPDF` class creates a `RooRealProxy` for each parameter in `fParSet`, such that the map given to `FitUtil::PdfEvaluate` contains an entry for each parameter `fParSet`. The nice thing in such a setup is that `FitPDF` does not need to be modified, parameters are added in this class and used in this class.
    -# extend the argument list of `FitUtil::RecoEvts` by the new parameter, extract the new parameter from the map in `FitUtil::PdfEvaluate` and give it to `FitUtil::RecoEvts`.
    -# Do something meaningful with the new parameter. In this case, the flux scale needs to be passed inside `FitUtil::RecoEvts` to `FitUtil::TrueEvts`, where the flux scale can multiply `atm_count_e` and `atm_count_m`.

*/
class FitUtil {

  //*********************************************************************************
  //public members and function of the `FitUtil` class
  //*********************************************************************************
  
 public:

  //------------------------------------------------------------------
  // constructors/destructors
  //------------------------------------------------------------------

  /** Default constructor. */
  FitUtil() {};

  FitUtil(Double_t op_time, TH3 *h_template,
	  Double_t emin, Double_t emax, Double_t ctmin, Double_t ctmax, Double_t bymin, Double_t bymax,
	  TString meff_file);
  ~FitUtil();

  //------------------------------------------------------------------
  // public functions that are called in `FitPDF`
  //------------------------------------------------------------------  
  std::pair<Double_t, Double_t> RecoEvts(Double_t E_reco, Double_t Ct_reco, Double_t By_reco, DetResponse *resp, const proxymap_t &proxymap);
  TH3D* Expectation(DetResponse *resp, const proxymap_t &proxymap, const char* rangeName);

  //------------------------------------------------------------------
  // other public functions
  //------------------------------------------------------------------  
  std::pair<Double_t, Double_t> TrueEvts(const TrueB &tb, const proxymap_t &proxymap);
  Double_t GetCachedFlux(UInt_t flav, Bool_t isnb, const TrueB &tb);
  Double_t GetCachedOsc(UInt_t flav_in, const TrueB &tb, const proxymap_t& proxymap);
  Double_t GetCachedXsec(const TrueB &tb);
  void SetNOlims();
  void SetNOcentvals();
  void SetIOlims();
  void SetIOcentvals();
  void FreeParLims();
  
  //------------------------------------------------------------------
  // setters/getters
  //------------------------------------------------------------------
  
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

  /** Get pointer to the member `RooRealVar` for the reco energy observable
      \return A pointer to the reco energy observable
   */
  RooRealVar* GetEobs()  { return fE_reco; }

  /** Get pointer to the member `RooRealVar` for the reco cos-theta observable
      \return A pointer to the reco cos-theta observable
   */
  RooRealVar* GetCTobs() { return fCt_reco; }

  /** Get pointer to the member `RooRealVar` for the reco bjorken-y observable
      \return A pointer to the reco bjorken-y observable
   */
  RooRealVar* GetBYobs() { return fBy_reco; }

  /** Get a pointer to the member histogram with binning configuration
      \return Pointer to a `TH3` with binning configuration
   */
  TH3D*       GetBinningHist() { return fHB; }

  /** Get the a pointer to `AtmFlux` member of `FitUtil` that is used to fill caches.
      \return pointer to `AtmFlux` instance.
   */
  AtmFlux* GetFluxCalculator() { return fFlux; }
  /** Get the a pointer to `NuXsec` member of `FitUtil` that is used to fill caches.
      \return pointer to `NuXsec` instance.
   */
  NuXsec*  GetXsecCalculator() { return fXsec; }
  /** Get the a pointer to `OscProb::PMNS_Fast` member of `FitUtil` that is used to fill caches.
      \return pointer to `OscProb::PMNS_Fast` instance.
   */
  OscProb::PMNS_Fast* GetOscCalculator() { return fProb; }
  /** Get the a pointer to `OscProb::PremModel` member of `FitUtil` that is used to fill caches.
      \return pointer to `OscProb::PremModel` instance.
   */
  OscProb::PremModel* GetEarthModel() { return fPrem; }

  RooRealVar* GetVar(TString varname);

  //*********************************************************************************
  //protected members and function of the `FitUtil` class, accessible to all inheriting classes
  //*********************************************************************************
  
 protected:

  //------------------------------------------------------------------
  // protected constants for detected neutrino count calculation
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
  // protected members for detected neutrino count calculation
  //------------------------------------------------------------------

  Double_t            fOpTime;                      //!< operation time in years
  TH3D               *fHB;                          //!< a template histogram that defines the binning
  AtmFlux            *fFlux;                        //!< atm flux calculator
  NuXsec             *fXsec;                        //!< xsec calculator
  EffMass            *fMeff;                        //!< effective mass calculator
  OscProb::PMNS_Fast *fProb;                        //!< oscillation probability calculator
  OscProb::PremModel *fPrem;                        //!< earth model

  //------------------------------------------------------------------
  // protected members for `RooFit` observable and parameter access
  //------------------------------------------------------------------

  RooArgSet   fParSet;   //!< set that includes all variables (observables+parameters)
  RooArgList  fObsList;  //!< list (ordered!) that includes only observables (Energy, cos-theta, bjorken-y)

  //------------------------------------------------------------------
  // protected internal structure to initialise private osc parameter limits below
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
  
  //*********************************************************************************
  //private members and function of the `FitUtil` class, accessible only inside this class
  //*********************************************************************************
  
 private:

  //------------------------------------------------------------------
  // private functions
  //------------------------------------------------------------------
  void ProbCacher(const proxymap_t& proxymap);  
  void InitFitVars(Double_t emin, Double_t emax, Double_t ctmin, Double_t ctmax, Double_t bymin, Double_t bymax);
  void InitCacheHists(TH3D *h_template);
  void FillFluxAndXsecCache(AtmFlux *flux, NuXsec *xsec, Double_t op_time);
  std::tuple<Double_t, Double_t, Int_t, Int_t> GetRange(Double_t min, Double_t max, TAxis *axis);
  enum rangeret { MIN=0, MAX, MINBIN, MAXBIN };           //!< enum for function `GetRange` return
    
  //------------------------------------------------------------------
  // central values and limits for oscillation parameters, used throughout the class
  //------------------------------------------------------------------
  
  const fpar f_NO_sinsqth12 = fpar(0.297, 0.25, 0.354);    //!< central value and limits for \f$ sin^2\theta_{12} \f$, normal ordering
  const fpar f_IO_sinsqth12 = fpar(0.297, 0.25, 0.354);    //!< central value and limits for \f$ sin^2\theta_{12} \f$, inverse ordering

  const fpar f_NO_sinsqth13 = fpar(0.0215, 0.019, 0.024);  //!< central value and limits for \f$ sin^2\theta_{13} \f$, normal ordering
  const fpar f_IO_sinsqth13 = fpar(0.0216, 0.019, 0.0242); //!< central value and limits for \f$ sin^2\theta_{13} \f$, inverse ordering

  const fpar f_NO_sinsqth23 = fpar(0.425, 0.381, 0.615);   //!< central value and limits for \f$ sin^2\theta_{23} \f$, normal ordering
  const fpar f_IO_sinsqth23 = fpar(0.589, 0.384, 0.636);   //!< central value and limits for \f$ sin^2\theta_{23} \f$, inverse ordering

  const fpar f_NO_dcp = fpar(1.38, -0.5, 2.5 ); //!< central value and limits for \f$ \delta_{CP} \f$, normal ordering; the fit parameter range needs to be wider than 2pi, otherwise some value will always coincide with parameter limit during minimization
  const fpar f_IO_dcp = fpar(1.31, f_NO_dcp.min, f_NO_dcp.max); //!< central value and limits for \f$ \delta_{CP} \f$, inverse ordering; the fit parameter range needs to be wider than 2pi, otherwise some value will always coincide with parameter limit during minimization

  const fpar f_NO_dm21 = fpar(7.37e-5, 6.93e-5, 7.96e-5);  //!< central value and limits for \f$ \Delta m_{21}^2 \f$, normal ordering
  const fpar f_IO_dm21 = fpar(7.37e-5, 6.93e-5, 7.96e-5);  //!< central value and limits for \f$ \Delta m_{21}^2 \f$, inverse ordering

  const fpar f_NO_dm31 = fpar(2.56e-3, 2.45e-3, 2.69e-3);  //!< central value and limits for \f$ \Delta m_{31}^2 \f$, normal ordering
  /// central value and limits for \f$ \Delta m_{31}^2 \f$, inverse ordering
  const fpar f_IO_dm31 = fpar(-2.54e-3+f_IO_dm21.cv, -2.66e-3+f_IO_dm21.cv, -2.42e-3+f_IO_dm21.cv);

  const fpar f_free_dm21 = fpar(7.37e-5,  5e-5, 1e-4); //!< free limits for \f$ \Delta m_{21}^2 \f$
  const fpar f_free_dm31 = fpar(2.56e-3, -5e-3, 5e-3); //!< free limits for \f$ \Delta m_{31}^2 \f$
    
  //------------------------------------------------------------------
  // private members for caching
  //------------------------------------------------------------------

  TH2D *fhFluxCache[fFlavs][fPols];         //!< atm flux cache with struct. [flav][is_nub]
  TH2D *fhOscCache [fFlavs][fFlavs][fPols]; //!< osc prob cache with struct. [flav_in][flav_out][isnb]
  TH1D *fhXsecCache[fFlavs][fInts][fPols];  //!< xsec cache with struct. [flav][is_cc][isnb]

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
  // private members for monitoring performance, used in `FitUtil::ProbCacher`
  //------------------------------------------------------------------
  Long64_t    fOscCalls;         //!< counts the number of oscillator calls
  TStopwatch *fOscCalcTime;      //!< counts the accumulated time for oscillation calculation

  //------------------------------------------------------------------
  // variables that need to be recognized by the RooFit minimizer and define
  // observables (E, ct, by) and fit parameters (oscillation parameters)
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

};

#endif
