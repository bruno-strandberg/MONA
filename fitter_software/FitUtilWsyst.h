#ifndef FitUtilWsyst_h
#define FitUtilWsyst_h

#include "FitUtil.h"
#include "TH3.h"
#include "RooRealVar.h"

/** This class inherits from `FitUtil` and overwrites the virtual function `TrueEvts` to include systematic effects for the fitting process. All added flux systematic parameters preserve the overall normalisation.*/
class FitUtilWsyst : protected FitUtil {

  //-------------------------------------------------------------
  // public functions and members
  //-------------------------------------------------------------
  
public:

  FitUtilWsyst(Double_t op_time, TH3 *h_template,
	       Double_t emin, Double_t emax, Double_t ctmin, Double_t ctmax,
	       Double_t bymin, Double_t bymax, TString meff_file);

  /** Constructor. All `RooRealVar` parameters of this class are added to `FitUtil::fParSet` and deleted in that constructor */
  virtual ~FitUtilWsyst() {};

  virtual std::pair<Double_t, Double_t> TrueEvts(const TrueB &tb, const proxymap_t &proxymap);
  Double_t GetFluxWsyst(UInt_t flav, Bool_t isnb, Int_t true_ebin, Int_t true_ctbin, const proxymap_t& proxymap);

  //-------------------------------------------------------------
  // protected functions and members
  //-------------------------------------------------------------
  
protected:
  
  Double_t FluxTiltCoeff(Double_t energy, Double_t costheta, Double_t e_tilt, Double_t ct_tilt);  
  void     CalcTiltedFluxNorms(Double_t e_tilt, Double_t ct_tilt);  
  Double_t GetTiltedFlux(UInt_t flav, Bool_t isnb, Int_t true_ebin, Int_t true_ctbin,
			 Double_t e_tilt, Double_t ct_tilt);  

  // flux systematics parameters
  RooRealVar* fE_tilt;        //!< parameter for atm. flux E tilt by multiplying with \f$energy^{E_{tilt}}$\f; 0 means no tilt
  RooRealVar* fCt_tilt;       //!< parameter for atm. flux ct tilt multiplying with \f$(1+ct_{tilt} * ct)\f$; 0 means no tilt
  RooRealVar* fSkew_mu_amu;   //!< parameter to skew muon to anti-muon flux, preserves mu+amu norm; 1 means no skew
  RooRealVar* fSkew_e_ae;     //!< parameter to skew elec to anti-elec flux, preserves e+ae norm; 1 means no skew
  RooRealVar* fSkew_mu_e;     //!< parameter to skew mu to e flux, preserves mu+amu+e+ae norm; 1 means no skew
  
  // cache for flux tilt
  Double_t f_cache_e_tilt;    //!< internal cache for `fE_tilt` value to reduce calls for normalisation calculation
  Double_t f_cache_ct_tilt;   //!< internal cache for `fCt_tilt` value to reduce calls for normalisation calculation
  Double_t fTiltedFluxNorms[fFlavs][fPols]; //!< cache for normalisation of atm. flux shape systematics (energy, ct tilt)
   
};

#endif
