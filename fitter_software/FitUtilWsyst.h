#ifndef FitUtilWsyst_h
#define FitUtilWsyst_h

#include "FitUtil.h"
#include "TH3.h"
#include "RooRealVar.h"

/** This class inherits from `FitUtil` and overwrites the virtual function `TrueEvts` and `RecoEvts` to include systematic effects for the fitting process.*/
class FitUtilWsyst : public FitUtil {

  //**************************************************************************
  // public functions and members
  //**************************************************************************
  
public:

  FitUtilWsyst(Double_t op_time, TH3 *h_temp_T, TH3 *h_temp_R,
	       Double_t emin, Double_t emax, Double_t ctmin, Double_t ctmax,
	       Double_t bymin, Double_t bymax, TString meff_file);

  FitUtilWsyst(Double_t op_time, TH3 *h_template,
	       Double_t emin, Double_t emax, Double_t ctmin, Double_t ctmax,
	       Double_t bymin, Double_t bymax, TString meff_file);

  /** Destructor */
  virtual ~FitUtilWsyst() {};

  virtual std::pair<Double_t, Double_t> TrueEvts(const TrueB &tb, const proxymap_t &proxymap);
  Double_t GetFluxWsyst(UInt_t flav, Bool_t isnb, Int_t true_ebin, Int_t true_ctbin, const proxymap_t& proxymap);
  
  //**************************************************************************
  // protected functions and members
  //**************************************************************************
  
protected:

  //--------------------------------------------------------------------
  // protected functions
  //--------------------------------------------------------------------

  virtual std::pair< Double_t, Double_t > RecoEvtsER(Double_t E_reco, Double_t Ct_reco, Double_t By_reco,
						     EvtResponse *resp, const proxymap_t &proxymap,
						     Bool_t AddMuonsNoise);  
  Double_t FluxTiltCoeff(Double_t energy, Double_t costheta, Double_t e_tilt, Double_t ct_tilt);  
  void     CalcTiltedFluxNorms(Double_t e_tilt, Double_t ct_tilt);  
  Double_t GetTiltedFlux(UInt_t flav, Bool_t isnb, Int_t true_ebin, Int_t true_ctbin,
			 Double_t e_tilt, Double_t ct_tilt);
  vector< std::pair<Int_t, Double_t> > GetBinFractions(Double_t lo, Double_t hi, TAxis* axis);

  //--------------------------------------------------------------------
  // protected new systematic parameters
  //--------------------------------------------------------------------
  
  RooRealVar* fE_tilt;      //!< parameter for atm. flux E tilt by multiplying with \f$E^{E_{tilt}}\f$; 0 means no tilt
  RooRealVar* fCt_tilt;     //!< parameter for atm. flux ct tilt multiplying with \f$(1+ct_{tilt} * ct)\f$; 0 means no tilt
  RooRealVar* fSkew_mu_amu; //!< parameter to skew muon to anti-muon flux, preserves mu+amu norm; 0 means no skew
  RooRealVar* fSkew_e_ae;   //!< parameter to skew elec to anti-elec flux, preserves e+ae norm; 0 means no skew
  RooRealVar* fSkew_mu_e;   //!< parameter to skew mu to e flux, preserves mu+amu+e+ae norm; 0 means no skew

  // detector systematics
  RooRealVar* fE_scale;   //!< energy scale parameter, such that \f$ E_{true} \rightarrow E_{true}E_{scale}\f$, 1 means no energy scaling
  Int_t       fE_tb_low;  //!< true bins below this bin nr are excluded from the analysis, this is necessary such that `fE_scale` cannot enter a region outside the MC range
  Int_t       fE_tb_high; //!< true bins above this bin nr are excluded from the analysis, this is necessary such that `fE_scale` cannot enter a region outside the MC range
  
  // xsec systematic parameters
  RooRealVar* fNC_norm;     //!< parameter for NC xsec normalisation; 1 means no scaling
  RooRealVar* fTau_norm;    //!< parameter for tau xsec normalisation; 1 means no scaling
  
  //--------------------------------------------------------------------
  // protected cache members for flux tilt cache
  //--------------------------------------------------------------------
  Double_t fTiltedFluxNorms[fFlavs][fPols]; //!< cache for normalisation of atm. flux shape systematics (energy, ct tilt)
   
};

#endif
