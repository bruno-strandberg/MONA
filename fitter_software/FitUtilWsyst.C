#include "FitUtilWsyst.h"

#include "TMath.h"

/** Constructor - see `FitUtil::FitUtil` constructor for parameter list.
    The added fit parameters are initialised here.
*/
FitUtilWsyst::FitUtilWsyst(Double_t op_time, TH3 *h_template,
			   Double_t emin, Double_t emax, Double_t ctmin, Double_t ctmax,
			   Double_t bymin, Double_t bymax, TString meff_file) :
  FitUtil(op_time, h_template, emin, emax, ctmin, ctmax, bymin, bymax, meff_file) {

  // init the new systematic parameters
  fE_tilt      = new RooRealVar("E_tilt"     ,"E_tilt"      , 0, -1, 1);
  fCt_tilt     = new RooRealVar("ct_tilt"    , "ct_tilt"    , 0, -1, 1);
  fSkew_mu_amu = new RooRealVar("skew_mu_amu", "skew_mu_amu", 1,  0, 2);
  fSkew_e_ae   = new RooRealVar("skew_e_ae"  , "skew_e_ae"  , 1,  0, 2);
  fSkew_mu_e   = new RooRealVar("skew_mu_e"  , "skew_mu_e"  , 1,  0, 2);

  // add the new parameters to `FitUtil::fParSet` to make them known to `FitPDF` and accessible in `FitUtil::proxymap_t`.
  fParSet.add( RooArgSet( *fE_tilt, *fCt_tilt, *fSkew_mu_amu, *fSkew_e_ae, *fSkew_mu_e) );

  // init flux shape systematics cache variables
  Double_t f_cache_e_tilt = 0;
  Double_t f_cache_ct_tilt = 0;
  CalcTiltedFluxNorms( f_cache_e_tilt, f_cache_ct_tilt );
    
};

//***************************************************************************************

/** Function that defines and returns the tilt coefficient for atmospheric flux
    \param energy   true neutrino energy
    \param costheta true neutrino cos-theta 
    \param e_tilt   energy tilt parameter value
    \param ct_tilt  cos-theta tilt parameter value
    \return \f$ energy^{e_{tilt}} * (1 + ct_{tilt}*\cos\theta) $\f
*/
Double_t FitUtilWsyst::FluxTiltCoeff( Double_t energy, Double_t costheta, Double_t e_tilt, Double_t ct_tilt ) {

  return TMath::Power(energy, e_tilt) * (1 + ct_tilt * costheta);
    
}

//***************************************************************************************

/** Protected function that re-calculates the normalisation coefficients each time one of the tilt parameters change.

    The normalisation coefficients are necessary to ensure that the tilt parameters do not affect the over-all normalisation of the atmospheric flux. I.e., the integral over energy and cos-theta of the flux without tilt systematics equals the integral with tilt systematics.

    \param e_tilt    energy tilt parameter value
    \param ct_tilt   cos-theta tilt parameter value

*/
void FitUtilWsyst::CalcTiltedFluxNorms(Double_t e_tilt, Double_t ct_tilt) {

  // calculate the normalisation factors for the tilted flux to preserve
  // the over-all normalisation of the flux

  for (UInt_t f = 0; f < fFlavs; f++) {
    for (UInt_t isnb = 0; isnb < fPols; isnb++) {

      Double_t numerator   = 0.;
      Double_t denominator = 0.;

      for (Int_t ebin = 1; ebin <= fHB->GetXaxis()->GetNbins(); ebin++) {
	for (Int_t ctbin = 1; ctbin <= fHB->GetYaxis()->GetNbins(); ctbin++) {

	  // get energy and cos-theta at bin center
	  Double_t E   = fHB->GetXaxis()->GetBinCenter(ebin);
	  Double_t ct  = fHB->GetYaxis()->GetBinCenter(ctbin);

	  // get the cached flux, add to numerator and denominator
	  Double_t flux = GetCachedFlux( f, isnb, ebin, ctbin );

	  numerator   += flux;
	  denominator += FluxTiltCoeff(E, ct, e_tilt, ct_tilt) * flux;

	}
      }

      // store the normalisation for the specific neutrino type
      fTiltedFluxNorms[f][isnb] = numerator/denominator;
	
    } // end loop over nu, anu
  } // end loop over flavors
    
}
  
//***************************************************************************************

/** Function to retrieve the atm. flux of a neutrino type that is modified by the tilt systematics.
    \param flav        neutrino flavor (0 - elec, 1 - muon, 2 - tau)
    \param isnb        flag for is-nubar (0 - nu, 1 - nubar)
    \param true_ebin   true energy bin of the flux
    \param true_ctbin  true cos-theta bin of the flux
    \param e_tilt      energy tilt parameter
    \param ct_tilt     cos-theta tilt parameter
    \return atmospheric neutrino flux in the specified bin for the specified neutrino type
 */
Double_t FitUtilWsyst::GetTiltedFlux(UInt_t flav, Bool_t isnb, Int_t true_ebin, Int_t true_ctbin,
				     Double_t e_tilt, Double_t ct_tilt) {

  if (e_tilt != f_cache_e_tilt || ct_tilt != f_cache_ct_tilt ) {
    f_cache_e_tilt  = e_tilt;
    f_cache_ct_tilt = ct_tilt;
    CalcTiltedFluxNorms( f_cache_e_tilt, f_cache_ct_tilt );
  }

  Double_t E  = fHB->GetXaxis()->GetBinCenter( true_ebin );
  Double_t ct = fHB->GetYaxis()->GetBinCenter( true_ctbin );
    
  return GetCachedFlux(flav, isnb, true_ebin, true_ctbin) * FluxTiltCoeff(E, ct, e_tilt, ct_tilt) * 
    fTiltedFluxNorms[flav][isnb];

}

//***************************************************************************************

/** Function to retrieve the atm. flux of a neutrino type that is modified by tilt systematics and skew systematics.
    \param flav        neutrino flavor (0 - elec, 1 - muon, 2 - tau)
    \param isnb        flag for is-nubar (0 - nu, 1 - nubar)
    \param true_ebin   true energy bin of the flux
    \param true_ctbin  true cos-theta bin of the flux
    \param proxymap    A proxymap from `FitPDF` that includes all of the fit parameters
    \return atmospheric neutrino flux in specified bin for the specified neutrino type
*/
Double_t FitUtilWsyst::GetFluxWsyst(UInt_t flav, Bool_t isnb, Int_t true_ebin, Int_t true_ctbin,
				    const proxymap_t& proxymap) {

  // get the skew parameters and tilt parameters
  Double_t e_tilt       = *( proxymap.at( (TString)fE_tilt->GetName() ) );
  Double_t ct_tilt      = *( proxymap.at( (TString)fCt_tilt->GetName() ) );
  Double_t skew_mu_amu  = *( proxymap.at( (TString)fSkew_mu_amu->GetName() ) );
  Double_t skew_e_ae    = *( proxymap.at( (TString)fSkew_e_ae->GetName() ) );
  Double_t skew_mu_e    = *( proxymap.at( (TString)fSkew_mu_e->GetName() ) );

  // get the tilted fluxes - overall normalisation is preserved
  Double_t atm_count_e   = GetTiltedFlux(ELEC, false, true_ebin, true_ctbin, e_tilt, ct_tilt);
  Double_t atm_count_mu  = GetTiltedFlux(MUON, false, true_ebin, true_ctbin, e_tilt, ct_tilt);
  Double_t atm_count_ae  = GetTiltedFlux(ELEC, true , true_ebin, true_ctbin, e_tilt, ct_tilt);
  Double_t atm_count_amu = GetTiltedFlux(MUON, true , true_ebin, true_ctbin, e_tilt, ct_tilt);

  
  // apply nu-antinu skew's that preserves nu+antinu flux
  Double_t mu_1  = atm_count_mu * skew_mu_amu;
  Double_t amu_1 = atm_count_amu + (1 - skew_mu_amu) * atm_count_mu;

  Double_t e_1  = atm_count_e * skew_e_ae;
  Double_t ae_1 = atm_count_ae + (1 - skew_e_ae) * atm_count_e;

  
  // apply muon-elec skew that preserves the over-all flux
  Double_t mu_2  =  mu_1 * skew_mu_e;
  Double_t amu_2 = amu_1 * skew_mu_e;

  Double_t e_2  =  e_1 +  e_1/(e_1 + ae_1) * (1 - skew_mu_e) * (mu_1 + amu_1);
  Double_t ae_2 = ae_1 + ae_1/(e_1 + ae_1) * (1 - skew_mu_e) * (mu_1 + amu_1);

  // create an array for return
  enum pols { NU = 0, ANU };
    
  Double_t fluxes[fFlavs][fPols];
  fluxes[MUON][NU]  = mu_2;
  fluxes[MUON][ANU] = amu_2;
  fluxes[ELEC][NU]  = e_2;
  fluxes[ELEC][ANU] = ae_2;
  fluxes[TAU][NU]   = 0.;
  fluxes[TAU][ANU]  = 0.;

  return fluxes[flav][isnb];
    
}

//***************************************************************************************
  
/** Overload of the virtual function `FitUtil::TrueEvts` to get the expected number of events in a true bin, given the set of parameters as specified in `proxymap`
    \param tb       `TrueB` with the true bin and neutrino type data, see `DetResponse.h`
    \param proxymap Container with all fit parameters known to `RooFit` from `FitPDF`.
*/
std::pair<Double_t, Double_t> FitUtilWsyst::TrueEvts(const TrueB &tb, const proxymap_t &proxymap) {

  // get the atm nu count
  Double_t atm_count_e = GetFluxWsyst(ELEC, tb.fIsNB, tb.fE_true_bin, tb.fCt_true_bin, proxymap);
  Double_t atm_count_m = GetFluxWsyst(MUON, tb.fIsNB, tb.fE_true_bin, tb.fCt_true_bin, proxymap);

  // get the oscillation probabilities
  Double_t prob_elec = GetCachedOsc(ELEC, tb, proxymap);
  Double_t prob_muon = GetCachedOsc(MUON, tb, proxymap);

  // get the oscillated nu count in operation time (in units 1/m2)
  Double_t osc_count = ( atm_count_e * prob_elec + atm_count_m * prob_muon );

  // get the interacted neutrino count in operation time (in units 1/MTon)
  Double_t int_count = osc_count * GetCachedXsec(tb)/fMN * fKg_per_MTon;

  // find true observable values necassary for the effective mass and to split events to bjorken-y bins
  Double_t e_true  = fHB->GetXaxis()->GetBinCenter( tb.fE_true_bin );
  Double_t ct_true = fHB->GetYaxis()->GetBinCenter( tb.fCt_true_bin );
  Double_t by_true = fHB->GetZaxis()->GetBinCenter( tb.fBy_true_bin );

  // get the effective mass
  Double_t meff  = fMeff->GetMeff( tb.fFlav, tb.fIsCC, tb.fIsNB, e_true, ct_true, by_true );

  // calculate the number of detected events (unitless)
  Double_t det_count = int_count * fXsec->GetBYfrac(tb.fFlav, tb.fIsCC, tb.fIsNB, e_true, by_true) * meff * 1e-6;

  // MC stat err, coming from eff mass and BY distribution (0 for now)
  Double_t det_err = 0.;

  return std::make_pair(det_count, det_err);

}
