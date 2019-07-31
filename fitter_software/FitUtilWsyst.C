#include "FitUtilWsyst.h"

#include "TMath.h"

/** Constructor - see `FitUtil::FitUtil` constructors for parameter list.
    The added fit parameters are initialised here.
*/
FitUtilWsyst::FitUtilWsyst(Double_t op_time, TH3 *h_temp_T, TH3* h_temp_R,
			   Double_t emin, Double_t emax, Double_t ctmin, Double_t ctmax,
			   Double_t bymin, Double_t bymax, TString meff_file) :
  FitUtil(op_time, h_temp_T, h_temp_R, emin, emax, ctmin, ctmax, bymin, bymax, meff_file) {

  /*----------------------------------------------------------------------------------
    Init the systematic parameters. It is important to add all common parameters to `FitUtil::fParSet` to make them known to `FitPDF` and accessible in `FitUtil::proxymap_t` that is passed to `RecoEvts` and subsequent functions.
  ----------------------------------------------------------------------------------*/

  // flux
  fE_tilt      = new RooRealVar("E_tilt"     ,"E_tilt"      , 0, -0.5, 0.5);
  fCt_tilt     = new RooRealVar("ct_tilt"    , "ct_tilt"    , 0, -0.5, 0.5);
  fSkew_mu_amu = new RooRealVar("skew_mu_amu", "skew_mu_amu", 0, -0.5, 0.5);
  fSkew_e_ae   = new RooRealVar("skew_e_ae"  , "skew_e_ae"  , 0, -0.5, 0.5);
  fSkew_mu_e   = new RooRealVar("skew_mu_e"  , "skew_mu_e"  , 0, -0.5, 0.5);
  fParSet.add( RooArgSet(*fE_tilt, *fCt_tilt, *fSkew_mu_amu, *fSkew_e_ae, *fSkew_mu_e) );
  
  // xsec
  fNC_norm     = new RooRealVar("NC_norm", "NC_norm", 1, 0, 2);
  fTau_norm    = new RooRealVar("Tau_norm", "Tau_norm", 1, 0, 2);
  fParSet.add( RooArgSet(*fNC_norm, *fTau_norm) );
  
  // detector
  fE_scale = new RooRealVar("E_scale", "E_scale", 0., -0.3, 0.3);
  fParSet.add( *fE_scale );
  
  // init flux shape systematics cache variables
  CalcTiltedFluxNorms( 0., 0. );

  //----------------------------------------------------------------------------------
  // the limits of `fE_reco` are set to bin edges in `FitUtil`. Check that the energy scale cannot
  // take the bin edges outside of the MC range
  //----------------------------------------------------------------------------------
  
  Double_t emin_scaled = fE_reco->getMin() * ( 1. + fE_scale->getMin() );
  Double_t emax_scaled = fE_reco->getMax() * ( 1. + fE_scale->getMax() );

  if ( emin_scaled < fHBR->GetXaxis()->GetXmin() ) {
    throw std::invalid_argument("ERROR! FitUtilWsyst::FitUtilWsyst() the input minimum energy " + to_string(emin) + " is adjusted to the bin edge " + to_string(fE_reco->getMin()) + ", which can be taken to the value " + to_string(emin_scaled) + " by the minimum of the energy scale parameter " + to_string(fE_scale->getMin()) + ". This is outside of the binning minimum " + to_string(fHBR->GetXaxis()->GetXmin()) + ". Choose a higher emin to avoid this problem." );
  }

  if ( emax_scaled >= fHBR->GetXaxis()->GetXmax() ) {
    throw std::invalid_argument("ERROR! FitUtilWsyst::FitUtilWsyst() the input maximum energy " + to_string(emax) + " is adjusted to the bin edge " + to_string(fE_reco->getMax()) + ", which can be taken to the value " + to_string(emax_scaled) + " by the maximum of the energy scale parameter " + to_string(fE_scale->getMax()) + ". This is outside of the binning maximum " + to_string(fHBR->GetXaxis()->GetXmax()) + ". Choose a lower emax to avoid this problem." );
  }
    
};

//***************************************************************************************

/** Same as the default constructor, but in this case the same binning is used for true bins and reco bins.*/
FitUtilWsyst::FitUtilWsyst(Double_t op_time, TH3 *h_template,
			   Double_t emin, Double_t emax, Double_t ctmin, Double_t ctmax,
			   Double_t bymin, Double_t bymax, TString meff_file) : 
  FitUtilWsyst(op_time, h_template, h_template, emin, emax, ctmin, ctmax, bymin, bymax, meff_file) {}

//***************************************************************************************

/** Function that defines and returns the tilt coefficient for atmospheric flux
    \param energy   true neutrino energy
    \param costheta true neutrino cos-theta 
    \param e_tilt   energy tilt parameter value
    \param ct_tilt  cos-theta tilt parameter value
    \return \f$ energy^{e_{tilt}} * (1 + ct_{tilt}*\cos\theta) \f$
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

      for (Int_t ebin = 1; ebin <= fHBT->GetXaxis()->GetNbins(); ebin++) {
	for (Int_t ctbin = 1; ctbin <= fHBT->GetYaxis()->GetNbins(); ctbin++) {

	  // get energy and cos-theta at bin center
	  Double_t E   = fHBT->GetXaxis()->GetBinCenter(ebin);
	  Double_t ct  = fHBT->GetYaxis()->GetBinCenter(ctbin);

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

/** Function to retrieve the atm flux of a neutrino type that is modified by the tilt systematics.
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

  // get cached parameter values
  Double_t cache_e_tilt  = GetCachedVar( (TString)fE_tilt->GetName() );
  Double_t cache_ct_tilt = GetCachedVar( (TString)fCt_tilt->GetName() );
  
  if (e_tilt != cache_e_tilt || ct_tilt != cache_ct_tilt ) {
    fParCache[ (TString)fE_tilt->GetName() ]  = e_tilt;
    fParCache[ (TString)fCt_tilt->GetName() ] = ct_tilt;
    CalcTiltedFluxNorms( e_tilt, ct_tilt );
  }

  Double_t E  = fHBT->GetXaxis()->GetBinCenter( true_ebin );
  Double_t ct = fHBT->GetYaxis()->GetBinCenter( true_ctbin );
    
  return GetCachedFlux(flav, isnb, true_ebin, true_ctbin) * FluxTiltCoeff(E, ct, e_tilt, ct_tilt) * 
    fTiltedFluxNorms[flav][isnb];

}

//***************************************************************************************

/** Function to retrieve the atm flux of a neutrino type that is modified by tilt systematics and skew systematics.
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
  Double_t e_tilt       = *( proxymap.at( (TString)fE_tilt->GetName() ) ) ;
  Double_t ct_tilt      = *( proxymap.at( (TString)fCt_tilt->GetName() ) );
  Double_t skew_mu_amu  = *( proxymap.at( (TString)fSkew_mu_amu->GetName() ) );
  Double_t skew_e_ae    = *( proxymap.at( (TString)fSkew_e_ae->GetName() ) )  ;
  Double_t skew_mu_e    = *( proxymap.at( (TString)fSkew_mu_e->GetName() ) )  ;

  // get the tilted fluxes - overall normalisation is preserved
  Double_t atm_count_e   = GetTiltedFlux(ELEC, false, true_ebin, true_ctbin, e_tilt, ct_tilt);
  Double_t atm_count_mu  = GetTiltedFlux(MUON, false, true_ebin, true_ctbin, e_tilt, ct_tilt);
  Double_t atm_count_ae  = GetTiltedFlux(ELEC, true , true_ebin, true_ctbin, e_tilt, ct_tilt);
  Double_t atm_count_amu = GetTiltedFlux(MUON, true , true_ebin, true_ctbin, e_tilt, ct_tilt);

  
  // apply nu-antinu skew's that preserves nu+antinu flux
  Double_t N_mu_amu = ( atm_count_mu + atm_count_amu )/( atm_count_mu * (1 + skew_mu_amu) + atm_count_amu );
  Double_t mu_1  = atm_count_mu * (1 + skew_mu_amu) * N_mu_amu;
  Double_t amu_1 = atm_count_amu * N_mu_amu;

  Double_t N_e_ae = ( atm_count_e + atm_count_ae )/( atm_count_e * (1 + skew_e_ae) + atm_count_ae );
  Double_t e_1  = atm_count_e * (1 + skew_e_ae) * N_e_ae;
  Double_t ae_1 = atm_count_ae * N_e_ae;

  
  // apply muon-elec skew that preserves the over-all flux
  Double_t N_mu_e = (mu_1 + amu_1 + e_1 + ae_1)/( (mu_1 + amu_1) * (1 + skew_mu_e) + e_1 + ae_1 );
  Double_t mu_2  =  mu_1 * (1 + skew_mu_e) * N_mu_e;
  Double_t amu_2 = amu_1 * (1 + skew_mu_e) * N_mu_e;

  Double_t e_2  =  e_1 * N_mu_e;
  Double_t ae_2 = ae_1 * N_mu_e;

  // create an array for return
  enum pols { NU = 0, ANU };
    
  Double_t fluxes[fFlavs][fPols];
  fluxes[MUON][NU]  = mu_2;
  fluxes[MUON][ANU] = amu_2;
  fluxes[ELEC][NU]  = e_2;
  fluxes[ELEC][ANU] = ae_2;
  fluxes[TAU][NU]   = 0.;
  fluxes[TAU][ANU]  = 0.;

  if ( fluxes[flav][isnb] < 0. ) {
    cout << "NOTICE FitUtilWsyst::GetFluxWsyst() parameter values at error:" << endl;
    for (auto kv: proxymap) {cout << "Parameter: " << kv.first << "\tValue: " << (Double_t)(*kv.second) << endl;}
    throw std::logic_error("ERROR! FitUtilWsyst::GetFluxWsyst() has calculated a negative flux, probably an issue with the implementation of systematic effects.");
  }

  return fluxes[flav][isnb];
    
}

//***************************************************************************************
  
/** Overload of the virtual function `FitUtil::TrueEvts` (see that for parameter definitions) to get the expected number of events in a true bin.
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

  // calculate xsec normalisation; compound of NC normalisation and tau normalisation
  Double_t xsec_norm = 1.;
  if ( !tb.fIsCC )       { xsec_norm *= *( proxymap.at( (TString)fNC_norm ->GetName() ) ); }
  if ( tb.fFlav == TAU ) { xsec_norm *= *( proxymap.at( (TString)fTau_norm->GetName() ) ); }

  // get the interacted neutrino count in operation time (in units 1/MTon), apply xsec systematic
  Double_t int_count = osc_count * GetCachedXsec(tb)/fMN * fKg_per_MTon * xsec_norm;

  // get the effective mass
  Double_t meff  = GetCachedMeff( tb );

  // calculate the number of detected events (unitless)
  Double_t det_count = int_count * GetCachedBYfrac(tb) * meff * 1e-6;

  // MC stat err, coming from eff mass and BY distribution (0 for now)
  Double_t det_err = 0.;

  return std::make_pair(det_count, det_err);

}

//***************************************************************************************
  
/** Overload of the virtual function `FitUtil::RecoEvts` (see that for parameter definitions) to get the expected number of events in a reco bin.
*/
std::pair<Double_t, Double_t> FitUtilWsyst::RecoEvts(Double_t E_reco, Double_t Ct_reco, Double_t By_reco, AbsResponse *resp, const proxymap_t &proxymap) {

  // get the energy scale
  Double_t e_scale = *( proxymap.at( (TString)fE_scale->GetName() ) );
  
  // scale the bin low and high edges and use these new limits to calculate the fractions that specify
  //how much of the reco event density should come from which bin
  Double_t ereco_bin = fHBR->GetXaxis()->FindBin( E_reco );
  Double_t ereco_min = fHBR->GetXaxis()->GetBinLowEdge( ereco_bin ) * (1. + e_scale);
  Double_t ereco_max = fHBR->GetXaxis()->GetBinUpEdge ( ereco_bin ) * (1. + e_scale);
  auto binfracs = GetBinFractions( ereco_min, ereco_max, fHBR->GetXaxis());  

  // loop over the bin fractions and calculate the event density from the relative fractions
  Double_t RE    = 0.;
  Double_t REerr = 0.;
  
  for (auto bf: binfracs) {

    // get the bin center and calculate the number 
    Double_t ereco = fHBR->GetXaxis()->GetBinCenter( bf.first );
    auto RE_bin = FitUtil::RecoEvts( ereco, Ct_reco, By_reco, resp, proxymap );
	    
    RE    += RE_bin.first * bf.second;
    REerr += TMath::Power(RE_bin.second * bf.second, 2);
          
  } // end loop over bin fractions
  
  // take sqrt of the error and return
  return std::make_pair( RE, TMath::Sqrt(REerr) );
   
}

//***************************************************************************************

/** Overload the virtual function `FitUtil::RecoEvtsER` to accommodate calculations with systematics, see that function for the parameter interface.

    The structuring of the virtual functions gets a bit complicated here, but that is to avoid duplicating code. As one can see, `FitUtilWsyst::RecoEvts` applies the energy scale systematic, but otherwise calls `FitUtil::RecoEvts`. The latter, in turn, calls two virtual functions, `FitUtil::RecoEvtsDR` and `FitUtil::RecoEvtsER`. The first of those performs the calculation using the binned response `DetResponse`, the latter uses the `EvtResponse`. `RecoEvtsDR` does not need to be modified here, because most of the systematics are taken care of in the `FitUtilWsyst::TrueEvts` function, that is used by `RecoEvtsDR`. On the other hand, `RecoEvtsER` needs to be modified, because it does not use `TrueEvts` and this one needs to account for the systematics here.
    
*/
std::pair< Double_t, Double_t > FitUtilWsyst::RecoEvtsER(Double_t E_reco, Double_t Ct_reco, Double_t By_reco,
							 EvtResponse *resp, const proxymap_t &proxymap,
							 Bool_t AddMuonsNoise) {

  Double_t det_count = 0.;
  Double_t det_err   = 0.;

  auto true_evts = resp->GetBinEvts( E_reco, Ct_reco, By_reco );

  for (const auto &te: true_evts) {

    // find bin coordinates for retrieving info for the cache
    Int_t e_bin  = fHBT->GetXaxis()->FindBin( te.GetTrueE() );
    Int_t ct_bin = fHBT->GetYaxis()->FindBin( te.GetTrueCt() );
    Int_t by_bin = fHBT->GetZaxis()->FindBin( te.GetTrueBy() );

    // get bin widths to calculate differential flux from the cache
    Double_t ew  = fHBT->GetXaxis()->GetBinWidth(e_bin);
    Double_t ctw = fHBT->GetYaxis()->GetBinWidth(ct_bin);

    // calculate constant to convert back to differential flux and calculate the fluxes
    Double_t flux_conv = ew * ctw * fSec_per_y * TMath::Pi() * 2;

    // use `GetFluxWsyst` instead of GetCachedFlux, this takes care of the flux systematics
    Double_t atm_count_e = GetFluxWsyst(ELEC, te.GetIsNB(), e_bin, ct_bin, proxymap) / flux_conv;
    Double_t atm_count_m = GetFluxWsyst(MUON, te.GetIsNB(), e_bin, ct_bin, proxymap) / flux_conv;

    if ( te.GetIsCC() ) {

      //create trueB object to access cached oscillation probabilities
      TrueB TB( te.GetFlav(), te.GetIsCC(), te.GetIsNB(), e_bin, ct_bin, by_bin);
	
      Double_t prob_elec = GetCachedOsc(ELEC, TB, proxymap);
      Double_t prob_muon = GetCachedOsc(MUON, TB, proxymap);
		   
      // calculate the oscillated differential flux for the neutrino type
      Double_t oscf = atm_count_e*prob_elec + atm_count_m*prob_muon;

      // apply tau normalisation
      if ( te.GetFlav() == TAU ) oscf *= *( proxymap.at( (TString)fTau_norm->GetName() ) );

      det_count += te.GetW1y() * oscf;
      det_err++;
		      
    }
    else {

      // for NC, the xsec and meff for all flavors are the ~same. For this reason, only elec-NC are simulated.
      // For this calculation, however, the NC flux of muons and tau's also needs to be taken into account.
      // The calculation for CC events after oscillation is osc_flux * xsec_{flavor} * meff_{flavor}. As
      // xsec_{flavor} and meff_{flavor} are the same for elec, muon, tau for NC events, the calculation
      // can be simplified and the oscillated flux is just equal to the un-oscillated elec+muon flux

      if ( te.GetFlav() != ELEC ) {
	throw std::invalid_argument("ERROR! FitUtilWsyst::RecoEvtsER() expects only elec events for NC. If the response is filled also with muon and tau NC data, this function requires modification, current implementation leads to double-counting.");
      }
      
      Double_t oscf = atm_count_e + atm_count_m;

      // calculate the fraction of tau events in the total flux to apply tau normalisation also to NC events
      TrueB TB_tau( TAU, kFALSE, te.GetIsNB(), e_bin, ct_bin, by_bin );	
      Double_t prob_e_tau  = GetCachedOsc(ELEC, TB_tau, proxymap); // probability for elec->tau
      Double_t prob_mu_tau = GetCachedOsc(MUON, TB_tau, proxymap); // probability for muon->tau
      Double_t oscf_tau    = atm_count_e*prob_e_tau + atm_count_m*prob_mu_tau; // oscillated tau flux

      oscf = oscf + oscf_tau * ( *( proxymap.at( (TString)fTau_norm ->GetName() ) ) - 1 );

      // apply NC normalisation
      oscf *= *( proxymap.at( (TString)fNC_norm ->GetName() ) );

      det_count += te.GetW1y() * oscf;
      det_err++;
	
    }

  } // end loop over true events

  // convert the error to the fraction of the neutrino events
  det_err = TMath::Sqrt(det_err)/det_err * det_count;

  if ( AddMuonsNoise ) {

    auto muons = resp->GetAtmMuCount1y(E_reco, Ct_reco, By_reco);
    auto noise = resp->GetNoiseCount1y(E_reco, Ct_reco, By_reco);
    det_count += muons.first * fOpTime;
    det_count += noise.first * fOpTime;
    Double_t errsq = det_err * det_err + 
      TMath::Power(muons.second, 2) * fOpTime + TMath::Power(noise.second, 2) * fOpTime;

    det_err = TMath::Sqrt( errsq );
  }

  return std::make_pair( det_count, det_err );

}

//***************************************************************************************

/** This function looks up the bin fractions from the axis in the low to high range.
    
    The function looks up all the bins on the axis in the range low to high. Bins that are fully enclosed in the range are assigned weight 1. For bins that are partially enclosed (these can only be the first bin and the last bin), the weight is calculated as a fraction of the covered range. For example, the fraction for the low bin is calculated as \f$ ( upEdge(bin_{lo}) - lo) )/binWidth(bin_{lo}) \f$.

    If the low or high value is outside of the axis range, an error is raised.
    
    \param lo             low value
    \param hi             high value
    \param axis           pointer to `TAxis` where the bins are searched for
    \return      A vector of pairs, where the first is the bin number and the second is the bin weight.
*/
vector< std::pair<Int_t, Double_t> > FitUtilWsyst::GetBinFractions(Double_t lo, Double_t hi, TAxis* axis) {

  if ( hi <= lo ) {
    throw std::invalid_argument("ERROR! FitUtilWsyst::GetBinFractions() high value is smaller than low value");
  }

  if ( lo < axis->GetXmin() ) {
    throw std::invalid_argument("ERROR! FitUtilWsyst::GetBinFractions() low value " + to_string(lo) + " is outside the axis range.");
  }

  if ( hi > axis->GetXmax() ) {
    throw std::invalid_argument("ERROR! FitUtilWsyst::GetBinFractions() high value " + to_string(hi) + " is outside the axis range.");
  }

  // find the bins corresponding to low and high and calculate the relative
  // fraction of the events from the corresponding bins
  //---------------------------------------------------------------------
  
  vector< std::pair<Int_t, Double_t> > bins;
  
  // find the bins where low value and high value reside
  Int_t bin_lo = axis->FindBin( lo );
  Int_t bin_hi = axis->FindBin( hi );

  Double_t binW_lo = ( axis->GetBinUpEdge( bin_lo ) - lo )/axis->GetBinWidth( bin_lo );
  Double_t binW_hi = ( hi - axis->GetBinLowEdge(bin_hi) )/axis->GetBinWidth( bin_hi );

  // add the bins to the return vector
  //---------------------------------------------------------------------

  // add the lower bin to the return vector, if bin weight is larger than 0
  if (binW_lo > 0.0) bins.push_back( std::make_pair(bin_lo, binW_lo) );

  // add all the intermediate bins with weight 1; this happens if the range encloses several bins
  for (Int_t binnr = bin_lo+1; binnr < bin_hi; binnr++) {
    bins.push_back( std::make_pair(binnr, 1.) );
  }
  
  // add the higher bin to the return vector, if binw is larger than 0
  if (binW_hi > 0.0) bins.push_back( std::make_pair(bin_hi, binW_hi) );

  return bins;
  
}
