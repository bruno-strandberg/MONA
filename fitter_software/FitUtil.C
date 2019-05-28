#include "FitUtil.h"
#include <iostream>
#include <typeinfo>
#include "TFile.h"
#include "TRandom.h"

using namespace std;

/** Constructor
    \param op_time         Operation time in years
    \param h_temp_T        A 3D histogram that represents the binning settings in true space
    \param h_temp_R        A 3D histogram that represents the binning settings in reco space
    \param emin            Minimum energy included in the fit range
    \param emax            Maximum energy included in the fit range
    \param ctmin           Minimum cos-theta included in the fit range
    \param ctmax           Maximum cos-theta included in the fit range
    \param bymin           Minimum bjorken-y included in the fit range
    \param bymax           Maximum bjorken-y included in the fit range
    \param meff_file       File with effective mass data (created with `EffMass` class)
*/
FitUtil::FitUtil(Double_t op_time, TH3 *h_temp_T, TH3* h_temp_R,
		 Double_t emin, Double_t emax, Double_t ctmin, Double_t ctmax, Double_t bymin, Double_t bymax,
		 TString meff_file) {

  fOpTime = op_time;

  // histogram with true binning, used for flux->osc->xsec->meff
  fHBT = (TH3D*)h_temp_T->Clone("HBT");
  fHBT->SetDirectory(0);
  fHBT->Reset();
  fHBT->SetNameTitle("HBT", "HBT");

  // histogram with reco binning
  fHBR = (TH3D*)h_temp_R->Clone("HBR");
  fHBR->SetDirectory(0);
  fHBR->Reset();
  fHBR->SetNameTitle("HBR", "HBR");
  
  fFlux = new AtmFlux;
  fXsec = new NuXsec( fHBT->GetZaxis()->GetNbins() );
  fMeff = new EffMass( meff_file, fHBT->GetXaxis()->GetNbins(), fHBT->GetYaxis()->GetNbins(), fHBT->GetZaxis()->GetNbins() );
  fProb = new OscProb::PMNS_Fast;
  fPrem = new OscProb::PremModel;

  // calculate the observable ranges to match exactly bin edges, set ranges for integration
  auto ENr = GetRange( emin , emax , fHBR->GetXaxis() );
  auto CTr = GetRange( ctmin, ctmax, fHBR->GetYaxis() );
  auto BYr = GetRange( bymin, bymax, fHBR->GetZaxis() );

  fEbin_min  = get<MINBIN>(ENr);
  fEbin_max  = get<MAXBIN>(ENr);
  fCtbin_min = get<MINBIN>(CTr);
  fCtbin_max = get<MAXBIN>(CTr);
  fBybin_min = get<MINBIN>(BYr);
  fBybin_max = get<MAXBIN>(BYr);
  
  InitFitVars(get<MIN>(ENr), get<MAX>(ENr), get<MIN>(CTr), get<MAX>(CTr), get<MIN>(BYr), get<MAX>(BYr));
  InitCacheHists(fHBT);

  fFluxSamplesN = 1;
  FillFluxCache(fFlux, fOpTime, fFluxSamplesN);
  FillXsecMeffCache(fXsec, fMeff);

  fOscSamplesN        = 1;
  f_cache_oscsamplesn = 0;
  
  fOscCalls    = 0;
  fOscCalcTime = new TStopwatch();
  fOscCalcTime->Stop(); fOscCalcTime->Reset();

}

//***************************************************************************

/** Constructor
    \param op_time         Operation time in years
    \param h_temp          A 3D histogram that represents the binning settings in both true and reco space
    \param emin            Minimum energy included in the fit range
    \param emax            Maximum energy included in the fit range
    \param ctmin           Minimum cos-theta included in the fit range
    \param ctmax           Maximum cos-theta included in the fit range
    \param bymin           Minimum bjorken-y included in the fit range
    \param bymax           Maximum bjorken-y included in the fit range
    \param meff_file       File with effective mass data (created with `EffMass` class)
*/
FitUtil::FitUtil(Double_t op_time, TH3 *h_temp,
		 Double_t emin, Double_t emax, Double_t ctmin, Double_t ctmax, Double_t bymin, Double_t bymax,
		 TString meff_file) : FitUtil(op_time, h_temp, h_temp, emin, emax, ctmin, ctmax, bymin, bymax, meff_file) {}

//***************************************************************************

/** Destructor */
FitUtil::~FitUtil() {

  cout << "FitUtil::~FitUtil() total oscillator calls: " << fOscCalls << endl;
  cout << "FitUtil::~FitUtil() duration of oscillation calculations [s]: "
       << (Double_t)fOscCalcTime->RealTime() << endl;

  if (fOscCalcTime) delete fOscCalcTime;

  TAxis* eaxis  = fHBT->GetXaxis();
  TAxis* ctaxis = fHBT->GetYaxis();

  // remove cache hists
  for (UInt_t f = 0; f < fFlavs; f++) {
    for (UInt_t isnb = 0; isnb < fPols; isnb++) {

      // osc
      for (UInt_t fout = 0; fout < fFlavs; fout++) {
	if ( fOscCache[f][fout][isnb] ) ClearCache2D( fOscCache[f][fout][isnb], eaxis );
      }

      // flux
      if ( fFluxCache[f][isnb] ) ClearCache2D(fFluxCache[f][isnb], eaxis );

      // xsec and meff
      for (UInt_t iscc = 0; iscc < fInts; iscc++) {
	
	if ( fXsecCache[f][iscc][isnb] ) ClearCache1D( fXsecCache[f][iscc][isnb] );
	if ( fBYfracCache[f][iscc][isnb] ) ClearCache2D( fBYfracCache[f][iscc][isnb], eaxis );
	if ( fMeffCache[f][iscc][isnb] ) ClearCache3D( fMeffCache[f][iscc][isnb], eaxis, ctaxis );
	if ( fTECache[f][iscc][isnb] ) ClearCache3D( fTECache[f][iscc][isnb], eaxis, ctaxis );
	
      }

    }
  }

  // remove binning hist and calculators
  if (fHBT)  delete fHBT;
  if (fHBR)  delete fHBR;
  if (fFlux) delete fFlux;
  if (fXsec) delete fXsec;
  if (fMeff) delete fMeff;
  if (fProb) delete fProb;
  if (fPrem) delete fPrem;
  
  // remove the RooFit variables
  RooArgList l( fParSet );
  for (Int_t i = 0; i < l.getSize(); i++) {
    RooRealVar *v = (RooRealVar*)l.at(i);
    if (v) delete v;
  }

}

//***************************************************************************

/** Private function to find a range in an axis that exactly matches bin edges.

    For example, consider 40 logarithmic energy bins from 1-100 GeV. The user can ask the fit to be performed in the range from 2 - 82 GeV. This function find the bin edges closest to the dialled range. For example, the range 2 - 82 would become something like 2.42 - 79.2. This is necessary to make sure that data import to `RooFit`, fit ranges and the detector response are consistent.

    \param _min   Minimum of the range
    \param _max   Maximum of the range
    \param axis   Pointer to an axis that defines the binning
    \return       A std::tuple with the minimum, the maximum, the minimum bin number and the maximum bin number
*/
std::tuple<Double_t, Double_t, Int_t, Int_t> FitUtil::GetRange(Double_t _min, Double_t _max, TAxis *axis) {

  Int_t bin_min = 1;
  Int_t bin_max = axis->GetNbins();
  Double_t min  = axis->GetBinLowEdge( bin_min );
  Double_t max  = axis->GetBinUpEdge( bin_max );
  
  if ( _min > min ) {
    bin_min = axis->FindBin( _min );
    min = axis->GetBinLowEdge( bin_min );
  }

  if ( _max < max ) {
    bin_max = axis->FindBin(_max); 
    max = axis->GetBinUpEdge( bin_max );
  }
  
  return std::make_tuple(min, max, bin_min, bin_max);
  
}

//***************************************************************************

/** Private function to initialise `RooFit` variables.
    \param emin   Minimum energy in the fit range
    \param emax   Maximum energy in the fit range
    \param ctmin  Minimum cos-theta in the fit range
    \param ctmax  Maximum cos-theta in the fit range
    \param bymin  Minimum bjorken-y in the fit range
    \param bymax  Maximum bjorken-y in the fit range
*/
void FitUtil::InitFitVars(Double_t emin, Double_t emax, Double_t ctmin, Double_t ctmax, 
			  Double_t bymin, Double_t bymax) {

  // observables
  fE_reco    = new RooRealVar("E_reco" , "E_reco" , emin, emax);
  fCt_reco   = new RooRealVar("ct_reco", "ct_reco", ctmin, ctmax);
  fBy_reco   = new RooRealVar("by_reco", "by_reco", bymin, bymax);

  // fit parameters, initialised at NO central values, free ranges
  fSinsqTh12 = new RooRealVar("SinsqTh12", "sin^2(theta12)", f_NO_sinsqth12.cv, 0., 1.);
  fSinsqTh13 = new RooRealVar("SinsqTh13", "sin^2(theta13)", f_NO_sinsqth13.cv, 0., 1.);
  fSinsqTh23 = new RooRealVar("SinsqTh23", "sin^2(theta23)", f_NO_sinsqth23.cv, 0., 1.);
  fDcp       = new RooRealVar(      "dcp",       "delta-cp",       f_NO_dcp.cv, f_NO_dcp.min   , f_NO_dcp.max);
  fDm21      = new RooRealVar(     "Dm21",         "dm21^2",      f_NO_dm21.cv, f_free_dm21.min, f_free_dm21.max);
  fDm31      = new RooRealVar(     "Dm31",         "dm31^2",      f_NO_dm31.cv, f_free_dm31.min, f_free_dm31.max);

  // add observables to the observables set
  fObsList.add( RooArgList(*fE_reco, *fCt_reco, *fBy_reco) );

  // add all variables to the parameter set
  fParSet.add( fObsList );
  fParSet.add( RooArgSet( *fSinsqTh12, *fSinsqTh13, *fSinsqTh23, *fDcp, *fDm21, *fDm31) );
  
}

//***************************************************************************

/** Function to set the oscillation parameter limits corresponding to normal mass ordering 
 */
void FitUtil::SetNOlims() {

  fSinsqTh12->setMin(f_NO_sinsqth12.min);
  fSinsqTh12->setMax(f_NO_sinsqth12.max);
  
  fSinsqTh13->setMin(f_NO_sinsqth13.min);
  fSinsqTh13->setMax(f_NO_sinsqth13.max);

  fSinsqTh23->setMin(f_NO_sinsqth23.min);
  fSinsqTh23->setMax(f_NO_sinsqth23.max);

  fDcp->setMin(f_NO_dcp.min);
  fDcp->setMax(f_NO_dcp.max);

  fDm21->setMin(f_NO_dm21.min);
  fDm21->setMax(f_NO_dm21.max);
  
  fDm31->setMax(f_NO_dm31.max);
  fDm31->setMin(f_NO_dm31.min);
    
}

//***************************************************************************

/** Function to set the oscillation parameters to central values that correspond to normal mass ordering */
void FitUtil::SetNOcentvals() {
  
  fSinsqTh12->setVal(f_NO_sinsqth12.cv);
  fSinsqTh13->setVal(f_NO_sinsqth13.cv);
  fSinsqTh23->setVal(f_NO_sinsqth23.cv);
  fDcp->setVal(f_NO_dcp.cv);
  fDm21->setVal(f_NO_dm21.cv);
  fDm31->setVal(f_NO_dm31.cv);
  
}

//***************************************************************************

/** Function to set the oscillation parameter limits corresponding to inverted mass ordering 
 */
void FitUtil::SetIOlims() {

  fSinsqTh12->setMin(f_IO_sinsqth12.min);
  fSinsqTh12->setMax(f_IO_sinsqth12.max);
  
  fSinsqTh13->setMin(f_IO_sinsqth13.min);
  fSinsqTh13->setMax(f_IO_sinsqth13.max);

  fSinsqTh23->setMin(f_IO_sinsqth23.min);
  fSinsqTh23->setMax(f_IO_sinsqth23.max);

  fDcp->setMin(f_IO_dcp.min);
  fDcp->setMax(f_IO_dcp.max);

  fDm21->setMin(f_IO_dm21.min);
  fDm21->setMax(f_IO_dm21.max);
  
  fDm31->setMin(f_IO_dm31.min);
  fDm31->setMax(f_IO_dm31.max);
    
}

//***************************************************************************

/** Function to set the oscillation parameters to central values that correspond to inverted mass ordering */
void FitUtil::SetIOcentvals() {
  
  fSinsqTh12->setVal(f_IO_sinsqth12.cv);
  fSinsqTh13->setVal(f_IO_sinsqth13.cv);
  fSinsqTh23->setVal(f_IO_sinsqth23.cv);
  fDcp->setVal(f_IO_dcp.cv);
  fDm21->setVal(f_IO_dm21.cv);
  fDm31->setVal(f_IO_dm31.cv);
  
}

//***************************************************************************

/** This function sets free parameter limits to the oscillation parameters.

    This means that the sin^2 of the angles are between 0 and 1, delta-cp is between 0-2 pi. dm21^2 is limited from 5e-5 to 1e-4, dm31^2 is limited from -5e-3 to 5e-3.

*/
void FitUtil::FreeParLims() {

  fSinsqTh12->setMin(0.);
  fSinsqTh12->setMax(1.);
  
  fSinsqTh13->setMin(0.);
  fSinsqTh13->setMax(1.);

  fSinsqTh23->setMin(0.);
  fSinsqTh23->setMax(1.);

  fDcp->setMin(f_NO_dcp.min);
  fDcp->setMax(f_NO_dcp.max);

  fDm21->setMin(f_free_dm21.min);
  fDm21->setMax(f_free_dm21.max);
  
  fDm31->setMin(f_free_dm31.min);
  fDm31->setMax(f_free_dm31.max);

}

//***************************************************************************

/** Private function to initialize histograms for caching variabes 
    \param h_template  a TH3 template histogram that stores the binning information
 */
void FitUtil::InitCacheHists(TH3D *h_template) {

  TAxis* eaxis  = h_template->GetXaxis();
  TAxis* ctaxis = h_template->GetYaxis();
  TAxis* byaxis = h_template->GetZaxis();

  // init the oscillation probability cache histograms
  for (UInt_t f_in = 0; f_in < fFlavs; f_in++) {
    for (UInt_t f_out = 0; f_out < fFlavs; f_out++) {
      for (UInt_t isnb = 0; isnb < fPols; isnb++) {
	fOscCache[f_in][f_out][isnb] = InitCache2D( eaxis, ctaxis );
      }
    }
  }

  // init the flux histograms
  for (UInt_t f = 0; f < fFlavs; f++) {
    for (UInt_t isnb = 0; isnb < fPols; isnb++) {
      fFluxCache[f][isnb] = InitCache2D( eaxis, ctaxis );
    }
  }
  
  // init the xsec histograms, meff histograms and TrueEvts cache histograms
  for (UInt_t f = 0; f < fFlavs; f++) {    
    for (UInt_t iscc = 0; iscc < fInts; iscc++) {
      for (UInt_t isnb = 0; isnb < fPols; isnb++) {

	fXsecCache[f][iscc][isnb] = InitCache1D( eaxis );
	fBYfracCache[f][iscc][isnb] = InitCache2D( eaxis, byaxis );
	fMeffCache[f][iscc][isnb] = InitCache3D( eaxis, ctaxis, byaxis );
	fTECache[f][iscc][isnb] = InitCache3D( eaxis, ctaxis, byaxis );
	
      }
    }
  }

}

//***************************************************************************

/** A private function to fill the flux cache hists
    \param flux     Pointer to the `AtmFlux` member instance
    \param op_time  Operation time in years.
    \param nsamples Parameter to increase flux calculation accuracy. If 1 bin center is used, otherwise each bin is sub-divided to nsamples along energy and cos-theta, e.g. nsamples=2 means that each bin is calculated as an average of 2^2=4 points.
*/
void FitUtil::FillFluxCache(AtmFlux *flux, Double_t op_time, UInt_t nsamples) {

  // fill the atm flux cache. Note that (currently) AtmFlux returns 0 for tau flux
  for (UInt_t f = 0; f < fFlavs; f++) {
    for (UInt_t isnb = 0; isnb < fPols; isnb++) {

      TH3D* hflux = fHBT; // get the true bin space template histogram

      for (Int_t ebin = 1; ebin <= hflux->GetXaxis()->GetNbins(); ebin++) {
	for (Int_t ctbin = 1; ctbin <= hflux->GetYaxis()->GetNbins(); ctbin++) {

	  // create a subdivision to improve the flux calculation accuracy
	  Double_t emin  = hflux->GetXaxis()->GetBinLowEdge(ebin);
	  Double_t emax  = hflux->GetXaxis()->GetBinUpEdge(ebin);
	  Double_t ctmin = hflux->GetYaxis()->GetBinLowEdge(ctbin);
	  Double_t ctmax = hflux->GetYaxis()->GetBinUpEdge(ctbin);
	  TH2D hb("hb", "hb", nsamples, emin, emax, nsamples, ctmin, ctmax);

	  // loop over the subdivision bins to calculate the average flux
	  Double_t atmflux = 0.;

	  for (Int_t eb = 1; eb <= hb.GetXaxis()->GetNbins(); eb++) {
	    for (Int_t cb = 1; cb <= hb.GetYaxis()->GetNbins(); cb++) {

	      Double_t E   = hb.GetXaxis()->GetBinCenter(eb);
	      Double_t ew  = hb.GetXaxis()->GetBinWidth(eb);
	      Double_t ct  = hb.GetYaxis()->GetBinCenter(cb);
	      Double_t ctw = hb.GetYaxis()->GetBinWidth(cb);

	      Double_t atm_flux_factor = ew * ctw * fSec_per_y * op_time;
	      atmflux += flux->Flux_dE_dcosz(f, isnb, E, ct) * atm_flux_factor;
	  	      
	    }
	  }

	  fFluxCache[f][isnb][ebin][ctbin] = std::make_pair( atmflux/(nsamples*nsamples), 0. );

	}
      }

    }
  }

}

//***************************************************************************

/** A private function to fill the xsec and eff mass cache hists
    \param xsec     Pointer to the `NuXsec` member instance
    \param meff     Pointer to the `EffMass` member instance
*/
void FitUtil::FillXsecMeffCache(NuXsec *xsec, EffMass *meff) {

  // fill xsec cache and meff cache
  for (UInt_t f = 0; f < fFlavs; f++) {
    for (UInt_t iscc = 0; iscc < fInts; iscc++) {
      for (UInt_t isnb = 0; isnb < fPols; isnb++) {
	
	TH3D* hb = fHBT;
	
	for (Int_t ebin = 1; ebin <= hb->GetXaxis()->GetNbins(); ebin++) {

	  Double_t E = hb->GetXaxis()->GetBinCenter(ebin);
	  fXsecCache[f][iscc][isnb][ebin] = std::make_pair( xsec->GetXsec(f, iscc, isnb, E), 0. );

	  for (Int_t bybin = 1; bybin <= hb->GetZaxis()->GetNbins(); bybin++) {

	    Double_t by = hb->GetZaxis()->GetBinCenter(bybin);
	    fBYfracCache[f][iscc][isnb][ebin][bybin] = std::make_pair( xsec->GetBYfrac(f, iscc, isnb, E, by), 0. );
	    
	    for (Int_t ctbin = 1; ctbin <= hb->GetYaxis()->GetNbins(); ctbin++) {

	      Double_t ct = hb->GetYaxis()->GetBinCenter(ctbin);
	      fMeffCache[f][iscc][isnb][ebin][ctbin][bybin] = std::make_pair( meff->GetMeff( f, iscc, isnb, E, ct, by ), 0.);
	      
	    } // end loop over by
	  } // end loop over ct
	} // end loop over energy
	
      } // end loop over isnb
    } // end loop over iscc
  } // end loop over flavors

}

//***************************************************************************

/** 
    Function that returns the cached atmospheric flux value for a certain true bin.
    \param flav        Neutrino flavor, 0 - elec, 1 - muon, 2 - tau
    \param isnb        Flag for anti-neutrino, 0 - neutrino, 1 - anti-neutrino
    \param true_ebin   True energy bin
    \param true_ctbin  True cos-theta bin
    \return      Cached atmoshperic neutrino flux in specified bin.
*/
Double_t FitUtil::GetCachedFlux(UInt_t flav, Bool_t isnb, Int_t true_ebin, Int_t true_ctbin) {

  if (flav > TAU) {
    throw std::invalid_argument("ERROR! FitUtil::GetCachedFlux() unknown flavor " + to_string(flav));
  }

  return fFluxCache[flav][(UInt_t)isnb][true_ebin][true_ctbin].first;
}

//***************************************************************************

/** 
    Function that returns the cached oscillation probability for a certain true bin.

    The function handles re-calculation of the cached values by calling `FitUtil::ProbCacher` whenever any of the oscillation parameters change.

    For example, to retrieve the oscillation probability of elec-->muon for a certain energy and cos-theta combination, the param `flav_in` should be `ELEC`, and `tb` should specify a muon-neutrino with a certain energy and cos-theta bin coordinates.

    \param flav_in     Incoming atmospheric neutrino flavor (0-elec, 1 - muon). Tau (2) not supported, as atmospheric tau flux in negligible.
    \param tb          A `TrueB` object (see `DetResponse.h`) that stores the neutrino type and true bin coordinate info. The flavor stored in `tb` is the oscillated neutrino flavor.
    \param proxymap    A proxy map with `RooFit` variables from `FitPDF` that contains the oscillation parameters. 
    \return            Cached neutrino oscillation probability in specified bin.
*/
Double_t FitUtil::GetCachedOsc(UInt_t flav_in, const TrueB &tb, const proxymap_t &proxymap) {

  if (flav_in > MUON) {
    throw std::invalid_argument("ERROR! FitUtil::GetCachedOsc() tau-->flavor oscillations are not calculated by FitUtil::ProbCacher(), as the atmospheric flux of tau's is negligible compared to muon and elec flux at ORCA energies");
  }

  // recalculate the oscillation parameter cache if anything has changed.
  ProbCacher(proxymap, fOscSamplesN);
  
  return fOscCache[flav_in][tb.fFlav][tb.fIsNB][tb.fE_true_bin][tb.fCt_true_bin].first;
}



//***************************************************************************

/** 
    Function that returns the cached neutrino cross-section value for a certain true bin.
    \param tb    A `TrueB` object (see `DetResponse.h`) that stores the neutrino type and true bin coordinate info.
    \return      Cached cross-section value in a specified bin.
*/
Double_t FitUtil::GetCachedXsec(const TrueB &tb) {
  return fXsecCache[tb.fFlav][tb.fIsCC][tb.fIsNB][tb.fE_true_bin].first;
}

//***************************************************************************

/**
   Function that returns the cached effective mass value for a certain true bin
   \param tb    A `TrueB` object (see `DetResponse.h`) that stores the neutrino type and true bin coordinate info.
   \return      Cached effective mass value in a specified bin in Ton
 */
Double_t FitUtil::GetCachedMeff(const TrueB &tb) {
  return fMeffCache[tb.fFlav][tb.fIsCC][tb.fIsNB][tb.fE_true_bin][tb.fCt_true_bin][tb.fBy_true_bin].first;
}

//***************************************************************************

/** 
    Function that returns the cached fraction of events expected in a certain BY bin at a certain energy.
   \param tb    A `TrueB` object (see `DetResponse.h`) that stores the neutrino type and true bin coordinate info.
   \return      Cached bjorken-y fraction
*/
Double_t FitUtil::GetCachedBYfrac(const TrueB &tb) {
  return fBYfracCache[tb.fFlav][tb.fIsCC][tb.fIsNB][tb.fE_true_bin][tb.fBy_true_bin].first;
}

//***************************************************************************

/** Function that configures the `OscProb` instance of the class for given parameter values.

    By default, the OscProb calculator `FitUtil::fProb` is initiated to `OscProb::PMNS_Fast`. If a different calculator is required, e.g. `OscProb::PMNS_NSI`, the class inheriting from `FitUtil` should do in the constructor

    \code{.cpp}
    if (fProb) delete fProb;
    fProb = new OscProb::PMNS_NSI()
    \endcode

    and the class should re-implement this function, such that correct fit parameters are passed to the `FitUtil::fProb` instance.

    \param proxymap A proxy map with `RooFit` variables from `FitPDF` that contain all shared observables and parameters, including the oscillation parameters. 
    \return True - osc parameters have changed, re-calculate the cache; False - osc parameters have not changed, re-calculation not required
*/
Bool_t FitUtil::ConfigOscProb(const proxymap_t& proxymap) {

  // check that the type is PMNS_Fast
  TString type = (TString)typeid( *fProb ).name();

  if ( !type.Contains("PMNS_Fast") ) {
    throw std::logic_error("ERROR! FitUtil::ConfigOscProb() is called with fProb instance type " + (string)type + ", whereas OscProb::PMNS_Fast is expected. If a class inheriting from FitUtil has initiated fProb with a different type, the virtual function `FitUtil::ConfigOscProb` needs to be over-loaded to account for the new/different parameters.");
  }

  Bool_t reCalc = kFALSE;

  // get the parameter values from the proxies
  Double_t SinsqTh12 = *( proxymap.at( (TString)fSinsqTh12->GetName() ) );
  Double_t SinsqTh13 = *( proxymap.at( (TString)fSinsqTh13->GetName() ) );
  Double_t SinsqTh23 = *( proxymap.at( (TString)fSinsqTh23->GetName() ) );
  Double_t Dcp       = *( proxymap.at( (TString)fDcp->GetName() ) );
  Double_t Dm21      = *( proxymap.at( (TString)fDm21->GetName() ) );
  Double_t Dm31      = *( proxymap.at( (TString)fDm31->GetName() ) );

  // get the pointers to the cache variables in the cache map
  Double_t cache_sinsqth12 = GetCachedVar( (TString)fSinsqTh12->GetName() );
  Double_t cache_sinsqth13 = GetCachedVar( (TString)fSinsqTh13->GetName() );
  Double_t cache_sinsqth23 = GetCachedVar( (TString)fSinsqTh23->GetName() );
  Double_t cache_dcp       = GetCachedVar( (TString)fDcp->GetName() );
  Double_t cache_dm21      = GetCachedVar( (TString)fDm21->GetName() );
  Double_t cache_dm31      = GetCachedVar( (TString)fDm31->GetName() );

  if (SinsqTh12 != cache_sinsqth12 || SinsqTh13 != cache_sinsqth13 || SinsqTh23 != cache_sinsqth23 || 
      Dcp != cache_dcp || Dm21 != cache_dm21 || Dm31 != cache_dm31) {

    // update the cache values
    fParCache[ (TString)fSinsqTh12->GetName() ] = SinsqTh12;
    fParCache[ (TString)fSinsqTh13->GetName() ] = SinsqTh13;
    fParCache[ (TString)fSinsqTh23->GetName() ] = SinsqTh23;
    fParCache[ (TString)fDcp->GetName() ]  = Dcp;
    fParCache[ (TString)fDm21->GetName() ] = Dm21;
    fParCache[ (TString)fDm31->GetName() ] = Dm31;
    
    // give variables to the oscillator
    fProb->SetAngle(1, 2, TMath::ASin( TMath::Sqrt( SinsqTh12 ) ) );
    fProb->SetAngle(1, 3, TMath::ASin( TMath::Sqrt( SinsqTh13 ) ) );
    fProb->SetAngle(2, 3, TMath::ASin( TMath::Sqrt( SinsqTh23 ) ) );
    fProb->SetDelta(1, 3, Dcp * TMath::Pi()  );
    fProb->SetDm(2, Dm21);
    fProb->SetDm(3, Dm31);

    // set flag that re-calculation should be performed
    reCalc = kTRUE;

  }
  else { reCalc = kFALSE; }

  return reCalc;

}

//***************************************************************************

/**
   Private function that handles the caching of the oscillation probabilities.

   This function re-calculates the oscillation probabilites and stores them in member histograms `fOscCache` whenever any of the oscillation parameters change. If none of them change, the calculation is not performed. The way the `DetResponse` class is set up and the way fitting works, this enables to save a large amount of time.

   NB! The oscillations tau-> (elec, muon, tau) are not calculated to save time; tau flux is assumed 0.

   \param proxymap    A proxy map with `RooFit` variables from `FitPDF` that contain all shared observables and parameters, including the oscillation parameters. 
   \param nsamples    Number of samples per bin
 */
void FitUtil::ProbCacher(const proxymap_t &proxymap, UInt_t nsamples) {

  if (nsamples < 1) {
    throw std::invalid_argument("ERROR! FitUtil::ProbCacher() cannot have below 1 samples per bin.");
  }
  
  // call function that configures OscProb and returns flag to indicate re-calculation
  Bool_t reCalc = ConfigOscProb(proxymap);

  // also recalculate if number of requested samples changes
  reCalc = reCalc || (nsamples != f_cache_oscsamplesn);

  // osc parameters have not changed, no re-calculation required, return
  if ( !reCalc ) return;

  //=======================================================================
  // recalculation of the cache 
  //=======================================================================

  fOscCalcTime->Start(kFALSE);

  // update the variable for samples; osc. par. cache handled by ConfigOscProb
  f_cache_oscsamplesn = nsamples;

  // calculate all oscillation probabilities and cache them  - NB! TAU-IN IS NOT CALCULATED!
  for (UInt_t f_in = 0; f_in < fFlavs-1; f_in++) {
    for (UInt_t f_out = 0; f_out < fFlavs; f_out++) {
      for (UInt_t isnb = 0; isnb < fInts; isnb++) {

	fProb->SetIsNuBar( isnb );

	Int_t N_ebins  = fHBT->GetXaxis()->GetNbins();
	Int_t N_ctbins = fHBT->GetYaxis()->GetNbins();

	for (Int_t ebin = 1; ebin <= N_ebins; ebin++) {
	  for (Int_t ctbin = 1; ctbin <= N_ctbins; ctbin++) {

	    Double_t prob_ave = 0;

	    Double_t emin    = fHBT->GetXaxis()->GetBinLowEdge(ebin);
	    Double_t emax    = fHBT->GetXaxis()->GetBinUpEdge(ebin);
	    Double_t e_step  = (emax - emin)/(nsamples+1);

	    Double_t ctmin   = fHBT->GetYaxis()->GetBinLowEdge(ctbin);
	    Double_t ctmax   = fHBT->GetYaxis()->GetBinUpEdge(ctbin);
	    Double_t ct_step = (ctmax - ctmin)/(nsamples+1);

	    // calculate average of N samples inside the bin
	    for (UInt_t n_e = 0; n_e < nsamples; n_e++) {

	      Double_t E  = emin + (n_e + 1) * e_step;

	      for (UInt_t n_ct = 0; n_ct < nsamples; n_ct++) {

		Double_t ct = ctmin + (n_ct + 1) * ct_step;

		fPrem->FillPath( ct );	      
		fProb->SetPath ( fPrem->GetNuPath() );

		prob_ave += fProb->Prob(f_in, f_out, E);

		fOscCalls++;
		
	      } // end loop over cos-theta samples
	    } // end loop over energy samples

	    fOscCache[f_in][f_out][isnb][ebin][ctbin] = std::make_pair( prob_ave/(nsamples*nsamples), 0. );

	  } // end loop over cos-theta bins
	} //end loop over energy bins

      } // end loop over isnb
    } // end loop over flavor_out
  } // end loop over flavor_in

  fOscCalcTime->Stop();
      
}

//***************************************************************************

/** Function that calculates the number of events in a \f$ (E_{\rm true}, cos\theta_{\rm true}, by_{\rm true}) \f$ bin.

    \param tb          A `TrueB` object (see `DetResponse.h`) that stores the neutrino type and true bin coordinate info.
    \param proxymap    A proxy map with `RooFit` variables from `FitPDF` that contain all shared observables and parameters, including the oscillation parameters.
    \return            A pair; first is the number of expected events in the true bin, second is the MC statistical uncertainty
*/
std::pair<Double_t, Double_t> FitUtil::TrueEvts(const TrueB &tb, const proxymap_t &proxymap) {
  
  // get the atm nu count
  Double_t atm_count_e = GetCachedFlux(ELEC, tb.fIsNB, tb.fE_true_bin, tb.fCt_true_bin);
  Double_t atm_count_m = GetCachedFlux(MUON, tb.fIsNB, tb.fE_true_bin, tb.fCt_true_bin);
  
  // get the oscillation probabilities
  Double_t prob_elec = GetCachedOsc(ELEC, tb, proxymap);
  Double_t prob_muon = GetCachedOsc(MUON, tb, proxymap);
  
  // get the oscillated nu count in operation time (in units 1/m2)
  Double_t osc_count = ( atm_count_e * prob_elec + atm_count_m * prob_muon );

  // get the interacted neutrino count in operation time (in units 1/MTon)
  Double_t int_count = osc_count * GetCachedXsec(tb)/fMN * fKg_per_MTon;
  
  // get the effective mass
  Double_t meff = GetCachedMeff(tb);

  // calculate the number of detected events (unitless)
  Double_t det_count = int_count * GetCachedBYfrac(tb) * meff * 1e-6;

  // MC stat err, coming from eff mass and BY distribution (0 for now)
  Double_t det_err = 0.;

  return std::make_pair(det_count, det_err);

}

//***************************************************************************

/** Function to retrieve the cached `TrueEvts` return value.
    \param tb  A `TrueB` object (see `DetResponse.h`) that stores the neurino type and true bin coordinate infp
    \return    A pair, first is the number of expected events in true bin, second is MC error
 */
std::pair<Double_t, Double_t> FitUtil::GetCachedTE(const TrueB &tb) {
  
  return fTECache[tb.fFlav][tb.fIsCC][tb.fIsNB][tb.fE_true_bin][tb.fCt_true_bin][tb.fBy_true_bin];
  
}

//***************************************************************************

/** Function that re-calculates the cached values of `TrueEvts`.

    This caching saves a lot of time. Firstly, for a single response (e.g. tracks), one true bin contributes to several reco bins. Without this cache, the function `RecoEvts` would call `TrueEvts` tens of times (for each reco bin that the true bin contributed to) when the fitter is looping over reco bins, while the other fit parameters are unchanged. Secondly, the true bin data is also the same for several event selections. This means that once this cache is filled for some parameter values, it is shared between several selections (e.g. tracks and showers) that share the `FitUtil` class exactly for such central caching functionality.

    \param proxymap A structure with all fit parameters known to RooFit.
*/
void FitUtil::FillTECache(const proxymap_t& proxymap) {

  // loop over flavors, interaction types and polarisations
  for (UInt_t f = 0; f < fFlavs; f++) {
    for (UInt_t i = 0; i < fInts; i++) {
      for (UInt_t p = 0; p < fPols; p++) {

	TH3D* hc = fHBT;

	// loop over true bins
	for (Int_t xbin = 1; xbin <= hc->GetXaxis()->GetNbins(); xbin++) {
	  for (Int_t ybin = 1; ybin <= hc->GetYaxis()->GetNbins(); ybin++) {
	    for (Int_t zbin = 1; zbin <= hc->GetZaxis()->GetNbins(); zbin++) {

	      TrueB _tb(f, i, p, xbin, ybin, zbin);
	      fTECache[f][i][p][xbin][ybin][zbin] = TrueEvts(_tb, proxymap);
	      
	    }
	  }
	}
	  
      }
    }
  }
  
}

//***************************************************************************

/** This function is called inside `FitPDF::analyticalIntegral` and `FitPDF::GetExpValHist` and fills a histogram with expectation values in reco bins.

    The argument map is created in `FitPDF` and contains names and corresponding proxies for all parameters in `fParSet`. The detector response is also part of the `FitPDF` class and configures what kind of an event selection the pdf is used to fit.
    
    \param resp       Pointer to `AbsResponse` instance used with the `FitPDF` class.
    \param proxymap   A proxy map with `RooFit` variables from `FitPDF` that contain all shared observables and parameters, including the oscillation parameters. 
    \param rangeName  Range string as used in `RooFit`, it will apply the defined rangeName to the E/Ct/By variables. If the range does not exists, the default ranges are used.
    \return           a 3D histogram in reco variables with exepectation values as bin contents

*/
TH3D* FitUtil::Expectation(AbsResponse *resp, const proxymap_t &proxymap, const char* rangeName) {

  // Default behaviour of getMin/getMax if there is no Range present: return the absolute min/max of 
  // the data type. For doubles: +-1E30
  // RooFit handles ranges on different RooRealVars with the same rangeName properly: they do not collide.
  // So using the same name for a collection of ranges on different variables is valid.
  Double_t E_min  = ((RooRealVar*)fParSet.find("E_reco")) ->getMin(rangeName);
  Double_t E_max  = ((RooRealVar*)fParSet.find("E_reco")) ->getMax(rangeName);
  Double_t ct_min = ((RooRealVar*)fParSet.find("ct_reco"))->getMin(rangeName);
  Double_t ct_max = ((RooRealVar*)fParSet.find("ct_reco"))->getMax(rangeName);
  Double_t by_min = ((RooRealVar*)fParSet.find("by_reco"))->getMin(rangeName);
  Double_t by_max = ((RooRealVar*)fParSet.find("by_reco"))->getMax(rangeName);

  // create the histogram with expectation values
  TH3D   *hexp  = (TH3D*)resp->GetHist3DReco()->Clone();
  TString hname = resp->GetRespName() + "_expct";
  hexp->SetDirectory(0);
  hexp->Reset();
  hexp->SetNameTitle(hname, hname);
  
  auto XFitRange = GetRange( E_min,  E_max,  hexp->GetXaxis() );
  auto YFitRange = GetRange( ct_min, ct_max, hexp->GetYaxis() );
  auto ZFitRange = GetRange( by_min, by_max, hexp->GetZaxis() );

  // Do not overwrite the default fEbin_min, etc. since they are used by all instances that call FitUtil!
  // Instead define new Ints to use in the for loop downstairs
  // Take the max to move the minimum value up and the min to move the maximum value down, this 
  // constrains the range.
  Int_t ebin_min  = max( fEbin_min,  std::get<MINBIN>(XFitRange) );
  Int_t ebin_max  = min( fEbin_max,  std::get<MAXBIN>(XFitRange) );
  Int_t ctbin_min = max( fCtbin_min, std::get<MINBIN>(YFitRange) );
  Int_t ctbin_max = min( fCtbin_max, std::get<MAXBIN>(YFitRange) );
  Int_t bybin_min = max( fBybin_min, std::get<MINBIN>(ZFitRange) );
  Int_t bybin_max = min( fBybin_max, std::get<MAXBIN>(ZFitRange) );

  // loop over bins and fill the expectation value histogram
  for (Int_t ebin = ebin_min; ebin <= ebin_max; ebin++) {
    for (Int_t ctbin = ctbin_min; ctbin <= ctbin_max; ctbin++) {
      for (Int_t bybin = bybin_min; bybin <= bybin_max; bybin++) {

        Double_t E  = hexp->GetXaxis()->GetBinCenter( ebin );
        Double_t ct = hexp->GetYaxis()->GetBinCenter( ctbin );
        Double_t by = hexp->GetZaxis()->GetBinCenter( bybin );

        // calculate the bin width, necessary to convert event density to event count
        Double_t e_w  = hexp->GetXaxis()->GetBinWidth( ebin  );
        Double_t ct_w = hexp->GetYaxis()->GetBinWidth( ctbin );
        Double_t by_w = hexp->GetZaxis()->GetBinWidth( bybin );
        Double_t binw = e_w * ct_w * by_w;

        auto recoevts = RecoEvts(E, ct, by, resp, proxymap);

	//cout<<"Expectation reco value in bin: "<<recoevts.first*binw<<endl;

        hexp->SetBinContent(ebin, ctbin, bybin, recoevts.first*binw  );
        hexp->SetBinError  (ebin, ctbin, bybin, recoevts.second*binw );

      }
    }
  }
  
  return hexp;
    
}

//***************************************************************************

/** This function is called from inside `FitPDF::evaluate()` and returns the expected event density at a \f$ (E_{\rm reco}, cos\theta_{\rm reco}, by_{\rm reco}) \f$ value.

    Note that the return of this function is normalised by RooFit by dividing with \f$ \int f(\vec{x})d\vec{x} \f$, where f(x) is the value returned by this function and \f$ \vec{x} \f$ are energy, cos-theta and bjorken-y.

    The argument map is created in `FitPDF` and contains names and corresponding proxies for all parameters in `fParSet`. The detector response is also part of the `FitPDF` class and configures what kind of an event selection the pdf is used to fit.

    Whereas the fit parameters (e.g. oscillation parameters) are successfully passed through the contents of the parameter `proxymap`, the observables reco energy, cos-theta and bjorken-y need to be passed explicitly. This is necessary, as otherwise the function `FitUtil::Expectation` would need to set the values for `fE_reco`, `fCt_reco` and `fBy_reco` manually in the loop over the bins, and this interferes in some way with `RooFit` internal cache for these variables, such that the fitting procedure gets interrupted and does not lead to convergence.

    \param E_reco      Reco energy
    \param Ct_reco     Reco cos-theta
    \param By_reco     Reco bjorken-y
    \param resp        Pointer to `AbsResponse` instance used with the `FitPDF` class.
    \param proxymap    A proxy map with `RooFit` variables from `FitPDF` that contains the shared fit parameters, including oscillation parameters.
    \return            a pair with the un-normalised event density (calculated by dividing the expected number of events in a bin by bin width) and the associated statistical uncertainty */
std::pair<Double_t, Double_t> FitUtil::RecoEvts(Double_t E_reco, Double_t Ct_reco, Double_t By_reco, AbsResponse *resp, const proxymap_t &proxymap) {

  //----------------------------------------------------------------------------------
  // check that observables are within the range of the detector response
  //----------------------------------------------------------------------------------

  Bool_t throwException = kFALSE;
  std::invalid_argument exception("");

  if ( E_reco < fHBR->GetXaxis()->GetXmin() || E_reco >= fHBR->GetXaxis()->GetXmax() ) {
    throwException = kTRUE;
    exception = std::invalid_argument("ERROR! FitUtil::RecoEvts energy " + to_string(E_reco) + " outside the binning range.");
  }

  if ( Ct_reco < fHBR->GetYaxis()->GetXmin() || Ct_reco >= fHBR->GetYaxis()->GetXmax() ) {
    throwException = kTRUE;
    exception = std::invalid_argument("ERROR! FitUtil::RecoEvts cos-theta " + to_string(Ct_reco) + " outside the binninb range.");
  }

  if ( By_reco < fHBR->GetZaxis()->GetXmin() || By_reco >= fHBR->GetZaxis()->GetXmax() ) {
    throwException = kTRUE;
    exception = std::invalid_argument("ERROR! FitUtil::RecoEvts bjorken-y " + to_string(By_reco) + " outside the binning range.");
  }

  if (throwException) {
    
    cout << "=====================================================================================" << endl;
    cout << "NOTICE FitUtil::RecoEvts() throwing exception, parameter value dump at exception: " << endl;
    for (auto p: proxymap) { cout << "Parameter: " << p.first << "\t Value: " << (Double_t)(*p.second) << endl; }
    cout << "=====================================================================================" << endl;

    throw exception;
  }  

  Double_t det_count = 0.;
  Double_t det_err   = 0.;

  //----------------------------------------------------------------------------------
  // perform the BinnedResponse calculation
  //----------------------------------------------------------------------------------
    
  if ( resp->GetResponseType() == AbsResponse::BinnedResponse ) {

    //----------------------------------------------------------------------------------
    // if any of the fit parameters has changed re-calculate the TrueEvts cache
    //----------------------------------------------------------------------------------

    if ( CheckVarCache(proxymap) ) {
      FillTECache(proxymap);
      UpdateVarCache(proxymap); // this function makes sure that the cache values of all variables is updated
    }
    
    auto true_bins = ((DetResponse*)resp)->GetBinWeights( E_reco, Ct_reco, By_reco );
	  
    for (const auto &tb: true_bins) {
	    
      if (tb.fIsCC) {
	      
	Double_t TE = GetCachedTE(tb).first;

	det_count += tb.fW * TE;
	det_err   += TMath::Power(tb.fWE * TE, 2);
	      
      }
      else {

	Double_t TE = 0.;

	/* for NC events the response only says how the elec-NC events from the considered true bin contribute to the reco-bin. However, in the flux chain we have we also have a contribution from muon-NC and tau-NC true bin to the reco bin, which are assumed to look identical to the detector as elec-NC (this is why we only simulate elec-NC). Hence for NC events I need to add `TrueEvts` contributions from the three NC flavours */

	TrueB elecTB( tb );
	TrueB muonTB( tb );
	TrueB tauTB ( tb );
	elecTB.fFlav = ELEC;
	muonTB.fFlav = MUON;
	tauTB.fFlav  = TAU;
	      
	TE += GetCachedTE(elecTB).first;
	TE += GetCachedTE(muonTB).first;
	TE += GetCachedTE(tauTB).first;
	      
	det_count += tb.fW * TE;
	det_err   += TMath::Power(tb.fWE * TE, 2);

      }

    } // end loop over true bins
    
    det_err = TMath::Sqrt(det_err);
    
  }

  //----------------------------------------------------------------------------------
  // perform the EvtResponse calculation
  //----------------------------------------------------------------------------------

  else if ( resp->GetResponseType() == AbsResponse::EvtResponse ) {
  
    auto true_evts = ((EvtResponse*)resp)->GetBinEvts( E_reco, Ct_reco, By_reco );

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
      Double_t atm_count_e = GetCachedFlux(ELEC, te.GetIsNB(), e_bin, ct_bin) / flux_conv;
      Double_t atm_count_m = GetCachedFlux(MUON, te.GetIsNB(), e_bin, ct_bin) / flux_conv;

      if ( te.GetIsCC() ) {

	//create trueB object to access cached oscillation probabilities
	TrueB TB( te.GetFlav(), te.GetIsCC(), te.GetIsNB(), e_bin, ct_bin, by_bin);
	
	Double_t prob_elec = GetCachedOsc(ELEC, TB, proxymap);
	Double_t prob_muon = GetCachedOsc(MUON, TB, proxymap);
		   
	// calculate the oscillated differential flux for the neutrino type
	Double_t oscf = atm_count_e*prob_elec + atm_count_m*prob_muon;

	det_count += te.GetW1y() * oscf;
	det_err++;
		      
      }
      else {

	// for NC, the xsec and meff for all flavors are the ~same. For this reason, only elec-NC are simulated.
	// For this calculation, however, the NC flux of muons and tau's also needs to be taken into account.
	// The calculation for CC events after oscillation is osc_flux * xsec_{flavor} * meff_{flavor}. As
	// xsec_{flavor} and meff_{flavor} are the same for elec, muon, tau for NC events, the calculation
	// can be simplified and the oscillated flux is just equal to the un-oscillated elec+muon flux

	Double_t oscf = atm_count_e + atm_count_m;
	det_count += te.GetW1y() * oscf;
	det_err++;
	
      }

    } // end loop over true events

    det_err = TMath::Sqrt(det_err)/det_err * det_count;

  //----------------------------------------------------------------------------------
  // unknown response type, throw error
  //----------------------------------------------------------------------------------
    
  }
  else {
    throw std::invalid_argument("ERROR! FitUtil::RecoEvts() unknown response type " + (string)typeid( *resp ).name() );
  }
  
  // should you wish to add atm muons and noise at some point
  //det_count += resp->GetAtmMuCount1y(E_reco, ct_reco, by_reco) * fOpTime;
  //det_count += resp->GetNoiseCount1y(E_reco, ct_reco, by_reco) * fOpTime;

  // finally, convert the event counts in reco bin to event density
  Double_t e_w  = fHBR->GetXaxis()->GetBinWidth( fHBR->GetXaxis()->FindBin( E_reco )  );
  Double_t ct_w = fHBR->GetYaxis()->GetBinWidth( fHBR->GetYaxis()->FindBin( Ct_reco ) );
  Double_t by_w = fHBR->GetZaxis()->GetBinWidth( fHBR->GetZaxis()->FindBin( By_reco ) );
  Double_t binw = e_w * ct_w * by_w;
  
  return std::make_pair(det_count/binw, det_err/binw);

}

//***************************************************************************

/** Public function to fetch a pointer to a variable in the member `fParSet`.

    This function knows of all the variables known to `RooFit`. If the variable is not found set an exception is thrown.

    \param varname   Name of the variable
    \return          Pointer to the variable
 */
RooRealVar* FitUtil::GetVar(TString varname) {

  RooRealVar *var = (RooRealVar*)fParSet.find(varname);

  if (var == NULL) {
    throw std::invalid_argument("ERROR! FitUtil::GetVar() cannot find variable " + (string)varname);
  }
    
  return var;

}

//***************************************************************************
 
/** Function to fetch the cached value of a variable, which is stored in `FitUtil::fParCache`.

    If no cache for the variable exists, the cache is created for the input variable name with value 0.

    \param varname Name of the variable
    \return        const reference to the cached variable
*/
const Double_t& FitUtil::GetCachedVar(TString varname) {

   // search for the element in the map
   auto it = fParCache.find( varname );

   // if element is not in the map create an entry for it in the cache
   // otherwise return the pointer to the cached variable
  
   if ( it == fParCache.end() ) {

     cout << "NOTICE FitUtil::GetCachedVar() started caching variable " << varname << endl;
     fParCache.insert( std::make_pair( varname, 0. ) );
     return fParCache[ varname ];
    
   }
   else {
    
     return it->second;

   }
  
 }
 
//***************************************************************************

/** Function that checks whether any of the parameters in the input proxymap has changed with respect to the values stored in the cache `FitUtil::fParCache`.
    \param proxymap  structure with all fit parameters known to `RooFit`.
    \return          true if one or more parameter(s) in the cache is/are different than in the proxymap.
*/
Bool_t FitUtil::CheckVarCache(const proxymap_t& proxymap) {

  // loop over all parameters in the proxymap (all parameters known to RooFit)
  for (auto p: proxymap) {

    // ignore observables (E, ct, by), we don't wish to cash these
    if ( fObsList.find( p.first ) != NULL ) { continue; }

    Double_t cache = GetCachedVar( p.first );  // value from the cache
    Double_t  val  = (Double_t)(*p.second);    // current value of the parameter
    
    if ( val != cache ) return kTRUE;          // when any parameter mismatches the cache return true
    
  }

  return kFALSE;
  
}

//***************************************************************************

/** Function to update all of the cached parameter values in `FitUtil::fParCache`.

    If `FitUtil::fParCache` does not contain an entry for one/some of the parameters in proxymap, the cache is created. As an exception, observables (E,ct,by) are not cached.

    \param proxymap structure with all fit parameters known to `RooFit`.
 */
void FitUtil::UpdateVarCache(const proxymap_t& proxymap) {

  // loop over all parameters in the proxymap (all parameters known to RooFit)
  for (auto p: proxymap) {

    // ignore observables (E, ct, by), we don't wish to cache these
    if ( fObsList.find( p.first ) != NULL ) { continue; }

    Double_t val = (Double_t)(*p.second);    // current value of the parameter

    // find the parameter in the par cache map and update it's value. If the parameter is not (yet)
    // in the map insert it. 
    auto it = fParCache.find( p.first );
    if ( it == fParCache.end() ) { 
      cout << "NOTICE FitUtil::UpdateVarCache() started caching variable " << p.first << endl;
      fParCache.insert( std::make_pair( p.first, val ) ); 
    }
    else { 
      it->second = val; 
    }
    
  }
    
}

//***************************************************************************

/** Function to initialise a 1D cache structure.
    
    The indexing of the cache bins is such that it matches the indices of the input axis bins, i.e. bin[0] is undeflow and bin[1] is the first counting bin. 

    \param  axis Pointer to a `TAxis` object from which the cache array dimesions are calculated
    \return Pointer to a 1D array of pairs where values can be cached, first is meant for value, second is meant for error
*/
cache1D_t FitUtil::InitCache1D(TAxis *axis) {

  // bin[0] is underflow, bin[axis->GetNbins+1] is overflow, hence +2
  Int_t nbins = axis->GetNbins() + 2;
  
  cache1D_t cache = new std::pair<Double_t, Double_t> [ nbins ];

  for (Int_t bin = 0; bin < nbins; bin++) {
    cache[bin].first  = 0.;
    cache[bin].second = 0.;
  }

  return cache;
  
}

//***************************************************************************

/** Function to initialise a 2D cache structure.
    
    The indexing of the cache bins is such that it matches the indices of the input axis bins, i.e. bin[0] is undeflow and bin[1] is the first counting bin. 

    \param  xaxis Pointer to a `TAxis` object from which the cache array dimesions are calculated
    \param  yaxis Pointer to a `TAxis` object from which the cache array dimesions are calculated
    \return Pointer to a 2D array of pairs where values can be cached, first is meant for value, second is meant for error
*/
cache2D_t FitUtil::InitCache2D(TAxis *xaxis, TAxis *yaxis) {

  // bin[0] is underflow, bin[axis->GetNbins+1] is overflow, hence +2
  Int_t nxbins = xaxis->GetNbins() + 2;
  
  cache2D_t cache = new std::pair<Double_t, Double_t>* [ nxbins ];

  for (Int_t xbin = 0; xbin < nxbins; xbin++) {
    cache[xbin] = InitCache1D( yaxis );
  }

  return cache;
  
}

//***************************************************************************

/** Function to initialise a 3D cache structure.
    
    The indexing of the cache bins is such that it matches the indices of the input axis bins, i.e. bin[0] is undeflow and bin[1] is the first counting bin. 

    \param  xaxis Pointer to a `TAxis` object from which the cache array dimesions are calculated
    \param  yaxis Pointer to a `TAxis` object from which the cache array dimesions are calculated
    \param  zaxis Pointer to a `TAxis` object from which the cache array dimesions are calculated
    \return Pointer to a 3D array of pairs where values can be cached, first is meant for value, second is meant for error
*/
cache3D_t FitUtil::InitCache3D(TAxis *xaxis, TAxis *yaxis, TAxis *zaxis) {

  // bin[0] is underflow, bin[axis->GetNbins+1] is overflow, hence +2
  Int_t nxbins = xaxis->GetNbins() + 2;
  
  cache3D_t cache = new std::pair<Double_t, Double_t>** [ nxbins ];
      
  for (Int_t xbin = 0; xbin < nxbins; xbin++) { 
    cache[xbin] = InitCache2D( yaxis, zaxis );
  }

  return cache;
}

//***************************************************************************

/** Function to clear a dynamically allocated 1D cache array
    \param cache  Pointer to the dynamically allocated 1D cache array
*/
void FitUtil::ClearCache1D(cache1D_t cache) {
  
  delete [] cache;
  
}

//***************************************************************************

/** Function to clear a dynamically allocated 2D cache array
    \param cache  Pointer to the dynamically allocated 2D cache array
    \param xaxis  Pointer to a `TAxis` object that was used to initialise the xaxis of the 2D cache array
*/
void FitUtil::ClearCache2D(cache2D_t cache, TAxis *xaxis) {

  Int_t nxbins = xaxis->GetNbins() + 2;
  
  for (Int_t xbin = 0; xbin < nxbins; xbin++) {
    delete [] cache[xbin];
  }
  delete [] cache;
  
}

//***************************************************************************

/** Function to clear a dynamically allocated 3D cache array
    \param cache  Pointer to the dynamically allocated 3D cache array
    \param xaxis  Pointer to a `TAxis` object that was used to initialise the xaxis of the 3D cache array
    \param yaxis  Pointer to a `TAxis` object that was used to initialise the yaxis of the 3D cache array
*/
void FitUtil::ClearCache3D(cache3D_t cache, TAxis *xaxis, TAxis *yaxis) {

  Int_t nxbins = xaxis->GetNbins() + 2;
  Int_t nybins = yaxis->GetNbins() + 2;
  
  for (Int_t xbin = 0; xbin < nxbins; xbin++) {
    for (Int_t ybin = 0; ybin < nybins; ybin++) {
      delete[] cache[xbin][ybin];
    }
    delete[] cache[xbin];
  }
  delete [] cache;
	
}
