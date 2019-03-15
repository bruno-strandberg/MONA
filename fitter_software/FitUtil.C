#include "FitUtil.h"
#include <iostream>

#include "TFile.h"

using namespace std;

/** Constructor
    \param op_time         Operation time in years
    \param h_template      A 3D histogram that represents the binning settings
    \param emin            Minimum energy included in the fit range
    \param emax            Maximum energy included in the fit range
    \param ctmin           Minimum cos-theta included in the fit range
    \param ctmax           Maximum cos-theta included in the fit range
    \param bymin           Minimum bjorken-y included in the fit range
    \param bymax           Maximum bjorken-y included in the fit range
    \param meff_file       File with effective mass data (created with `EffMass` class)
*/
FitUtil::FitUtil(Double_t op_time, TH3 *h_template,
		 Double_t emin, Double_t emax, Double_t ctmin, Double_t ctmax, Double_t bymin, Double_t bymax,
		 TString meff_file) {

  fOpTime = op_time;

  fHB = (TH3D*)h_template->Clone("HB");
  fHB->SetDirectory(0);
  fHB->Reset();
  fHB->SetNameTitle("HB", "HB");

  fFlux = new AtmFlux;
  fXsec = new NuXsec( fHB->GetZaxis()->GetNbins() );
  fMeff = new EffMass( meff_file, fHB->GetXaxis()->GetNbins(), fHB->GetYaxis()->GetNbins(), fHB->GetZaxis()->GetNbins() );
  fProb = new OscProb::PMNS_Fast;
  fPrem = new OscProb::PremModel;

  // calculate the observable ranges to match exactly bin edges, set ranges for integration
  auto ENr = GetRange( emin , emax , fHB->GetXaxis() );
  auto CTr = GetRange( ctmin, ctmax, fHB->GetYaxis() );
  auto BYr = GetRange( bymin, bymax, fHB->GetZaxis() );

  fEbin_min  = get<MINBIN>(ENr);
  fEbin_max  = get<MAXBIN>(ENr);
  fCtbin_min = get<MINBIN>(CTr);
  fCtbin_max = get<MAXBIN>(CTr);
  fBybin_min = get<MINBIN>(BYr);
  fBybin_max = get<MAXBIN>(BYr);
  
  InitFitVars(get<MIN>(ENr), get<MAX>(ENr), get<MIN>(CTr), get<MAX>(CTr), get<MIN>(BYr), get<MAX>(BYr));
  InitCacheHists(fHB);
  Fill_Flux_Xsec_Meff_cache(fFlux, fXsec, fMeff, fOpTime);

  f_cache_sinsqth12 = 0;
  f_cache_sinsqth13 = 0;
  f_cache_sinsqth23 = 0;
  f_cache_dcp       = 0;
  f_cache_dm21      = 0;
  f_cache_dm31      = 0;

  fOscCalls    = 0;
  fOscCalcTime = new TStopwatch();
  fOscCalcTime->Stop(); fOscCalcTime->Reset();

}

//***************************************************************************

/** Destructor */
FitUtil::~FitUtil() {

  cout << "FitUtil::~FitUtil() total oscillator calls: " << fOscCalls << endl;
  cout << "FitUtil::~FitUtil() duration of oscillation calculations [s]: "
       << (Double_t)fOscCalcTime->RealTime() << endl;

  if (fOscCalcTime) delete fOscCalcTime;

  if (fHB)   delete fHB;
  if (fFlux) delete fFlux;
  if (fXsec) delete fXsec;
  if (fMeff) delete fMeff;
  if (fProb) delete fProb;
  if (fPrem) delete fPrem;

  // remove cache hists
  for (UInt_t f = 0; f < fFlavs; f++) {
    for (UInt_t isnb = 0; isnb < fPols; isnb++) {

      // osc
      for (UInt_t fout = 0; fout < fFlavs; fout++) {
	if ( fhOscCache[f][fout][isnb] ) delete fhOscCache[f][fout][isnb];
      }

      // flux
      if ( fhFluxCache[f][isnb] ) delete fhFluxCache[f][isnb];

      // xsec and meff
      for (UInt_t iscc = 0; iscc < fInts; iscc++) {
	if ( fhXsecCache[f][iscc][isnb] ) delete fhXsecCache[f][iscc][isnb];
	if ( fhMeffCache[f][iscc][isnb] ) delete fhMeffCache[f][iscc][isnb];
	if ( fhBYfracCache[f][iscc][isnb] ) delete fhBYfracCache[f][iscc][isnb];
      }

    }
  }

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

    \param min    Minimum of the range
    \param max    Maximum of the range
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
    bin_max = axis->FindBin(_max) - 1;
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
    \param bymin  Maximum bjorken-y in the fit range
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

  // maps to create histogram names
  std::map<UInt_t, TString> pol_map  = { {0, "nu"}, {1, "nub"} }; 
  std::map<UInt_t, TString> iscc_map = { {0, "NC"}, {1, "CC"} }; 
  std::map<UInt_t, TString> flav_map = { {0, "elec"}, {1, "muon"}, {2, "tau"} };

  // project the 3D template to 1D and 2D templates
  TH2D *h_template_EvsCT = (TH2D*)h_template->Project3D("yx");
  TH2D *h_template_EvsBY = (TH2D*)h_template->Project3D("zx");
  TH1D *h_template_E     = (TH1D*)h_template->Project3D("x");

  // init the oscillation probability cache histograms
  for (UInt_t f_in = 0; f_in < fFlavs; f_in++) {
    for (UInt_t f_out = 0; f_out < fFlavs; f_out++) {
      for (UInt_t isnb = 0; isnb < fPols; isnb++) {

	TString name = "oscprob_" + flav_map[f_in] + "_to_" + flav_map[f_out] + "_" + pol_map[isnb];

        fhOscCache[f_in][f_out][isnb] = (TH2D*)h_template_EvsCT->Clone();
        fhOscCache[f_in][f_out][isnb]->SetDirectory(0);
        fhOscCache[f_in][f_out][isnb]->Reset();
        fhOscCache[f_in][f_out][isnb]->SetNameTitle(name, name);

      }
    }
  }

  // init the flux histograms
  for (UInt_t f = 0; f < fFlavs; f++) {
    for (UInt_t isnb = 0; isnb < fPols; isnb++) {

      TString name = "flux_" + flav_map[f] + "_" + pol_map[isnb];

      fhFluxCache[f][isnb] = (TH2D*)h_template_EvsCT->Clone();
      fhFluxCache[f][isnb]->SetDirectory(0);
      fhFluxCache[f][isnb]->Reset();
      fhFluxCache[f][isnb]->SetNameTitle(name, name);

    }
  }

  // init the xsec histograms and meff histograms
  for (UInt_t f = 0; f < fFlavs; f++) {
    for (UInt_t iscc = 0; iscc < fInts; iscc++) {
      for (UInt_t isnb = 0; isnb < fPols; isnb++) {

        TString name_xsec   = "xsec_"   + flav_map[f] + "_" + iscc_map[iscc] + "_" + pol_map[isnb];
	TString name_byfrac = "byfrac_" + flav_map[f] + "_" + iscc_map[iscc] + "_" + pol_map[isnb];
	TString name_meff   = "meff_"   + flav_map[f] + "_" + iscc_map[iscc] + "_" + pol_map[isnb];

        fhXsecCache[f][iscc][isnb] = (TH1D*)h_template_E->Clone();
        fhXsecCache[f][iscc][isnb]->SetDirectory(0);
        fhXsecCache[f][iscc][isnb]->Reset();
        fhXsecCache[f][iscc][isnb]->SetNameTitle(name_xsec, name_xsec);

	fhBYfracCache[f][iscc][isnb] = (TH2D*)h_template_EvsBY->Clone();
	fhBYfracCache[f][iscc][isnb]->SetDirectory(0);
	fhBYfracCache[f][iscc][isnb]->Reset();
	fhBYfracCache[f][iscc][isnb]->SetNameTitle(name_byfrac, name_byfrac);
	
	fhMeffCache[f][iscc][isnb] = (TH3D*)h_template->Clone();
	fhMeffCache[f][iscc][isnb]->SetDirectory(0);
	fhMeffCache[f][iscc][isnb]->Reset();
	fhMeffCache[f][iscc][isnb]->SetNameTitle(name_meff, name_meff);
	
      }
    }
  }

}

//***************************************************************************

/** A private function to fill the flux and xsec cache hists
    \param flux     Pointer to the `AtmFlux` member instance
    \param xsec     Pointer to the `NuXsec` member instance
    \param meff     Pointer to the `EffMass` member instance
    \param op_time  Operation time in years.
*/
void FitUtil::Fill_Flux_Xsec_Meff_cache(AtmFlux *flux, NuXsec *xsec, EffMass *meff, Double_t op_time) {

  // fill the atm flux cache. Note that (currently) AtmFlux returns 0 for tau flux
  for (UInt_t f = 0; f < fFlavs; f++) {
    for (UInt_t isnb = 0; isnb < fPols; isnb++) {

      TH2D* hflux = fhFluxCache[f][isnb];

      for (Int_t ebin = 1; ebin <= hflux->GetXaxis()->GetNbins(); ebin++) {
	for (Int_t ctbin = 1; ctbin <= hflux->GetYaxis()->GetNbins(); ctbin++) {

	  Double_t E   = hflux->GetXaxis()->GetBinCenter(ebin);
	  Double_t ew  = hflux->GetXaxis()->GetBinWidth(ebin);
	  Double_t ct  = hflux->GetYaxis()->GetBinCenter(ctbin);
	  Double_t ctw = hflux->GetYaxis()->GetBinWidth(ctbin);

	  Double_t atm_flux_factor = ew * ctw * fSec_per_y * op_time;
	  Double_t atmflux = flux->Flux_dE_dcosz(f, isnb, E, ct) * atm_flux_factor;
	  
	  hflux->SetBinContent(ebin, ctbin, atmflux);

	}
      }

    }
  }

  // fill xsec cache and meff cache
  for (UInt_t f = 0; f < fFlavs; f++) {
    for (UInt_t iscc = 0; iscc < fInts; iscc++) {
      for (UInt_t isnb = 0; isnb < fPols; isnb++) {
	
	TH1D* hxsec   = fhXsecCache[f][iscc][isnb];
	TH3D* hmeff   = fhMeffCache[f][iscc][isnb];
	TH2D* hbyfrac = fhBYfracCache[f][iscc][isnb];
	
	for (Int_t ebin = 1; ebin <= hxsec->GetXaxis()->GetNbins(); ebin++) {

	  Double_t E = hxsec->GetXaxis()->GetBinCenter(ebin);
	  hxsec->SetBinContent(ebin, xsec->GetXsec(f, iscc, isnb, E) );

	  for (Int_t bybin = 1; bybin <= hmeff->GetZaxis()->GetNbins(); bybin++) {

	    Double_t by = hmeff->GetZaxis()->GetBinCenter(bybin);
	    hbyfrac->SetBinContent( ebin, bybin, xsec->GetBYfrac(f, iscc, isnb, E, by) );
	    
	    for (Int_t ctbin = 1; ctbin <= hmeff->GetYaxis()->GetNbins(); ctbin++) {

	      Double_t ct = hmeff->GetYaxis()->GetBinCenter(ctbin);
	      hmeff->SetBinContent( ebin, ctbin, bybin, meff->GetMeff( f, iscc, isnb, E, ct, by ) );
	      
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
    \param flav  Neutrino flavor, 0 - elec, 1 - muon, 2 - tau
    \param isnb  Flag for anti-neutrino, 0 - neutrino, 1 - anti-neutrino
    \param tb    A `TrueB` object (see `DetResponse.h`) that stores the true bin coordinate info.
    \return      Cached atmoshperic neutrino flux in specified bin.
*/
Double_t FitUtil::GetCachedFlux(UInt_t flav, Bool_t isnb, Int_t true_ebin, Int_t true_ctbin) {

  if (flav > TAU) {
    throw std::invalid_argument("ERROR! FitUtil::GetCachedFlux() unknown flavor " + to_string(flav));
  }

  return fhFluxCache[flav][(UInt_t)isnb]->GetBinContent(true_ebin, true_ctbin);
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
  ProbCacher(proxymap);
  
  return fhOscCache[flav_in][tb.fFlav][tb.fIsNB]->GetBinContent(tb.fE_true_bin, tb.fCt_true_bin);
}

//***************************************************************************

/** 
    Function that returns the cached neutrino cross-section value for a certain true bin.
    \param tb    A `TrueB` object (see `DetResponse.h`) that stores the neutrino type and true bin coordinate info.
    \return      Cached cross-section value in a specified bin.
*/
Double_t FitUtil::GetCachedXsec(const TrueB &tb) {
  return fhXsecCache[tb.fFlav][tb.fIsCC][tb.fIsNB]->GetBinContent(tb.fE_true_bin);
}

//***************************************************************************

/**
   Function that returns the cached effective mass value for a certain true bin
   \param tb    A `TrueB` object (see `DetResponse.h`) that stores the neutrino type and true bin coordinate info.
   \return      Cached effective mass value in a specified bin in Ton
 */
Double_t FitUtil::GetCachedMeff(const TrueB &tb) {
  return fhMeffCache[tb.fFlav][tb.fIsCC][tb.fIsNB]->GetBinContent( tb.fE_true_bin, tb.fCt_true_bin, tb.fBy_true_bin );
}

//***************************************************************************

/** 
    Function that returns the cached fraction of events expected in a certain BY bin at a certain energy.
   \param tb    A `TrueB` object (see `DetResponse.h`) that stores the neutrino type and true bin coordinate info.
   \return      Cached bjorken-y fraction
*/
Double_t FitUtil::GetCachedBYfrac(const TrueB &tb) {
  return fhBYfracCache[tb.fFlav][tb.fIsCC][tb.fIsNB]->GetBinContent( tb.fE_true_bin, tb.fBy_true_bin );
}

//***************************************************************************

/**
   Private function that handles the caching of the oscillation probabilities.

   This function re-calculates the oscillation probabilites and stores them in member histograms `fhOscCache` whenever any of the oscillation parameters change. If none of them change, the calculation is not performed. The way the `DetResponse` class is set up and the way fitting works, this enables to save a large amount of time.

   NB! The oscillations tau-> (elec, muon, tau) are not calculated to save time; tau flux is assumed 0.

   \param proxymap    A proxy map with `RooFit` variables from `FitPDF` that contain all shared observables and parameters, including the oscillation parameters. 
 */
void FitUtil::ProbCacher(const proxymap_t &proxymap) {

  // get the parameter values from the proxies
  Double_t SinsqTh12 = *( proxymap.at( (TString)fSinsqTh12->GetName() ) );
  Double_t SinsqTh13 = *( proxymap.at( (TString)fSinsqTh13->GetName() ) );
  Double_t SinsqTh23 = *( proxymap.at( (TString)fSinsqTh23->GetName() ) );
  Double_t Dcp       = *( proxymap.at( (TString)fDcp->GetName() ) );
  Double_t Dm21      = *( proxymap.at( (TString)fDm21->GetName() ) );
  Double_t Dm31      = *( proxymap.at( (TString)fDm31->GetName() ) );

  if (SinsqTh12 != f_cache_sinsqth12 || SinsqTh13 != f_cache_sinsqth13 || SinsqTh23 != f_cache_sinsqth23 || 
      Dcp != f_cache_dcp || Dm21 != f_cache_dm21 || Dm31 != f_cache_dm31) {

    fOscCalcTime->Start(kFALSE);

    // update the cache variables
    f_cache_sinsqth12 = SinsqTh12;
    f_cache_sinsqth13 = SinsqTh13;
    f_cache_sinsqth23 = SinsqTh23;
    f_cache_dcp       = Dcp;
    f_cache_dm21      = Dm21;
    f_cache_dm31      = Dm31;

    // give variables to the oscillator
    fProb->SetAngle(1, 2, TMath::ASin( TMath::Sqrt( f_cache_sinsqth12 ) ) );
    fProb->SetAngle(1, 3, TMath::ASin( TMath::Sqrt( f_cache_sinsqth13 ) ) );
    fProb->SetAngle(2, 3, TMath::ASin( TMath::Sqrt( f_cache_sinsqth23 ) ) );
    fProb->SetDelta(1, 3, f_cache_dcp * TMath::Pi()  );
    fProb->SetDm(2, f_cache_dm21);
    fProb->SetDm(3, f_cache_dm31);

    // calculate all oscillation probabilities and cache them  - NB! TAU-IN IS NOT CALCULATED!
    for (UInt_t f_in = 0; f_in < fFlavs-1; f_in++) {
      for (UInt_t f_out = 0; f_out < fFlavs; f_out++) {
	for (UInt_t isnb = 0; isnb < fInts; isnb++) {

	  Int_t N_ebins  = fhOscCache[f_in][f_out][isnb]->GetXaxis()->GetNbins();
	  Int_t N_ctbins = fhOscCache[f_in][f_out][isnb]->GetYaxis()->GetNbins();

	  for (Int_t ebin = 1; ebin <= N_ebins; ebin++) {
	    for (Int_t ctbin = 1; ctbin <= N_ctbins; ctbin++) {

	      Double_t E  = fhOscCache[f_in][f_out][isnb]->GetXaxis()->GetBinCenter(ebin);
	      Double_t ct = fhOscCache[f_in][f_out][isnb]->GetYaxis()->GetBinCenter(ctbin);

	      // set oscillation path and is-nu-bar flag
	      fPrem->FillPath( ct );	      
	      fProb->SetPath ( fPrem->GetNuPath() );
	      fProb->SetIsNuBar( isnb );

	      fhOscCache[f_in][f_out][isnb]->SetBinContent(ebin, ctbin, fProb->Prob(f_in, f_out, E) );

	      fOscCalls++;

	    } // end loop over cos-theta bins
	  } //end loop over energy bins

	} // end loop over isnb
      } // end loop over flavor_out
    } // end loop over flavor_in

    fOscCalcTime->Stop();

  }
  // osc pars have not changed, return
  else { return; }
      
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

/** This function is called inside `FitPDF::analyticalIntegral` and `FitPDF::GetExpValHist` and fills a histogram with expectation values in reco bins.

    The argument map is created in `FitPDF` and contains names and corresponding proxies for all parameters in `fParSet`. The detector response is also part of the `FitPDF` class and configures what kind of an event selection the pdf is used to fit.
    
    \param resp       Pointer to `DetResponse` instance used with the `FitPDF` class.
    \param proxymap   A proxy map with `RooFit` variables from `FitPDF` that contain all shared observables and parameters, including the oscillation parameters. 
    \param norm       An overall normalisation parameter
    \param rangeName  Range string as used in `RooFit`, currently dummy.
    \return           a 3D histogram in reco variables with exepectation values as bin contents

*/
TH3D* FitUtil::Expectation(DetResponse *resp, const proxymap_t &proxymap, Double_t norm, const char* rangeName) {

  // create the histogram with expectation values
  TH3D   *hexp  = (TH3D*)resp->GetHist3D()->Clone();
  TString hname = resp->Get_RespName() + "_expct";
  hexp->SetDirectory(0);
  hexp->Reset();
  hexp->SetNameTitle(hname, hname);
  
  // loop over bins and fill the expectation value histogram
  for (Int_t ebin = fEbin_min; ebin <= fEbin_max; ebin++) {
    for (Int_t ctbin = fCtbin_min; ctbin <= fCtbin_max; ctbin++) {
      for (Int_t bybin = fBybin_min; bybin <= fBybin_max; bybin++) {

        Double_t E  = hexp->GetXaxis()->GetBinCenter( ebin );
        Double_t ct = hexp->GetYaxis()->GetBinCenter( ctbin );
        Double_t by = hexp->GetZaxis()->GetBinCenter( bybin );

	// calculate the bin width, necessary to convert event density to event count
	Double_t e_w  = hexp->GetXaxis()->GetBinWidth( ebin  );
	Double_t ct_w = hexp->GetYaxis()->GetBinWidth( ctbin );
	Double_t by_w = hexp->GetZaxis()->GetBinWidth( bybin );
	Double_t binw = e_w * ct_w * by_w;

        auto recoevts = RecoEvts(E, ct, by, resp, proxymap, norm);
	
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

    Also, an overall normalisation parameter is specific for each `FitPDF` instance, is not part of `proxymap_t` and is thus passed explicitly to the function as argument (see `FitPDF::evaluate`).

    \param E_reco      Reco energy
    \param Ct_reco     Reco cos-theta
    \param By_reco     Reco bjorken-y
    \param resp        Pointer to `DetResponse` instance used with the `FitPDF` class.
    \param proxymap    A proxy map with `RooFit` variables from `FitPDF` that contains the shared fit parameters, including oscillation parameters.
    \param norm        An overall normalisation parameter.
    \return            a pair with the un-normalised event density (calculated by dividing the expected number of events in a bin by bin width) and the associated statistical uncertainty */
std::pair<Double_t, Double_t> FitUtil::RecoEvts(Double_t E_reco, Double_t Ct_reco, Double_t By_reco, DetResponse *resp, const proxymap_t &proxymap, Double_t norm) {

  if ( E_reco < fHB->GetXaxis()->GetXmin() || E_reco >= fHB->GetXaxis()->GetXmax() ) {
    throw std::invalid_argument("ERROR! FitUtil::RecoEvts energy " + to_string(E_reco) + " outside the binning range.");
  }

  if ( Ct_reco < fHB->GetYaxis()->GetXmin() || Ct_reco >= fHB->GetYaxis()->GetXmax() ) {
    throw std::invalid_argument("ERROR! FitUtil::RecoEvts cos-theta " + to_string(Ct_reco) + " outside the binninb range.");
  }

  if ( By_reco < fHB->GetZaxis()->GetXmin() || By_reco >= fHB->GetZaxis()->GetXmax() ) {
    throw std::invalid_argument("ERROR! FitUtil::RecoEvts bjorken-y " + to_string(By_reco) + " outside the binning range.");
  }
  
  auto true_bins = resp->GetBinWeights( E_reco, Ct_reco, By_reco );

  Double_t det_count = 0.;
  Double_t det_err   = 0.;
  
  for (auto &tb: true_bins) {
    
    if (tb.fIsCC) {
      
      Double_t TE = TrueEvts(tb, proxymap).first;

      det_count += tb.fW * TE;
      det_err   += TMath::Power(tb.fWE * TE, 2);
      
    }
    else {

      Double_t TE = 0.;

      /* for NC events the response only says how the elec-NC events from the considered true bin contribute
	 to the reco-bin. However, in the flux chain we have we also have a contribution from muon-NC and tau-NC true bin
	 to the reco bin, which are assumed to look identical to the detector as elec-NC (this is why we only simulate
	 elec-NC). Hence for NC events I need to add `TrueEvts` contributions from the three NC flavours
      */

      TrueB elecTB( tb );
      TrueB muonTB( tb );
      TrueB tauTB ( tb );
      elecTB.fFlav = ELEC;
      muonTB.fFlav = MUON;
      tauTB.fFlav  = TAU;
      
      TE += TrueEvts(elecTB, proxymap).first;
      TE += TrueEvts(muonTB, proxymap).first;
      TE += TrueEvts(tauTB , proxymap).first;
      
      det_count += tb.fW * TE;
      det_err   += TMath::Power(tb.fWE * TE, 2);
      
    }

  } // end loop over true bins

  // should you wish to add atm muons and noise at some point
  //det_count += resp->GetAtmMuCount1y(E_reco, ct_reco, by_reco) * fOpTime;
  //det_count += resp->GetNoiseCount1y(E_reco, ct_reco, by_reco) * fOpTime;

  det_err = TMath::Sqrt(det_err);

  // finally, convert the event counts in reco bin to event density
  Double_t e_w  = fHB->GetXaxis()->GetBinWidth( fHB->GetXaxis()->FindBin( E_reco )  );
  Double_t ct_w = fHB->GetYaxis()->GetBinWidth( fHB->GetYaxis()->FindBin( Ct_reco ) );
  Double_t by_w = fHB->GetZaxis()->GetBinWidth( fHB->GetZaxis()->FindBin( By_reco ) );
  Double_t binw = e_w * ct_w * by_w;
  
  return std::make_pair(det_count/binw * norm, det_err/binw * norm);

}

//***************************************************************************

/** Public function to fetch a pointer to a variable in the member `fParSet` or `fNormSet`.

    This function knows of all the variables known to `RooFit`. If the variable is not found set an exception is thrown.

    \param varname   Name of the variable
    \return          Pointer to the variable
 */
RooRealVar* FitUtil::GetVar(TString varname) {

  RooRealVar *var1 = (RooRealVar*)fParSet.find(varname);
  RooRealVar *var2 = (RooRealVar*)fNormSet.find(varname);

  if (var1 == NULL && var2 == NULL) {
    throw std::invalid_argument("ERROR! FitUtil::GetVar() cannot find variable " + (string)varname);
  }

  if (var1 != NULL && var2 != NULL) {
    throw std::logic_error("ERROR! FitUtil::GetVar() variable " + (string)varname + " is found in both fParSet and fNormSet members.");
  }
  
  RooRealVar *var = ( (var1 != NULL) ? var1 : var2 );
  
  return var;
}

//***************************************************************************

/** Public function to add a normalisation constant for a specific `FitPDF` instance.
    
    Unlike other fit parameters, the overall normalisation parameters are not shared between different event classes, e.g. track normalisation does not affect the shower channel. This function is used in `FitPDF` to automatically add a normalisation constant for each `FitPDF` instance.

    \param name   constant name
    \param title  constant title
    \param val    value
    \param min    minimum
    \param max    maximum						
*/
void FitUtil::AddNormConst(TString name, TString title, Double_t val, Double_t min, Double_t max, Bool_t setConstant) {

  // check whether such a variable exists already, throw error if that is the case
  RooRealVar *var = (RooRealVar*)fNormSet.find( name );
  if ( var != NULL ) {
    throw std::invalid_argument("ERROR! FitUtil::AddNormConst() a normalisation constant with name " + (string)name + " already exists." );
  }
  
  fNorms.push_back( new RooRealVar(name, title, val, min, max) );
  if ( setConstant ) fNorms.back()->setConstant(kTRUE);
  fNormSet.add( *fNorms.back() );
  
}
