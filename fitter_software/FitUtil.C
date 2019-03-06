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
  FillFluxAndXsecCache(fFlux, fXsec, fOpTime);

  f_cache_sinsqth12 = 0;
  f_cache_sinsqth13 = 0;
  f_cache_sinsqth23 = 0;
  f_cache_dcp       = 0;
  f_cache_dm21      = 0;
  f_cache_dm31      = 0;

  fOscCalls    = 0;
  fTrueEvtsCalls = 0;
  fOscCalcTime = new TStopwatch();
  fOscCalcTime->Stop(); fOscCalcTime->Reset();

}

//***************************************************************************

/** Destructor */
FitUtil::~FitUtil() {

  cout << "FitUtil::~FitUtil() total oscillator calls: " << fOscCalls << endl;
  cout << "FitUtil::~FitUtil() duration of oscillation calculations [s]: "
       << (Double_t)fOscCalcTime->RealTime() << endl;
  cout << "FitUtil::~FitUtil() total `FitUtil::TrueEvts` calls: " << fTrueEvtsCalls << endl;

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
  TH2D *h_template_2D = (TH2D*)h_template->Project3D("yx");
  TH1D *h_template_1D = (TH1D*)h_template->Project3D("x");

  // init the oscillation probability cache histograms
  for (UInt_t f_in = 0; f_in < fFlavs; f_in++) {
    for (UInt_t f_out = 0; f_out < fFlavs; f_out++) {
      for (UInt_t isnb = 0; isnb < fPols; isnb++) {

	TString name = "oscprob_" + flav_map[f_in] + "_to_" + flav_map[f_out] + "_" + pol_map[isnb];

        fhOscCache[f_in][f_out][isnb] = (TH2D*)h_template_2D->Clone();
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

      fhFluxCache[f][isnb] = (TH2D*)h_template_2D->Clone();
      fhFluxCache[f][isnb]->SetDirectory(0);
      fhFluxCache[f][isnb]->Reset();
      fhFluxCache[f][isnb]->SetNameTitle(name, name);

    }
  }

  // init the xsec histograms
  for (UInt_t f = 0; f < fFlavs; f++) {
    for (UInt_t iscc = 0; iscc < fInts; iscc++) {
      for (UInt_t isnb = 0; isnb < fPols; isnb++) {

        TString name = "xsec_" + flav_map[f] + "_" + iscc_map[iscc] + "_" + pol_map[isnb];

        fhXsecCache[f][iscc][isnb] = (TH1D*)h_template_1D->Clone();
        fhXsecCache[f][iscc][isnb]->SetDirectory(0);
        fhXsecCache[f][iscc][isnb]->Reset();
        fhXsecCache[f][iscc][isnb]->SetNameTitle(name, name);

      }
    }
  }

}

//***************************************************************************

/** A private function to fill the flux and xsec cache hists
    \param flux     Pointer to the `AtmFlux` member instance
    \param xsec     Pointer to the `NuXsec` member instance
    \param op_time  Operation time in years.
*/
void FitUtil::FillFluxAndXsecCache(AtmFlux *flux, NuXsec *xsec, Double_t op_time) {

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

  // fill xsec cache
  for (UInt_t f = 0; f < fFlavs; f++) {
    for (UInt_t iscc = 0; iscc < fInts; iscc++) {
      for (UInt_t isnb = 0; isnb < fPols; isnb++) {
	
	TH1D* hxsec = fhXsecCache[f][iscc][isnb];

	for (Int_t ebin = 1; ebin <= hxsec->GetXaxis()->GetNbins(); ebin++) {

	  Double_t E = hxsec->GetXaxis()->GetBinCenter(ebin);
	  hxsec->SetBinContent(ebin, xsec->GetXsec(f, iscc, isnb, E) );

	}
	
      }
    }
  }

}

//***************************************************************************

Double_t FitUtil::GetCachedFlux(TrueB &tb) {
  return fhFluxCache[tb.fFlav][tb.fIsNB]->GetBinContent(tb.fE_true_bin, tb.fCt_true_bin);
}

//***************************************************************************

Double_t FitUtil::GetCachedOsc(UInt_t flav_in, TrueB &tb, proxymap_t& proxymap) {

  if (flav_in > MUON) {
    throw std::invalid_argument("ERROR! FitUtil::GetCachedOsc() tau-->flavor oscillations are not calculated by FitUtil::ProbCacher(), as the atmospheric flux of tau's is negligible compared to muon and elec flux at ORCA energies");
  }

  // recalculate the oscillation parameter cache if anything has changed.
  ProbCacher(proxymap);
  
  return fhOscCache[flav_in][tb.fFlav][tb.fIsNB]->GetBinContent(tb.fE_true_bin, tb.fCt_true_bin);
}

//***************************************************************************

Double_t FitUtil::GetCachedXsec(TrueB &tb) {
  return fhXsecCache[tb.fFlav][tb.fIsCC][tb.fIsNB]->GetBinContent(tb.fE_true_bin);
}

//***************************************************************************

void FitUtil::ProbCacher(proxymap_t& proxymap) {

  // get the parameter values from the proxies
  Double_t SinsqTh12 = *( proxymap.at( (TString)fSinsqTh12->GetName() ) );
  Double_t SinsqTh13 = *( proxymap.at( (TString)fSinsqTh13->GetName() ) );
  Double_t SinsqTh23 = *( proxymap.at( (TString)fSinsqTh23->GetName() ) );
  Double_t Dcp       = *( proxymap.at( (TString)fDcp->GetName() ) );
  Double_t Dm21      = *( proxymap.at( (TString)fDm21->GetName() ) );
  Double_t Dm31      = *( proxymap.at( (TString)fDm31->GetName() ) );

  ProbCacher(SinsqTh12, SinsqTh13, SinsqTh23, Dcp, Dm21, Dm31);
    
}

//***************************************************************************

/**
   Private function that handles the caching of the oscillation probabilities.

   This function re-calculates the oscillation probabilites and stores them in member histograms `fhOscCache` whenever any of the oscillation parameters change. If none of them change, the calculation is not performed. The way the `DetResponse` class is set up and the way fitting works, this enables to save a large amount of time.

   NB! The oscillations tau-> (elec, muon, tau) are not calculated to save time; tau flux is assumed 0.

   \param SinsqTh12   \f$ sin^2\theta_{12} \f$ value
   \param SinsqTh13   \f$ sin^2\theta_{13} \f$ value
   \param SinsqTh23   \f$ sin^2\theta_{23} \f$ value
   \param Dcp         \f$ \delta_{CP} \f$ value
   \param Dm21        \f$ \Delta m_{21}^2 \f$ value
   \param Dm31        \f$ \Delta m_{31}^2 \f$ value
 */
void FitUtil::ProbCacher(Double_t SinsqTh12, Double_t SinsqTh13, Double_t SinsqTh23, 
			 Double_t Dcp, Double_t Dm21, Double_t Dm31) {

  if (SinsqTh12 != f_cache_sinsqth12 || SinsqTh13 != f_cache_sinsqth13 || SinsqTh23 != f_cache_sinsqth23 || 
      Dcp != f_cache_dcp || Dm21 != f_cache_dm21 || Dm31 != f_cache_dm31) {

    fOscCalcTime->Start(kFALSE);

    // set the cache variables
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

/** Private function that calculates the number of events in a \f$ (E_{\rm true}, cos\theta_{\rm true}, by_{\rm true}) \f$ bin.
    \param ebin_true   Number of true energy bin
    \param ctbin_true  Number of true cos-theta bin
    \param bybin_true  Number of true bjorken-y bin
    \param flav        neutrino flavor (0 - elec, 1 - tau, 2 - muon)
    \param iscc        0 - NC, 1 - CC
    \param isnb        0 - nu, 1 - nub
    \param SinsqTh12   \f$ sin^2\theta_{12} \f$ value
    \param SinsqTh13   \f$ sin^2\theta_{13} \f$ value
    \param SinsqTh23   \f$ sin^2\theta_{23} \f$ value
    \param Dcp         \f$ \delta_{CP} \f$ value
    \param Dm21        \f$ \Delta m_{21}^2 \f$ value
    \param Dm31        \f$ \Delta m_{31}^2 \f$ value
    \return            A pair; first is the number of expected events in the true bin, second is the MC statistical uncertainty
*/
std::pair<Double_t, Double_t> FitUtil::TrueEvts(Int_t ebin_true, Int_t ctbin_true, Int_t bybin_true, 
						UInt_t flav, UInt_t iscc, UInt_t isnb, 
						Double_t SinsqTh12, Double_t SinsqTh13, Double_t SinsqTh23, 
						Double_t Dcp, Double_t Dm21, Double_t Dm31) {

  fTrueEvtsCalls++;

  // check that the bins are in the range
  if ( ebin_true  < 1 || ebin_true  > fHB->GetXaxis()->GetNbins() || 
       ctbin_true < 1 || ctbin_true > fHB->GetYaxis()->GetNbins() || 
       bybin_true < 1 || bybin_true > fHB->GetZaxis()->GetNbins() ) {
    throw std::invalid_argument("ERROR! FitUtil::TrueEvts() event out of binning bounds (ebin, ctbin, bybin): " + 
				to_string(ebin_true) + " " + to_string(ctbin_true) + " " + to_string(bybin_true));
  }

  if (flav >= fFlavs) {
    throw std::invalid_argument("ERROR! FitUtil::TrueEvts() unknown flavor " + to_string(flav) );
  }

  if (iscc >= fInts) {
    throw std::invalid_argument("ERROR! FitUtil::TrueEvts() unknown interaction " + to_string(iscc) );
  }

  if (isnb >= fPols) {
    throw std::invalid_argument("ERROR! FitUtil::TrueEvts() unknown polarisation " + to_string(isnb) );
  }

  // pass the oscillation parameters to the oscillation cache handler
  ProbCacher(SinsqTh12, SinsqTh13, SinsqTh23, Dcp, Dm21, Dm31);

  // get the atm nu count
  Double_t atm_count_e = fhFluxCache[ELEC][isnb]->GetBinContent(ebin_true, ctbin_true);
  Double_t atm_count_m = fhFluxCache[MUON][isnb]->GetBinContent(ebin_true, ctbin_true);
  
  // get the oscillation probabilities
  Double_t prob_elec = fhOscCache[ELEC][flav][isnb]->GetBinContent(ebin_true, ctbin_true);
  Double_t prob_muon = fhOscCache[MUON][flav][isnb]->GetBinContent(ebin_true, ctbin_true);

  // get the oscillated nu count in operation time (in units 1/m2)
  Double_t osc_count = ( atm_count_e * prob_elec + atm_count_m * prob_muon );

  // get the interacted neutrino count in operation time (in units 1/MTon)
  Double_t int_count = osc_count * fhXsecCache[flav][iscc][isnb]->GetBinContent(ebin_true)/fMN * fKg_per_MTon;
  
  // find true observable values necassary for the effective mass and to split events to bjorken-y bins
  Double_t e_true  = fHB->GetXaxis()->GetBinCenter( ebin_true );
  Double_t ct_true = fHB->GetYaxis()->GetBinCenter( ctbin_true );
  Double_t by_true = fHB->GetZaxis()->GetBinCenter( bybin_true );

  // get the effective mass
  Double_t meff  = fMeff->GetMeff( flav, iscc, isnb, e_true, ct_true, by_true );

  // calculate the number of detected events (unitless)
  Double_t det_count = int_count * fXsec->GetBYfrac(flav, iscc, isnb, e_true, by_true) * meff * 1e-6;

  // MC stat err, coming from eff mass and BY distribution (0 for now)
  Double_t det_err   = 0.;

  return std::make_pair(det_count, det_err);

}

//***************************************************************************

/** This function is called from inside `FitPDF::evaluate()` and returns the expected event density at a \f$ (E_{\rm reco}, cos\theta_{\rm reco}, by_{\rm reco}) \f$ value.

    Note that the return of this function is normalised by RooFit by dividing with \f$ \int f(\vec{x})d\vec{x} \f$, where f(x) is the value returned by this function and \f$ \vec{x} \f$ are energy, cos-theta and bjorken-y.

    The argument map is created in `FitPDF` and contains names and corresponding proxies for all parameters in `fParSet`. The detector response is also part of the `FitPDF` class and configures what kind of an event selection the pdf is used to fit.

    \param parmap   Reference to a map with parameter names and corresponding `RooRealProxy`'s.
    \param resp     Pointer to `DetResponse` instance used with the `FitPDF` class.
    \return         a pair with the un-normalised event density (calculated by dividing the expected number of events in a bin by bin width) and the associated statistical uncertainty

*/
std::pair<Double_t, Double_t> FitUtil::PdfEvaluate(const proxymap_t &parmap,
						   DetResponse *resp) {

  // get the parameter values from the proxies

  Double_t E_reco    = *( parmap.at( (TString)fE_reco->GetName() ) );
  Double_t Ct_reco   = *( parmap.at( (TString)fCt_reco->GetName() ) );
  Double_t By_reco   = *( parmap.at( (TString)fBy_reco->GetName() ) );
  Double_t SinsqTh12 = *( parmap.at( (TString)fSinsqTh12->GetName() ) );
  Double_t SinsqTh13 = *( parmap.at( (TString)fSinsqTh13->GetName() ) );
  Double_t SinsqTh23 = *( parmap.at( (TString)fSinsqTh23->GetName() ) );
  Double_t Dcp       = *( parmap.at( (TString)fDcp->GetName() ) );
  Double_t Dm21      = *( parmap.at( (TString)fDm21->GetName() ) );
  Double_t Dm31      = *( parmap.at( (TString)fDm31->GetName() ) );

  // calculate the bin width to convert the number of events to event density

  Double_t e_w  = fHB->GetXaxis()->GetBinWidth( fHB->GetXaxis()->FindBin( E_reco )  );
  Double_t ct_w = fHB->GetYaxis()->GetBinWidth( fHB->GetYaxis()->FindBin( Ct_reco ) );
  Double_t by_w = fHB->GetZaxis()->GetBinWidth( fHB->GetZaxis()->FindBin( By_reco ) );
  Double_t binw = e_w * ct_w * by_w;

  auto recoevts = RecoEvts(resp, E_reco, Ct_reco, By_reco, SinsqTh12, SinsqTh13, SinsqTh23, Dcp, Dm21, Dm31);
  
  return std::make_pair(recoevts.first/binw, recoevts.second/binw);

}

//***************************************************************************

/** This function is called inside `FitPDF::analyticalIntegral` and `FitPDF::GetExpValHist` and fills a histogram with expectation values in reco bins.

    The argument map is created in `FitPDF` and contains names and corresponding proxies for all parameters in `fParSet`. The detector response is also part of the `FitPDF` class and configures what kind of an event selection the pdf is used to fit.
    
    \param parmap     Reference to a map with parameter names and corresponding `RooRealProxy`'s.
    \param resp       Pointer to `DetResponse` instance used with the `FitPDF` class.
    \param rangeName  Range string as used in `RooFit`, currently dummy.
    \return           a 3D histogram in reco variables with exepectation values as bin contents

*/
TH3D* FitUtil::PdfExpectation(const proxymap_t &parmap,
			      DetResponse *resp, const char* rangeName) {

  // get the parameter values from the proxies
  Double_t SinsqTh12 = *( parmap.at( (TString)fSinsqTh12->GetName() ) );
  Double_t SinsqTh13 = *( parmap.at( (TString)fSinsqTh13->GetName() ) );
  Double_t SinsqTh23 = *( parmap.at( (TString)fSinsqTh23->GetName() ) );
  Double_t Dcp       = *( parmap.at( (TString)fDcp->GetName() ) );
  Double_t Dm21      = *( parmap.at( (TString)fDm21->GetName() ) );
  Double_t Dm31      = *( parmap.at( (TString)fDm31->GetName() ) );

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
	
        auto recoevts = RecoEvts(resp, E, ct, by, SinsqTh12, SinsqTh13, SinsqTh23, Dcp, Dm21, Dm31);
	
	hexp->SetBinContent(ebin, ctbin, bybin, recoevts.first  );
	hexp->SetBinError  (ebin, ctbin, bybin, recoevts.second );
	
      }
    }
  }
  
  return hexp;
    
}

//***************************************************************************

/** Private function that calculates the number of events in a \f$ (E_{\rm reco}, cos\theta_{\rm reco}, by_{\rm reco}) \f$ bin.
    \param E_reco      Reconstructed energy
    \param Ct_reco     Reconstructed cos-theta
    \param By_reco     Reconstructed bjorken-y
    \param SinsqTh12   \f$ sin^2\theta_{12} \f$ value
    \param SinsqTh13   \f$ sin^2\theta_{13} \f$ value
    \param SinsqTh23   \f$ sin^2\theta_{23} \f$ value
    \param Dcp         \f$ \delta_{CP} \f$ value
    \param Dm21        \f$ \Delta m_{21}^2 \f$ value
    \param Dm31        \f$ \Delta m_{31}^2 \f$ value
    \return            A pair; first is the number of expected events in the reco bin, second is the MC statistical uncertainty
*/
std::pair<Double_t, Double_t> FitUtil::RecoEvts(DetResponse *resp, 
						Double_t E_reco, Double_t Ct_reco, Double_t By_reco,
						Double_t SinsqTh12, Double_t SinsqTh13, Double_t SinsqTh23,
						Double_t Dcp, Double_t Dm21, Double_t Dm31) {

  auto true_bins = resp->GetBinWeights( E_reco, Ct_reco, By_reco );

  Double_t det_count = 0.;
  Double_t det_err   = 0.;
  
  for (auto &tb: true_bins) {
    
    if (tb.fIsCC) {
      
      Double_t TE = TrueEvts(tb.fE_true_bin, tb.fCt_true_bin, tb.fBy_true_bin, 
			     tb.fFlav, tb.fIsCC, tb.fIsNB,
			     SinsqTh12, SinsqTh13, SinsqTh23, Dcp, Dm21, Dm31).first;

      det_count += tb.fW * TE;
      det_err   += TMath::Power(tb.fWE * TE, 2);
      
    }
    else {

      Double_t TE = 0.;

      TE += TrueEvts(tb.fE_true_bin, tb.fCt_true_bin, tb.fBy_true_bin, ELEC, tb.fIsCC, tb.fIsNB,
		     SinsqTh12, SinsqTh13, SinsqTh23, Dcp, Dm21, Dm31).first;
      
      TE += TrueEvts(tb.fE_true_bin, tb.fCt_true_bin, tb.fBy_true_bin, MUON, tb.fIsCC, tb.fIsNB,
		     SinsqTh12, SinsqTh13, SinsqTh23, Dcp, Dm21, Dm31).first;

      TE += TrueEvts(tb.fE_true_bin, tb.fCt_true_bin, tb.fBy_true_bin,  TAU, tb.fIsCC, tb.fIsNB,
		     SinsqTh12, SinsqTh13, SinsqTh23, Dcp, Dm21, Dm31).first;

      det_count += tb.fW * TE;
      det_err   += TMath::Power(tb.fWE * TE, 2);
      
    }

  } // end loop over true bins

  // should you wish to add atm muons and noise at some point
  //det_count += resp->GetAtmMuCount1y(E_reco, ct_reco, by_reco) * fOpTime;
  //det_count += resp->GetNoiseCount1y(E_reco, ct_reco, by_reco) * fOpTime;

  det_err = TMath::Sqrt(det_err);
  
  return std::make_pair(det_count, det_err);

}

//***************************************************************************

/** Public function to fetch a pointer to a variable in the member `fParSet`.

    If the variable is not in the parameter set an exception is thrown.

    \param varname   Name of the variable
    \return          Pointer to the variable
 */
RooRealVar* FitUtil::GetVar(TString varname) {

  RooRealVar * var = (RooRealVar*)fParSet.find(varname);

  if (var == NULL) {
    throw std::invalid_argument("ERROR! FitUtil::GetVar() cannot find variable " + (string)varname);
  }
  
  return var;
}
