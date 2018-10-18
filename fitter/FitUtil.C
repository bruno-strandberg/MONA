#include "FitUtil.h"
#include <iostream>

#include "TFile.h"

using namespace std;

/** Constructor*/
FitUtil::FitUtil(Double_t op_time, TH3 *h_template,
		 Double_t emin, Double_t emax, Double_t ctmin, Double_t ctmax, Double_t bymin, Double_t bymax,
		 TString meffh_elec_cc, TString meffh_muon_cc, TString meffh_tau_cc, TString meffh_elec_nc) {

  fOpTime = op_time;

  fHB = (TH3D*)h_template->Clone("HB");
  fHB->SetDirectory(0);
  fHB->Reset();
  fHB->SetNameTitle("HB", "HB");

  fFlux = new AtmFlux;
  fXsec = new NuXsec( fHB->GetZaxis()->GetNbins() );
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
  ReadMeffHists(fHB, meffh_elec_cc, meffh_muon_cc, meffh_tau_cc, meffh_elec_nc);
  FillFluxAndXsecCache(fFlux, fXsec, fOpTime);

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

/** Destructor.*/
FitUtil::~FitUtil() {

  cout << "FitUtil::~FitUtil() total oscillator calls: " << fOscCalls << endl;
  cout << "FitUtil::~FitUtil() duration of oscillation calculations [s]: "
       << (Double_t)fOscCalcTime->RealTime() << endl;

  if (fOscCalcTime) delete fOscCalcTime;

  if (fHB)   delete fHB;
  if (fFlux) delete fFlux;
  if (fXsec) delete fXsec;
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
	if ( fhMeff[f][iscc][isnb] ) delete fhMeff[f][iscc][isnb];
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

void FitUtil::InitFitVars(Double_t emin, Double_t emax, Double_t ctmin, Double_t ctmax, 
			  Double_t bymin, Double_t bymax) {

  // observables
  fE_reco    = new RooRealVar("E_reco" , "E_reco" , emin, emax);
  fCt_reco   = new RooRealVar("ct_reco", "ct_reco", ctmin, ctmax);
  fBy_reco   = new RooRealVar("by_reco", "by_reco", bymin, bymax);

  // fit parameters, this will probably require some sort of an interface...
  fSinsqTh12 = new RooRealVar("SinsqTh12", "sin^2(theta12)",     0.297,      0.25,     0.354);
  fSinsqTh13 = new RooRealVar("SinsqTh13", "sin^2(theta13)",    0.0215,     0.019,    0.0242);
  fSinsqTh23 = new RooRealVar("SinsqTh23", "sin^2(theta23)",     0.425,     0.381,     0.636);
  fDcp       = new RooRealVar(      "dcp",       "delta-cp",      1.38,         0,      3.14);
  fDm21      = new RooRealVar(     "Dm21",         "dm21^2", 7.37*1e-5, 6.93*1e-5, 7.96*1e-5);
  fDm31      = new RooRealVar(     "Dm31",         "dm31^2", 2.56*1e-3, 2.42*1e-3, 2.69*1e-3);;   

  // add observables to the observables set
  fObsList.add( RooArgList(*fE_reco, *fCt_reco, *fBy_reco) );

  // add all variables to the parameter set
  fParSet.add( fObsList );
  fParSet.add( RooArgSet( *fSinsqTh12, *fSinsqTh13, *fSinsqTh23, *fDcp, *fDm21, *fDm31) );
  
}

//***************************************************************************

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
	  xsec->SelectInteraction(f, iscc, isnb);
	  hxsec->SetBinContent(ebin, xsec->GetXsec(E) );

	}
	
      }
    }
  }

}

//***************************************************************************

/**
   Function that handles the caching of the oscillation probabilities.

   NB! The oscillations tau-> (elec, muon, tau) are not calculated to save time; tau flux is assumed 0.
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

/**
   Private function to read effective mass histograms, copied from `evt_sampler/FluxChain.C`

   \param h_template     3D histogram template that carries the binning information
   \param meffh_elec_cc  Name of the effective mass histograms file for elec CC
   \param meffh_muon_cc  Name of the effective mass histograms file for muon CC
   \param meffh_tau_cc   Name of the effective mass histograms file for tau CC
   \param meffh_elec_nc  Name of the effective mass histograms file for elec NC

 */
void FitUtil::ReadMeffHists(TH3D* h_template, TString meffh_elec_cc, TString meffh_muon_cc, 
			    TString meffh_tau_cc, TString meffh_elec_nc) {

  //  map of flavor numbers and strings for histogram names from file
  map < Int_t, TString > flavs  = { {0, "elec"},   
				    {1, "muon"},
				    {2, "tau" } };

  // map of nc/cc numbers and strings for histogram names from file
  map < Int_t, TString > itypes = { {0, "nc"},
				    {1, "cc"} };
  
  vector<TString> meff_filenames = {meffh_elec_cc, meffh_muon_cc, meffh_tau_cc};

  for (UInt_t f = 0; f < fFlavs; f++) {
    for (UInt_t cc = 0; cc < fInts; cc++) {

      // for nc events the effective mass is identical for all flavors, only elec_NC simulated
      TString meff_fname = meff_filenames[f];
      if (cc == 0) meff_fname = meffh_elec_nc;
    
      TFile meff_file(meff_fname, "READ");
      if ( !meff_file.IsOpen() ) {
	throw std::invalid_argument( "ERROR! FitUtil::ReadMeffHists() could not find file " + (string)meff_fname );
      }

      //---------------------------------------------------------------
      // get the histograms, rebin and divide to get effective mass histos
      //---------------------------------------------------------------

      TString hname = "Meff_" + flavs[f] + "_" + itypes[cc];

      // get the histograms from file
      TH3D *h_gen_nu   = (TH3D*)meff_file.Get("Generated_scaled_nu");
      TH3D *h_gen_nub  = (TH3D*)meff_file.Get("Generated_scaled_nub");
      fhMeff[f][cc][0] = (TH3D*)meff_file.Get("Detected_nu")->Clone(hname + "_nu");
      fhMeff[f][cc][1] = (TH3D*)meff_file.Get("Detected_nub")->Clone(hname + "_nub");

      // detach from meff_file
      fhMeff[f][cc][0]->SetDirectory(0);
      fhMeff[f][cc][1]->SetDirectory(0);

      // calculate rebinning
      Int_t existing_ebins  = fhMeff[f][cc][0]->GetXaxis()->GetNbins();
      Int_t existing_ctbins = fhMeff[f][cc][0]->GetYaxis()->GetNbins();
      Int_t existing_bybins = fhMeff[f][cc][0]->GetZaxis()->GetNbins();

      Int_t rebinning_ebins  = existing_ebins /h_template->GetXaxis()->GetNbins();
      Int_t rebinning_ctbins = existing_ctbins/h_template->GetYaxis()->GetNbins();
      Int_t rebinning_bybins = existing_bybins/h_template->GetZaxis()->GetNbins();

      // check that rebinning is valid
      if ( existing_ebins % h_template->GetXaxis()->GetNbins() != 0 ) {
	throw std::invalid_argument( "ERROR! FitUtil::ReadMeffHists() energy axis nbins=" + to_string(existing_ebins) + " of file " + (string)meff_fname + " cannot be rebinned to " + to_string( h_template->GetXaxis()->GetNbins() ) + ". Change the binning template or change and re-run effective_mass/EffMhists.C with suitable binning." );
      }

      if ( existing_ctbins % h_template->GetYaxis()->GetNbins() != 0 ) {
	throw std::invalid_argument( "ERROR! FitUtil::ReadMeffHists() costheta axis nbins=" + to_string(existing_ctbins) + " of file " + (string)meff_fname + " cannot be rebinned to " + to_string( h_template->GetYaxis()->GetNbins() ) + ". Change the binning template or change and re-run effective_mass/EffMhists.C with suitable binning." );
      }
      
      if ( existing_bybins % h_template->GetZaxis()->GetNbins() != 0 ) {
	throw std::invalid_argument( "ERROR! FitUtil::ReadMeffHists() bjorken-y axis nbins=" + to_string(existing_bybins) + " of file " + (string)meff_fname + " cannot be rebinned to " + to_string( h_template->GetZaxis()->GetNbins() ) + ". Change the binning template or change and re-run effective_mass/EffMhists.C with suitable binning." );
      }

      // perform rebinning
      h_gen_nu ->Rebin3D(rebinning_ebins, rebinning_ctbins, rebinning_bybins);
      h_gen_nub->Rebin3D(rebinning_ebins, rebinning_ctbins, rebinning_bybins);
      fhMeff[f][cc][0]->Rebin3D(rebinning_ebins, rebinning_ctbins, rebinning_bybins);
      fhMeff[f][cc][1]->Rebin3D(rebinning_ebins, rebinning_ctbins, rebinning_bybins);

      // divide to get effective mass in Ton
      fhMeff[f][cc][0]->Divide(h_gen_nu);
      fhMeff[f][cc][1]->Divide(h_gen_nub);
            
      meff_file.Close();

    }
  }

}

//***************************************************************************

Double_t FitUtil::TrueEvts(Int_t ebin_true, Int_t ctbin_true, Int_t bybin_true, 
			   UInt_t flav, UInt_t iscc, UInt_t isnb, 
			   Double_t SinsqTh12, Double_t SinsqTh13, Double_t SinsqTh23, 
			   Double_t Dcp, Double_t Dm21, Double_t Dm31) {

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
  
  // get the effective mass
  Double_t meff = fhMeff[flav][iscc][isnb]->GetBinContent(ebin_true, ctbin_true, bybin_true);

  // find true observable values necassary to split events to bjorken-y bins
  Double_t e_true  = fHB->GetXaxis()->GetBinCenter( ebin_true );
  Double_t by_true = fHB->GetZaxis()->GetBinCenter( bybin_true );
  fXsec->SelectInteraction(flav, iscc, isnb);

  // calculate the number of detected events (unitless)
  Double_t det_count = int_count * fXsec->GetBYfrac(e_true, by_true) * meff * 1e-6;

  return det_count;

}

//***************************************************************************

Double_t FitUtil::PdfEvaluate(const std::map<TString, RooRealProxy*> &parmap, DetResponse *resp) {

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

  return RecoEvts(resp, E_reco, Ct_reco, By_reco, SinsqTh12, SinsqTh13, SinsqTh23, Dcp, Dm21, Dm31);

}

//***************************************************************************

/** Development: will need to think about the range here!*/
Double_t FitUtil::PdfIntegrate(const std::map<TString, RooRealProxy*> &parmap, DetResponse *resp,
			       const char* rangeName) {

  // get the parameter values from the proxies

  Double_t SinsqTh12 = *( parmap.at( (TString)fSinsqTh12->GetName() ) );
  Double_t SinsqTh13 = *( parmap.at( (TString)fSinsqTh13->GetName() ) );
  Double_t SinsqTh23 = *( parmap.at( (TString)fSinsqTh23->GetName() ) );
  Double_t Dcp       = *( parmap.at( (TString)fDcp->GetName() ) );
  Double_t Dm21      = *( parmap.at( (TString)fDm21->GetName() ) );
  Double_t Dm31      = *( parmap.at( (TString)fDm31->GetName() ) );

  return GetIntegral(resp, SinsqTh12, SinsqTh13, SinsqTh23, Dcp, Dm21, Dm31);
  
}

//***************************************************************************

Double_t FitUtil::RecoEvts(DetResponse *resp, 
			   Double_t E_reco, Double_t Ct_reco, Double_t By_reco,
			   Double_t SinsqTh12, Double_t SinsqTh13, Double_t SinsqTh23,
			   Double_t Dcp, Double_t Dm21, Double_t Dm31) {

  auto true_bins = resp->GetBinWeights( E_reco, Ct_reco, By_reco );

  Double_t det_count = 0;

  for (auto &tb: true_bins) {
    
    if (tb.fIsCC) {
      
      Double_t TE = TrueEvts(tb.fE_true_bin, tb.fCt_true_bin, tb.fBy_true_bin, 
			     tb.fFlav, tb.fIsCC, tb.fIsNB,
			     SinsqTh12, SinsqTh13, SinsqTh23, Dcp, Dm21, Dm31);

      det_count += tb.fW * TE;
      
    }
    else {

      Double_t TE = 0.;

      TE += TrueEvts(tb.fE_true_bin, tb.fCt_true_bin, tb.fBy_true_bin, ELEC, tb.fIsCC, tb.fIsNB,
		     SinsqTh12, SinsqTh13, SinsqTh23, Dcp, Dm21, Dm31);
      
      TE += TrueEvts(tb.fE_true_bin, tb.fCt_true_bin, tb.fBy_true_bin, MUON, tb.fIsCC, tb.fIsNB,
		     SinsqTh12, SinsqTh13, SinsqTh23, Dcp, Dm21, Dm31);

      TE += TrueEvts(tb.fE_true_bin, tb.fCt_true_bin, tb.fBy_true_bin,  TAU, tb.fIsCC, tb.fIsNB,
		     SinsqTh12, SinsqTh13, SinsqTh23, Dcp, Dm21, Dm31);

      det_count += tb.fW * TE;
      
    }

  } // end loop over true bins

  // should you wish to add atm muons and noise at some point
  //det_count += resp->GetAtmMuCount1y(E_reco, ct_reco, by_reco) * fOpTime;
  //det_count += resp->GetNoiseCount1y(E_reco, ct_reco, by_reco) * fOpTime;

  return det_count;

}

//***************************************************************************

/**
   Development - will also need to think about the integration range!
 */
Double_t FitUtil::GetIntegral(DetResponse *resp,
			      Double_t SinsqTh12, Double_t SinsqTh13, Double_t SinsqTh23,
			      Double_t Dcp, Double_t Dm21, Double_t Dm31) {
  
  Double_t integral = 0.;
  TH3D *hb = resp->GetHist3D();
  
  for (Int_t ebin = fEbin_min; ebin <= fEbin_max; ebin++) {
    for (Int_t ctbin = fCtbin_min; ctbin <= fCtbin_max; ctbin++) {
      for (Int_t bybin = fBybin_min; bybin <= fBybin_max; bybin++) {

        Double_t E  = hb->GetXaxis()->GetBinCenter( ebin );
        Double_t ct = hb->GetYaxis()->GetBinCenter( ctbin );
        Double_t by = hb->GetZaxis()->GetBinCenter( bybin );

	integral += RecoEvts(resp, E, ct, by, SinsqTh12, SinsqTh13, SinsqTh23, Dcp, Dm21, Dm31);
	
      }
    }
  }
  
  return integral;
  
}

//***************************************************************************
