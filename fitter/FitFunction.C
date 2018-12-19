#include "FitFunction.h"

#include "TFile.h"

#include <stdexcept>

/**
   Constructor
 */
FitFunction::FitFunction(DetResponse *DR, Double_t op_time,
			 TString meff_elec_cc, TString meff_muon_cc, 
			 TString meff_tau_cc, TString meff_elec_nc) {

  fResponse = DR;
  fOpTime   = op_time;
  fHB       = fResponse->GetHist3D();
  fFlux     = new AtmFlux;
  fXsec     = new NuXsec( fHB->GetZaxis()->GetNbins() );
  fProb     = new OscProb::PMNS_Fast;
  fPrem     = new OscProb::PremModel;
  ReadMeffHists(meff_elec_cc, meff_muon_cc, meff_tau_cc, meff_elec_nc);

  //timers for development
  fPathCalc = new TStopwatch();
  fOscCalc = new TStopwatch();
  fAtmCalc = new TStopwatch();
  fRestCalc = new TStopwatch();

  fPathCalc->Stop();
  fOscCalc->Stop();
  fAtmCalc->Stop();
  fRestCalc->Stop();

  fPathCalc->Reset();
  fOscCalc->Reset();
  fAtmCalc->Reset();
  fRestCalc->Reset();

  fOscCalls = 0;

  InitCacheHists(fHB, fFlux, fXsec, fOpTime);
  f_cache_th12 = 0.;
  f_cache_th13 = 0.;
  f_cache_th23 = 0.;
  f_cache_dcp  = 0.;
  f_cache_dm21 = 0.;
  f_cache_dm31 = 0.;

}

//*******************************************************************************

/**
   Destructor
 */
FitFunction::~FitFunction() {

  //cout << "NOTICE FitFuction() timer data*********************************" << endl;
  //cout << "NOTICE FitFunction() path calculation: "  << (Double_t)fPathCalc->RealTime() << " seconds" << endl;
  //cout << "NOTICE FitFunction() osc calculation : "  << (Double_t)fOscCalc->RealTime() << " seconds" << endl;
  //cout << "NOTICE FitFunction() flux calculation: "  << (Double_t)fAtmCalc->RealTime() << " seconds" << endl;
  //cout << "NOTICE FitFunction() rest calculation: "  << (Double_t)fRestCalc->RealTime() << " seconds" << endl;
  //cout << "NOTICE Probability calculation calls : "  << fOscCalls << endl;
  //cout << "NOTICE FitFuction() timer data*********************************" << endl;

  // delete fFlux;
  // delete fXsec;
  // delete fProb;
  // delete fPrem;

  // for (Int_t f = 0; f < 3; f++) {
  //   for (Int_t iscc = 0; iscc < 2; iscc++) {
  //     for (Int_t isnb = 0; isnb < 2; isnb++) {
  // 	if ( fhMeff[f][iscc][isnb] ) delete fhMeff[f][iscc][isnb];
  //     }
  //   }
  // }

}

//*******************************************************************************

/**
   Function that returns the expected number of events in a reco bin
 */
std::tuple<double, double> FitFunction::operator() (double *x, double *p) {

  // input variables are reco (E, ct, by)
  Double_t E_reco  = x[0];
  Double_t ct_reco = x[1];
  Double_t by_reco = x[2];

  // get the 'true' bins that contribute to this reco bin and loop over them
  auto true_bins = fResponse->GetBinWeights( E_reco, ct_reco, by_reco );

  Double_t det_count = 0;
  Double_t det_err   = 0;
  for (auto &tb: true_bins) {

    if (tb.fIsCC) {

      // increase the count in reco_bin by the fraction this bin counts towards it
      Double_t td = TrueDetected (tb.fE_true_bin, tb.fCt_true_bin, tb.fBy_true_bin, tb.fFlav, tb.fIsCC, tb.fIsNB, p);
      det_count += tb.fW  * td;
      det_err   += std::pow(tb.fWE * td, 2.);
      
    }
    else {

      Double_t elec_nc = TrueDetected(tb.fE_true_bin, tb.fCt_true_bin, tb.fBy_true_bin, 0, tb.fIsCC, tb.fIsNB, p);
      Double_t muon_nc = TrueDetected(tb.fE_true_bin, tb.fCt_true_bin, tb.fBy_true_bin, 1, tb.fIsCC, tb.fIsNB, p);
      Double_t tau_nc  = TrueDetected(tb.fE_true_bin, tb.fCt_true_bin, tb.fBy_true_bin, 2, tb.fIsCC, tb.fIsNB, p);

      det_count += tb.fW  * (elec_nc + muon_nc + tau_nc);
      det_err   += std::pow(tb.fWE * (elec_nc + muon_nc + tau_nc), 2.);

    }

  }

  det_err = std::sqrt(det_err);
  return std::make_tuple(det_count, det_err);

}

//*******************************************************************************

/**
   Function that returns a vector containing the weights used to create a reco bin
  */
std::vector<Double_t> FitFunction::Weights (double *x) {

  // input variables are reco (E, ct, by)
  Double_t E_reco  = x[0];
  Double_t ct_reco = x[1];
  Double_t by_reco = x[2];

  // get the 'true' bins that contribute to this reco bin and loop over them
  auto true_bins = fResponse->GetBinWeights( E_reco, ct_reco, by_reco );

  std::vector<Double_t> weights;
  for (auto &tb: true_bins) {
    weights.push_back(tb.fW);
  }

  return weights;
}

//*******************************************************************************

/**
   Function that returns the expected number of events in true bin
   \param e_true    True energy bin
   \param ct_true   True cos-theta bin
   \param by_true   True bjorken-y bin
   \param flav      Flavor (0 elec, 1 muon, 2 tau)
   \param iscc      Is charged current flag
   \param isnb      Is anti-neutrino flag
   \param p         Input array with oscillation parameters

   \return          The number of true detected neutrinos for the given input
 */
Double_t FitFunction::TrueDetected (Int_t ebin_true, Int_t ctbin_true, Int_t bybin_true, 
				    Int_t flav, Int_t iscc, Int_t isnb, double *p) {

  // input parameters
  Double_t th12 = TMath::ASin( TMath::Sqrt(p[0]) );
  Double_t th13 = TMath::ASin( TMath::Sqrt(p[1]) );
  Double_t th23 = TMath::ASin( TMath::Sqrt(p[2]) );
  Double_t dcp  = p[3] * TMath::Pi();
  Double_t dm21 = p[4];
  Double_t dm31 = p[5];

  // cacher dev
  ProbCacher(th12, th13, th23, dcp, dm21, dm31);

  Double_t e_true  = fHB->GetXaxis()->GetBinCenter( ebin_true  );
  Double_t ct_true = fHB->GetYaxis()->GetBinCenter( ctbin_true );
  Double_t by_true = fHB->GetZaxis()->GetBinCenter( bybin_true );

  Double_t ew      = fHB->GetXaxis()->GetBinWidth ( ebin_true  );
  Double_t ctw     = fHB->GetYaxis()->GetBinWidth ( ctbin_true );

  // configure xsec
  fXsec->SelectInteraction(flav, iscc, isnb);

  // to convert to atm nu count in op time
  Double_t atm_flux_factor = ew * ctw * fSec_per_y * fOpTime;
  
  // get the atmospheric counts in operation time (in units 1/m2)
  fAtmCalc->Start(kFALSE);
  // Double_t atm_count_e = fFlux->Flux_dE_dcosz(0, isnb, e_true, ct_true) * atm_flux_factor;
  // Double_t atm_count_m = fFlux->Flux_dE_dcosz(1, isnb, e_true, ct_true) * atm_flux_factor;
  Double_t atm_count_e = GetCachedFlux(0, isnb, ebin_true, ctbin_true);
  Double_t atm_count_m = GetCachedFlux(1, isnb, ebin_true, ctbin_true);
  fAtmCalc->Stop();

  // get oscillation probabilities
  Double_t prob_elec = GetCachedProb(isnb, 0, flav, ebin_true, ctbin_true);
  Double_t prob_muon = GetCachedProb(isnb, 1, flav, ebin_true, ctbin_true);

  fRestCalc->Start(kFALSE);
  // get the oscillated counts in operation time (in units 1/m2)
  Double_t osc_count = ( atm_count_e * prob_elec + atm_count_m * prob_muon );

  // get the interacted neutrino count in operation time (units 1/MTon)
  //Double_t int_count = osc_count * fXsec->GetXsec(e_true)/fMN * fKg_per_Mton;
  Double_t int_count = osc_count * GetCachedXsec(flav, iscc, isnb, ebin_true)/fMN * fKg_per_Mton;

  // get the effective mass (in units Ton)
  Double_t meff = fhMeff[flav][iscc][isnb]->GetBinContent(ebin_true, ctbin_true, bybin_true);

  // get the detected count in the true bin in operation time (unitless)
  Double_t det_count = int_count * fXsec->GetBYfrac(e_true, by_true) * meff * 1e-6;
  fRestCalc->Stop();

  return det_count;

}

//*******************************************************************************

/**
   Function that initialises the cached hists
 */
void FitFunction::InitCacheHists(TH3D *h_template, AtmFlux *flux, NuXsec *xsec, Double_t op_time) {

  std::map<Int_t, TString> pol_map  = { {0, "nu"}, {1, "nub"} }; 
  std::map<Int_t, TString> iscc_map = { {0, "NC"}, {1, "CC"} }; 
  std::map<Int_t, TString> flav_map = { {0, "elec"}, {1, "muon"}, {2, "tau"} };

  TH2D *h_template_2D = (TH2D*)h_template->Project3D("yx");
  TH1D *h_template_1D = (TH1D*)h_template->Project3D("x");

  //--------------------------------------------------------------
  // init the oscillation probability cache histograms
  //--------------------------------------------------------------
  for (Int_t isnb = 0; isnb < 2; isnb++) {
    for (Int_t f_in = 0; f_in < 2; f_in++) {
      for (Int_t f_out = 0; f_out < 3; f_out++) {
	
	TString name = "oscprob_" + pol_map[isnb] + "_" + flav_map[f_in] + "_to_" + flav_map[f_out];

	fhOscProb[isnb][f_in][f_out] = (TH2D*)h_template_2D->Clone();
	fhOscProb[isnb][f_in][f_out]->SetDirectory(0);
	fhOscProb[isnb][f_in][f_out]->Reset();
	fhOscProb[isnb][f_in][f_out]->SetNameTitle(name, name);
	
      }
    }
  }

  //--------------------------------------------------------------
  // init the atm flux cache and fill
  //--------------------------------------------------------------
  for (Int_t f = 0; f < 2; f++) {
    for (Int_t isnb = 0; isnb < 2; isnb++) {

      TString name = "flux_" + flav_map[f] + "_" + pol_map[isnb];

      fhAtmFluxCache[f][isnb] = (TH2D*)h_template_2D->Clone();
      fhAtmFluxCache[f][isnb]->SetDirectory(0);
      fhAtmFluxCache[f][isnb]->Reset();
      fhAtmFluxCache[f][isnb]->SetNameTitle(name, name);

      for (Int_t ebin = 1; ebin <= fhAtmFluxCache[f][isnb]->GetXaxis()->GetNbins(); ebin++) {
	for (Int_t ctbin = 1; ctbin <= fhAtmFluxCache[f][isnb]->GetYaxis()->GetNbins(); ctbin++) {

	  Double_t E   = fhAtmFluxCache[f][isnb]->GetXaxis()->GetBinCenter(ebin);
	  Double_t ct  = fhAtmFluxCache[f][isnb]->GetYaxis()->GetBinCenter(ctbin);
	  Double_t ew  = fhAtmFluxCache[f][isnb]->GetXaxis()->GetBinWidth( ebin );
	  Double_t ctw = fhAtmFluxCache[f][isnb]->GetYaxis()->GetBinWidth( ctbin );

	  // convert to atm nu count in op time
	  Double_t atm_flux_factor = ew * ctw * fSec_per_y * op_time;
	  Double_t atmflux         = flux->Flux_dE_dcosz(f, isnb, E, ct) * atm_flux_factor;

	  // set the bin content
	  fhAtmFluxCache[f][isnb]->SetBinContent(ebin, ctbin, atmflux);

	}
      }

    }
  }

  //--------------------------------------------------------------
  // init the xsec cache and fill
  //--------------------------------------------------------------

  for (Int_t f = 0; f < 3; f++) {
    for (Int_t iscc = 0; iscc < 2; iscc++) {
      for (Int_t isnb = 0; isnb < 2; isnb++) {
	
	TString name = "xsec_" + flav_map[f] + "_" + iscc_map[iscc] + "_" + pol_map[isnb];

	fhXsecCache[f][iscc][isnb] = (TH1D*)h_template_1D->Clone();
	fhXsecCache[f][iscc][isnb]->SetDirectory(0);
	fhXsecCache[f][iscc][isnb]->Reset();
	fhXsecCache[f][iscc][isnb]->SetNameTitle(name, name);
	
	for (Int_t ebin = 1; ebin <= fhXsecCache[f][iscc][isnb]->GetXaxis()->GetNbins(); ebin++) {

	  Double_t E = fhXsecCache[f][iscc][isnb]->GetXaxis()->GetBinCenter(ebin);
	  xsec->SelectInteraction(f, iscc, isnb);
	  fhXsecCache[f][iscc][isnb]->SetBinContent(ebin, xsec->GetXsec(E) );

	}

      }
    }
  }

}

//*******************************************************************************

/**
   Function that handles caching the oscillation probabilities
 */
void FitFunction::ProbCacher(Double_t th12, Double_t th13, Double_t th23, Double_t dcp, 
			     Double_t dm21, Double_t dm31) {

  if (th12 != f_cache_th12 || th13 != f_cache_th13 || th23 != f_cache_th23 || dcp != f_cache_dcp ||
      dm21 != f_cache_dm21 || dm31 != f_cache_dm31) {

    // set the cache variables
    f_cache_th12 = th12;
    f_cache_th13 = th13;
    f_cache_th23 = th23;
    f_cache_dcp  = dcp;
    f_cache_dm21 = dm21;
    f_cache_dm31 = dm31;

    // give variables to the oscillator
    fProb->SetAngle(1, 2, f_cache_th12 );
    fProb->SetAngle(1, 3, f_cache_th13 );
    fProb->SetAngle(2, 3, f_cache_th23 );
    fProb->SetDelta(1, 3, f_cache_dcp  );
    fProb->SetDm(2, f_cache_dm21);
    fProb->SetDm(3, f_cache_dm31);

    // calculate all oscillation probabilities and cache them
    for (Int_t isnb = 0; isnb < 2; isnb++) {
      for (Int_t f_in = 0; f_in < 2; f_in++) {
	for (Int_t f_out = 0; f_out < 3; f_out++) {

	  Int_t N_ebins  = fhOscProb[isnb][f_in][f_out]->GetXaxis()->GetNbins();
	  Int_t N_ctbins = fhOscProb[isnb][f_in][f_out]->GetYaxis()->GetNbins();

	  for (Int_t ebin = 1; ebin <= N_ebins; ebin++) {
	    for (Int_t ctbin = 1; ctbin <= N_ctbins; ctbin++) {

	      Double_t E  = fhOscProb[isnb][f_in][f_out]->GetXaxis()->GetBinCenter(ebin);
	      Double_t ct = fhOscProb[isnb][f_in][f_out]->GetYaxis()->GetBinCenter(ctbin);

	      // set oscillation path and is-nu-bar flag
	      fPathCalc->Start(kFALSE);
	      fPrem->FillPath( ct );
	      fPathCalc->Stop();

	      
	      fOscCalc->Start(kFALSE);
	      fProb->SetPath ( fPrem->GetNuPath() );
	      fProb->SetIsNuBar( isnb );
	      fhOscProb[isnb][f_in][f_out]->SetBinContent(ebin, ctbin, fProb->Prob(f_in, f_out, E) );
	      fOscCalls++;
	      fOscCalc->Stop();

	    } // end loop over cos-theta bins
	  } //end loop over energy bins

	} // end loop over flavor_out
      } // end loop over flavor_in
    } // end loop over isnb

  }
  // osc pars have not changed, return
  else { return; }

}

//*******************************************************************************

/**
   Function that returns the cached oscillation probability
 */
Double_t FitFunction::GetCachedProb(Int_t isnb, Int_t flav_in, Int_t flav_out, Int_t ebin, Int_t ctbin) {

  if ( ( ebin < 1 ) || ( ebin > fhOscProb[isnb][flav_in][flav_out]->GetXaxis()->GetNbins() ) ||
       ( ctbin < 1 ) || ( ctbin > fhOscProb[isnb][flav_in][flav_out]->GetYaxis()->GetNbins() ) ) {

    Double_t ct = fHB->GetYaxis()->GetBinCenter(ctbin);
    Double_t E  = fHB->GetXaxis()->GetBinCenter(ebin);

    cout << "WARNING! FitFunction::GetCachedProb() Energy, costheta " << E << ", " << ct 
	 << " is outside the cache range, calculating individual point." << endl;

    fPrem->FillPath( ct );
    fProb->SetPath ( fPrem->GetNuPath() );
    fProb->SetIsNuBar( isnb );
    return fProb->Prob(flav_in, flav_out, E);

  }
  else {
    return fhOscProb[isnb][flav_in][flav_out]->GetBinContent(ebin, ctbin);
  }

}

//*******************************************************************************

/**
   Function that returns the cached atmospheric flux
 */
Double_t FitFunction::GetCachedFlux(Int_t flav, Int_t isnb, Int_t ebin, Int_t ctbin) {

  return fhAtmFluxCache[flav][isnb]->GetBinContent(ebin, ctbin);

}

//*******************************************************************************

/**
   Function that returns the cached xsec
 */
Double_t FitFunction::GetCachedXsec(Int_t flav, Int_t iscc, Int_t isnb, Int_t ebin) {

  return fhXsecCache[flav][iscc][isnb]->GetBinContent(ebin);

}

//*******************************************************************************

/**
   Inline function to read effective mass histograms, copied from `evt_sampler/FluxChain.C`

   \param meff_elec_cc  Name of the effective mass file for elec CC
   \param meff_muon_cc  Name of the effective mass file for muon CC
   \param meff_tau_cc   Name of the effective mass file for tau CC
   \param meff_elec_nc  Name of the effective mass file for elec NC

 */
void FitFunction::ReadMeffHists(TString meff_elec_cc, TString meff_muon_cc, 
				TString meff_tau_cc, TString meff_elec_nc) {

  //  map of flavor numbers and strings for histogram names from file
  map < Int_t, TString > flavs  = { {0, "elec" },   
				     {1, "muon" },
				     {2, "tau"  } };

  // map of nc/cc numbers and strings for histogram names from file
  map < Int_t, TString > itypes = { {0, "nc"  },
				    {1, "cc" } };
  
  vector<TString> meff_filenames = {meff_elec_cc, meff_muon_cc, meff_tau_cc};

  for (Int_t f = 0; f < 3; f++) {
    for (Int_t cc = 0; cc < 2; cc++) {

      // for nc events the effective mass is identical for all flavors, only elec_NC simulated
      TString meff_fname = meff_filenames[f];
      if (cc == 0) meff_fname = meff_elec_nc;
    
      TFile meff_file(meff_fname, "READ");
      if ( !meff_file.IsOpen() ) {
	throw std::invalid_argument( "ERROR! FitFunction::ReadMeffHists() could not find file " + (string)meff_fname );
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

      Int_t rebinning_ebins  = existing_ebins /fHB->GetXaxis()->GetNbins();
      Int_t rebinning_ctbins = existing_ctbins/fHB->GetYaxis()->GetNbins();
      Int_t rebinning_bybins = existing_bybins/fHB->GetZaxis()->GetNbins();

      // check that rebinning is valid
      if ( existing_ebins % fHB->GetXaxis()->GetNbins() != 0 ) {
	throw std::invalid_argument( "ERROR! FitFunction::ReadMeffHists() energy axis nbins=" + to_string(existing_ebins) + " of file " + (string)meff_fname + " cannot be rebinned to " + to_string( fHB->GetXaxis()->GetNbins() ) + ". Change the binning of detector response or change and re-run NMHDIR/effective_mass/EffMhists.C with suitable binning." );
      }

      if ( existing_ctbins % fHB->GetYaxis()->GetNbins() != 0 ) {
	throw std::invalid_argument( "ERROR! FitFunction::ReadMeffHists() costheta axis nbins=" + to_string(existing_ctbins) + " of file " + (string)meff_fname + " cannot be rebinned to " + to_string( fHB->GetYaxis()->GetNbins() ) + ". Change the binning of the detector response or change and re-run NMHDIR/effective_mass/EffMhists.C with suitable binning." );
      }
      
      if ( existing_bybins % fHB->GetZaxis()->GetNbins() != 0 ) {
	throw std::invalid_argument( "ERROR! FitFunction::ReadMeffHists() bjorken-y axis nbins=" + to_string(existing_bybins) + " of file " + (string)meff_fname + " cannot be rebinned to " + to_string( fHB->GetZaxis()->GetNbins() ) + ". Change the binning of the detector response or change and re-run NMHDIR/effective_mass/EffMhists.C with suitable binning." );
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
