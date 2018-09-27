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

}

//*******************************************************************************

/**
   Destructor
 */
FitFunction::~FitFunction() {

  cout << "NOTICE FitFuction() timer data*********************************" << endl;
  cout << "NOTICE FitFunction() path calculation: "  << (Double_t)fPathCalc->RealTime() << " seconds" << endl;
  cout << "NOTICE FitFunction() osc calculation : "  << (Double_t)fOscCalc->RealTime() << " seconds" << endl;
  cout << "NOTICE FitFunction() flux calculation: "  << (Double_t)fAtmCalc->RealTime() << " seconds" << endl;
  cout << "NOTICE FitFunction() rest calculation: "  << (Double_t)fRestCalc->RealTime() << " seconds" << endl;
  cout << "NOTICE Probability calculation calls : "  << fOscCalls << endl;
  cout << "NOTICE FitFuction() timer data*********************************" << endl;

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
double FitFunction::operator() (double *x, double *p) {

  // input variables are reco (E, ct, by)
  Double_t E_reco  = x[0];
  Double_t ct_reco = x[1];
  Double_t by_reco = x[2];

  // get the 'true' bins that contribute to this reco bin and loop over them
  auto true_bins = fResponse->GetBinWeights( E_reco, ct_reco, by_reco );

  Double_t det_count = 0;
  for (auto &tb: true_bins) {

    // get true e, ct, by as bin centers
    Double_t e_true  = fHB->GetXaxis()->GetBinCenter( tb.fE_true_bin );
    Double_t ct_true = fHB->GetYaxis()->GetBinCenter( tb.fCt_true_bin );
    Double_t by_true = fHB->GetZaxis()->GetBinCenter( tb.fBy_true_bin );

    if (tb.fIsCC) {

      // increase the count in reco_bin by the fraction this bin counts towards it
      det_count += tb.fW * 
	TrueDetected (e_true, ct_true, by_true, tb.fFlav, tb.fIsCC, tb.fIsNB, p);
      
    }
    else {

      Double_t elec_nc = TrueDetected(e_true, ct_true, by_true, 0, tb.fIsCC, tb.fIsNB, p);
      Double_t muon_nc = TrueDetected(e_true, ct_true, by_true, 1, tb.fIsCC, tb.fIsNB, p);
      Double_t tau_nc  = TrueDetected(e_true, ct_true, by_true, 2, tb.fIsCC, tb.fIsNB, p);

      det_count += tb.fW * (elec_nc + muon_nc + tau_nc);

    }

  }

  return det_count;

}

//*******************************************************************************

/**
   Function that returns the expected number of events in true bin
   \param e_true    True energy
   \param ct_true   True cos-theta
   \param by_true   True bjorken-y
   \param flav      Flavor (0 elec, 1 muon, 2 tau)
   \param iscc      Is charged current flag
   \param isnb      Is anti-neutrino flag
   \param p         Input array with oscillation parameters

   \return          The number of true detected neutrinos for the given input
 */
Double_t FitFunction::TrueDetected (Double_t e_true, Double_t ct_true, Double_t by_true, 
				    Int_t flav, Int_t iscc, Int_t isnb, double *p) {

  // input parameters
  Double_t th12 = TMath::ASin( TMath::Sqrt(p[0]) );
  Double_t th13 = TMath::ASin( TMath::Sqrt(p[1]) );
  Double_t th23 = TMath::ASin( TMath::Sqrt(p[2]) );
  Double_t dcp  = p[3] * TMath::Pi();
  Double_t dm21 = p[4];
  Double_t dm31 = p[5];

  // give them to the oscillator
  fProb->SetAngle(1, 2, th12 );
  fProb->SetAngle(1, 3, th13 );
  fProb->SetAngle(2, 3, th23 );
  fProb->SetDelta(1, 3, dcp  );
  fProb->SetDm(2, dm21);
  fProb->SetDm(3, dm31);

  Int_t ebin  = fHB->GetXaxis()->FindBin( e_true );
  Int_t ctbin = fHB->GetYaxis()->FindBin( ct_true );
  Int_t bybin = fHB->GetZaxis()->FindBin( by_true );

  Double_t ew  = fHB->GetXaxis()->GetBinWidth( ebin );
  Double_t ctw = fHB->GetYaxis()->GetBinWidth( ctbin );

  // set oscillation path
  fPathCalc->Start(kFALSE);
  fPrem->FillPath( ct_true );
  fPathCalc->Stop();
  fProb->SetPath ( fPrem->GetNuPath() );

  // configure oscillator and xsec
  fProb->SetIsNuBar( isnb );
  fXsec->SelectInteraction(flav, iscc, isnb);

  // to convert to atm nu count in op time
  Double_t atm_flux_factor = ew * ctw * fSec_per_y * fOpTime;
  
  // get the atmospheric counts in operation time (in units 1/m2)
  fAtmCalc->Start(kFALSE);
  Double_t atm_count_e = fFlux->Flux_dE_dcosz(0, isnb, e_true, ct_true) * atm_flux_factor;
  Double_t atm_count_m = fFlux->Flux_dE_dcosz(1, isnb, e_true, ct_true) * atm_flux_factor;
  fAtmCalc->Stop();

  // get oscillation probabilities
  fOscCalc->Start(kFALSE);
  Double_t prob_elec = fProb->Prob(0, flav, e_true);
  Double_t prob_muon = fProb->Prob(1, flav, e_true);
  fOscCalls += 2;
  fOscCalc->Stop();

  fRestCalc->Start(kFALSE);
  // get the oscillated counts in operation time (in units 1/m2)
  Double_t osc_count = ( atm_count_e * prob_elec + atm_count_m * prob_muon );

  // get the interacted neutrino count in operation time (units 1/MTon)
  Double_t int_count = osc_count * fXsec->GetXsec(e_true)/fMN * fKg_per_Mton;

  // get the effective mass (in units Ton)
  Double_t meff = fhMeff[flav][iscc][isnb]->GetBinContent(ebin, ctbin, bybin);

  // get the detected count in the true bin in operation time (unitless)
  Double_t det_count = int_count * fXsec->GetBYfrac(e_true, by_true) * meff * 1e-6;
  fRestCalc->Stop();

  return det_count;

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
