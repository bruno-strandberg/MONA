
#include "PMNS_Fast.h"
#include "PremModel.h"
#include "AtmFlux.h"
#include "NuXsec.h"
#include "NMHUtils.h"
#include "TString.h"
#include "TSystem.h"
#include "TMath.h"
#include "TH2.h"
#include "TRandom3.h"
#include "TFile.h"
#include <iostream>

using namespace std;

//*****************************************************************
// functions
//*****************************************************************
void   InitClasses();
void   InitOscPars(Bool_t NH);
void   InitHists();
Bool_t ReadMeffHists();
void   FillHists(Double_t op_time, Int_t flav, Int_t is_cc, Int_t is_nb, Int_t nsamples);
void   CleanUp();
void   WriteToFile(TString output_name);

//*****************************************************************
// globally used variables in this script
//*****************************************************************
OscProb::PMNS_Fast *fProb;     //!< probability Calculator
OscProb::PremModel *fPrem;     //!< PREM Model
AtmFlux  *fFlux;               //!< flux calculator
NuXsec   *fXsec;               //!< cross-section calculator
TRandom3 *fRand;               //!< random number generator for over-sampling bins

TH2D    *fhAtmFlux[3][2][2];   //!< atmospheric flux histograms, indices [flavor][is_cc][is_nb]
TH2D    *fhOscFlux[3][2][2];   //!< oscillated flux histograms
TH2D    *fhIntFlux[3][2][2];   //!< interacted flux histograms
TH2D    *fhDetFlux[3][2][2];   //!< detected flux histograms
TH2D       *fhMeff[3][2][2];   //!< effective mass histograms

map < Int_t, TString > fFlavs  = { {0, "elec" },   //!< map of flavor numbers and strings
				   {1, "muon" },
				   {2, "tau"  } };

map < Int_t, TString > fItypes = { {0, "nc"  },    //!< map of nc/cc numbers and strings
				   {1, "cc" } };

map < Int_t, TString > fPtypes = { {0, "nu"  },    //!< map of particle/antip. numbers and strings
				   {1, "nub" } };

// constants used in FillHists()
Double_t fMp = 1.672621898e-27;                     //!< proton mass in kg
Double_t fMn = 1.674927471e-27;                     //!< neutron mass in kg
Double_t fMN = (fMp+fMn)/2;                         //!< nucleon-average mass
Double_t fSec_per_y   = 365.2421897 * 24 * 60 * 60; //!< seconds in a tropical year
Double_t fKg_per_Mton = 1e9;                        //!< kg per MTon (MTon = 1e6 Ton; Ton = 1e3 kg)

//*****************************************************************

/**
 *  This routine creates a file with E vs costheta histograms with expected numbers of nu events.
 *
 *  The output root file will have directories atmflux, oscflux, intflux, detflux and meff. Eech
 *  of the flux directories contains E vs costheta histograms for {flavor}_{nc/cc}_{nu/nub} (12 in
 *  total), each bin will indicate the number of neutrinos of this type after op_time years.
 *  Directory atmflux/ has histograms with atmospheric neutrinos. oscflux/ shows atmflux/ histograms
 *  after oscillation through Earth, intflux/ shows oscflux/ histograms after nu + H2O xsec is taken
 *  into account, detflux/ shows intflux/ histograms after effective mass is taken into account.
 *  meff/ directory displays the effective mass histograms used for the generation of detflux/
 *  histograms. In the region of fast oscillations it is necessary to sample each bin several times,
 *  parameter nsamples determines how many samples per bin are calculated.
 *
 * \param  op_time       Operation time in years.
 * \param  output_name   Name of the file where histograms are written.
 * \param  NH            True - normal nu mass hierarchy, False - inverted nu mass hierarchy.
 * \param  nsamples      Number of samples per filling one bin (recommended > 10).
 *
 */
void FluxChain(Double_t op_time     = 3.,
	       TString  output_name = "flux_chain_out.root",
	       Bool_t   NH          = true,
	       Int_t    nsamples    = 50) {

  gSystem->Load("$NMHDIR/common_software/libnmhsoft.so");
  gSystem->Load("$OSCPROBDIR/libOscProb.so");
  
  InitClasses();
  InitOscPars(NH);
  InitHists();
  if ( !ReadMeffHists() ) {
    cout << "ERROR! FluxChain() Reading of effective mass histogram failed, quitting." << endl;
    return;
  }
 
  for (Int_t f = 0; f < 3; f++) {
    for (Int_t cc = 0; cc < 2; cc++) {
      for (Int_t nb = 0; nb < 2; nb++) {
	cout << "FluxChain() Filling flux for: " << fFlavs[f] << "\t" << fItypes[cc]
	     << "\t" << fPtypes[nb] << endl;
	FillHists(op_time, f, cc, nb, nsamples);
      }
    }
  }

  WriteToFile(output_name);
  CleanUp();
  
}

//*****************************************************************

/**
 *  Inline function to initialise classes.
 */
void InitClasses() {

  fProb = new OscProb::PMNS_Fast;
  fPrem = new OscProb::PremModel;
  fFlux = new AtmFlux;
  fXsec = new NuXsec;
  fRand = new TRandom3(0);
}

//*****************************************************************

/**
 *  Inline function to initialise osc parameters (PDG 2016) and give them to the osc calculator.
 */
void InitOscPars(Bool_t NH) {

  //------------------------------------------------
  // Oscillation parameters PDG 2016
  //------------------------------------------------
  Double_t sinsq_th12  = 0.304;
  Double_t sinsq_th23  = 0.51;
  Double_t sinsq_th13  = 2.19e-2;
  Double_t dcp         = 0.;
  Double_t dm21        = 7.53e-5;
  Double_t dm32        = 2.44e-3;
  Double_t dm31        = dm32 + dm21;
  
  if (!NH) {
    sinsq_th23 = 0.50;
    dm32       = -2.51e-3;
    dm31       = dm32 + dm21;
  }

  //------------------------------------------------
  // Pass to the oscillation probability calculator
  //------------------------------------------------

  fProb->SetAngle(1, 2, TMath::ASin( TMath::Sqrt(sinsq_th12) ) );
  fProb->SetAngle(1, 3, TMath::ASin( TMath::Sqrt(sinsq_th13) ) );
  fProb->SetAngle(2, 3, TMath::ASin( TMath::Sqrt(sinsq_th23) ) );

  fProb->SetDm(2, dm21);
  fProb->SetDm(3, dm31);
  
}

//*****************************************************************

/**
 *  Inline function to init E vs costheta histograms.
 */
void InitHists() {

  Int_t    ebins  = 40;
  Double_t elow   = 1;
  Double_t ehigh  = 100;
  vector<Double_t> e_edges = NMHUtils::GetLogBins(ebins, elow, ehigh);
  
  Int_t    ctbins = 40;
  Double_t ctlow  = -1;
  Double_t cthigh =  1;
  
  for (Int_t f = 0; f < 3; f++) {
    for (Int_t cc = 0; cc < 2; cc++) {
      for (Int_t nb = 0; nb < 2; nb++) {
      
	TString atmname = "atmflux_" + fFlavs[f] + "_" + fItypes[cc] + "_" + fPtypes[nb];
	TString oscname = "oscflux_" + fFlavs[f] + "_" + fItypes[cc] + "_" + fPtypes[nb];
	TString intname = "intflux_" + fFlavs[f] + "_" + fItypes[cc] + "_" + fPtypes[nb];
	TString detname = "detflux_" + fFlavs[f] + "_" + fItypes[cc] + "_" + fPtypes[nb];

	fhAtmFlux[f][cc][nb] = new TH2D(atmname, atmname, ebins, &e_edges[0], ctbins, ctlow, cthigh);
	fhOscFlux[f][cc][nb] = new TH2D(oscname, oscname, ebins, &e_edges[0], ctbins, ctlow, cthigh);
	fhIntFlux[f][cc][nb] = new TH2D(intname, intname, ebins, &e_edges[0], ctbins, ctlow, cthigh);
	fhDetFlux[f][cc][nb] = new TH2D(detname, detname, ebins, &e_edges[0], ctbins, ctlow, cthigh);

	fhAtmFlux[f][cc][nb]->Sumw2();
	fhOscFlux[f][cc][nb]->Sumw2();
	fhIntFlux[f][cc][nb]->Sumw2();
	fhDetFlux[f][cc][nb]->Sumw2();
       
      }
    }
  }
  
}

//*****************************************************************

/**
 *  Inline function to read effective mass histograms from $NMHDIR/data/eff_mass/.
 */
Bool_t ReadMeffHists() {

  for (Int_t f = 0; f < 3; f++) {
    for (Int_t cc = 0; cc < 2; cc++) {

      // for nc events the effective mass is identical for all flavors, only elec_NC simulated
      TString meff_fname = "$NMHDIR/data/eff_mass/EffMass_" + fFlavs[f] + "_CC.root";
      if (cc == 0) meff_fname = "$NMHDIR/data/eff_mass/EffMass_elec_NC.root";
    
      TFile meff_file(meff_fname, "READ");
      if ( !meff_file.IsOpen() ) {
	cout << "ERROR! InitHists() could not find file " << meff_fname << endl;
	return false;
      }

      TString hname = "Meff_" + fFlavs[f] + "_" + fItypes[cc];
      fhMeff[f][cc][0] = (TH2D*)meff_file.Get("Meff_nu")->Clone(hname + "_nu");
      fhMeff[f][cc][1] = (TH2D*)meff_file.Get("Meff_nu")->Clone(hname + "_nub");
      fhMeff[f][cc][0]->SetDirectory(0); //detach these histograms from meff_file
      fhMeff[f][cc][1]->SetDirectory(0);
    
      meff_file.Close();

    }
  }

  return true;
  
}

//*****************************************************************

/**
 *  Inline function to fill E vs costheta histograms.
 */
void FillHists(Double_t op_time, Int_t flav, Int_t is_cc, Int_t is_nb, Int_t nsamples) {
  
  fProb->SetIsNuBar(is_nb);
  fXsec->SelectInteraction(flav, is_cc, is_nb);
  
  TH2D *h = fhOscFlux[flav][is_cc][is_nb]; //pointer for convenient access to bin data
  
  for (Int_t ybin = 1; ybin <= h->GetYaxis()->GetNbins(); ybin++) {

    Double_t ct_low  = h->GetYaxis()->GetBinLowEdge(ybin);
    Double_t ct_high = h->GetYaxis()->GetBinUpEdge(ybin);
    Double_t ct_w    = h->GetYaxis()->GetBinWidth(ybin);
    
    for (Int_t xbin = 1; xbin <= h->GetXaxis()->GetNbins(); xbin++) {

      Double_t en_low  = h->GetXaxis()->GetBinLowEdge(xbin);
      Double_t en_high = h->GetXaxis()->GetBinUpEdge(xbin);
      Double_t en_w    = h->GetXaxis()->GetBinWidth(xbin);

      // this serves to calculate averages per bin, necessary in regions of fast oscillation
      // sumw is necessary because the flux can have a steep decrease in the range [E, E+dt]
      // thus, the sum of muon and electron atm. nu flux is used as a weight for averaging osc_flux
      // this can be understood by considering that fRand should not sample from a uniform energy
      // distribution, but from something resembling 1/E
      
      Double_t sumw = 0;
      for (Int_t n = 0; n < nsamples; n++) {

	// get costheta and energy random values inside this bin
	Double_t ct = fRand->Uniform( ct_low, ct_high );
	Double_t en = fRand->Uniform( en_low, en_high );
	
	// calculate and set the neutrino path
	fPrem->FillPath(ct);
	fProb->SetPath( fPrem->GetNuPath() );

	// effective mass (in units Ton)
	Int_t meff_xbin = fhMeff[flav][is_cc][is_nb]->GetXaxis()->FindBin( en );
	Int_t meff_ybin = fhMeff[flav][is_cc][is_nb]->GetYaxis()->FindBin( ct );
	Double_t meff   = fhMeff[flav][is_cc][is_nb]->GetBinContent(meff_xbin, meff_ybin);
	
	// atm neutrino count in operation time (in units 1/m^2)
	Double_t atm_flux_factor = en_w * ct_w * fSec_per_y * op_time;
	Double_t atm_flux_e = fFlux->Flux_dE_dcosz(0, is_nb, en, ct) * atm_flux_factor;
	Double_t atm_flux_m = fFlux->Flux_dE_dcosz(1, is_nb, en, ct) * atm_flux_factor;
	Double_t weight = atm_flux_e + atm_flux_m;
	sumw += weight;
	
	// oscillated neutrino count in operation time (in units 1/m^2)
	Double_t osc_flux = ( atm_flux_e * fProb->Prob(0, flav, en) +
			      atm_flux_m * fProb->Prob(1, flav, en) );

	// interacted neutrino count in operation time (in units 1/MTon)
	Double_t int_flux = osc_flux * fXsec->GetXsec(en)/fMN * fKg_per_Mton;

	// detected neutrino count in operation time (unitless)
	Double_t det_flux = int_flux * meff * 1e-6;

	// fill histograms
	Double_t atm_fluxes[3] = {atm_flux_e, atm_flux_m, 0. }; //to help filling atmflux hist

	fhAtmFlux[flav][is_cc][is_nb]->Fill( en, ct, atm_fluxes[flav]/nsamples );
	fhOscFlux[flav][is_cc][is_nb]->Fill( en, ct, osc_flux * weight );
	fhIntFlux[flav][is_cc][is_nb]->Fill( en, ct, int_flux * weight );
	fhDetFlux[flav][is_cc][is_nb]->Fill( en, ct, det_flux * weight );

      } // end loop over samples

      // divide by the sum of weights
      Double_t oflux = fhOscFlux[flav][is_cc][is_nb]->GetBinContent( xbin, ybin );
      Double_t iflux = fhIntFlux[flav][is_cc][is_nb]->GetBinContent( xbin, ybin );
      Double_t dflux = fhDetFlux[flav][is_cc][is_nb]->GetBinContent( xbin, ybin );
      
      fhOscFlux[flav][is_cc][is_nb]->SetBinContent( xbin, ybin, oflux/sumw );
      fhIntFlux[flav][is_cc][is_nb]->SetBinContent( xbin, ybin, iflux/sumw );
      fhDetFlux[flav][is_cc][is_nb]->SetBinContent( xbin, ybin, dflux/sumw );

      // also do this for errors
      Double_t oflux_err = fhOscFlux[flav][is_cc][is_nb]->GetBinError( xbin, ybin );
      Double_t iflux_err = fhIntFlux[flav][is_cc][is_nb]->GetBinError( xbin, ybin );
      Double_t dflux_err = fhDetFlux[flav][is_cc][is_nb]->GetBinError( xbin, ybin );
      fhOscFlux[flav][is_cc][is_nb]->SetBinError( xbin, ybin, oflux_err/sumw );
      fhIntFlux[flav][is_cc][is_nb]->SetBinError( xbin, ybin, iflux_err/sumw );
      fhDetFlux[flav][is_cc][is_nb]->SetBinError( xbin, ybin, dflux_err/sumw );

    } // end loop over xbins
    
  }// end loop over ybins
  
}

//*****************************************************************

/**
 *  Inline function to write to file.
 */
void WriteToFile(TString output_name) {

  TFile fout(output_name, "RECREATE");
  TDirectory *atmflux = fout.mkdir("atmflux");
  TDirectory *oscflux = fout.mkdir("oscflux");
  TDirectory *intflux = fout.mkdir("intflux");
  TDirectory *detflux = fout.mkdir("detflux");
  TDirectory *meff    = fout.mkdir("meff");
  
  for (Int_t f = 0; f < 3; f++) {
    for (Int_t cc = 0; cc < 2; cc++) {
      for (Int_t nb = 0; nb < 2; nb++) {
	atmflux->cd();
	fhAtmFlux[f][cc][nb]->Write();
	oscflux->cd();
	fhOscFlux[f][cc][nb]->Write();
	intflux->cd();
	fhIntFlux[f][cc][nb]->Write();
	detflux->cd();
	fhDetFlux[f][cc][nb]->Write();
	meff->cd();
	fhMeff[f][cc][nb]->Write();
      }
    }
  }
  
  fout.Close();
  
}

//*****************************************************************

/**
 *  Inline function to clean up dynamic memory.
 */
void CleanUp() {

  delete fProb;
  delete fPrem;
  delete fFlux;
  delete fXsec;
  delete fRand;
  
  for (Int_t f = 0; f < 3; f++) {
    for (Int_t cc = 0; cc < 2; cc++) {
      for (Int_t nb = 0; nb < 2; nb++) {
	delete fhAtmFlux[f][cc][nb];
	delete fhOscFlux[f][cc][nb];
	delete fhIntFlux[f][cc][nb];
	delete fhDetFlux[f][cc][nb];
	delete fhMeff[f][cc][nb];
      }
    }
  }
  
}

