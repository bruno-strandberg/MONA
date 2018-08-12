
#include "PMNS_Fast.h"
#include "PremModel.h"
#include "AtmFlux.h"
#include "NuXsec.h"
#include "NMHUtils.h"
#include "FileHeader.h"
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
Bool_t InitOscPars(Bool_t NH, FileHeader &h,
		   Double_t _sinsq_th12=0., Double_t _sinsq_th23=0., Double_t _sinsq_th13=0., 
		   Double_t _dcp=-10., Double_t _dm21=0., Double_t _dm32=0.);
void   SetToPDG(Bool_t NH, 
		Double_t &sinsq_th12, Double_t &sinsq_th23, Double_t &sinsq_th13, 
		Double_t &dcp, Double_t &dm21, Double_t &dm32, Double_t &dm31);
void   InitHists();
Bool_t ReadMeffHists(FileHeader &h, TString meff_elec_cc, TString meff_muon_cc, 
		     TString meff_tau_cc, TString meff_elec_nc);
void   FillHists(Double_t op_time, Int_t flav, Int_t is_cc, Int_t is_nb, Int_t nsamples);
void   CleanUp();
void   WriteToFile(TString output_name, FileHeader &h);

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
TH2D    *fhOscProb[3][3][2][2];//!< oscillation probabilities in format[flavor_in][flavor_out][is_cc][is_nb]

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
   This routine creates a file with E vs costheta histograms with expected numbers of nu events.
 
   The output root file will have directories atmflux, oscflux, intflux, detflux and meff. Eech
   of the flux directories contains E vs costheta histograms for {flavor}_{nc/cc}_{nu/nub} (12 in
   total), each bin will indicate the number of neutrinos of this type after op_time years.
   Directory atmflux/ has histograms with atmospheric neutrinos. oscflux/ shows atmflux/ histograms
   after oscillation through Earth, intflux/ shows oscflux/ histograms after nu + H2O xsec is taken
   into account, detflux/ shows intflux/ histograms after effective mass is taken into account.
   meff/ directory displays the effective mass histograms used for the generation of detflux/
   histograms. In the region of fast oscillations it is necessary to sample each bin several times,
   parameter nsamples determines how many samples per bin are calculated.
 
   \param  op_time       Operation time in years.
   \param  output_name   Name of the file where histograms are written.
   \param  nsamples      Number of samples per filling one bin (recommended > 10). If set to 1 then bin central values are used.
   \param  NH            True - normal nu mass hierarchy, False - inverted nu mass hierarchy.
   \param  UseMeff       If true, program will look for effective mass histograms in the input files to calculate detected number of events for detflux/ directory. If one inputs the FluxChain output to GSGSampler, meff is not necessary. If false, the input effective mass strings can be empty.
   \param meff_elec_cc   ROOT file with effective mass hists for elec-CC (create with `EffMass.C`)
   \param meff_muon_cc   ROOT file with effective mass hists for muon-CC (create with `EffMass.C`)
   \param meff_tau_cc    ROOT file with effective mass hists for tau-CC (create with `EffMass.C`)
   \param meff_elec_nc   ROOT file with effective mass hists for NC (create with `EffMass.C`)
   \param sinsq_th12     \f$\sin^2\theta_{12}\f$ value. If 0., PDG value is used.
   \param sinsq_th23     \f$\sin^2\theta_{23}\f$ value. If 0., PDG value is used.
   \param sinsq_th13     \f$\sin^2\theta_{13}\f$ value. If 0., PDG value is used.
   \param dcp            \f$\delta_{CP}\f$ value in \f$\pi\f$'s, as given by PDG group (e.g 1.38). If -10, PDG value is used.
   \param dm21           \f$\Delta m^2_{21}\f$ value. If 0., PDG value is used.
   \param dm32           \f$\Delta m^2_{32}\f$ value; this must be > 0 for NH and < 0 for IH. If 0., PDG value is used.

 */
void FluxChain(Double_t op_time      = 3.,
	       TString  output_name  = "flux_chain_out.root",
	       Int_t    nsamples     = 20,
	       Bool_t   NH           = true,
	       Bool_t   UseMeff      = false,
	       TString  meff_elec_cc = "../data/eff_mass/EffMass_elec_CC.root",
	       TString  meff_muon_cc = "../data/eff_mass/EffMass_muon_CC.root",
	       TString  meff_tau_cc  = "../data/eff_mass/EffMass_tau_CC.root",
	       TString  meff_elec_nc = "../data/eff_mass/EffMass_elec_NC.root",
	       Double_t sinsq_th12   = 0., 
	       Double_t sinsq_th23   = 0., 
	       Double_t sinsq_th13   = 0., 
	       Double_t dcp          = -10., 
	       Double_t dm21         = 0., 
	       Double_t dm32         = 0.) {

  gSystem->Load("$OSCPROBDIR/libOscProb.so");

  // create header
  FileHeader h("FluxChain");
  h.AddParameter("op_time"    , (TString)to_string(op_time) );
  h.AddParameter("output_name", output_name );
  h.AddParameter("nsamples"   , (TString)to_string(nsamples) );
  h.AddParameter("NH"         , (TString)to_string(NH) );
  h.AddParameter("UseMeff"    , (TString)to_string(UseMeff) );
  
  InitClasses();
  if ( !InitOscPars(NH, h, sinsq_th12, sinsq_th23, sinsq_th13, dcp, dm21, dm32) ) return;
  InitHists();
  if ( UseMeff ) {
    if ( !ReadMeffHists(h, meff_elec_cc, meff_muon_cc, meff_tau_cc, meff_elec_nc) ) return;
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

  WriteToFile(output_name, h);
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
 *  Inline function to initialise osc parameters and give them to the osc calculator.
 *
 * \param NH           - true for normal hierarchy, false for inverted hierarchy
 * \param sinsq_th12   - \f$\sin^2\theta_{12}\f$ value
 * \param sinsq_th23   - \f$\sin^2\theta_{23}\f$ value
 * \param sinsq_th13   - \f$\sin^2\theta_{13}\f$ value
 * \param dcp          - \f$\delta_{CP}\f$ value in \f$\pi\f$'s, as given by PDG group (e.g 1.38)
 * \param dm21         - \f$\Delta m^2_{21}\f$ value
 * \param dm32         - \f$\Delta m^2_{32}\f$ value; this must be > 0 for NH and < 0 for IH
 *
 */
Bool_t InitOscPars(Bool_t NH, FileHeader &h,
		   Double_t _sinsq_th12 , Double_t _sinsq_th23, Double_t _sinsq_th13, 
		   Double_t _dcp, Double_t _dm21, Double_t _dm32) {

  // init oscillation parameters

  Double_t sinsq_th12, sinsq_th23, sinsq_th13, dcp, dm21, dm32, dm31;
  
  SetToPDG(NH, sinsq_th12, sinsq_th23, sinsq_th13, dcp, dm21, dm32, dm31);

  // modify only these parameters for which value has been set
  if (_sinsq_th12 != 0.) sinsq_th12 = _sinsq_th12;
  if (_sinsq_th23 != 0.) sinsq_th23 = _sinsq_th23;
  if (_sinsq_th13 != 0.) sinsq_th13 = _sinsq_th13;
  if (_dcp  != -10.)     dcp        = _dcp;
  if (_dm21 != 0.)       dm21       = _dm21;
  if (_dm32 != 0.)       dm32       = _dm32;
  if (_dm21 != 0. || _dm32 != 0.) dm31 = dm21 + dm32;

  if ( ( NH && (_dm32 < 0) ) || ( !NH && (_dm32 > 0) ) ) { 
    cout << "ERROR! InitOscPars() wrong sign of dm32 " << dm32 << " for hierarchy (0=IH, 1=NH) " 
	 << NH << ", should be dm32 > 0 for NH and dm32 < 0 for IH."<< endl;
    return false;
  }

  //------------------------------------------------
  // Pass to the oscillation probability calculator, save them to header
  //------------------------------------------------
  h.AddParameter("sinsq_th12", (TString)to_string(sinsq_th12) );
  h.AddParameter("sinsq_th23", (TString)to_string(sinsq_th23) );
  h.AddParameter("sinsq_th13", (TString)to_string(sinsq_th13) );
  h.AddParameter("dcp"       , (TString)to_string(dcp) );
  h.AddParameter("dm21_1e5"  , (TString)to_string(dm21*1e5) );
  h.AddParameter("dm31_1e5"  , (TString)to_string(dm31*1e5) );
  h.AddParameter("dm32_1e5"  , (TString)to_string(dm32*1e5) );

  fProb->SetAngle(1, 2, TMath::ASin( TMath::Sqrt(sinsq_th12) ) );
  fProb->SetAngle(1, 3, TMath::ASin( TMath::Sqrt(sinsq_th13) ) );
  fProb->SetAngle(2, 3, TMath::ASin( TMath::Sqrt(sinsq_th23) ) );
  fProb->SetDelta(1, 3, dcp * TMath::Pi() );

  fProb->SetDm(2, dm21);
  fProb->SetDm(3, dm31);
  
  return true;
}

//*****************************************************************

/**
 * Inline function to set oscillation parameters to PDG values, depending on hierarchy.
 *
 * Values from http://pdg.lbl.gov/2017/reviews/rpp2017-rev-neutrino-mixing.pdf
 *
 * \param NH           - true for normal hierarchy, false for inverted hierarchy
 * \param sinsq_th12   - reference to a variable where \f$\sin^2\theta_{12}\f$ is stored
 * \param sinsq_th23   - reference to a variable where \f$\sin^2\theta_{23}\f$ is stored
 * \param sinsq_th13   - reference to a variable where \f$\sin^2\theta_{13}\f$ is stored
 * \param dcp          - reference to a variable where \f$\delta_{CP}\f$ is stored (in \f$\pi\f$'s)
 * \param dm21         - reference to a variable where \f$\Delta m^2_{21}\f$ is stored
 * \param dm32         - reference to a variable where \f$\Delta m^2_{32}\f$ is stored
 * \param dm31         - reference to a variable where \f$\Delta m^2_{31}\f$ is stored
 */
void SetToPDG(Bool_t NH, 
	      Double_t &sinsq_th12, Double_t &sinsq_th23, Double_t &sinsq_th13, 
	      Double_t &dcp, Double_t &dm21, Double_t &dm32, Double_t &dm31) {

  sinsq_th12  = 0.297;
  sinsq_th23  = 0.425;
  sinsq_th13  = 0.0215;
  dcp         = 1.38;
  dm21        = 7.37e-5;
  dm31        = 2.56e-3;
  dm32        = dm31 - dm21;

  if (!NH) {
    sinsq_th23 = 0.589;
    sinsq_th13 = 0.0216;
    dcp        = 1.31;
    dm32       = -2.54e-3;
    dm31       = dm32 + dm21;
  }
  
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
	fhMeff[f][cc][nb] = NULL;

	fhAtmFlux[f][cc][nb]->Sumw2();
	fhOscFlux[f][cc][nb]->Sumw2();
	fhIntFlux[f][cc][nb]->Sumw2();
	fhDetFlux[f][cc][nb]->Sumw2();

	for (Int_t _f = 0; _f < 3; _f++) {
	  TString oscprobn = "oscprob_" + fFlavs[f] + "_to_" + fFlavs[_f] + "_" + fItypes[cc] + "_" + fPtypes[nb];
	  fhOscProb[f][_f][cc][nb] = new TH2D(oscprobn, oscprobn, ebins, &e_edges[0], ctbins, ctlow, cthigh);
	  fhOscProb[f][_f][cc][nb]->Sumw2();
	}
       
      }
    }
  }
  
}

//*****************************************************************

/**
 *  Inline function to read effective mass histograms from $NMHDIR/data/eff_mass/.
 */
Bool_t ReadMeffHists(FileHeader &h, TString meff_elec_cc, TString meff_muon_cc, 
		     TString meff_tau_cc, TString meff_elec_nc) {

  vector<TString> meff_filenames = {meff_elec_cc, meff_muon_cc, meff_tau_cc};

  for (Int_t f = 0; f < 3; f++) {
    for (Int_t cc = 0; cc < 2; cc++) {

      // for nc events the effective mass is identical for all flavors, only elec_NC simulated
      TString meff_fname = meff_filenames[f];
      if (cc == 0) meff_fname = meff_elec_nc;
    
      TFile meff_file(meff_fname, "READ");
      if ( !meff_file.IsOpen() ) {
	cout << "ERROR! InitHists() could not find file " << meff_fname << endl;
	return false;
      }

      // read the header and append to the header of this application
      h.ReadHeader(meff_fname);

      TString hname = "Meff_" + fFlavs[f] + "_" + fItypes[cc];
      fhMeff[f][cc][0] = (TH2D*)meff_file.Get("Meff_nu")->Clone(hname + "_nu");
      fhMeff[f][cc][1] = (TH2D*)meff_file.Get("Meff_nub")->Clone(hname + "_nub");
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
	
	// as an exception, if there is only one sample, use bin center
	if (nsamples == 1) {
	  ct = h->GetYaxis()->GetBinCenter( ybin );
	  en = h->GetXaxis()->GetBinCenter( xbin );
	}

	// calculate and set the neutrino path
	fPrem->FillPath(ct);
	fProb->SetPath( fPrem->GetNuPath() );

	// effective mass (in units Ton); if meff not used then hists not defined and detflux
	// histograms will be 0
	Double_t meff = 0.;
	if ( fhMeff[flav][is_cc][is_nb] ) {
	  Int_t meff_xbin = fhMeff[flav][is_cc][is_nb]->GetXaxis()->FindBin( en );
	  Int_t meff_ybin = fhMeff[flav][is_cc][is_nb]->GetYaxis()->FindBin( ct );
	  meff = fhMeff[flav][is_cc][is_nb]->GetBinContent(meff_xbin, meff_ybin);
	}
	
	// atm neutrino count in operation time (in units 1/m^2)
	Double_t atm_flux_factor = en_w * ct_w * fSec_per_y * op_time;
	Double_t atm_flux_e = fFlux->Flux_dE_dcosz(0, is_nb, en, ct) * atm_flux_factor;
	Double_t atm_flux_m = fFlux->Flux_dE_dcosz(1, is_nb, en, ct) * atm_flux_factor;
	Double_t weight = atm_flux_e + atm_flux_m;
	sumw += weight;
	
	// oscillated neutrino count in operation time (in units 1/m^2)
	Double_t prob_elec = fProb->Prob(0, flav, en);
	Double_t prob_muon = fProb->Prob(1, flav, en);
	Double_t osc_flux = ( atm_flux_e * prob_elec + atm_flux_m * prob_muon );

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
	fhOscProb[0][flav][is_cc][is_nb]->Fill( en, ct, prob_elec * weight);
	fhOscProb[1][flav][is_cc][is_nb]->Fill( en, ct, prob_muon * weight);

      } // end loop over samples

      // divide by the sum of weights
      Double_t oflux = fhOscFlux[flav][is_cc][is_nb]->GetBinContent( xbin, ybin );
      Double_t iflux = fhIntFlux[flav][is_cc][is_nb]->GetBinContent( xbin, ybin );
      Double_t dflux = fhDetFlux[flav][is_cc][is_nb]->GetBinContent( xbin, ybin );
      Double_t oprob_elec = fhOscProb[0][flav][is_cc][is_nb]->GetBinContent( xbin, ybin );
      Double_t oprob_muon = fhOscProb[1][flav][is_cc][is_nb]->GetBinContent( xbin, ybin );
      
      fhOscFlux[flav][is_cc][is_nb]->SetBinContent( xbin, ybin, oflux/sumw );
      fhIntFlux[flav][is_cc][is_nb]->SetBinContent( xbin, ybin, iflux/sumw );
      fhDetFlux[flav][is_cc][is_nb]->SetBinContent( xbin, ybin, dflux/sumw );
      fhOscProb[0][flav][is_cc][is_nb]->SetBinContent( xbin, ybin, oprob_elec/sumw );
      fhOscProb[1][flav][is_cc][is_nb]->SetBinContent( xbin, ybin, oprob_muon/sumw );

      // also do this for errors
      Double_t oflux_err = fhOscFlux[flav][is_cc][is_nb]->GetBinError( xbin, ybin );
      Double_t iflux_err = fhIntFlux[flav][is_cc][is_nb]->GetBinError( xbin, ybin );
      Double_t dflux_err = fhDetFlux[flav][is_cc][is_nb]->GetBinError( xbin, ybin );
      Double_t eprob_err = fhOscProb[0][flav][is_cc][is_nb]->GetBinError( xbin, ybin );
      Double_t mprob_err = fhOscProb[1][flav][is_cc][is_nb]->GetBinError( xbin, ybin );
      fhOscFlux[flav][is_cc][is_nb]->SetBinError( xbin, ybin, oflux_err/sumw );
      fhIntFlux[flav][is_cc][is_nb]->SetBinError( xbin, ybin, iflux_err/sumw );
      fhDetFlux[flav][is_cc][is_nb]->SetBinError( xbin, ybin, dflux_err/sumw );
      fhOscProb[0][flav][is_cc][is_nb]->SetBinError( xbin, ybin, eprob_err/sumw );
      fhOscProb[1][flav][is_cc][is_nb]->SetBinError( xbin, ybin, mprob_err/sumw );

    } // end loop over xbins
    
  }// end loop over ybins
  
}

//*****************************************************************

/**
 *  Inline function to write to file.
 */
 void WriteToFile(TString output_name, FileHeader &h) {

  TFile fout(output_name, "RECREATE");

  h.WriteHeader(&fout);

  TDirectory *atmflux = fout.mkdir("atmflux");
  TDirectory *oscflux = fout.mkdir("oscflux");
  TDirectory *intflux = fout.mkdir("intflux");
  TDirectory *detflux = fout.mkdir("detflux");
  TDirectory *meff    = fout.mkdir("meff");
  TDirectory *oscprob = fout.mkdir("oscprob");
  
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
	if (fhMeff[f][cc][nb]) fhMeff[f][cc][nb]->Write();
	oscprob->cd();
	for (Int_t _f = 0; _f < 3; _f++) fhOscProb[f][_f][cc][nb]->Write();
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
	if (fhMeff[f][cc][nb]) delete fhMeff[f][cc][nb];
	for (Int_t _f = 0; _f < 3; _f++) delete fhOscProb[f][_f][cc][nb];
      }
    }
  }
  
}

