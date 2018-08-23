
#include "PMNS_Fast.h"
#include "PremModel.h"

#include "AtmFlux.h"
#include "NuXsec.h"
#include "NMHUtils.h"
#include "FileHeader.h"

//#include "TString.h"
#include "TSystem.h"
#include "TMath.h"
#include "TH2.h"
#include "TH3.h"
#include "TRandom3.h"
#include "TFile.h"

#include <iostream>
#include <stdexcept>

/**
   Namespace that stores functions and globally used variables/members for `FluxChain.C`.
 */
namespace FLUXCHAIN {

  //*****************************************************************
  // functions
  //*****************************************************************
  void   InitClasses();
  void   InitOscPars(Bool_t NH, FileHeader &h,
		     Double_t _sinsq_th12=0., Double_t _sinsq_th23=0., Double_t _sinsq_th13=0., 
		     Double_t _dcp=-10., Double_t _dm21=0., Double_t _dm32=0.);
  void   SetToPDG(Bool_t NH, 
		  Double_t &sinsq_th12, Double_t &sinsq_th23, Double_t &sinsq_th13, 
		  Double_t &dcp, Double_t &dm21, Double_t &dm32, Double_t &dm31);
  void   InitHists();
  void   ReadMeffHists(FileHeader &h, TString meff_elec_cc, TString meff_muon_cc, 
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
  TH3D    *fhDetFlux[3][2][2];   //!< detected flux histograms
  TH3D       *fhMeff[3][2][2];   //!< effective mass histograms
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

  // binning for histograms
  Int_t    fEbins  =  40;
  Double_t fEmin   =   1.;
  Double_t fEmax   = 100.;

  Int_t    fCtbins =  40;
  Double_t fCtmin  =  -1.;
  Double_t fCtmax  =   1.;

  Int_t    fBybins =   4;
  Double_t fBymin  =   0.;
  Double_t fBymax  =   1.;

};

using namespace std;
using namespace FLUXCHAIN;

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
   \param meff_elec_cc   ROOT file with 'detected' and 'generated' hists for elec-CC (create with `EffMhists.C`)
   \param meff_muon_cc   ROOT file with 'detected' and 'generated' hists for muon-CC (create with `EffMhists.C`)
   \param meff_tau_cc    ROOT file with 'detected' and 'generated' hists for tau-CC (create with `EffMhists.C`)
   \param meff_elec_nc   ROOT file with 'detected' and 'generated' hists for NC (create with `EffMhists.C`)
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
	       TString  meff_elec_cc = "../data/eff_mass/EffMhists_elec_CC.root",
	       TString  meff_muon_cc = "../data/eff_mass/EffMhists_muon_CC.root",
	       TString  meff_tau_cc  = "../data/eff_mass/EffMhists_tau_CC.root",
	       TString  meff_elec_nc = "../data/eff_mass/EffMhists_elec_NC.root",
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
  InitOscPars(NH, h, sinsq_th12, sinsq_th23, sinsq_th13, dcp, dm21, dm32);
  InitHists();
  if ( UseMeff ) ReadMeffHists(h, meff_elec_cc, meff_muon_cc, meff_tau_cc, meff_elec_nc);
 
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
void FLUXCHAIN::InitClasses() {

  fProb = new OscProb::PMNS_Fast;
  fPrem = new OscProb::PremModel;
  fFlux = new AtmFlux;
  fXsec = new NuXsec(fBybins);
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
void FLUXCHAIN::InitOscPars(Bool_t NH, FileHeader &h,
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
    throw std::invalid_argument( "ERROR! FLUXCHAIN::InitOscPars() wrong sign of dm32 " + to_string(dm32) + " for hierarchy (0=IH, 1=NH) " 
				 + to_string(NH) + ", should be dm32 > 0 for NH and dm32 < 0 for IH." );
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
void FLUXCHAIN::SetToPDG(Bool_t NH, 
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
void FLUXCHAIN::InitHists() {

  vector<Double_t> e_edges  = NMHUtils::GetLogBins(fEbins, fEmin, fEmax);
  vector<Double_t> ct_edges = NMHUtils::GetBins(fCtbins, fCtmin, fCtmax);
  vector<Double_t> by_edges = NMHUtils::GetBins(fBybins, fBymin, fBymax);
  
  for (Int_t f = 0; f < 3; f++) {
    for (Int_t cc = 0; cc < 2; cc++) {
      for (Int_t nb = 0; nb < 2; nb++) {
      
	TString atmname = "atmflux_" + fFlavs[f] + "_" + fItypes[cc] + "_" + fPtypes[nb];
	TString oscname = "oscflux_" + fFlavs[f] + "_" + fItypes[cc] + "_" + fPtypes[nb];
	TString intname = "intflux_" + fFlavs[f] + "_" + fItypes[cc] + "_" + fPtypes[nb];
	TString detname = "detflux_" + fFlavs[f] + "_" + fItypes[cc] + "_" + fPtypes[nb];

	fhAtmFlux[f][cc][nb] = new TH2D(atmname, atmname, fEbins, &e_edges[0], fCtbins, &ct_edges[0]);
	fhOscFlux[f][cc][nb] = new TH2D(oscname, oscname, fEbins, &e_edges[0], fCtbins, &ct_edges[0]);
	fhIntFlux[f][cc][nb] = new TH2D(intname, intname, fEbins, &e_edges[0], fCtbins, &ct_edges[0]);
	fhDetFlux[f][cc][nb] = new TH3D(detname, detname, fEbins, &e_edges[0], fCtbins, &ct_edges[0], fBybins, &by_edges[0]);
	fhMeff[f][cc][nb] = NULL;

	fhAtmFlux[f][cc][nb]->Sumw2();
	fhOscFlux[f][cc][nb]->Sumw2();
	fhIntFlux[f][cc][nb]->Sumw2();
	fhDetFlux[f][cc][nb]->Sumw2();

	for (Int_t _f = 0; _f < 3; _f++) {
	  TString oscprobn = "oscprob_" + fFlavs[f] + "_to_" + fFlavs[_f] + "_" + fItypes[cc] + "_" + fPtypes[nb];
	  fhOscProb[f][_f][cc][nb] = new TH2D(oscprobn, oscprobn, fEbins, &e_edges[0], fCtbins, &ct_edges[0]);
	  fhOscProb[f][_f][cc][nb]->Sumw2();
	}
       
      }
    }
  }
  
}

//*****************************************************************

/**
   Inline function to read effective mass histograms.

   \param h             Address of a FileHeader instance for information storage
   \param meff_elec_cc  Name of the effective mass file for elec CC
   \param meff_muon_cc  Name of the effective mass file for muon CC
   \param meff_tau_cc   Name of the effective mass file for tau CC
   \param meff_elec_nc  Name of the effective mass file for elec NC

   \return True if effective mass histograms successfully found, False otherwise.
 */
void FLUXCHAIN::ReadMeffHists(FileHeader &h, TString meff_elec_cc, TString meff_muon_cc, 
				TString meff_tau_cc, TString meff_elec_nc) {
  
  vector<TString> meff_filenames = {meff_elec_cc, meff_muon_cc, meff_tau_cc};

  for (Int_t f = 0; f < 3; f++) {
    for (Int_t cc = 0; cc < 2; cc++) {

      // for nc events the effective mass is identical for all flavors, only elec_NC simulated
      TString meff_fname = meff_filenames[f];
      if (cc == 0) meff_fname = meff_elec_nc;
    
      TFile meff_file(meff_fname, "READ");
      if ( !meff_file.IsOpen() ) {
	throw std::invalid_argument( "ERROR! FLUXCHAIN::ReadMeffHists() could not find file " + (string)meff_fname );
      }

      //---------------------------------------------------------------
      // read the header from EffMhists output and manually add some parameters to the header
      // of this app. Typically EffMhists outputs are hadd'ed, adding header fields manually
      // removes duplicate entries
      //---------------------------------------------------------------
      
      FileHeader emh("emhheader");     // read the header from effmhists output
      emh.ReadHeader(meff_fname);

      h.AddParameter( emh, "Rvol" );   // add parameter Rvol from emh to h
      h.AddParameter( emh, "Zmin" );
      h.AddParameter( emh, "Zmax" );
      h.AddParameter( emh, "atmmu_cut" );
      h.AddParameter( emh, "noise_cut" );

      //---------------------------------------------------------------
      // get the histograms, rebin and divide to get effective mass histos
      //---------------------------------------------------------------

      TString hname = "Meff_" + fFlavs[f] + "_" + fItypes[cc];

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

      Int_t rebinning_ebins  = existing_ebins /fEbins;
      Int_t rebinning_ctbins = existing_ctbins/fCtbins;
      Int_t rebinning_bybins = existing_bybins/fBybins;

      // check that rebinning is valid
      if ( existing_ebins % fEbins != 0 ) {
	throw std::invalid_argument( "ERROR! FLUXCHAIN::ReadMeffHists() energy axis nbins=" + to_string(existing_ebins) + " of file " + (string)meff_fname + " cannot be rebinned to " + to_string(fEbins) + ". Change fEbins in FluxChain.C or change and re-run NMHDIR/effective_mass/EffMhists.C with suitable binning." );
      }

      if ( existing_ctbins % fCtbins != 0 ) {
	throw std::invalid_argument( "ERROR! FLUXCHAIN::ReadMeffHists() costheta axis nbins=" + to_string(existing_ctbins) + " of file " + (string)meff_fname + " cannot be rebinned to " + to_string(fCtbins) + ". Change fCtbins in FluxChain.C or change and re-run NMHDIR/effective_mass/EffMhists.C with suitable binning." );
      }
      
      if ( existing_bybins % fBybins != 0 ) {
	throw std::invalid_argument( "ERROR! FLUXCHAIN::ReadMeffHists() bjorken-y axis nbins=" + to_string(existing_bybins) + " of file " + (string)meff_fname + " cannot be rebinned to " + to_string(fBybins) + ". Change fBybins in FluxChain.C or change and re-run NMHDIR/effective_mass/EffMhists.C with suitable binning." );
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

//*****************************************************************

/**
   Inline function to fill E vs costheta histograms of this macro.

   \param op_time   Operation time
   \param flav      Nu flavor
   \param is_cc     CC or NC event
   \param is_nb     nu or nubar
   \param nsamples  Number of over-samples per bin
 */
void FLUXCHAIN::FillHists(Double_t op_time, Int_t flav, Int_t is_cc, Int_t is_nb, Int_t nsamples) {
  
  fProb->SetIsNuBar(is_nb);
  fXsec->SelectInteraction(flav, is_cc, is_nb);
  
  TH3D *h = fhDetFlux[flav][is_cc][is_nb]; //pointer for convenient access to bin data
  
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

	// fill hists from atm. flux to int. flux and osc. probs
	// define array to help filling atmflux hist
	Double_t atm_fluxes[3] = {atm_flux_e, atm_flux_m, 0. };

	fhAtmFlux[flav][is_cc][is_nb]->Fill( en, ct, atm_fluxes[flav]/nsamples );
	fhOscFlux[flav][is_cc][is_nb]->Fill( en, ct, osc_flux * weight );
	fhIntFlux[flav][is_cc][is_nb]->Fill( en, ct, int_flux * weight );
	fhOscProb[0][flav][is_cc][is_nb]->Fill( en, ct, prob_elec * weight);
	fhOscProb[1][flav][is_cc][is_nb]->Fill( en, ct, prob_muon * weight);

	// distribute the interacted flux in each E, ct bin to bjorken-y bins
	for (Int_t zbin = 1; zbin <= h->GetZaxis()->GetNbins(); zbin++) {

	  Double_t by = h->GetZaxis()->GetBinCenter(zbin);
	  
	  // effective mass (in units Ton); if meff not used then hists not defined and detflux histograms will be 0
	  Double_t meff = 0.;
	  if ( fhMeff[flav][is_cc][is_nb] ) {
	    Int_t meff_xbin = fhMeff[flav][is_cc][is_nb]->GetXaxis()->FindBin( en );
	    Int_t meff_ybin = fhMeff[flav][is_cc][is_nb]->GetYaxis()->FindBin( ct );
	    Int_t meff_zbin = fhMeff[flav][is_cc][is_nb]->GetZaxis()->FindBin( by );
	    meff = fhMeff[flav][is_cc][is_nb]->GetBinContent(meff_xbin, meff_ybin, meff_zbin);
	  }

	  // detected neutrino count in operation time (unitless)
	  Double_t det_flux = int_flux * fXsec->GetBYfrac(en, by) * meff * 1e-6;
	  fhDetFlux[flav][is_cc][is_nb]->Fill( en, ct, by, det_flux * weight );

	} // end loop over bjorken-y bins
	
      } // end loop over samples
      
      // divide by the sum of weights
      Double_t oflux      = fhOscFlux[flav][is_cc][is_nb]   ->GetBinContent( xbin, ybin );
      Double_t iflux      = fhIntFlux[flav][is_cc][is_nb]   ->GetBinContent( xbin, ybin );
      Double_t oprob_elec = fhOscProb[0][flav][is_cc][is_nb]->GetBinContent( xbin, ybin );
      Double_t oprob_muon = fhOscProb[1][flav][is_cc][is_nb]->GetBinContent( xbin, ybin );
      
      fhOscFlux[flav][is_cc][is_nb]   ->SetBinContent( xbin, ybin, oflux/sumw );
      fhIntFlux[flav][is_cc][is_nb]   ->SetBinContent( xbin, ybin, iflux/sumw );
      fhOscProb[0][flav][is_cc][is_nb]->SetBinContent( xbin, ybin, oprob_elec/sumw );
      fhOscProb[1][flav][is_cc][is_nb]->SetBinContent( xbin, ybin, oprob_muon/sumw );

      // also do this for errors
      Double_t oflux_err = fhOscFlux[flav][is_cc][is_nb]->GetBinError( xbin, ybin );
      Double_t iflux_err = fhIntFlux[flav][is_cc][is_nb]->GetBinError( xbin, ybin );
      Double_t eprob_err = fhOscProb[0][flav][is_cc][is_nb]->GetBinError( xbin, ybin );
      Double_t mprob_err = fhOscProb[1][flav][is_cc][is_nb]->GetBinError( xbin, ybin );
      fhOscFlux[flav][is_cc][is_nb]->SetBinError( xbin, ybin, oflux_err/sumw );
      fhIntFlux[flav][is_cc][is_nb]->SetBinError( xbin, ybin, iflux_err/sumw );
      fhOscProb[0][flav][is_cc][is_nb]->SetBinError( xbin, ybin, eprob_err/sumw );
      fhOscProb[1][flav][is_cc][is_nb]->SetBinError( xbin, ybin, mprob_err/sumw );

      for (Int_t zbin = 1; zbin <= h->GetZaxis()->GetNbins(); zbin++) {
	Double_t dflux = fhDetFlux[flav][is_cc][is_nb]->GetBinContent( xbin, ybin, zbin );
	Double_t dflux_err = fhDetFlux[flav][is_cc][is_nb]->GetBinError( xbin, ybin, zbin );
	fhDetFlux[flav][is_cc][is_nb]->SetBinContent( xbin, ybin, zbin, dflux/sumw );
	fhDetFlux[flav][is_cc][is_nb]->SetBinError( xbin, ybin, zbin, dflux_err/sumw );
      }
      
    } // end loop over xbins
    
  }// end loop over ybins
  
}

//*****************************************************************

/**
   Inline function to write histograms to file.

   \param output_name   Name of the output file
   \param h             Reference to a FileHeader instance that is written to output
 */
void FLUXCHAIN::WriteToFile(TString output_name, FileHeader &h) {

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
   Inline function to clean up dynamic memory.
 */
void FLUXCHAIN::CleanUp() {

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

