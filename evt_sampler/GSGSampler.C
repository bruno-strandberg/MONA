#include "SummaryParser.h"
#include "GSGParser.h"
#include "TSystem.h"
#include "TH2.h"
#include "TRandom3.h"

//*****************************************************************
// functions
//*****************************************************************
Bool_t GetIntHists(TString flux_chain_file, Int_t flavor, Int_t is_cc);
void   InitVars(Int_t flavor, Int_t is_cc);
void   CleanUp();

//*****************************************************************
// globally used variables in this script
//*****************************************************************
TH2D *fhInt_nu;        //!< histogram with expected number of interacted nu events
TH2D *fhInt_nub;       //!< histogram with expected number of interacted nubar events
TH2D *fhGSG_nu;        //!< histogram with all available MC nu events in GSG files
TH2D *fhGSG_nub;       //!< histogram with all available MC nubar events in GSG files
SummaryParser *fSp;    //!< class to parse summary events
TRandom3      *fRand;  //!< random number generator
Int_t fEbins;          //!< number of energy bins in the event vectors
Int_t fCtbins;         //!< number of costheta bins in the event vectors
//! 2D array of vectors, each vector holds summary nu event numbers of (E, costheta) bin
vector<Double_t> **fGSGEvts_nu;
//! 2D array of vectors, each vector holds summary nub event numbers of (E, costheta) bin
vector<Double_t> **fGSGEvts_nub;
//! 2D array of vectors to hold a sub-sample of events in fGSGEvts_nu
vector<Double_t> **fSampleEvts_nu;
//! 2D array of vectors to hold a sub-sample of events in fGSGEvts_nub
vector<Double_t> **fSampleEvts_nub;

//! map of flavor numbers and strings
map < Int_t, TString > fFlavs  = { {0, "elec" },
				   {1, "muon" },
				   {2, "tau"  } };

//*****************************************************************
// main function
//*****************************************************************

void GSGSampler(TString flux_chain_file, TString gsg_file_list, Int_t flavor, Int_t is_cc, Int_t nsamples = 1) {

  gSystem->Load("$NMHDIR/common_software/libnmhsoft.so");

  if ( !GetDetHists(flux_chain_file, flavor, is_cc) ) {
    cout << "ERROR! GSGSampler() problem opening flux histograms." << endl;
    return;
  }

  InitVars(flavor, is_cc);
  
  
}

//*****************************************************************

/**
 *  Inline function to fetch the histograms with # of interactions in operation time per MTon
 *  from FluxChain.C output.
 *
 *  The histograms in the intflux/ directory of the FluxChain.C output are in the units
 *  [number of events in operation time/MTon].
 *
 * \param  flux_chain_file  Output of FluxChain.C.
 * \param  flavor           Neutrino flavor.
 * \param  is_cc            0 - NC, 1 -CC.
 * \return                  True if histograms found, False if not found.
 *
 */
Bool_t GetIntHists(TString flux_chain_file, Int_t flavor, Int_t is_cc) {

  TFile f(flux_chain_file, "READ");
  if ( !f.IsOpen() ) {
    cout << "ERROR! GetIntHists() cannot open file " << flux_chain_file << endl;
    return false;
  }

  if (is_cc > 0) {

    TString hname_nu  = "intflux/intflux_" + fFlavs[flavor] + "_cc_nu";
    TString hname_nub = "intflux/intflux_" + fFlavs[flavor] + "_cc_nub";
    TH2D   *h_nu  = (TH2D*)f.Get(hname_nu);
    TH2D   *h_nub = (TH2D*)f.Get(hname_nub);

    if ( !h_nu || !h_nub ) {
      cout << "ERROR! GetIntHists() cannot find hists " << hname_nu << "\t" << hname_nub << endl;
      return false;
    }

    fhInt_nu  = (TH2D*)h_nu->Clone();
    fhInt_nub = (TH2D*)h_nub->Clone();
    fhInt_nu->SetDirectory(0);
    fhInt_nub->SetDirectory(0);

  }
  else {

    //for nc events we only have MC events for elec_NC, hence need to sample them together
    TString hname_elec_nu  = "intflux/intflux_elec_nc_nu";
    TString hname_elec_nub = "intflux/intflux_elec_nc_nub";
    TString hname_muon_nu  = "intflux/intflux_muon_nc_nu";
    TString hname_muon_nub = "intflux/intflux_muon_nc_nub";
    TString hname_tau_nu   = "intflux/intflux_tau_nc_nu";
    TString hname_tau_nub  = "intflux/intflux_tau_nc_nub";

    TH2D   *h_elec_nu  = (TH2D*)f.Get(hname_elec_nu);
    TH2D   *h_elec_nub = (TH2D*)f.Get(hname_elec_nub);
    TH2D   *h_muon_nu  = (TH2D*)f.Get(hname_muon_nu);
    TH2D   *h_muon_nub = (TH2D*)f.Get(hname_muon_nub);
    TH2D   *h_tau_nu   = (TH2D*)f.Get(hname_tau_nu);
    TH2D   *h_tau_nub  = (TH2D*)f.Get(hname_tau_nub);

    if (!h_elec_nu || !h_elec_nub || !h_muon_nu || !h_muon_nub || !h_tau_nu || !h_tau_nub) {
      cout << "ERROR! GetIntHists() cannot find NC hists" << endl;
      return false;
    }

    fhInt_nu  = (TH2D*)h_elec_nu->Clone("intflux_allflav_nc_nu");
    fhInt_nu->Add(h_muon_nu);
    fhInt_nu->Add(h_tau_nu);

    fhInt_nub = (TH2D*)h_elec_nub->Clone("intflux_allflav_nc_nub");
    fhInt_nub->Add(h_muon_nub);
    fhInt_nub->Add(h_tau_nub);

    fhInt_nu->SetTitle("intflux_allflav_nc_nu");
    fhInt_nub->SetTitle("intflux_allflav_nc_nub");

    fhInt_nu->SetDirectory(0);
    fhInt_nub->SetDirectory(0);
    
  }

  f.Close();
  return true;
  
}

//*****************************************************************

/**
 *  Inline function to initialize globally used classes, histograms and vectors.
 *
 * \param  flavor           Neutrino flavor.
 * \param  is_cc            0 - NC, 1 -CC.
 *
 */
void InitVars(Int_t flavor, Int_t is_cc) {

  // create hist names; we only simulate NC for elec, so there is an exception for that
  
  TString hname_nu  = "gsgevts_"  + fFlavs[flavor] + "_cc_nu";
  TString hname_nub = "gsgevts_"  + fFlavs[flavor] + "_cc_nub";
  
  if (is_cc == 0) {
    hname_nu  = "gsgevts_allflav_nc_nu";
    hname_nub = "gsgevts_allflav_nc_nub";
  }

  // histograms to display distributions of all available gsg events
  
  fhGSG_nu  = (TH2D*)fhInt_nu->Clone(hname_nu);
  fhGSG_nub = (TH2D*)fhInt_nub->Clone(hname_nub); 
  fhGSG_nu ->Reset();
  fhGSG_nub->Reset();
  fhGSG_nu ->SetTitle(hname_nu);
  fhGSG_nub->SetTitle(hname_nub);

  // allocate vectors, depending on the input hist binning, to store gsg event IDs
  // bin [0] is underflow, bin[ X/Yaxis->GetNbins() ] is the last counting bin (hence array
  // length +1), bin[ X/Yaxis->GetNbins()+1 ] is overflow (hence array length + 2)
  fEbins  = fhGSG_nu->GetXaxis()->GetNbins() + 2;
  fCtbins = fhGSG_nu->GetYaxis()->GetNbins() + 2;
  
  fGSGEvts_nu     = new vector<Double_t>* [fEbins];
  fGSGEvts_nub    = new vector<Double_t>* [fEbins];
  fSampleEvts_nu  = new vector<Double_t>* [fEbins];
  fSampleEvts_nub = new vector<Double_t>* [fEbins];
  for (Int_t eb = 0; eb < fEbins; eb++) {
    fGSGEvts_nu[eb]     = new vector<Double_t> [fCtbins];
    fGSGEvts_nub[eb]    = new vector<Double_t> [fCtbins];
    fSampleEvts_nu[eb]  = new vector<Double_t> [fCtbins];
    fSampleEvts_nub[eb] = new vector<Double_t> [fCtbins];
  }

  // initialize the summary parser and random generator
  fSp = new SummaryParser;
  fRand = new TRandom3(0);
}

//*****************************************************************

void ReadGSGData(TString gsg_file_list) {

  vector<TString> fnames = NMHUtils::ReadLines(gsg_file_list);

  for (auto fname: fnames) {

    GSGParser gp(fname);

    while ( gp.NextEvent() ) {

    }
    
  }
  
}

//*****************************************************************

/**
 *  Inline function to clear dynamically allocated memory.
 */
void CleanUp() {

  //clean up the dynamically allocated vectors/arrays
  
  if (fGSGEvts_nu && fGSGEvts_nub && fhSum_nu) {

    for (Int_t xbin = 0; xbin < fEbins; xbin++) {
      if ( fGSGEvts_nu[xbin]     ) delete[] fGSGEvts_nu[xbin];
      if ( fGSGEvts_nub[xbin]    ) delete[] fGSGEvts_nub[xbin];
      if ( fSampleEvts_nu[xbin]  ) delete[] fSampleEvts_nu[xbin];
      if ( fSampleEvts_nub[xbin] ) delete[] fSampleEvts_nub[xbin];

    }

    if ( fGSGEvts_nu     ) delete[] fGSGEvts_nu;
    if ( fGSGEvts_nub    ) delete[] fGSGEvts_nub;
    if ( fSampleEvts_nu  ) delete[] fSampleEvts_nu;
    if ( fSampleEvts_nub ) delete[] fSampleEvts_nub;
    
  }

  //clean up hists/classes/etc
  
  if (fhInt_nu)  delete fhInt_nu; 
  if (fhInt_nub) delete fhInt_nub;
  if (fhGSG_nu)  delete fhGSG_nu; 
  if (fhGSG_nub) delete fhGSG_nub;

  if (fSp)   delete fSp;
  if (fRand) delete fRand;
  
}
