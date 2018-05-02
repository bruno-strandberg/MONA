#include "SummaryParser.h"
#include "GSGParser.h"
#include "TSystem.h"
#include "TH2.h"
#include "TRandom3.h"
#include "NMHUtils.h"

//*****************************************************************
// functions
//*****************************************************************
Bool_t GetIntHists(TString flux_chain_file, Int_t flavor, Int_t is_cc);
void   InitVars(Int_t flavor, Int_t is_cc);
void   CleanUp();
void   ReadGSGData(TString gsg_file_list, Int_t flavor, Int_t is_cc);
Bool_t SampleEvents(TH2D *h_expected, TH2D *h_smeared,
		    vector<Double_t> **store, vector<Double_t> **sample,
		    Int_t low_sample_lim);

//*****************************************************************
// structure to store run nr and event nr in a vector and sort
//*****************************************************************

/**
 *  Struct to store the Monte Carlo run number and event number pair.
 *
 */
struct evtid {
  Int_t run_nr;
  Int_t evt_nr;

  //! Default constructor
  evtid(): run_nr(0), evt_nr(0) {}

  /**
   * Constructor
   * \param _run_nr Monte-Carlo file number
   * \param _evt_nr Monte-Carlo event number in file
   */
  evtid(Int_t _run_nr, Int_t _evt_nr) {
    run_nr = _run_nr;
    evt_nr = _evt_nr;
  }

  /**
   * Operator to sort events by increasing file number and event id.
   * \param i first evtid structure
   * \param j second evtid structure
   * \return  true if second run nr is larger; if equal run numbers true if second event nr is larger.
   */
  bool operator < (const evtid& rhs) {
    
    if      ( this->run_nr < rhs.run_nr )  { return true;                  }
    else if ( this->run_nr == rhs.run_nr ) { return (this->evt_nr < rhs.evt_nr); }
    else                                   { return false;                 }
  }


};

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
Double_t fVcan = 0.;   //!< can size, set in ReadGSGData()
Double_t fRhoSW = 0.;  //!< sea water density, set in ReadGSGData()
//! 2D array of vectors, each vector holds summary nu event numbers of (E, costheta) bin
vector<evtid> **fGSGEvts_nu;
//! 2D array of vectors, each vector holds summary nub event numbers of (E, costheta) bin
vector<evtid> **fGSGEvts_nub;
//! 2D array of vectors to hold a sub-sample of events in fGSGEvts_nu
vector<evtid> **fSampleEvts_nu;
//! 2D array of vectors to hold a sub-sample of events in fGSGEvts_nub
vector<evtid> **fSampleEvts_nub;

//! map of flavor numbers and strings
map < Int_t, TString > fFlavs  = { {0, "elec" },
				   {1, "muon" },
				   {2, "tau"  } };

//! map of interaction numbers and strings
map < Int_t, TString > fInts  = { {0, "NC" },
				  {1, "CC" } };

//*****************************************************************
// main function
//*****************************************************************

void GSGSampler(TString flux_chain_file, TString gsg_file_list, Int_t flavor, Int_t is_cc, Int_t nsamples = 1) {

  gSystem->Load("$NMHDIR/common_software/libnmhsoft.so");

  if ( !GetIntHists(flux_chain_file, flavor, is_cc) ) {
    cout << "ERROR! GSGSampler() problem opening flux histograms." << endl;
    return;
  }

  InitVars(flavor, is_cc);
  ReadGSGData(gsg_file_list, flavor, is_cc);

  // convert from unit 1/MTon to unitless by 1e-6 [MTon/Ton] * V [m3] * rho [Ton/m3]
  // after this, each bin will have expected nr of interactions per E,ct bin in
  // certain operation time, as defined when FluxChain is run.
  
  fhInt_nu->Scale(1e-6 * fVcan * fRhoSw);
  fhInt_nub->Scale(1e-6 * fVcan * fRhoSw);
  
  for (Int_t N = 0; N < nsamples; N++) {

    TString out_name = "output/EvtSample_" + fFlavs[flavor] + "-CC_" +
      (TString)to_string(N) + ".root";
    if (is_cc == 0.) out_name = "output/EvtSample_allflavs-NC_" + (TString)to_string(N) + ".root";

    TH2D *smeared_nu  = (TH2D*)fhInt_nu->Clone("sample_nu");
    TH2D *smeared_nub = (TH2D*)fhInt_nu->Clone("sample_nub");
    smeared_nu->Reset();
    smeared_nub->Reset();

    Bool_t SampleOK_nu  = SampleEvents(fhInt_nu,  smeared_nu , fGSGEvts_nu , fSampleEvts_nu );
    Bool_t SampleOK_nub = SampleEvents(fhInt_nub, smeared_nub, fGSGEvts_nub, fSampleEvts_nub);

    
  }
  
  //std::sort( fGSGEvts_nu[10][5].begin(), fGSGEvts_nu[10][5].end() );


  CleanUp();
  
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
  
  fGSGEvts_nu     = new vector<evtid>* [fEbins];
  fGSGEvts_nub    = new vector<evtid>* [fEbins];
  fSampleEvts_nu  = new vector<evtid>* [fEbins];
  fSampleEvts_nub = new vector<evtid>* [fEbins];
  for (Int_t eb = 0; eb < fEbins; eb++) {
    fGSGEvts_nu[eb]     = new vector<evtid> [fCtbins];
    fGSGEvts_nub[eb]    = new vector<evtid> [fCtbins];
    fSampleEvts_nu[eb]  = new vector<evtid> [fCtbins];
    fSampleEvts_nub[eb] = new vector<evtid> [fCtbins];
  }

  // initialize the summary parser and random generator
  fSp = new SummaryParser;
  fRand = new TRandom3(0);
}

//*****************************************************************

void ReadGSGData(TString gsg_file_list, Int_t flavor, Int_t is_cc) {

  vector<TString> fnames = NMHUtils::ReadLines(gsg_file_list);

  // loop over file names

  for (auto fname: fnames) {

    // check that the files correspond to the expected flavor and interaction

    if ( !fname.Contains( fFlavs[flavor] ) || !fname.Contains( fInts[is_cc] ) ) {
      cout << "WARNING! ReadGSGData() specified flavor " << fFlavs[flavor] << " and interaction " << fInts[is_cc]
	   << " seem to mismatch flavor and interaction in filename, skipping file " << fname << endl;
      continue;
    }

    cout << "NOTICE ReadGSGData() Reading from file: " << fname << endl;

    // initalise parser, set detector can size and sea water density

    GSGParser gp(fname);

    if (fVcan == 0.) fVcan = gp.fVcan;
    else if ( fVcan != gp.fVcan) {
      cout << "ERROR! ReadGSGData() detector can change from  " << fVcan << " to " << gp.fVcan
    	   << ", exiting." << endl;
      return;
    }

    if (fRhoSW == 0.) fRhoSW = gp.fRho_seawater;
    else if ( fRhoSW != gp.fRho_seawater) {
      cout << "ERROR! ReadGSGData() sea water density change from  " << fRhoSW << " to " << gp.fRho_seawater
    	   << ", exiting." << endl;
      return;
    }

    //loop over events

    while ( gp.NextEvent() ) {

      // skip events with vertices outside the can
      if ( !gp.VertexInCan() ) continue;

      Double_t energy =  gp.Neutrino_E;
      Double_t ct     = -gp.Neutrino_D3;
    
      Int_t xbin = fhGSG_nu->GetXaxis()->FindBin( energy );
      Int_t ybin = fhGSG_nu->GetYaxis()->FindBin( ct );

      if (gp.Neutrino_PdgCode > 0)  {
      	fhGSG_nu->Fill( energy, ct );
      	fGSGEvts_nu[xbin][ybin].push_back( evtid(gp.fRunNr, gp.iEvt) );
      }
      else {
      	fhGSG_nub->Fill( energy, ct );
      	fGSGEvts_nub[xbin][ybin].push_back( evtid(gp.fRunNr, gp.iEvt) );
      }

    } // end loop over events
    
  } //end loop over files
  
  cout << "NOTICE ReadGSGData() finished reading GSG data" << endl;

}
//*****************************************************************

/**
 *  Inline function to create an event sample from the GSG events.
 *
 * \param  h_expected     Pointer to a histogram with expected event distribution in (E,ct).
 * \param  h_smeared      Pointer to a histogram where the Poisson-smeared event distribution is saved.
 * \param  store          Pointer to a 2D array of vectors, each store[ebin][ctbin] vector stores the
 *                        MC event IDs available in this (E,ct) bin.
 * \param  sample         Pointer to a 2D array of vectors, the MC event IDs of sampled events will be
 *                        stored into each vector[ebin][ctbin]
 * \param low_sample_lim  In some bins there are very few events expected and available, which may lead
 *                        to situations where e.g. 6 events should be sampled, but only 5 are in store.
 *                        If sampled < low_sample_lim, this problem is ignored and all events in store
 *                        will be used for this bin.
 *
 */
Bool_t SampleEvents(TH2D *h_expected, TH2D *h_smeared,
		    vector<Double_t> **store, vector<Double_t> **sample,
		    Int_t low_sample_lim) {

  // make sure sample is empty

  for (Int_t ebin = 0; ebin < fEbins; ebin++) {
    for (Int_t ctbin = 0; ctbin < fCtbins; ctbin++) {
      sample[ebin][ctbin].clear();
    }
  }

  // loop over bins
  
  for (Int_t xbin = 1; xbin <= h_expected->GetXaxis()->GetNbins(); xbin++) {
    for (Int_t ybin = 1; ybin <= h_expected->GetYaxis()->GetNbins(); ybin++) {

      // get the expected number of events in this bin, perform poisson smearing
      
      Double_t expected = h_expected->GetBinContent(xbin, ybin);
      Double_t smeared  = fRand->Poisson(expected);
      
      // check that enough events are available in store
      
      if ( smeared > store[xbin][ybin].size() ) {

	Double_t en = h_expected->GetXaxis()->GetBinCenter(xbin);
	Double_t ct = h_expected->GetYaxis()->GetBinCenter(ybin);
	
	if (expected <= low_sample_lim) {
	  cout << "WARNING! SampleEvents() for E, ct " << en << ", " << ct <<  " wanted " << smeared
	       << " events (expectation " << expected << "), but only " << store[xbin][ybin].size()
	       << " available, using all available events." << endl;
	  smeared = store[xbin][ybin].size();
	}
	else {
	  cout << "ERROR! SampleEvents() for E, ct " << en << ", " << ct <<  " wanted " << smeared
	       << " events (expectation " << expected << "), but only " << store[xbin][ybin].size()
	       << " available, exiting." << endl;
	  return false;
	}
      }

      h_smeared->SetBinContent(xbin, ybin, smeared);
      
      // randomly draw smeared number of events from store, insert them to sample

      Int_t counter = 0;
      Int_t limit   = 1e6;
      
      while ( sample[xbin][ybin].size() != smeared ) {

	// random event index in store
	Int_t idx = fRand->Integer( store[xbin][ybin].size() );

	// check that this event has not already been drawn
	if ( std::find( sample[xbin][ybin].begin(), sample[xbin][ybin].end(), store[xbin][ybin][idx] )
	     == sample[xbin][ybin].end() ) {

	  // insert the event id to sample
	  sample[xbin][ybin].push_back( store[xbin][ybin][idx] );
	  
	} // end if

	if (counter++ > limit) {
	  cout << "ERROR! SampleEvents() while loop exceeds sampling limit of " << limit
	       << " per xbin, ybin " << xbin << "\t" << ybin << endl;
	  cout << "Smeared, store size, sample size: " << smeared << "\t" << store[xbin][ybin].size()
	       << "\t" << sample[xbin][ybin].size() << endl;
	  return false;
	}
	
      } //end while
      
    } //end loop over ybins
  } //end loop over xbins

  return true;
}

//*****************************************************************

/**
 *  Inline function to clear dynamically allocated memory.
 */
void CleanUp() {

  //clean up the dynamically allocated vectors/arrays
  
  if (fGSGEvts_nu && fGSGEvts_nub) {

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
