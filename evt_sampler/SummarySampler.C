#include "SummaryParser.h"
#include "TSystem.h"
#include "TH2.h"
#include "TRandom3.h"
#include "TCanvas.h"

//*****************************************************************
// functions
//*****************************************************************
Bool_t GetDetHists(TString flux_chain_file, Int_t flavor, Int_t is_cc);
void   InitVars(Int_t flavor, Int_t is_cc);
void   ReadSummaryData(Int_t flavor, Int_t is_cc);
Bool_t SampleEvents(TH2D *h_expected, TH2D *h_smeared,
		    vector<Double_t> **store, vector<Double_t> **sample,
		    Int_t low_sample_lim = 10);
void WriteToFile(TString out_name, TH2D *h_sample_nu, TH2D *h_sample_nub,
		 vector<Double_t>** sample_nu, vector<Double_t>** sample_nub);
void   TestPrint();
void   CleanUp();

//*****************************************************************
// globally used variables in this script
//*****************************************************************
TH2D *fhDet_nu;        //!< histogram with expected number of detected nu events
TH2D *fhDet_nub;       //!< histogram with expected number of detected nubar events
TH2D *fhSum_nu;        //!< histogram with all available MC nu events
TH2D *fhSum_nub;       //!< histogram with all available MC nubar events
SummaryParser *fSp;    //!< class to parse summary events
TRandom3      *fRand;  //!< random number generator
Int_t fEbins;          //!< number of energy bins in the event vectors
Int_t fCtbins;         //!< number of costheta bins in the event vectors
//! 2D array of vectors, each vector holds summary nu event numbers of (E, costheta) bin
vector<Double_t> **fSumEvts_nu;
//! 2D array of vectors, each vector holds summary nub event numbers of (E, costheta) bin
vector<Double_t> **fSumEvts_nub;
//! 2D array of vectors to hold a sub-sample of events in fSumEvts_nu
vector<Double_t> **fSampleEvts_nu;
//! 2D array of vectors to hold a sub-sample of events in fSumEvts_nub
vector<Double_t> **fSampleEvts_nub;

map < Int_t, TString > fFlavs  = { {0, "elec" },   //!< map of flavor numbers and strings
				   {1, "muon" },
				   {2, "tau"  } };

//*****************************************************************

/**
 *  This routine creates a sample of simulated MC events that corresponds to a pattern expexted
 *  from running the detector for a certain period of time.
 *
 *  The output of this macro is in output/SummarySample_{flavor}_{NC/CC}_{sample#}.root. If NC, then
 *  the output will contain an event sample of muon-NC + tau-NC + elec-NC.
 *
 * \param  flux_chain_file  Output of FluxChain.C, provides histograms with expected number of events
 *                          in (E,ct) bins after a certain number of years. These histgrams are used
 *                          to create event samples.
 * \param  flavor           Neutrino flavor.
 * \param  is_cc            0 - NC, 1 -CC. If NC, the output will have samples over all flavors.
 * \param  nsamples         Number of event samples (=experiments) to draw.
 *
 */
void SummarySampler(TString flux_chain_file, Int_t flavor, Int_t is_cc, Int_t nsamples = 1) {

  gSystem->Load("$NMHDIR/common_software/libnmhsoft.so");
  
  if ( !GetDetHists(flux_chain_file, flavor, is_cc) ) {
    cout << "ERROR! SummarySampler() problem opening flux histograms." << endl;
    return;
  }

  InitVars(flavor, is_cc);
  ReadSummaryData(flavor, is_cc);

  for (Int_t N = 0; N < nsamples; N++) {
    
    TString out_name = "output/SummarySample_" + fFlavs[flavor] + "-CC_" +
      (TString)to_string(N) + ".root";
    if (is_cc == 0.) out_name = "output/SummarySample_allflavs-NC_" + (TString)to_string(N) + ".root";

    TH2D *smeared_nu  = (TH2D*)fhDet_nu->Clone("sample_nu");
    TH2D *smeared_nub = (TH2D*)fhDet_nu->Clone("sample_nub");
    smeared_nu->Reset();
    smeared_nub->Reset();

    Bool_t SampleOK_nu  = SampleEvents(fhDet_nu,  smeared_nu , fSumEvts_nu , fSampleEvts_nu );
    Bool_t SampleOK_nub = SampleEvents(fhDet_nub, smeared_nub, fSumEvts_nub, fSampleEvts_nub);
    if (SampleOK_nu && SampleOK_nub) {
      WriteToFile(out_name, smeared_nu, smeared_nub, fSampleEvts_nu, fSampleEvts_nub);
    }

    if (smeared_nu) delete smeared_nu;
    if (smeared_nub) delete smeared_nub;

  }
  
  //TestPrint();
  CleanUp();
  
}

//*****************************************************************

/**
 *  Inline function to fetch the histograms with expected event distributions from FluxChain.C output.
 *
 * \param  flux_chain_file  Output of FluxChain.C.
 * \param  flavor           Neutrino flavor.
 * \param  is_cc            0 - NC, 1 -CC.
 * \return                  True if histograms found, False if not found.
 *
 */
Bool_t GetDetHists(TString flux_chain_file, Int_t flavor, Int_t is_cc) {

  TFile f(flux_chain_file, "READ");
  if ( !f.IsOpen() ) {
    cout << "ERROR! GetDetHists() cannot open file " << flux_chain_file << endl;
    return false;
  }

  if (is_cc > 0) {

    TString hname_nu  = "detflux/detflux_" + fFlavs[flavor] + "_cc_nu";
    TString hname_nub = "detflux/detflux_" + fFlavs[flavor] + "_cc_nub";
    TH2D   *h_nu  = (TH2D*)f.Get(hname_nu);
    TH2D   *h_nub = (TH2D*)f.Get(hname_nub);

    if ( !h_nu || !h_nub ) {
      cout << "ERROR! GetDetHists() cannot find hists " << hname_nu << "\t" << hname_nub << endl;
      return false;
    }

    fhDet_nu  = (TH2D*)h_nu->Clone();
    fhDet_nub = (TH2D*)h_nub->Clone();
    fhDet_nu->SetDirectory(0);
    fhDet_nub->SetDirectory(0);

  }
  else {

    //for nc events we only have MC events for elec_NC, hence need to sample them together
    TString hname_elec_nu  = "detflux/detflux_elec_nc_nu";
    TString hname_elec_nub = "detflux/detflux_elec_nc_nub";
    TString hname_muon_nu  = "detflux/detflux_muon_nc_nu";
    TString hname_muon_nub = "detflux/detflux_muon_nc_nub";
    TString hname_tau_nu   = "detflux/detflux_tau_nc_nu";
    TString hname_tau_nub  = "detflux/detflux_tau_nc_nub";

    TH2D   *h_elec_nu  = (TH2D*)f.Get(hname_elec_nu);
    TH2D   *h_elec_nub = (TH2D*)f.Get(hname_elec_nub);
    TH2D   *h_muon_nu  = (TH2D*)f.Get(hname_muon_nu);
    TH2D   *h_muon_nub = (TH2D*)f.Get(hname_muon_nub);
    TH2D   *h_tau_nu   = (TH2D*)f.Get(hname_tau_nu);
    TH2D   *h_tau_nub  = (TH2D*)f.Get(hname_tau_nub);

    if (!h_elec_nu || !h_elec_nub || !h_muon_nu || !h_muon_nub || !h_tau_nu || !h_tau_nub) {
      cout << "ERROR! GetDetHists() cannot find NC hists" << endl;
      return false;
    }

    fhDet_nu  = (TH2D*)h_elec_nu->Clone("detflux_allflav_nc_nu");
    fhDet_nu->Add(h_muon_nu);
    fhDet_nu->Add(h_tau_nu);

    fhDet_nub = (TH2D*)h_elec_nub->Clone("detflux_allflav_nc_nub");
    fhDet_nub->Add(h_muon_nub);
    fhDet_nub->Add(h_tau_nub);

    fhDet_nu->SetTitle("detflux_allflav_nc_nu");
    fhDet_nub->SetTitle("detflux_allflav_nc_nub");

    fhDet_nu->SetDirectory(0);
    fhDet_nub->SetDirectory(0);
    
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

  // determine summary files where events should be read, create hist names
  // we only simulate NC for elec, so there is an exception for that
  
  TString hname_nu  = "sumevts_"  + fFlavs[flavor] + "_cc_nu";
  TString hname_nub = "sumevts_"  + fFlavs[flavor] + "_cc_nub";
  
  if (is_cc == 0) {
    hname_nu  = "sumevts_allflav_nc_nu";
    hname_nub = "sumevts_allflav_nc_nub";
  }

  // histograms to display distributions of all available summary events
  
  fhSum_nu  = (TH2D*)fhDet_nu->Clone(hname_nu);
  fhSum_nub = (TH2D*)fhDet_nub->Clone(hname_nub); 
  fhSum_nu ->Reset();
  fhSum_nub->Reset();
  fhSum_nu ->SetTitle(hname_nu);
  fhSum_nub->SetTitle(hname_nub);

  // allocate vectors, depending on the input hist binning, to store summary event IDs
  // bin [0] is underflow, bin[ X/Yaxis->GetNbins() ] is the last counting bin (hence array
  // length +1), bin[ X/Yaxis->GetNbins()+1 ] is overflow (hence array length + 2)
  fEbins  = fhSum_nu->GetXaxis()->GetNbins() + 2;
  fCtbins = fhSum_nu->GetYaxis()->GetNbins() + 2;
  
  fSumEvts_nu     = new vector<Double_t>* [fEbins];
  fSumEvts_nub    = new vector<Double_t>* [fEbins];
  fSampleEvts_nu  = new vector<Double_t>* [fEbins];
  fSampleEvts_nub = new vector<Double_t>* [fEbins];
  for (Int_t eb = 0; eb < fEbins; eb++) {
    fSumEvts_nu[eb]     = new vector<Double_t> [fCtbins];
    fSumEvts_nub[eb]    = new vector<Double_t> [fCtbins];
    fSampleEvts_nu[eb]  = new vector<Double_t> [fCtbins];
    fSampleEvts_nub[eb] = new vector<Double_t> [fCtbins];
  }

  // initialize the summary parser and random generator
  fSp = new SummaryParser;
  fRand = new TRandom3(0);
}


//*****************************************************************

/**
 *  Inline function to read in summary data from NMH/data/mc_end/.
 *
 *  Given the flavor and is_cc flag, the function will look for filenames in NMH/data/mc_end/ in
 *  specific format (e.g. for muon cc it will look for NMH/data/mc_end/data_atmnu/summary_muon-CC*).
 *  It will read the data from all of the matched files and sort the event IDs to the globally used 
 *  vectors fSumEvts_nu[ebin][ctbin], fSumEvts_nub[ebin][ctbin]. After this routine, e.g.
 *  fSumEvts_nu[ebin][ctbin] will store all of the event IDs (=event entry number it fSp->fChain)
 *  simulated for the specified bin.
 *
 * \param  flavor           Neutrino flavor.
 * \param  is_cc            0 - NC, 1 -CC.
 *
 */
void ReadSummaryData(Int_t flavor, Int_t is_cc) {
  
  TString fnames    = "$NMHDIR/data/mc_end/data_atmnu/summary_" + fFlavs[flavor] + "-CC*.root";
  if (is_cc == 0) fnames = "$NMHDIR/data/mc_end/data_atmnu/summary_elec-NC*.root"; 
  
  // add files to the summary parser
  fSp->fChain->Add(fnames);

  cout << "ReadSummaryData() reading in " << fSp->fChain->GetEntries() << " summary events" << endl;
  
  for (Int_t evt = 0; evt < fSp->fChain->GetEntries(); evt++ ) {

    fSp->fChain->GetEntry(evt);

    Double_t energy =  fSp->MC_energy;
    Double_t ct     = -fSp->MC_dir_z ;
    
    Int_t xbin = fhSum_nu->GetXaxis()->FindBin( energy );
    Int_t ybin = fhSum_nu->GetYaxis()->FindBin( ct );
    
    if (fSp->MC_type > 0)  {
      fhSum_nu->Fill( energy, ct );
      fSumEvts_nu[xbin][ybin].push_back( evt );
    }
    else {
      fhSum_nub->Fill( energy, ct );
      fSumEvts_nub[xbin][ybin].push_back( evt );
    }
  }

  cout << "ReadSummaryData() finished reading in events" << endl;
}

//*****************************************************************

/**
 *  Inline function to create an event sample from the summary events.
 *
 * \param  h_expected     Pointer to a histogram with expected event distribution in (E,ct).
 * \param  h_smeared      Pointer to a histogram where the Poisson-smeared event distribution is saved.
 * \param  store          Pointer to a 2D array of vectors, each store[ebin][ctbin] vector stores the
 *                        MC event IDs available in this (E,ct) bin.
 * \param  sample         Pointer to a 2D array of vectors, the MC event IDs of sampled events will be
 *                        stored into each vector[ebin][ctbin]
 * \param low_sample_lim  In some bins there are very few events expected and available, which may
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
	
	if (smeared <= low_sample_lim) {
	  cout << "WARNING! SampleEvents() for E, ct " << en << ", " << ct <<  " wanted " << smeared
	       << " events, but only " << store[xbin][ybin].size()
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
      Int_t limit   = 1e4;
      
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
 *  Inline function to write the created event samples to a root file.
 *
 * \param  out_name     Name of the file where events are written.
 * \param  h_sample_nu  2D histogram that displays the E,ct distribution of sampled nu events.
 * \param  h_sample_nub 2D histogram that displays the E,ct distribution of sampled nub events.
 * \param  sample_nu    Pointer to a 2D array of vectors, sample_nu[ebin][ctbin] holds the nu MC event
 *                      IDs that were sampled for this and are to be written out.
 * \param  sample_nub   Pointer to a 2D array of vectors, sample_nub[ebin][ctbin] holds the nub 
 *                      MC event IDs that were sampled for this and are to be written out.
 *
 */
void WriteToFile(TString out_name, TH2D *h_sample_nu, TH2D *h_sample_nub,
		 vector<Double_t>** sample_nu, vector<Double_t>** sample_nub) {

  //--------------------------------------------------------------
  // The way fChain accesses entries, it is extremely faster to sort event numbers
  // before looping and writing
  //--------------------------------------------------------------
  
  vector<Double_t> all_evts;

  for (Int_t ebin = 0; ebin < fEbins; ebin++) {
    for (Int_t ctbin = 0; ctbin < fCtbins; ctbin++) {
      
      for ( auto evt: sample_nu[ebin][ctbin] ) {
	all_evts.push_back( evt );
      }

      for ( auto evt: sample_nub[ebin][ctbin] ) {
	all_evts.push_back( evt );
      }

    }
  }

  std::sort( all_evts.begin(), all_evts.end() );

  //--------------------------------------------------------------
  // After sorting write to file
  //--------------------------------------------------------------
  
  TFile fout(out_name, "RECREATE");
  TTree *tout = fSp->fChain->CloneTree(0);

  Int_t counter = 0;
  for (auto evt: all_evts) {
    fSp->fChain->GetEntry(evt);
    tout->Fill();

    if (counter % 25000 == 0) cout << "WriteToFile() written " << counter << " events" << endl;
    counter++;
  }
  
  h_sample_nu->Write();
  h_sample_nub->Write();
  tout->Write();
  fout.Close();
  
}

//*****************************************************************

/**
 *  Inline function to test that summary data was read to global vectors correctly.
 */
void TestPrint() {

  vector<TH2D*>                hists = {fhSum_nu, fhSum_nub};
  vector< vector<Double_t>** > vecs  = {fSumEvts_nu, fSumEvts_nub};

  for (Int_t i = 0; i < (Int_t)hists.size(); i++) {
  
    TH2D              *h = hists[i];
    vector<Double_t> **v = vecs[i];
    
    for (Int_t xbin = 0; xbin <= h->GetXaxis()->GetNbins(); xbin++) {

      for (Int_t ybin = 0; ybin <= h->GetYaxis()->GetNbins(); ybin++) {
      
	if ( h->GetBinContent(xbin, ybin) != v[xbin][ybin].size() ) {

	  cout << "TestPrint() ERROR! in xbin, ybin " << xbin << "\t" << ybin << endl;
	  cout << "Histogram " << h->GetName() << " has " << h->GetBinContent(xbin, ybin)
	       << " events in this bin, the vector has " << v[xbin][ybin].size() << endl;
	
	}
      }
    }
  }
  
}

//*****************************************************************

/**
 *  Inline function to clear dynamically allocated memory.
 */
void CleanUp() {

  //clean up the dynamically allocated vectors/arrays
  
  if (fSumEvts_nu && fSumEvts_nub && fhSum_nu) {

    for (Int_t xbin = 0; xbin < fEbins; xbin++) {
      if ( fSumEvts_nu[xbin]     ) delete[] fSumEvts_nu[xbin];
      if ( fSumEvts_nub[xbin]    ) delete[] fSumEvts_nub[xbin];
      if ( fSampleEvts_nu[xbin]  ) delete[] fSampleEvts_nu[xbin];
      if ( fSampleEvts_nub[xbin] ) delete[] fSampleEvts_nub[xbin];

    }

    if ( fSumEvts_nu     ) delete[] fSumEvts_nu;
    if ( fSumEvts_nub    ) delete[] fSumEvts_nub;
    if ( fSampleEvts_nu  ) delete[] fSampleEvts_nu;
    if ( fSampleEvts_nub ) delete[] fSampleEvts_nub;
    
  }

  //clean up hists/classes/etc
  
  if (fhDet_nu)  delete fhDet_nu; 
  if (fhDet_nub) delete fhDet_nub;
  if (fhSum_nu)  delete fhSum_nu; 
  if (fhSum_nub) delete fhSum_nub;

  if (fSp)   delete fSp;
  if (fRand) delete fRand;
  
}
