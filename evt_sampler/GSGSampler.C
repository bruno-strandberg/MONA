//header
#include "GSGSampler.h"

//root
#include "TSystem.h"
#include "TStopwatch.h"

//NMH
#include "SummaryParser.h"
#include "GSGParser.h"
#include "NMHUtils.h"

//generic cpp

using namespace GSGS;

/**
 *  Function to create a sample of Monte-Carlo events.
 *
 *  This function creates a sample of Monte-Carlo events for a given flavor and interaction type.
 *  Each FluxChain.C output file contains histograms in TDirectory intflux/, each histogram holds
 *  the number of interacted neutrino events per specific flavor and interaction type in a certain
 *  operation time per MTon. The histograms are scaled by the can size and water density, after that
 *  each (E, ct) bin will hold the number of interacted events inside the detector can. These numbers
 *  are in turn Poisson-smeared to simulate statistical fluctuations.
 *
 *  In ORCA Monte-Carlo chain, interacting events are simulated at gSeaGen level. The gSeaGen data
 *  from the files in gsg_flist are read into vectors, each (E, ct) bin will have an associated
 *  vector that will hold the Monte-Carlo event IDs (see struct evtid) available in that bin. Then,
 *  N number of random event ID's are sampled from the vector, where N is the Poisson-smeared
 *  number of expected interacted events from the intflux histogram.
 *
 *  After sampling, a loop over all summary files for this flavor and interaction type are performed.
 *  The purpose of the loop is to check which of the sampled interacted (=gSeaGen) events made it 
 *  to the end of the Monte Carlo chain (i.e. into the summary file). The events that made it to the
 *  end are the events that will be stored in the output.
 * 
 *  The code will output a file output/EvtSample_{flavor}-{interaction}_{Nf}_{Ns}.root, where
 *  Nf is the sequence number of the flux chain file and Ns is the sample number.
 *
 * \param flux_chain_list a text file with a list of FluxChain.C macro outputs
 * \param gsg_flist       a text file with a list of gSeaGen files
 * \param flavor          neutrino flavor
 * \param is_cc           0 - neutral current, 1 - charged current
 * \param nsamples        Number of Monte-Carlo samples to be generated.
 * \param memory_lim      If more RAM than memory_lim used to store data, trigger write-out
 * \param exp_lim         If more experiments than exp_lim in memory, trigger write-out
 * \param dbg_scale       Option to scale down the number of interacted events, 
 *                        such that the program can be run with a reduced nr of gSeaGen files
 *                        without getting errors related to not having enough Monte-Carlo events.
 */
void GSGSampler(TString flux_chain_flist, 
		TString gsg_flist, 
		Int_t flavor, Int_t is_cc, Int_t nsamples = 1, 
		Double_t memory_lim = 2,
		Double_t exp_lim    = 10,
		Double_t dbg_scale  = 1.) {

  // timers
  TStopwatch OverallTimer, DataReadTimer, DataWriteTimer, SamplerTimer;
  DataReadTimer.Stop();
  DataWriteTimer.Stop();
  SamplerTimer.Stop();
  DataReadTimer.Reset();
  DataWriteTimer.Reset();
  SamplerTimer.Reset();

  cout << "NOTICE GSGSampler() running for flavor, nc/cc: " << flavor << "\t" << is_cc << endl;

  gSystem->Load("$NMHDIR/common_software/libnmhsoft.so");

  Bool_t initialized = false; //variable to control that init and gSeaGen reading is only done once
  Double_t gsg_data_size = 0.;

  vector<TString> flux_files = NMHUtils::ReadLines(flux_chain_flist);

  //------------------------------------------------------
  // loop over flux files
  //------------------------------------------------------

  for (Int_t F = 0; F < (Int_t)flux_files.size(); F++) {

    TString fc_file = flux_files[F];

    if ( !GetIntHists(fc_file, flavor, is_cc) ) {
      cout << "ERROR! GSGSampler() problem opening flux histograms from " << fc_file << endl;
      return;
    }

    if ( !initialized ) {
      initialized = true;
      InitVars(flavor, is_cc);

      DataReadTimer.Start(kFALSE);
      gsg_data_size = ReadGSGData(gsg_flist, flavor, is_cc);
      if ( gsg_data_size == 0. ) {
	cout << "ERROR! GSGSampler() problem reading gSeaGen files." << endl;
	return;
      }
      DataReadTimer.Stop();
    }

    // convert from unit 1/MTon to unitless by 1e-6 [MTon/Ton] * V [m3] * rho [Ton/m3]
    // after this, each bin will have expected nr of interactions per E,ct bin in
    // certain operation time, as defined when FluxChain is run.
  
    fhInt_nu ->Scale(1e-6 * fVcan * fRhoSW * dbg_scale);
    fhInt_nub->Scale(1e-6 * fVcan * fRhoSW * dbg_scale);

    //------------------------------------------------------
    // loop over samples
    //------------------------------------------------------
    
    SamplerTimer.Start(kFALSE);
    for (Int_t N = 0; N < nsamples; N++) {

      cout << "NOTICE GSGSampler() flux file " << fc_file << " sampling experiment " << N << endl;

      TH2D *smeared_nu  = (TH2D*)fhInt_nu->Clone("sample_nu");
      TH2D *smeared_nub = (TH2D*)fhInt_nu->Clone("sample_nub");
      smeared_nu->Reset();
      smeared_nub->Reset();

      Bool_t SampleOK_nu  = SampleEvents(fhInt_nu,  smeared_nu , fGSGEvts_nu , fSampleEvts_nu , 10);
      Bool_t SampleOK_nub = SampleEvents(fhInt_nub, smeared_nub, fGSGEvts_nub, fSampleEvts_nub, 10);

      StoreForWriting( (SampleOK_nu && SampleOK_nub), smeared_nu, smeared_nub, F, N);

      if (smeared_nu) delete smeared_nu;
      if (smeared_nub) delete smeared_nub;
    
    } //end loop over samples
    SamplerTimer.Stop();

    // calculate the RAM used by gsg_data and sample_data; write out data if
    // too much ram in use or too many experiments in memory or last flux file has been sampled

    Double_t sample_data_size = 0.;
    for (auto E: fExps) { sample_data_size += E.size() * sizeof(evtid); }
    cout << "NOTICE GSGSampler() RAM under data: " << sample_data_size/1e9 + gsg_data_size << endl;

    if ( ( (sample_data_size/1e9 + gsg_data_size) > memory_lim ) || 
	 ( (Double_t)fExps.size() > exp_lim ) || ( F == Int_t( flux_files.size()-1 ) ) ) {
      DataWriteTimer.Start(kFALSE);
      WriteToFiles(flavor, is_cc);
      DataWriteTimer.Stop();
    }

  } // end loop over flux files

  CleanUp();
  
  cout << "NOTICE GSGSampler() data reading time: "  << (Double_t)DataReadTimer.RealTime()/3600. << " hours" << endl;
  cout << "NOTICE GSGSampler() total sampling time: "  << (Double_t)SamplerTimer.RealTime()/3600. << " hours" << endl;
  cout << "NOTICE GSGSampler() data writing time: " << (Double_t)DataWriteTimer.RealTime()/3600. << " hours" << endl;
  cout << "NOTICE GSGSampler() Overall time consumed: " << (Double_t)OverallTimer.RealTime()/3600. << " hours" << endl;
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
Bool_t GSGS::GetIntHists(TString flux_chain_file, Int_t flavor, Int_t is_cc) {

  if (fhInt_nu)  delete fhInt_nu;
  if (fhInt_nub) delete fhInt_nub;

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
void GSGS::InitVars(Int_t flavor, Int_t is_cc) {

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

  // initialize random generator
  fRand = new TRandom3(0);
}

//*****************************************************************

/**
 * Inline function to read gSeaGen data into globally used variables.
 *
 * The data are read into vectors fGSGEvts_nu and fGSGEvts_nub. Additionally, the can size and
 * sea water density fVcan and fRho_seawater are set.
 *
 * \param gsg_file_list  list of gSeaGen files
 * \param flavor         nu flavor
 * \param is_cc          0 - NC, 1 -CC
 * \return               Size of data in RAM if reading successful, 0. if unsuccessful
 *
*/
Double_t GSGS::ReadGSGData(TString gsg_file_list, Int_t flavor, Int_t is_cc) {

  cout << "NOTICE ReadGSGData() started reading GSG data" << endl;

  Double_t gsg_data_size = 0.;
  vector<TString> fnames = NMHUtils::ReadLines(gsg_file_list);

  // create a hash from the input gsg filenames; if a cache for this input has
  // been created, read data from cache instead of from GSGfiles

  std::string hash_str = "";
  for (auto f: fnames) hash_str += (std::string)f;

  TString hash       = (TString)to_string( std::hash<std::string>{}(hash_str) );
  TString cache_name = "cache/" + hash + "_" + fFlavs[flavor] + "_" + fInts[is_cc] + ".root";

  if ( NMHUtils::FileExists(cache_name) ) {
    gsg_data_size = ReadFromCache(cache_name);
    return gsg_data_size;
  }

  // no cache --> loop over file names and read data from GSG files

  for (auto fname: fnames) {

    // check that the files correspond to the expected flavor and interaction

    if ( !fname.Contains( fFlavs[flavor] ) || !fname.Contains( fInts[is_cc] ) ) {
      cout << "WARNING! ReadGSGData() specified flavor " << fFlavs[flavor] << " and interaction " 
	   << fInts[is_cc] << " seem to mismatch flavor and interaction in filename, skipping file "
	   << fname << endl;
      continue;
    }

    // initalise parser, set detector can size and sea water density

    cout << "NOTICE ReadGSGData() Reading from file: " << fname << endl;

    GSGParser gp(fname);

    if (fVcan == 0.) fVcan = gp.fVcan;
    else if ( fVcan != gp.fVcan) {
      cout << "ERROR! ReadGSGData() detector can change from  " << fVcan << " to " << gp.fVcan
    	   << ", exiting." << endl;
      return 0.;
    }

    if (fRhoSW == 0.) fRhoSW = gp.fRho_seawater;
    else if ( fRhoSW != gp.fRho_seawater) {
      cout << "ERROR! ReadGSGData() sea water density change from  " << fRhoSW
	   << " to " << gp.fRho_seawater << ", exiting." << endl;
      return 0.;
    }

    // check that the corresponding summary file exists

    TString summary_file = GetSummaryName(flavor, is_cc, gp.fE_min, gp.fE_max, gp.fRunNr);
    if ( !NMHUtils::FileExists( summary_file ) ) {
      cout << "WARNING! ReadGSGData() cannot find summary file " << summary_file
	   << ", skipping gSeaGen file " << fname << endl;
      continue;
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
      	fGSGEvts_nu[xbin][ybin].push_back( evtid(gp.fRunNr, gp.iEvt, gp.fE_min) );
      }
      else {
      	fhGSG_nub->Fill( energy, ct );
      	fGSGEvts_nub[xbin][ybin].push_back( evtid(gp.fRunNr, gp.iEvt, gp.fE_min) );
      }

      gsg_data_size += sizeof(evtid);

    } // end loop over events
    
  } //end loop over files
  
  cout << "NOTICE ReadGSGData() finished reading GSG data, " 
       << gsg_data_size/1e9 << " gb of RAM used." << endl;

  CacheGSGdata(cache_name);

  return gsg_data_size/1e9;

}

//*****************************************************************

/**
 * Function to cache GSG data to vectors.
 *\param fname  Name of the root file where data is written.
 */
void GSGS::CacheGSGdata(TString fname) {

  if ( NMHUtils::FileExists(fname) ) {
    cout << "NOTICE CacheGSGdata() file " << fname << " exists, not overwriting." << endl;
    return;
  }

  Double_t run_nr, evt_nr, e_min, is_nub, en_bin, ct_bin;

  TFile *fout = new TFile(fname, "RECREATE");
  TTree *tout = new TTree("gsgcache", "Tree with evtid");
  tout->Branch("run_nr" ,  &run_nr, "run_nr/D");
  tout->Branch("evt_nr" ,  &evt_nr, "evt_nr/D");
  tout->Branch("e_min"  ,   &e_min,  "e_min/D");
  tout->Branch("is_nub" ,  &is_nub, "is_nub/D");
  tout->Branch("en_bin" ,  &en_bin, "en_bin/D");
  tout->Branch("ct_bin" ,  &ct_bin, "ct_bin/D");

  for (Int_t ebin = 0; ebin < fEbins; ebin++) {
    for (Int_t ctbin = 0; ctbin < fCtbins; ctbin++) {

      en_bin = ebin;
      ct_bin = ctbin;

      is_nub = 0.;
      for (auto evt: fGSGEvts_nu[ebin][ctbin]) {
	run_nr = evt.run_nr;
	evt_nr = evt.evt_nr;
	e_min  = evt.e_min;
	tout->Fill();
      }

      is_nub = 1.;
      for (auto evt: fGSGEvts_nub[ebin][ctbin]) {
	run_nr = evt.run_nr;
	evt_nr = evt.evt_nr;
	e_min  = evt.e_min;
	tout->Fill();
      }

    }
  }
  
  TH1D h("can_and_rho","can_and_rho", 5, 0, 5);
  h.SetBinContent(1, fVcan);
  h.GetXaxis()->SetBinLabel(1, "Vcan");
  h.SetBinContent(2, fRhoSW);
  h.GetXaxis()->SetBinLabel(2, "RhoSW");

  tout->Write();
  h.Write();
  fout->Close();
  delete fout;

}

//*****************************************************************

/**
 * Function to read cached gSeaGen data to program.
 *
 * \param fname  Name of the file where the data is stored
 * \return       Size of the data in RAM if successful, 0. if read-in failed
 */
Double_t GSGS::ReadFromCache(TString fname) {

  Double_t gsg_data_size = 0.;

  //-------------------------------------------------------------
  // open file and get the tree and the histogram
  //-------------------------------------------------------------
  TFile *f = new TFile(fname, "READ");

  if ( !f->IsOpen() ) {
    cout << "ERROR! ReadFromCache() could open file " << fname  << endl;
    return 0.;
  }

  TH1D  *h = (TH1D*)f->Get("can_and_rho");
  TTree *t = (TTree*)f->Get("gsgcache");

  if (t == NULL || h == NULL) {
    cout << "ERROR! ReadFromCache() could not find data in file " << fname  << endl;
    return 0.;
  }

  fVcan  = h->GetBinContent(1);
  fRhoSW = h->GetBinContent(2);

  Double_t run_nr, evt_nr, e_min, is_nub, en_bin, ct_bin;

  t->SetBranchAddress("run_nr" ,  &run_nr);
  t->SetBranchAddress("evt_nr" ,  &evt_nr);
  t->SetBranchAddress("e_min"  ,   &e_min);
  t->SetBranchAddress("is_nub" ,  &is_nub);
  t->SetBranchAddress("en_bin" ,  &en_bin);
  t->SetBranchAddress("ct_bin" ,  &ct_bin);

  //-------------------------------------------------------------
  // loop over the tree and store the data to vectors
  //-------------------------------------------------------------

  for (Int_t i = 0; i < t->GetEntries(); i++) {

    t->GetEntry(i);

    if (is_nub < 0.5)  {
      fhGSG_nu->SetBinContent(en_bin, ct_bin, fhGSG_nu->GetBinContent(en_bin, ct_bin) + 1 );
      fGSGEvts_nu[(Int_t)en_bin][(Int_t)ct_bin].push_back( evtid(run_nr, evt_nr, e_min) );
    }
    else {
      fhGSG_nub->SetBinContent(en_bin, ct_bin, fhGSG_nub->GetBinContent(en_bin, ct_bin) + 1 );
      fGSGEvts_nub[(Int_t)en_bin][(Int_t)ct_bin].push_back( evtid(run_nr, evt_nr, e_min) );
    }

    gsg_data_size += sizeof(evtid);

  }

  f->Close();
  delete f;

  cout << "NOTICE ReadFromCache() finished reading GSG data, " 
       << gsg_data_size/1e9 << " gb of RAM used." << endl;

  return gsg_data_size/1e9;

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
 *                        stored into each vector sample[ebin][ctbin]
 * \param low_sample_lim  In some bins there are very few events expected and available, which may lead
 *                        to situations where e.g. 6 events should be sampled, but only 5 are in store.
 *                        If sampled < low_sample_lim, this problem is ignored and all events in store
 *                        will be used for this bin.
 *
 */
Bool_t GSGS::SampleEvents(TH2D *h_expected, TH2D *h_smeared,
				vector<evtid> **store, vector<evtid> **sample,
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
      Int_t limit   = 1e8;
      
      while ( sample[xbin][ybin].size() != (UInt_t)smeared ) {

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
 * Inline function to store data for later writing to output.
 *
 * \param SampleOK     Boolean to indicate whether the sample in fSampleEvts_nu(b) is OK.
 * \param smeared_nu   Pointer to 2D histogram with a Poisson-smeared interaction counts in E,ct bins
 * \param smeared_nub  Pointer to 2D histogram with a Poisson-smeared interaction counts in E,ct bins
 * \param F            Index of the flux file in the flux file list
 * \param N            Index of the sample
 *
 */
void GSGS::StoreForWriting(Bool_t SampleOK, TH2D *smeared_nu, TH2D *smeared_nub, Int_t F, Int_t N) {
  
  // create a new vector for this experiment; vector for hists of this experiment; name of experiment
  fExps.push_back( vector<evtid>() );
  fExpHists.push_back( vector<TH2D*>() );
  fExpNames.push_back( "flux" + (TString)to_string(F) + "_sample" + (TString)to_string(N) );

  // store the histograms
  fExpHists.back().push_back( (TH2D*)fhInt_nu   ->Clone() );
  fExpHists.back().push_back( (TH2D*)fhInt_nub  ->Clone() );
  fExpHists.back().push_back( (TH2D*)fhGSG_nu   ->Clone() );
  fExpHists.back().push_back( (TH2D*)fhGSG_nub  ->Clone() );
  fExpHists.back().push_back( (TH2D*)smeared_nu ->Clone() );
  fExpHists.back().push_back( (TH2D*)smeared_nub->Clone() );

  if (!SampleOK) return; //if bad sample leave the vector empty, deal with this when writing to file

  // push all sampled events into the vector
  for (Int_t ebin = 0; ebin < fEbins; ebin++) {
    for (Int_t ctbin = 0; ctbin < fCtbins; ctbin++) {

      for ( auto evt: fSampleEvts_nu[ebin][ctbin]  ) { fExps.back().push_back(evt); }
      for ( auto evt: fSampleEvts_nub[ebin][ctbin] ) { fExps.back().push_back(evt); }
      
    }
  }
 
  // sort the sampled events by run number, energy range and event ID
  sort( fExps.back().begin(), fExps.back().end() );

}

//*****************************************************************

/**
 * Inline function to write sampled data to files.
 *
 * \param flavor  nu flavor
 * \param is_cc   0 - nc, 1 - cc
 *
 */
void GSGS::WriteToFiles(Int_t flavor, Int_t is_cc) {

  cout << "NOTICE WriteToFiles() writing out sampled data" << endl;

  // select relevant summary files and add to the summary parser

  TString fnames = "$NMHDIR/data/mc_end/data_atmnu/summary_" + fFlavs[flavor] + "-CC*.root";
  if (is_cc == 0) fnames = "$NMHDIR/data/mc_end/data_atmnu/summary_elec-NC*.root"; 

  SummaryParser sp;
  sp.fChain->Add(fnames);

  // create output files and trees; create search limits that indicate which
  // range of the fExps[i] vector should be searched for the given run_nr and e_min

  vector< TFile* > files;
  vector< TTree* > trees;
  vector< std::pair<Int_t, Int_t> > search_lims;

  for (Int_t N = 0; N < (Int_t)fExps.size(); N++) {

    TString out_name = "output/EvtSample_" + fFlavs[flavor] + "-CC_" + fExpNames[N] + ".root";
    if (is_cc == 0) out_name = "output/EvtSample_allflavs-NC_" + fExpNames[N] + ".root";

    files.push_back( new TFile(out_name, "RECREATE") );            //create file
    trees.push_back( sp.fChain->CloneTree(0) );                    //add empty tree
    search_lims.push_back( std::make_pair( 0, fExps[N].size() ) ); //by default search in full range 
  }

  // loop over summary events

  Double_t run_nr = -1;
  Double_t e_min  = -1;

  for (Int_t i = 0; i < sp.fChain->GetEntries(); i++) {

    sp.fChain->GetEntry(i);

    // if file changes update the limits in the vectors in which the events are searched

    if ( run_nr != sp.MC_runID || e_min != sp.MC_erange_start ) {

      run_nr = sp.MC_runID;
      e_min  = sp.MC_erange_start;

      for (Int_t N = 0; N < (Int_t)fExps.size(); N++) {
	search_lims[N] = get_start_stop( fExps[N], run_nr, e_min );
      }

    }

    evtid this_evt(sp.MC_runID, sp.MC_evtID, sp.MC_erange_start);

    // see if this summary event is in any of the event samples, if so write it out to
    // relevant file

    for (Int_t N = 0; N < (Int_t)fExps.size(); N++) {

      if ( std::binary_search( fExps[N].begin() + search_lims[N].first, 
			       fExps[N].begin() + search_lims[N].second, 
			       this_evt ) ) {
	trees[N]->Fill();
      }

    }

  } //end loop over summary events

  // write trees and hists; close the files; remove files where sampling failed; reset fExp.. vecs
  
  for (Int_t N = 0; N < (Int_t)fExps.size(); N++) {

    files[N]->cd();
    Bool_t remove = ( trees[N]->GetEntries() == 0 );
    TString name  = files[N]->GetName();
    trees[N]->Write();

    for (auto h: fExpHists[N]) {
      h->Write();
      delete h;
    }

    files[N]->Close();
    delete files[N];

    if (remove) {
      cout << "NOTICE WriteToFiles() empty tree in " << name << ", removing file." << endl;
      Int_t sysret = system("rm " + name);
    }

  }

  fExps.clear();
  fExpHists.clear();
  fExpNames.clear();

}

//*****************************************************************

/**
 *  Inline function to clear dynamically allocated memory.
 */
void GSGS::CleanUp() {

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

  if (fRand) delete fRand;
  
}

//*****************************************************************

/**
 *  Inline function to generate summary file name.
 *
 * \param flavor   nu flavor
 * \param is_cc    0 - NC, 1 - CC
 * \param emin     gSeaGen energy range start
 * \param emax     gSeaGen energy range stop
 * \param runnr    gSeaGen run number
 * \return         summary file name as generated by NMHDIR/data_sorting/RestoreParity.C
 */
TString GSGS::GetSummaryName(Int_t flavor, Int_t is_cc, Int_t emin, Int_t emax, Int_t runnr) {

  TString suffix = fInts[is_cc] + "_" + (TString)to_string(emin) + "-" + (TString)to_string(emax) + 
    "GeV_" + (TString)to_string(runnr) + ".root";
   
  TString flav = fFlavs[flavor];
  if (is_cc == 0) flav = fFlavs[0];
  
  TString dir = getenv("NMHDIR");
  return dir + "/data/mc_end/data_atmnu/summary_" + flav + "-" + suffix;  
}
