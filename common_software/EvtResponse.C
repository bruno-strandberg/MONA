#include "EvtResponse.h"
#include "FileHeader.h"
#include <iostream>
#include <vector>
#include <stdexcept>

#include "TCanvas.h"
#include "TH2.h"
#include "TGraph.h"
#include "TFile.h"
#include "TTree.h"

using namespace std;

/** Constructor to initialise the response with the same binning configuration in true space and reco space.
    See the constructor with different binning configurations for more info.
*/
EvtResponse::EvtResponse(reco reco_type, TString resp_name,
			 Int_t ebins , Double_t emin , Double_t emax ,
			 Int_t ctbins, Double_t ctmin, Double_t ctmax,
			 Int_t bybins, Double_t bymin, Double_t bymax,
			 Double_t memlim) : EvtResponse(reco_type, resp_name, ebins, emin, emax, ctbins, ctmin, ctmax, bybins, bymin, bymax, ebins, emin, emax, ctbins, ctmin, ctmax, bybins, bymin, bymax, memlim) {};

//=====================================================================================================

/** Constructor to initialise the response with a binning configuration in true space and reco space.

    The true space binning can be used in the software that uses the `EvtResponse` class for the dimensioning of e.g. flux cache. For all parameters other than memlim, see `AbsResponse` corresponding constructor. The additional parameterer memlim is to limit the RAM the response can occupy; if exceeded, throw an error and let user to take action. If the user has access to a machine with e.g. 16gb of memory, this number can be increased.
*/
EvtResponse::EvtResponse(reco reco_type, TString resp_name,
			 Int_t t_ebins , Double_t t_emin , Double_t t_emax ,
			 Int_t t_ctbins, Double_t t_ctmin, Double_t t_ctmax,
			 Int_t t_bybins, Double_t t_bymin, Double_t t_bymax,
			 Int_t r_ebins , Double_t r_emin , Double_t r_emax ,
			 Int_t r_ctbins, Double_t r_ctmin, Double_t r_ctmax,
			 Int_t r_bybins, Double_t r_bymin, Double_t r_bymax,
			 Double_t memlim) :
  AbsResponse(reco_type, resp_name,
	      t_ebins, t_emin, t_emax, t_ctbins, t_ctmin, t_ctmax, t_bybins, t_bymin, t_bymax,
	      r_ebins, r_emin, r_emax, r_ctbins, r_ctmin, r_ctmax, r_bybins, r_bymin, r_bymax) {

  fNormalised = false;
  fNEvts = 0;
  fMemLim = memlim;
  
  //----------------------------------------------------------
  // calculate the binning for the fResp structure and init fResp. bin [0] is underflow, bin[ axis->GetNbins() ]
  // is the last counting bin (hence dimension length +1), bin[ axis->GetNbins()+1 ] is overflow
  // (hence dimension length + 2)
  //----------------------------------------------------------

  fEbins  = fhBinsReco->GetXaxis()->GetNbins() + 2;
  fCtbins = fhBinsReco->GetYaxis()->GetNbins() + 2;
  fBybins = fhBinsReco->GetZaxis()->GetNbins() + 2;

  InitResponse( fEbins, fCtbins, fBybins );

  fhAtmMuCount1y = CloneFromTemplate(fhBinsReco, "hAtmMuCount1y_" + fRespName);
  fhNoiseCount1y = CloneFromTemplate(fhBinsReco, "hNoiseCount1y_" + fRespName);

}

//=====================================================================================================

/** Destructor */
EvtResponse::~EvtResponse() {

  CleanResponse();

  delete fhAtmMuCount1y;
  delete fhNoiseCount1y;
  
}

//=====================================================================================================

/** Function that checks the input summary event against cuts configured for this response and fills the events that pass the cuts to the response.
    \param evt   A pointer to a `SummaryEvent` with MC data
 */
void EvtResponse::Fill(SummaryEvent *evt) {

  if (fNormalised) {
    throw std::logic_error("ERROR! EvtResponse::Fill() cannot fill an already normalised response!");
  }
  
  // check that the event is recognised
  UInt_t flav;

  try {
    flav  = fType_to_Supported.at( (UInt_t)TMath::Abs( evt->Get_MC_type() ) );
  }
  catch (const std::out_of_range& oor) {
    throw std::invalid_argument("ERROR! EvtResponse::Fill() unknown particle type " +
				to_string( TMath::Abs( evt->Get_MC_type() ) ) );
  }

  UInt_t is_cc = evt->Get_MC_is_CC();
  UInt_t is_nb = (UInt_t)(evt->Get_MC_type() < 0 );

  CountEvents(flav, is_cc, is_nb, evt);

  // set the reconstruction observables, implemented in EventFilter.C
  SetObservables(evt);
  
  // fill neutrino events
  if ( flav <= TAU ) {

    // exclude events that are outside of the true bins range
    Double_t nebins  = fhBinsTrue->GetXaxis()->GetNbins();
    Double_t emin    = fhBinsTrue->GetXaxis()->GetBinLowEdge(1);
    Double_t emax    = fhBinsTrue->GetXaxis()->GetBinUpEdge(nebins);

    Double_t nctbins = fhBinsTrue->GetYaxis()->GetNbins();
    Double_t ctmin   = fhBinsTrue->GetYaxis()->GetBinLowEdge(1);
    Double_t ctmax   = fhBinsTrue->GetYaxis()->GetBinUpEdge(nctbins);

    Double_t nbybins = fhBinsTrue->GetZaxis()->GetNbins();
    Double_t bymin   = fhBinsTrue->GetZaxis()->GetBinLowEdge(1);
    Double_t bymax   = fhBinsTrue->GetZaxis()->GetBinUpEdge(nbybins);

    if (  evt->Get_MC_energy()   < emin  ||  evt->Get_MC_energy()   >= emax  ||
	 -evt->Get_MC_dir_z()    < ctmin || -evt->Get_MC_dir_z()    >= ctmax ||
	  evt->Get_MC_bjorkeny() < bymin ||  evt->Get_MC_bjorkeny() >= bymax ) {      
      return;
    }
    
    // exclude events that do not pass the selection cuts
    if ( !PassesCuts(evt) ) return;

    // find the reco bin based on reco variables and fill the response
    Int_t  e_reco_bin = fhBinsReco->GetXaxis()->FindBin(  fEnergy   );
    Int_t ct_reco_bin = fhBinsReco->GetYaxis()->FindBin( -fDir.z()  );
    Int_t by_reco_bin = fhBinsReco->GetZaxis()->FindBin(  fBy       );

    fResp[e_reco_bin][ct_reco_bin][by_reco_bin].push_back( TrueEvt(flav, is_cc, is_nb, evt) );
    fNEvts++;
  }

  // fill noise and atm muon events
  else {

    // if event does not pass the cuts return
    if ( !PassesCuts(evt) ) return;
    
    if      ( flav == ATMMU ) { fhAtmMuCount1y->Fill( fEnergy, -fDir.z(), fBy ); }
    else if ( flav == NOISE ) { fhNoiseCount1y->Fill( fEnergy, -fDir.z(), fBy ); }
    else {
      throw std::invalid_argument( "ERROR! EvtResponse::Fill() unknown particle with flavor " + to_string(flav) );
    }
 
  }

  // throw an error if too much RAM is eaten up
  Double_t ram_used = (Double_t)sizeof( TrueEvt ) * fNEvts * 1e-9; // sizeof [byte] * # of TrueEvt * 1e-9 [gbyte/byte]
  if ( ram_used > fMemLim ) {
    throw std::invalid_argument("ERROR! EvtResponse::Fill() the " + to_string(fNEvts) + " filled events occupy " + to_string(ram_used) + " gb of RAM, which is above the limit " + to_string(fMemLim) + " specified in the constructor. Either reduce the MC sample filled to the response or increate the memory limit at construction." );
  }
  
}

//=====================================================================================================

/** Private function to count the data necessary for the calculation of weight 1 year.

    To calculate the event weight per year for neutrinos, the W2 from gSeaGen has to be divided by N_tot, where N_tot is the total number of generated events over several runs. For muons and noise, the weight is calculated as 1/total_livetime_sec * sec_per_year, where total_livetime_sec is accumulated from e.g. mupage runs. The book-keeping is further complicated by the fact that there can be overlaps, i.e. there can be muon-CC productions in the range 1-5 GeV and 3-100 GeV that are fed to this response.

    Continuing with the example for muon-CC 1-5 and 3-100, this function does the following. For each neutrino type (in this case [muon][cc=1][nu] and [muon][cc][nubar]) the `EvtResponse` stores a map with structure < rangeID, vector<runID> >. Range ID is just a pair with upper and lower limit, the vector of runID's contains all the different run numbers with their respective N_tot from gSeaGen. So, for muon-CC 1-5 and 3-100, the map [muon][cc][nu] will contain two elements, one pair < {1,5}, vector<runID>{all runs} > and < {3,100}, vector<runID>{all runs} >. With this data, after the response has been filled one can calculate the total number of generated events over all of the runs filled to the response to re-normalise W2 to weight-one-year.

    \param flav Flavor (0 - elec, 1 - muon, 2 - tau, 3 - atmmu, 4 - noise)
    \param iscc 0 - NC event, 1 - CC event ( for noise and muons ignored )
    \param isnb 0 - nu event, 1 - anti-nu event ( for noise and muons ignored )
    \param evt  Pointer to a `SummaryEvent` with event data

*/
void  EvtResponse::CountEvents(UInt_t flav, UInt_t iscc, UInt_t isnb, SummaryEvent *evt) {

  rangeID range( evt->Get_MC_erange_start(), evt->Get_MC_erange_stop() ); // emin and emax of the gSeaGen simulation range; some default value for noise and muons
  runID   run  ( evt->Get_MC_runID(), evt->Get_MC_w2denom() );            // run ID and N_tot(gSeaGen) or livetime (mupage, noise)

  bool _iscc = iscc;
  bool _isnb = isnb;

  // for muons and noise there is only 1 type of event
  if ( evt->Get_MC_is_neutrino() < 0.5 ) {
    _iscc = 0;
    _isnb = 0;
  }

  auto& M = fRuns[flav][_iscc][_isnb]; // map with < range,vec<run> > for this event type
  
  auto element = M.find( range ); // try to find by the range in the map

  //-----------------------------------------------------------------
  // no vector exists for this range in the map, insert it
  //-----------------------------------------------------------------
  if ( element == M.end() ) {
    M.insert( std::make_pair( range, std::vector<runID>{ run } ) );
  }

  //-----------------------------------------------------------------
  // vector exists, see if the run exists already in the vector; if not, enter it
  //-----------------------------------------------------------------
  else {

    auto& V = element->second; // get the vector
    
    // try to find the runID in the vector
    auto RID = std::find_if(V.begin(), V.end(), [run](const runID& other){ return other.first == run.first; });

    // if found, check that N_tot is also the same and continue;
    if ( RID != V.end() ) {

      if ( RID->second != run.second ) {
	throw std::invalid_argument("ERROR! EvtResponse::CountEvents() runID " + to_string(run.first) + " identified with N_tot/livetime " + to_string(run.second) + " and " + to_string( RID->second ) );
      }
      
    }
    else {

      // run is not found, add it to the vector
      V.push_back( run );
      
    }
    
  } // end if range exists
  
}

//=====================================================================================================

/** Private function that converts the W2 weights of each event to weight_per_year.

    The gSeaGen weight w2 can be converted to weight per year by dividing it by the total number N_tot of generated events by gSeaGen. If there are several runs, N_tot is the sum over N_tot in individual files.

    For muons & noise, the weight per 1 year can be calculated as 1/livetime * sec_per_year.
*/
void  EvtResponse::Normalise() {

  // first count together N_tot for each neutrino type
  std::map< rangeID, Double_t > N_TOT[TAU+1][2][2];

  for (Int_t F = 0; F <= TAU; F++) {
    for (Int_t iscc = 0; iscc <= 1; iscc++) {
      for (Int_t isnb = 0; isnb <= 1; isnb++) {

	// get the map with < rangeID, vector<runID> > for this neutrino type
	auto& M = fRuns[ F ][ iscc ][ isnb ];

	// loop over the ranges, for each range calculate the total number of simulated events
	for (auto kv: M) {

	  Double_t ntot = 0.;
	  auto& range   = kv.first;
	  auto& runs    = kv.second;
	  for (auto run: runs) { ntot += run.second; }

	  N_TOT[F][iscc][isnb].insert( std::make_pair(range, ntot) );
	}

      }
    }
  }
  
  // loop over all of the events in the response (neutrino events only)
  for (Int_t ebin = 0; ebin < fEbins; ebin++) {
    for (Int_t ctbin = 0; ctbin < fCtbins; ctbin++) {
      for (Int_t bybin = 0; bybin < fBybins; bybin++) {

	std::vector<TrueEvt>& binevts = fResp[ebin][ctbin][bybin];

	for (auto& evt: binevts) {

	  Double_t w2denom_tot = 0;

	  // get the map with < rangeID, ntot > for this neutrino type
	  auto& M = N_TOT[ evt.GetFlav() ][ evt.GetIsCC() ][ evt.GetIsNB() ];
	  
	  // loop over the energy ranges of this neutrino type
	  for (auto kv: M) {

	    rangeID  range = kv.first;
	    Double_t ntot  = kv.second;

	    // if the event is inside the range, accumuate Ntot  (saved in w2denom for neutrinos)
	    if ( evt.GetTrueE() >= range.first && evt.GetTrueE() < range.second ) {
	      w2denom_tot += ntot;
	    }
	    
	  } // end loop over ranges

	  // convert W2 (stored in w1y) to weight in one year
	  evt.SetW1y( evt.GetW1y()/w2denom_tot );
	  
	} // end loop over TrueEvt's
	
      }
    }
  }

  // normalise muons and noise
  vector <UInt_t> others = {ATMMU, NOISE};

  for (auto t: others) {

    auto M = fRuns[t][0][0]; // for atm muons and noise I only use 1 type

    TH3D *h;
    if      ( t == ATMMU ) h = fhAtmMuCount1y;
    else if ( t == NOISE ) h = fhNoiseCount1y;
    else {
      throw std::logic_error("ERROR! EvtResponse::Normalise() unexpected type in muon & noise normalisation");
    }
    
    if ( M.size() > 1 ) {
      throw std::logic_error("ERROR! EvtResponse::Normalise() expecting 1 E-range for atm. mu or noise");
    }
    else if ( M.size() == 1 ) {

      vector<runID> runs = M.begin()->second;
      Double_t tot_livetime = 0;
      for (auto RID: runs) { tot_livetime += RID.second; }
      
      h->Scale( 1./tot_livetime * fSec_per_y );
      
    }
    else {

      // no muons or noise runs, hence the histogram should be empty
      if (h->GetEntries() != 0) {
	std::map<Int_t, string> helper;
	helper.insert( std::make_pair<Int_t, string> (ATMMU, "atm. muon") );
	helper.insert( std::make_pair<Int_t, string> (NOISE, "noise") );
	throw std::logic_error("ERROR! EvtResponse::Normalise() no muons/noise runs identified, but relevant histograms not empty for " + helper[t] );
      }
      
    }
    
  }
  
  fNormalised = kTRUE;
  
}

//=====================================================================================================

/** Function to print the neutrino run data of the events that have been filled to the response.

    This data is used by the function `EvtResponse::Normalise()` for the calculation of weight_1_year.
*/
void EvtResponse::PrintRunData() {

  cout << "NOTICE EvtResponse::PrintRunData() printing data of the runs filled to the response" << endl;

  // print neutrino run data
  for (Int_t F = 0; F <= TAU; F++) {
    for (Int_t iscc = 0; iscc <= 1; iscc++) {
      for (Int_t isnb = 0; isnb <= 1; isnb++) {
	
	auto M = fRuns[F][iscc][isnb];

	for (auto kv: M) {

	  auto range = kv.first;
	  auto runs  = kv.second;
	  
	  cout << "NOTICE EvtResponse::PrintRunData() Neutrino (f, iscc, isnb): " << F << " " << iscc << " " << isnb << " (emin, emax): " << range.first << " " << range.second << "; runs " << runs.size() << "; (runnr, ntot) of first element: " << runs[0].first << " " << runs[0].second << endl;
	  
	}

      }
    }
  }

  // print atm. muon and noise data
  for (Int_t F = ATMMU; F <= NOISE; F++) {
    auto M = fRuns[F][0][0];

    for (auto kv: M) {

      auto range = kv.first;
      auto runs  = kv.second;
	  
      cout << "NOTICE EvtResponse::PrintRunData() atm. muon/noise data (emin, emax): " << range.first << " " << range.second << "; runs " << runs.size() << "; (runnr, livetime) of first element: " << runs[0].first << " " << runs[0].second << endl;
	  
    }

  }
  
}

//=====================================================================================================

/** Private function to faciliate easier cloning of histograms
    \param tmpl   Template `TH3D` histogram
    \param name   Name and title for the created histogram
    \return       Pointer to the new empty histogram with identical binning compared to the template
*/
TH3D* EvtResponse::CloneFromTemplate(TH3D* tmpl, TString name) {

  TH3D* ret = (TH3D*)tmpl->Clone(name);
  ret->SetNameTitle(name, name);
  ret->Reset();
  ret->SetDirectory(0);

  return ret;

}

//=====================================================================================================

/** Function to get a vector of `TrueEvt`'s that contributed to the reco bin located at the provided coordinates.
    \param   E_reco  reco energy
    \param   ct_reco reco cos-theta
    \param   by_reco reco bjorken-y
    \return  vector of `TrueEvt`'s that contributed to the reco bin
*/
std::vector<TrueEvt>& EvtResponse::GetBinEvts(Double_t E_reco, Double_t ct_reco, Double_t by_reco) {

  if (!fNormalised) Normalise();
  
  Int_t ebin  = fhBinsReco->GetXaxis()->FindBin(E_reco);
  Int_t ctbin = fhBinsReco->GetYaxis()->FindBin(ct_reco);
  Int_t bybin = fhBinsReco->GetZaxis()->FindBin(by_reco);

  return fResp[ebin][ctbin][bybin];
  
}

//=====================================================================================================

/**
   Returns the number of atmospheric muon events in 1 year in the bin specified by the reconstruction variables.
   
   \param E_reco    Reconstructed energy
   \param ct_reco   Reconstructed cos-theta
   \param by_reco   Reconstructed bjorken-y
   \return          a pair with the atmospheric muon count in 1 year (first) and the MC statistical error (second)
 */
std::pair<Double_t, Double_t> EvtResponse::GetAtmMuCount1y(Double_t E_reco, Double_t ct_reco, Double_t by_reco) {

  Int_t ebin  = fhAtmMuCount1y->GetXaxis()->FindBin(E_reco);
  Int_t ctbin = fhAtmMuCount1y->GetYaxis()->FindBin(ct_reco);
  Int_t bybin = fhAtmMuCount1y->GetZaxis()->FindBin(by_reco);
  
  return std::make_pair( fhAtmMuCount1y->GetBinContent(ebin, ctbin, bybin), fhAtmMuCount1y->GetBinError(ebin, ctbin, bybin) );

}

//=====================================================================================================

/**
   Returns the number of noise events in 1 year in the bin specified by the reconstruction variables.
   
   \param E_reco    Reconstructed energy
   \param ct_reco   Reconstructed cos-theta
   \param by_reco   Reconstructed bjorken-y
   \return          a pair with noise count in 1 year (first) and the MC statistical error (second)
 */
std::pair<Double_t, Double_t> EvtResponse::GetNoiseCount1y(Double_t E_reco, Double_t ct_reco, Double_t by_reco) {

  Int_t ebin  = fhNoiseCount1y->GetXaxis()->FindBin(E_reco);
  Int_t ctbin = fhNoiseCount1y->GetYaxis()->FindBin(ct_reco);
  Int_t bybin = fhNoiseCount1y->GetZaxis()->FindBin(by_reco);
  
  return std::make_pair( fhNoiseCount1y->GetBinContent(ebin, ctbin, bybin), fhNoiseCount1y->GetBinError(ebin, ctbin, bybin) );

}

//=====================================================================================================

/** Function to display the true events contributing to the reco bin at the input values
    
    A summation over the bjorken-y bins is performed, as visualisation in 3D is difficult.

    \param e_reco  Reco energy
    \param ct_reco Reco cos-theta
    \param outname If specified, graphs are written to the output file
    \return        Canvas with the visualisation

*/
TCanvas* EvtResponse::DisplayResponse(Double_t e_reco, Double_t ct_reco, TString outname) {

   // graphs with true events in the reco bin by nu type
  TGraph G[fFlavs.size()][fInts.size()][fPols.size()];
  vector<TGraph*> glist;
  
  TString suffix = "_Ereco=" + (TString)to_string(e_reco) + "_ctreco=" + (TString)to_string(ct_reco);
  
  for (auto f: fFlavs) {
    for (auto i: fInts) {
      for (auto p: fPols) {
	TGraph& gr = G[f.first][i.first][p.first];
	TString gname = "events_" + f.second + "_" + i.second + "_" + p.second + suffix;
	gr.SetNameTitle(gname, gname);
	gr.GetXaxis()->SetTitle("E^{true}_{nu} [GeV]");
	gr.GetYaxis()->SetTitle("cos#theta^{true}_{nu}");
	glist.push_back( &gr );
      }
    }
  }

  // sum over bjorken-y
  vector<TrueEvt> TE;
  for (Int_t byb = 1; byb <= fhBinsReco->GetZaxis()->GetNbins(); byb++) {
    Double_t bjorkeny = fhBinsReco->GetZaxis()->GetBinCenter( byb );
    auto te = GetBinEvts(e_reco, ct_reco, bjorkeny);
    for (auto t: te) TE.push_back( t );
  }

  //create the histogram with reco bin
  TH2D* hr = (TH2D*)fhBinsReco->Project3D("yx")->Clone("RecoBin_" + suffix);
  hr->SetDirectory(0);
  hr->Reset();
  
  // add points to the graphs and the histogram
  for (auto t: TE) {
    TGraph& gr = G[ t.GetFlav() ][ t.GetIsCC() ][ t.GetIsNB() ];
    gr.SetPoint( gr.GetN(), t.GetTrueE(), t.GetTrueCt() );
    hr->Fill(e_reco, ct_reco);
  }

  // set point styles and ranges; only draw graphs that are not empty
  vector<TGraph*> _glist;
  for (auto gr: glist) {
    gr->SetMinimum(-1);
    gr->SetMaximum(1);
    gr->GetXaxis()->SetRangeUser(1,100);
    gr->SetMarkerStyle(5);
    gr->SetMarkerColor(kBlue);
    if ( gr->GetN() > 0 ) _glist.push_back( gr );
  }
  
  // draw the graphs and the true bin
  TCanvas *c1 = new TCanvas("c1","c1",1);
  c1->DivideSquare( _glist.size() + 1 );

  c1->cd(1);
  hr->Draw("colz");
  
  for (Int_t i = 0; i < (Int_t)_glist.size(); i++) {
    c1->cd(i+2);
    ((TGraph*)_glist[i]->Clone())->Draw("AP");
  }

  if (outname != "") {
    TFile fout(outname, "RECREATE");
    for (auto g: _glist) g->Write();
    fout.Close();
  }
  
  return c1;
  
}

//=====================================================================================================

/**
   Function to write the response to file.

   \param filename  Name of the output file.

 */
void EvtResponse::WriteToFile(TString filename) {

  if (!fNormalised) Normalise();

  FileHeader h("EvtResponse");
  h.AddParameter( "fRespName", fRespName);
  h.AddParameter( "fEbins"   , (TString)to_string(fEbins) );
  h.AddParameter( "fCtbins"  , (TString)to_string(fCtbins) );
  h.AddParameter( "fBybins"  , (TString)to_string(fBybins) );
  
  TFile fout(filename, "RECREATE");
  TTree tout("evtresponse","Detector response data");

  TrueEvt *te = new TrueEvt();
  Int_t E_reco_bin, ct_reco_bin, by_reco_bin;

  tout.Branch("E_reco_bin" , &E_reco_bin , "E_reco_bin/I");
  tout.Branch("ct_reco_bin", &ct_reco_bin, "ct_reco_bin/I");
  tout.Branch("by_reco_bin", &by_reco_bin, "by_reco_bin/I");
  tout.Branch("TrueEvt", &te, 2);

  for (Int_t ebin = 0; ebin < fEbins; ebin++) {
    for (Int_t ctbin = 0; ctbin < fCtbins; ctbin++) {
      for (Int_t bybin = 0; bybin < fBybins; bybin++) {

	E_reco_bin  = ebin;
	ct_reco_bin = ctbin;
	by_reco_bin = bybin;

	for (auto true_evt: fResp[ebin][ctbin][bybin]) {
	  *te = true_evt;
	  tout.Fill();
	}

      }
    }
  }

  tout.Write();

  fhBinsTrue->Write("hbinstrue");
  fhBinsReco->Write("hbinsreco");
  fhAtmMuCount1y->Write("atmmucount");
  fhNoiseCount1y->Write("noisecount");

  h.WriteHeader(&fout);

  fout.Close();
  delete te;
  
}

//=====================================================================================================

/**
   Function to read the response from file.

   \param filename   Name of the file where the response is stored.
 */
void EvtResponse::ReadFromFile(TString filename) {

  // if the evtresponse is read from file, the binning is changed to match that
  // of the evtresponse that is being read in
  CleanResponse();

  FileHeader h("for_read_in");
  h.ReadHeader(filename);

  fRespName   = h.GetParameter("fRespName");
  fNormalised = true;
  fEbins      = std::stoi( (string)h.GetParameter("fEbins")  );
  fCtbins     = std::stoi( (string)h.GetParameter("fCtbins") );
  fBybins     = std::stoi( (string)h.GetParameter("fBybins") );

  InitResponse(fEbins, fCtbins, fBybins);

  // read in the true bin data to the response

  TFile fin(filename, "READ");
  TTree *tin = (TTree*)fin.Get("evtresponse");

  TrueEvt *te = new TrueEvt();
  Int_t E_reco_bin, ct_reco_bin, by_reco_bin;

  tin->SetBranchAddress("E_reco_bin", &E_reco_bin);
  tin->SetBranchAddress("ct_reco_bin", &ct_reco_bin);
  tin->SetBranchAddress("by_reco_bin", &by_reco_bin);
  tin->SetBranchAddress("TrueEvt", &te);

  for (Int_t i = 0; i < tin->GetEntries(); i++) {
    tin->GetEntry(i);
    fResp[E_reco_bin][ct_reco_bin][by_reco_bin].push_back( TrueEvt(*te) );
  }

  // get the histgrams

  fhBinsTrue = (TH3D*)fin.Get("hbinstrue")->Clone();
  fhBinsTrue->SetDirectory(0);
  fhBinsTrue->SetNameTitle("hbinstrue_" + fRespName, "hbinstrue_" + fRespName);

  fhBinsReco = (TH3D*)fin.Get("hbinsreco")->Clone();
  fhBinsReco->SetDirectory(0);
  fhBinsReco->SetNameTitle("hbinsreco_" + fRespName, "hbinsreco_" + fRespName);

  fhAtmMuCount1y = (TH3D*)fin.Get("atmmucount")->Clone();
  fhAtmMuCount1y->SetDirectory(0);
  fhAtmMuCount1y->SetNameTitle("hAtmMuCount1y_" + fRespName, "hAtmMuCount1y_" + fRespName);

  fhNoiseCount1y = (TH3D*)fin.Get("noisecount")->Clone();
  fhNoiseCount1y->SetDirectory(0);
  fhNoiseCount1y->SetNameTitle("hNoiseCount1y_" + fRespName, "hNoiseCount1y_" + fRespName);

  delete te;
  fin.Close();

}

//=====================================================================================================

/**
   Private function to de-allocate memory of the member `fResp` structure
 */
void EvtResponse::CleanResponse() {

  for (Int_t ebin = 0; ebin < fEbins; ebin++) {
    for (Int_t ctbin = 0; ctbin < fCtbins; ctbin++) {
      if (fResp[ebin][ctbin]) delete[] fResp[ebin][ctbin];
    }
    if (fResp[ebin]) delete[] fResp[ebin];
  }
  delete[] fResp;

}

//=====================================================================================================

/**
   Private function to initialise the member `fResp` structure.
   \param ebins  Number of energy bins
   \param ctbins Number of cos-theta bins
   \param bybins Number of bjorken-y bins
 */
void EvtResponse::InitResponse(Int_t ebins, Int_t ctbins, Int_t bybins) {

  //----------------------------------------------------------
  // initialise the response structure. This is kind-of like a TH3D, but instead of doubles, the array at
  // [xbin][ybin][zbin] stores a vector of TrueEvt's. Thus it allows to store all of the MC data on
  // event-by-event basis, sorted to reconstruction bins.
  //----------------------------------------------------------
  
  fResp = new vector<TrueEvt>** [ebins]();
  for (Int_t ebin = 0; ebin < ebins; ebin++) {
    fResp[ebin] = new vector<TrueEvt>* [ctbins]();
    for (Int_t ctbin = 0; ctbin < ctbins; ctbin++) {
      fResp[ebin][ctbin] = new vector<TrueEvt> [bybins]();
    }
  }

}
