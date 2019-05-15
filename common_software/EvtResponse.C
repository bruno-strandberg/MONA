#include "EvtResponse.h"
#include <iostream>
#include <vector>
#include <stdexcept>

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
  // is the last counting bin (hence dimension length +1), bin[ axis->GetNbins()+1 ] is overflow (hence dimension length + 2)
  //----------------------------------------------------------

  fEbins  = fhBinsReco->GetXaxis()->GetNbins() + 2;
  fCtbins = fhBinsReco->GetYaxis()->GetNbins() + 2;
  fBybins = fhBinsReco->GetZaxis()->GetNbins() + 2;

  //----------------------------------------------------------
  // initialise the response structure. This is kind-of like a TH3D, but instead of doubles, the array at [xbin][ybin][zbin] stores
  // a vector of TrueEvt's. Thus it allows to store all of the MC data on event-by-event basis, sorted to reconstruction bins.
  //----------------------------------------------------------

  fResp = new vector<TrueEvt>** [fEbins]();
  for (Int_t ebin = 0; ebin < fEbins; ebin++) {
    fResp[ebin] = new vector<TrueEvt>* [fCtbins]();
    for (Int_t ctbin = 0; ctbin < fCtbins; ctbin++) {
      fResp[ebin][ctbin] = new vector<TrueEvt> [fBybins]();
    }
  }

  fhAtmMuCount1y = CloneFromTemplate(fhBinsReco, "hAtmMuCount1y_" + fRespName);
  fhNoiseCount1y = CloneFromTemplate(fhBinsReco, "hNoiseCount1y_" + fRespName);
}

//=====================================================================================================

/** Destructor */
EvtResponse::~EvtResponse() {

  // de-allocate the memory used by fResp
  for (Int_t ebin = 0; ebin < fEbins; ebin++) {
    for (Int_t ctbin = 0; ctbin < fCtbins; ctbin++) {
      if (fResp[ebin][ctbin]) delete[] fResp[ebin][ctbin];
    }
    if (fResp[ebin]) delete[] fResp[ebin];
  }
  delete[] fResp;
  delete fhAtmMuCount1y;
  delete fhNoiseCount1y;
  
}

//=====================================================================================================

/** Function that checks the input summary event against cuts configured for this response and fills the events that pass the cuts to the response.
    \param evt   A pointer to a `SummaryEvent` with MC data
 */
void EvtResponse::Fill(SummaryEvent *evt) {

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
    
    if      ( flav == ATMMU ) fhAtmMuCount1y->Fill( fEnergy, -fDir.z(), fBy, evt->Get_MC_w1y() );
    else if ( flav == NOISE ) fhNoiseCount1y->Fill( fEnergy, -fDir.z(), fBy, evt->Get_MC_w1y() );
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

*/
void  EvtResponse::CountEvents(UInt_t flav, UInt_t iscc, UInt_t isnb, SummaryEvent *evt) {

  rangeID range( evt->Get_MC_erange_start(), evt->Get_MC_erange_stop() ); // emin and emax of the gSeaGen simulation range; some default value for noise and muons
  runID   run  ( evt->Get_MC_runID(), evt->Get_MC_w2denom() );            // run ID and N_tot(gSeaGen) or livetime (mupage, noise)

  auto M = fRuns[flav][bool(iscc)][bool(isnb)]; // map with < range,vec<run> > for this event type

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

    auto V = element->second; // get the vector
    
    // try to find the runID in the vector
    auto RID = std::find_if( V.begin(), V.end(), [run](const runID& other){ return other.first == run.first; } );

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

void  EvtResponse::Normalise() {

  //TBD
  
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
