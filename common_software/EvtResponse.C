#include "EvtResponse.h"
#include <vector>

using namespace std;

/** Constructor.
    For all parameters other than memlim, see `AbsResponse` corresponding constructor.
    \param memlim RAM limit the response can occupy; if exceeded, throw an error and let user to take action. If the user has access to a machine with e.g. 16gb of memory, this number can be increased.
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
  
}

//=====================================================================================================

/** Function that checks the input summary event against cuts and fills the events that pass the cuts to response
    \param evt summary event
 */
void EvtResponse::Fill(SummaryEvent *evt) {

  //---------------------------------------------------------------------------------------
  // check that type is supported
  //---------------------------------------------------------------------------------------

  try {
    fType_to_Supported.at( (UInt_t)TMath::Abs( evt->Get_MC_type() ) );
  }
  catch (const std::out_of_range& oor) {
    throw std::invalid_argument("ERROR! EvtResponse::Fill() unknown particle type " +
				to_string( TMath::Abs( evt->Get_MC_type() ) ) );
  }

  /*********************************************************************************
   Any logic here to do something different with atm muons and noise? 
  **********************************************************************************/
  //---------------------------------------------------------------------------------------
  // if the event does not pass the cuts return
  //---------------------------------------------------------------------------------------

  if ( !PassesCuts(evt) ) return;

  //---------------------------------------------------------------------------------------
  // set reco observables, locate the reco bin this event belongs to and add to response
  //---------------------------------------------------------------------------------------

  SetObservables(evt); //implemented in EventFilter.C

  Int_t  e_reco_bin = fhBinsReco->GetXaxis()->FindBin(  fEnergy   );
  Int_t ct_reco_bin = fhBinsReco->GetYaxis()->FindBin( -fDir.z()  );
  Int_t by_reco_bin = fhBinsReco->GetZaxis()->FindBin(  fBy       );

  fResp[e_reco_bin][ct_reco_bin][by_reco_bin].push_back( TrueEvt(evt) );

  /*********************************************************************************
   Add test here to check the RAM usage once TrueEvt class is finalised
  **********************************************************************************/

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