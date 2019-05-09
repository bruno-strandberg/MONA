#include "EvtResponse.h"
#include <iostream>
#include <vector>

using namespace std;

/** Constructor.
    For all parameters other than memlim, see `AbsResponse` corresponding constructor. The additional parameterer memlim is to limit the RAM the response can occupy; if exceeded, throw an error and let user to take action. If the user has access to a machine with e.g. 16gb of memory, this number can be increased.
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

/** Function that checks the input summary event against cuts and fills the events that pass the cuts to response
    \param evt summary event
 */
void EvtResponse::Fill(SummaryEvent *evt) {

  //---------------------------------------------------------------------------------------
  // check that type is supported
  //---------------------------------------------------------------------------------------

  UInt_t flav;
  try {
    flav  = fType_to_Supported.at( (UInt_t)TMath::Abs( evt->Get_MC_type() ) );
  }
  catch (const std::out_of_range& oor) {
    throw std::invalid_argument("ERROR! DetResponse::Fill() unknown particle type " +
				to_string( TMath::Abs( evt->Get_MC_type() ) ) );
  }

  UInt_t is_cc = evt->Get_MC_is_CC();
  UInt_t is_nb = (UInt_t)(evt->Get_MC_type() < 0 );

  /*********************************************************************************
   Any logic here to do something different with atm muons and noise? 
  **********************************************************************************/
  

  // set the reconstruction observables
  SetObservables(evt); //implemented in EventFilter.C


  if	  ( flav == ELEC || flav == MUON || flav == TAU){
	
	  if ( !PassesCuts(evt) ) return;

	  //---------------------------------------------------------------------------------------
	  // set reco observables, locate the reco bin this event belongs to and add to response
	  //---------------------------------------------------------------------------------------

	  //SetObservables(evt); //implemented in EventFilter.C

	  Int_t  e_reco_bin = fhBinsReco->GetXaxis()->FindBin(  fEnergy   );
	  Int_t ct_reco_bin = fhBinsReco->GetYaxis()->FindBin( -fDir.z()  );
	  Int_t by_reco_bin = fhBinsReco->GetZaxis()->FindBin(  fBy       );


	  fResp[e_reco_bin][ct_reco_bin][by_reco_bin].push_back( TrueEvt(flav, is_cc, is_nb, evt) );
  }
  else if ( flav == ATMMU ) fhAtmMuCount1y->Fill( fEnergy, -fDir.z(), fBy, evt->Get_MC_w1y() );
  else if ( flav == NOISE ) fhNoiseCount1y->Fill( fEnergy, -fDir.z(), fBy, evt->Get_MC_w1y() );
  else {
      throw std::invalid_argument( "ERROR! DetResponse::Fill() unknown particle with flavor " + to_string(flav) );
  }

  //---------------------------------------------------------------------------------------
  // if the event does not pass the cuts return
  //---------------------------------------------------------------------------------------



  /*********************************************************************************
   Add test here to check the RAM usage once TrueEvt class is finalised
  **********************************************************************************/
  

}


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
