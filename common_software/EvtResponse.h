#ifndef EvtResponse_h
#define EvtResponse_h

#include "AbsResponse.h"
#include "SummaryEvent.h"

/** FILL DOCUMENTATION */
class TrueEvt {

 public:

  /** Default constructor */
 TrueEvt() : fNuType(0), fE_true(0), fCt_true(0), fBy_true(0) {};

  /** Constructor from summary event.
      Packs the necessary info from SummaryEvent to member variables.
      \param evt   Pointer to a summary event
   */
  TrueEvt(SummaryEvent *evt) {

    fNuType = (Short_t) evt->Get_MC_type();
    fE_true = (Float_t) evt->Get_MC_energy();
    fCt_true= (Float_t) -evt->Get_MC_dir_z();
    fBy_true = (Float_t) evt->Get_MC_bjorkeny();
    fIsCC = (UInt_t) evt->Get_MC_is_CC();
    f_W1y = (Float_t) evt->Get_MC_w1y();
    fIchan = 0; //it might be useful for the future cross section systematics, for now not found in PID summary files
    
  }

  /** Destructor */
  ~TrueEvt() {};

  // plus interface functions to "unpack" data from private members to double's
  
 //private:

  UInt_t fIsCC;
  Short_t fIchan;   //!< variable that hold the info about the interaction channel (dis, res, qe)
  Short_t fNuType;  //!< variable that holds all the info to determine flavor, iscc, isnb
  Float_t fE_true;  //!< variable that stores energy
  //Short_t fCt_true; //!< variable that stores cos-theta with 4-digit resolution
  Float_t fCt_true;
  //Short_t fBy_true; //!< variable that stores bjorken-y with 4-digit resolution
  Float_t fBy_true;

  Float_t f_W1y;    //!< weight 1 year, for testing, to be deleted in the final version

  // plus other members necessary for weight calculation...
  
};

/** This class implements an event-by-event detector response.
    
    It allows the configuration of filters on `SummaryEvent`'s through inheriting functionality from `EventFilter`. It stores each event that passes the filters in a response matrix. The response matrix can be used later to calculate the weight (which depends on osc. parameters and systematics) for each event contributing to a specific reco bin to arrive at a number of expected reconstructed events in a certrain reco bin.

    FILL DOCUMENTATION 
*/
class EvtResponse : public AbsResponse {

 public:

  EvtResponse(reco reco_type, TString resp_name,
	      Int_t t_ebins , Double_t t_emin , Double_t t_emax ,
  	      Int_t t_ctbins, Double_t t_ctmin, Double_t t_ctmax,
  	      Int_t t_bybins, Double_t t_bymin, Double_t t_bymax,
	      Int_t r_ebins , Double_t r_emin , Double_t r_emax ,
  	      Int_t r_ctbins, Double_t r_ctmin, Double_t r_ctmax,
  	      Int_t r_bybins, Double_t r_bymin, Double_t r_bymax,
	      Double_t memlim = 4);

  ~EvtResponse();

  void Fill(SummaryEvent *evt);
  resp GetResponseType() { return AbsResponse::EvtResponse; }
  std::vector<TrueEvt>& GetBinEvts(Double_t E_reco, Double_t ct_reco, Double_t by_reco);

 private:

  Double_t fMemLim;               //!< user-defined limit to how much RAM the response can eat up
  Int_t fEbins;                   //!< number of reco energy bins in `fResp`
  Int_t fCtbins;                  //!< number of reco cos-theta bins in `fResp`
  Int_t fBybins;                  //!< number of reco bjorken-y bins in `fResp`
  TH3D    *fhAtmMuCount1y;      //!< atmospheric muon count in 1 year
  TH3D    *fhNoiseCount1y;      //!< noise event count in 1 year
  std::vector<TrueEvt> ***fResp;  //!< Response structure in 3D in [Ereco][CTreco][BYreco] = vector<TrueEvt> {contributing true events in reco bin}
  
};

#endif
