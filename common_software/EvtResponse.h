#ifndef EvtResponse_h
#define EvtResponse_h

#include "AbsResponse.h"
#include "SummaryEvent.h"

/** A structure to store event data in an event-by-event detector response implemented in the class `EvtResponse`.

    The function `EvtResponse::GetBinEvts()` returns a vector of `TrueEvt` objects associated with the bin in reco space. The `TrueEvt` object stores all of the information necessary to calculate the number of expected events in a reco bin (see e.g. the class `fitter_software/FitUtil::RecoEvts`). Some data conversions are applied in this object in order to save RAM space, as all of the monte carlo events need to be read into memory for `EvtResponse`.

 */
class TrueEvt {

 public:

  /** Default constructor */
 TrueEvt() : fFlav(0), fIsCC(0), fIsNB(0), fIchan(0), fE_true(0), fCt_true(0), fBy_true(0), f_W1y(0) {};

  /** Constructor from `SummaryEvent` that packs the necessary info to member variables.
      \param flav  Neutrino flavor (0 - elec, 1 - muon, 2 - tau)
      \param iscc  Interaction type 0 - NC, 1 - CC
      \param isnb  nu or antineutrino 0 - neutrino, 1 - antineutrino
      \param evt   Pointer to a summary event
   */
  TrueEvt(UInt_t flav, UInt_t iscc, UInt_t isnb, SummaryEvent *evt) {

    if (flav > 2 || iscc > 1 || isnb > 1) {
      throw std::invalid_argument( "ERROR! TrueEvt::TrueEvt() unknown neutrino type with flavor, is_cc (flag), is_nb (flag) "
				   + std::to_string(flav) + " " + std::to_string(iscc) + " " + std::to_string(isnb) );
    }
    
    fFlav    = (Short_t)flav;
    fIsCC    = (bool)iscc;
    fIsNB    = (bool)isnb;
    fIchan   = 0;                             //future placeholder for xsec systematics, for now not used
    fE_true  = (Float_t)evt->Get_MC_energy();
    fCt_true = (Short_t)( -evt->Get_MC_dir_z()    * 1e4 );
    fBy_true = (Short_t)(  evt->Get_MC_bjorkeny() * 1e4 );
    //f_W1y    = (Float_t)evt->Get_MC_w2(); // what it will be
    f_W1y    = (Float_t)evt->Get_MC_w1y(); // debug for testing
    
  }

  /** Destructor */
  ~TrueEvt() {};
  
  // interface: getter functions
  UInt_t   GetFlav()   const { return fFlav; }
  UInt_t   GetIsCC()   const { return (UInt_t)fIsCC; }
  UInt_t   GetIsNB()   const { return (UInt_t)fIsNB; }
  UInt_t   GetIchan()  const { return (UInt_t)fIchan; }
  Double_t GetTrueE()  const { return (Double_t)fE_true; }
  Double_t GetTrueCt() const { return (Double_t)fCt_true * 1e-4; }
  Double_t GetTrueBy() const { return (Double_t)fBy_true * 1e-4; }
  Double_t GetW1y()    const { return f_W1y; }

  // interface: setter functions
  void SetW1y(Float_t w1y) { f_W1y = w1y; }
  
 private:

  Short_t fFlav;    //!< neutrino flavor (0 - elec, 1 - muon, 2 - tau)
  Bool_t  fIsCC;    //!< interaction type (0 - NC, 1 - CC)
  Bool_t  fIsNB;    //!< neutrino polarisation (0 - nu, 1 - nubar)
  Short_t fIchan;   //!< variable that hold the info about the interaction channel (dis, res, qe)
  Float_t fE_true;  //!< variable that stores energy
  Short_t fCt_true; //!< variable that stores cos-theta with 4-digit resolution
  Short_t fBy_true; //!< variable that stores bjorken-y with 4-digit resolution
  Float_t f_W1y;    //!< event weight in 1 year
  
};

/** This class implements an event-by-event detector response.
    
    It allows the configuration of filters on `SummaryEvent`'s through inheriting functionality from `EventFilter`. It stores each event that passes the filters in a response matrix. The response matrix can be used later to calculate the weight (which depends on osc. parameters and systematics) for each event contributing to a specific reco bin to arrive at a number of expected reconstructed events in a certrain reco bin.

    Example usage:
    
   \code{.cpp}
   //---------------------------------------------------
   // setup the response
   //---------------------------------------------------
   EvtResponse dr(EvtResponse::track, "response_tracks"); // init response that uses track reco variables (see EventFilter::SetObservables)
   dr.AddCut( &SummaryEvent::Get_track_ql0       , std::greater<double>()   ,  0.5, true ); // select events with quality level 0
   dr.AddCut( &SummaryEvent::Get_RDF_track_score , std::greater<double>()   ,  0.6, true ); // select events with track score > 0.6
   dr.AddCut( &SummaryEvent::Get_RDF_muon_score  , std::less_equal<double>(), 0.05, true ); // select events with atm. muon score < 0.05

   //---------------------------------------------------
   // fill the response
   //---------------------------------------------------
   SummaryParser sp( getenv("MONADIR") + "/data/ORCA_MC_summary_all_10Apr2018.root");       // init summary parser
   for (Int_t i = 0; i < sp.GetTree()->GetEntries(); i++) {                                 // loop over events and fill
       dr.Fill( sp.GetEvt(i) );
   }

   //---------------------------------------------------
   // get the true events that contribute to one reco bin with E = 10, ct = -0.8, by = 0.2
   // then the total number of expected events in that reco bin can be calculated
   //---------------------------------------------------
   auto true_evts = dr.GetBinEvts(10, -0.8, 0.2);

   Double_t n_reco = 0;
   for (auto e: true_evts) {
       n_reco += operation_time * e.GetW1y() * GetOscillatedFlux( (TrueEvt)e );
   }

   \endcode

   Note that `GetOscillatedFlux` is pseudo-code, it is a function that returns the oscillated atmospheric neutrino flux at the detector cite for the neutrino type specified in the input argument `TrueEvt`. The usage of the class `EvtResponse` is illustrated in `fitter_software/FitUtil::RecoEvts`.

*/
class EvtResponse : public AbsResponse {

 public:

  EvtResponse(reco reco_type, TString resp_name,
	      Int_t ebins , Double_t emin , Double_t emax ,
  	      Int_t ctbins, Double_t ctmin, Double_t ctmax,
  	      Int_t bybins, Double_t bymin, Double_t bymax,
	      Double_t memlim = 4);
  
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
  std::vector<TrueEvt>&         GetBinEvts     (Double_t E_reco, Double_t ct_reco, Double_t by_reco);
  std::pair<Double_t, Double_t> GetAtmMuCount1y(Double_t E_reco, Double_t ct_reco, Double_t by_reco);
  std::pair<Double_t, Double_t> GetNoiseCount1y(Double_t E_reco, Double_t ct_reco, Double_t by_reco);
  
  /** Returns that this implementation is of the response type `EvtResponse` 
      \return `AbsResponse::EvtResponse`
  */
  resp                GetResponseType() { return AbsResponse::EvtResponse; }
  
  /** Get pointer to the histogram with atmospheric muon counts in 1y
      \return Pointer to a `TH3` with all the atm muon events expected in one year.
   */
  TH3D*               GetHistAtmMu1y() const { return fhAtmMuCount1y; }

  /** Get pointer to the histogram with noise counts in 1y
      \return Pointer to to `TH3` with all the noise events expected in one year.
   */
  TH3D*               GetHistNoise1y() const { return fhNoiseCount1y; }
  
 private:

  TH3D* CloneFromTemplate(TH3D* tmpl, TString name);
  void  CountEvents(UInt_t flav, UInt_t iscc, UInt_t isnb, SummaryEvent *evt);
  void  Normalise();
  
  typedef std::pair<Double_t, Double_t>   rangeID; //!< type definition for identifying a range
  typedef std::pair<Double_t, Double_t>   runID;   //!< type definition for identifying a run
  std::map< rangeID, std::vector<runID> > fRuns[NOISE+1][2][2]; //!< structure to count the total number of gSeaGen events (nu's) or total livetime for normalisation

  Bool_t   fNormalised;           //!< flag to indicate that the weights have been normalised to weight one year
  Double_t fNEvts;                //!< calculates the total number of events in a response
  Double_t fMemLim;               //!< user-defined limit to how much RAM the response can eat up
  Int_t    fEbins;                //!< number of reco energy bins in `fResp`
  Int_t    fCtbins;               //!< number of reco cos-theta bins in `fResp`
  Int_t    fBybins;               //!< number of reco bjorken-y bins in `fResp`
  TH3D    *fhAtmMuCount1y;        //!< atmospheric muon count in 1 year
  TH3D    *fhNoiseCount1y;        //!< noise event count in 1 year
  std::vector<TrueEvt> ***fResp;  //!< Response structure in 3D in [Ereco][CTreco][BYreco] = vector<TrueEvt> {contributing true events in reco bin}
  
};

#endif
