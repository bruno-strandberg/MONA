#ifndef DetResponse_h
#define DetResponse_h

#include "AbsResponse.h"
#include "SummaryEvent.h"

#include "TCanvas.h"
#include "TH3.h"
#include "TList.h"

#include <vector>
#include <map>
#include <tuple>

/**
   class to store E_true, ct_true, by_true bins that contribute to E_reco, ct_reco, by_reco bin for neutrino events.

   The function `DetResponse::GetBinWeights` returns a vector of `TrueB` objects that can be subsequently used to map 'detected' neutrino events from E_true, ct_true, by_true to E_reco, ct_reco, by_reco.
 */
struct TrueB : public TObject {
  
  UInt_t   fFlav;        //!< flavor index (0 - elec, 1 - muon, 2 - tau)
  UInt_t   fIsCC;        //!< is cc (0 or 1)
  UInt_t   fIsNB;        //!< is anti-neutrino (0 or 1)
  Int_t    fE_true_bin;  //!< true energy bin index
  Int_t    fCt_true_bin; //!< true costheta bin index
  Int_t    fBy_true_bin; //!< true bjorken-y bin index
  Double_t fW;           //!< fraction of events (weight) from true 'detected' bin that contribute to 'detected' reco bin
  Double_t fWE;          //!< MC statistical error of the weight


  /** Default constructor. */
 TrueB(): fFlav(0), fIsCC(0), fIsNB(0), fE_true_bin(0), fCt_true_bin(0), fBy_true_bin(0), fW(0), fWE(0) {};
  
  /** Constructor.
      \param flav   Flavor (0 - elec, 1 - muon, 2 - tau)
      \param iscc   Is CC (0 or 1)
      \param isnb   Is nubar (0 or 1)
      \param ebin   True energy bin
      \param ctbin  True costheta bin
      \param bybin  True bjorken-y bin
  */
  TrueB(UInt_t flav, UInt_t iscc, UInt_t isnb,
	Int_t ebin, Int_t ctbin, Int_t bybin) {

    if (flav > 2 || iscc > 1 || isnb > 1) {
      throw std::invalid_argument( "ERROR! TrueB::TrueB() unknown neutrino type with flavor, is_cc (flag), is_nb (flag) "
				   + std::to_string(flav) + " " + std::to_string(iscc) + " " + std::to_string(isnb) );
    }

    fFlav = flav;
    fIsCC = iscc;
    fIsNB = isnb;
    fE_true_bin  = ebin;
    fCt_true_bin = ctbin;
    fBy_true_bin = bybin;

    //counters are initialized to 1 and error to 0. These converted to weights/error in `DetResponse::Normalise`
    fW        = 1;
    fWE       = 0.;
    
  };
  
  /** Destructor. */
  ~TrueB() {};

  /** Function to increment member counters */
  void Increment() {
    fW++;
  };
  
  /** 
      Operator to find a bin in a vector by using std::find.

      \param rhs instance of TrueB that this object is compared to
      \return true if the bins and neutrino type match, false otherwise
  */
  bool operator == (const TrueB& rhs) const {
    return ( ( this->fFlav == rhs.fFlav ) && 
	     ( this->fIsCC == rhs.fIsCC ) && 
	     ( this->fIsNB == rhs.fIsNB ) &&
	     ( this->fE_true_bin  == rhs.fE_true_bin  ) &&
	     ( this->fCt_true_bin == rhs.fCt_true_bin ) &&
	     ( this->fBy_true_bin == rhs.fBy_true_bin ) );
  }

  /** 
      Not-equal operator
      \param rhs instance of TrueB that this object is compared to
      \return true if bins or neutrino type mismatch, false otherwise
   */
  bool operator != (const TrueB& rhs) const {
    return !(*this == rhs);
  }

  
  /**
     Stream operator for cout.
  */
  friend std::ostream &operator << ( std::ostream &output, const TrueB &tb ) { 
    output << "Flavor, is-cc, is-nb; true e, ct, by bin; weight, weight err: "
	   << tb.fFlav << ' ' << tb.fIsCC << ' ' << tb.fIsNB << ' '
	   << tb.fE_true_bin << ' ' << tb.fCt_true_bin << ' ' << tb.fBy_true_bin << ' '
	   << tb.fW << ' ' << tb.fWE;
    return output;
  }

  ClassDef(TrueB, 1)
};

//==========================================================================================

/**
   This class implements the detector response structure.

   The response to neutrino events is realised in the class member `fResp[E_reco_bin][ct_reco_bin][by_reco_bin] = vector<TrueB>{...}`. Each `TrueB` object in the vector stores a weight for each `E_true, ct_true, by_true` bin that contributes to the reconstruction bin. In this way, the expected number of 'detected' events in the reconstruction bin can be calculated as \f$ \rm Detected_{reco} = \sum\limits_{i=0}^{true\:bins} Detected_{true} * Weight_{true} \f$. The function `GetBinWeights` returns the vector of `TrueB` objects associated with the reconstruction bin to enable the calculation of the number of events in the reco bin.

   Additionally, the response provides counts for the expected number of atmoshperic muon and noise events in 1 year through the functions `GetAtmMuCount1y` and `GetNoiseCount1y`. Atmospheric muons and noise are treated separately from neutrino events because they are simulated differently and they do not depend on oscillation parameters. Hence, for them one can just use the Monte-Carlo event weights to predict the number of events in 1 year.

   This class implements a detector response for a selection of events. In some sense it is a counterpart of the `EventSelection` class. Like `EventSelection`, this class inherits from `EventFilter`. `EventFilter` allows the user to define a custom filter of `SummaryEvent`'s through the `EventFilter::AddCut` function. Envisaged usage is such that for each `EventSelection` object, one can define a `DetResponse` counterpart, using the same cuts and reconstruction variables. The `DetResponse` structure helps to map 'detected' events (as defined in `FluxChain.C`) from [E_true][ct_true][by_true] to [E_reco][ct_reco][by_reco].

   Example usage:
   ```
   //---------------------------------------------------
   // setup the response
   //---------------------------------------------------
   DetResponse dr(DetResponse::track, "response_tracks"); // init response that uses track reco variables (see EventFilter::SetObservables)
   dr.AddCut( &SummaryEvent::Get_track_ql1       , std::greater<double>()   ,  0.5, true ); // select events with quality level 1
   dr.AddCut( &SummaryEvent::Get_RDF_track_score , std::greater<double>()   ,  0.6, true ); // select events with track score > 0.6
   dr.AddCut( &SummaryEvent::Get_RDF_muon_score  , std::less_equal<double>(), 0.05, true ); // select events with atm. muon score < 0.05

   //---------------------------------------------------
   // fill the response
   //---------------------------------------------------
   SummaryParser sp( getenv("MONADIR") + "/data/ORCA_MC_summary_all_10Apr2018.root");       // init summary parser
   for (Int_t i = 0; i < sp.GetTree()->GetEntries(); i++) {                                 // loop over events and fill
       sp.GetTree()->GetEntry(i);
       SummaryEvent *evt = sp.GetEvt();
       dr.Fill(evt);
   }

   //---------------------------------------------------
   // get the true bins that contribute to one reco bin with E = 10, ct = -0.8, by = 0.2
   // then the total number of expected events in that reco bin can be calculated
   //---------------------------------------------------
   auto true_bins = dr.GetBinWeights(10, -0.8, 0.2);

   Double_t n_reco = 0;
   Double_t n_reco_err = 0;
   for (auto b: true_bins) {
   Double_t n_true = GetDetectedTrue(b.fFlav, b.fIsCC, b.fIsNB, b.fE_true_bin, b.fCt_true_bin, b.fBy_true_bin); 
   n_reco     += b.fW * n_true;
   n_reco_err += TMath::Power(b.fWE * n_true, 2);
   }

   //---------------------------------------------------
   // add the atmospheric muon count and calculate the overall stat err
   //---------------------------------------------------
   n_reco += dr.GetAtmMuCount1y(10, -0.8, 0.2).first;
   n_reco_err += dr.GetAtmMuCount1y(10, -0.8, 0.2).second;
   n_reco_err = TMath::Sqrt(n_reco_err);   

   ```

   Note that `GetDetectedTrue()` is pseudo-code - this function needs to be implemented in the code that uses the `DetResponse`. It returns the number of 'detected' events (interacted per Mton * effective mass) for a given (flavor, is_cc, is_nubar, true energy, cos-theta and bjorken-y) combination. For example, such info is stored in the output histograms of `FluxChain.C` in directory `detflux`.

   NB! The binning of the true detected events needs to match the binning used for the `DetResponse`. The detector response should always span the full width of the simulation, i.e., if events are simulated in the range from [1, 100] GeV, [-1, 1] cos theta and [0, 1] bjorken-y, then the detector response needs to cover the full range.

 */
class DetResponse : public AbsResponse {

 public:
  DetResponse(reco reco_type, TString resp_name="",
  	      Int_t ebins  = 40, Double_t emin  =  1., Double_t emax  = 100.,
  	      Int_t ctbins = 40, Double_t ctmin = -1., Double_t ctmax = 1.,
  	      Int_t bybins =  1, Double_t bymin =  0., Double_t bymax = 1.);
  DetResponse(TString name, const DetResponse &detresp);
  ~DetResponse();

  std::vector<TrueB>& GetBinWeights(Double_t E_reco, Double_t ct_reco, Double_t by_reco);
  std::vector<TrueB>& GetBinWeights(SummaryEvent *evt);
  std::pair<Double_t, Double_t> GetAtmMuCount1y(Double_t E_reco, Double_t ct_reco, Double_t by_reco);
  std::pair<Double_t, Double_t> GetNoiseCount1y(Double_t E_reco, Double_t ct_reco, Double_t by_reco);
  void                Fill(SummaryEvent *evt);
  void                WriteToFile(TString filename);
  void                ReadFromFile(TString filename);
  std::tuple<TCanvas*, TCanvas*, TCanvas*> DisplayResponse(Double_t e_reco, Double_t ct_reco, TString outname="");

  /** Get pointer to the 3D histogram with selected reco events
      \return Pointer to a `TH3` with all the neutrino events that passed the selection cuts
   */
  TH3D*               GetHist3D() { return fHResp; }

  /** Get pointer to the histogram with atmospheric muon counts in 1y
      \return Pointer to a `TH3` with all the atm muon events expected in one year.
   */
  TH3D*               GetHistAtmMu1y() { return fhAtmMuCount1y; }

  /** Get pointer to the histogram with noise counts in 1y
      \return Pointer to to `TH3` with all the noise events expected in one year.
   */
  TH3D*               GetHistNoise1y() { return fhNoiseCount1y; }

  /** Returns that this implementation is of the response type `BINNED` 
      \return `AbsResponse::BINNED`
   */
  resp                GetResponseType() { return BINNED; }
  
 private:

  TH3D* CloneFromTemplate(TH3D* tmpl, TString name);
  void  InitResponse(Int_t ebins, Int_t ctbins, Int_t bybins);
  void  CleanResponse();
  void  Normalise();
  void  FillNuEvents(UInt_t flav, UInt_t is_cc, UInt_t is_nb, SummaryEvent *evt);
  
  Bool_t  fNormalised; //!< boolean to check that normalisation is called
  
  Int_t fEbins;        //!< number of reco energy bins
  Int_t fCtbins;       //!< number of reco cos-theta bins
  Int_t fBybins;       //!< number of reco bjorken-y bins

  TH3D    *fhSim[3][2][2];      //!< total numbers of simulated events [flavor][nc/cc][nu/nub]
  TH3D    *fhAtmMuCount1y;      //!< atmospheric muon count in 1 year
  TH3D    *fhNoiseCount1y;      //!< noise event count in 1 year
  TH3D    *fHResp;              //!< 3D histogram to help with binning functionality; stores the events with reco observables that pass cuts
  std::vector<TrueB> ***fResp;  //!< Response structure in 3D in [Ereco][CTreco][BYreco] = vector<TrueB> {contributing true bins}

  TList fHeap;
  
};

#endif
