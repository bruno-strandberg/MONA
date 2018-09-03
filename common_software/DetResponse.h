#ifndef DetResponse_h
#define DetResponse_h

#include "EventFilter.h"
#include "SummaryEvent.h"

#include "TCanvas.h"
#include "TH3.h"
#include "TList.h"

#include <vector>
#include <map>

/**
   class to store E_true, ct_true, by_true bins that contribute to E_reco, ct_reco, by_reco bin.

   The function `DetResponse::GetBinWeights` returns a vector of `TrueB` objects that can be subsequently used to map 'detected' events from E_true, ct_true, by_true to E_reco, ct_reco, by_reco.
 */
struct TrueB : public TObject {
  
  Int_t    fFlav;        //!< flavor index
  Int_t    fIsCC;        //!< is cc
  Int_t    fIsNB;        //!< is anti-neutrino
  Int_t    fE_true_bin;  //!< true energy bin index
  Int_t    fCt_true_bin; //!< true costheta bin index
  Int_t    fBy_true_bin; //!< true bjorken-y bin index
  Double_t fW;           //!< fraction of events (weight) from true 'detected' bin that contribute to 'detected' reco bin
  Double_t fWE;          //!< MC statistical error of the weight
  Double_t fFracReco;    //!< relative signal contribution of true bin to RecoB, used for `DetResponse::DisplayRespnse`


  /** Default constructor. */
 TrueB(): fFlav(0), fIsCC(0), fIsNB(0), fE_true_bin(0), fCt_true_bin(0), fBy_true_bin(0), fW(0), fWE(0), fFracReco(0) {};
  
  /** Constructor.
      \param flav   Flavor (see DetResponse::fFlavs)
      \param iscc   Is CC?
      \param isnb   Is nubar?
      \param ebin   True energy bin
      \param ctbin  True costheta bin
      \param bybin  True bjorken-y bin
  */
  TrueB(Int_t flav, Int_t iscc, Int_t isnb,
	Int_t ebin, Int_t ctbin, Int_t bybin) {
    fFlav = flav;
    fIsCC = iscc;
    fIsNB = isnb;
    fE_true_bin  = ebin;
    fCt_true_bin = ctbin;
    fBy_true_bin = bybin;

    //counters are initialized to 1 and error to 0. These converted to weights/error in `DetResponse::Normalise`
    fW        = 1;
    fFracReco = 1;
    fWE       = 0.;
    
  };
  
  /** Destructor. */
  ~TrueB() {};

  /** Function to increment member counters */
  void Increment() {
    fW++;
    fFracReco++;
  };
  
  /** 
      Operator to find a bin in a vector by using std::find.

      \param rhs instance of TrueB that this object is compared to
      \return true if the bins match, false otherwise
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
     Stream operator for cout.
  */
  friend std::ostream &operator << ( std::ostream &output, const TrueB &tb ) { 
    output << "Flavor, is-cc, is-nb; true e, ct, by bin; weight, weight err, frac reco: "
	   << tb.fFlav << ' ' << tb.fIsCC << ' ' << tb.fIsNB << ' '
	   << tb.fE_true_bin << ' ' << tb.fCt_true_bin << ' ' << tb.fBy_true_bin << ' '
	   << tb.fW << ' ' << tb.fWE << ' ' << tb.fFracReco;
    return output;
  }

  ClassDef(TrueB, 1)
};

//==========================================================================================

/**
   This class implements the detector response structure.

   The response is realised in the class member `fResp[E_reco_bin][ct_reco_bin][by_reco_bin] = vector<TrueB>{...}`. Each `TrueB` object in the vector stores a weight for each `E_true, ct_true, by_true` bin that contributes to the reconstruction bin. In this way, the expected number of 'detected' events in the reconstruction bin can be calculated as \f$ \rm Detected_{reco} = \sum\limits_{i=0}^{true\:bins} Detected_{true} * Weight_{true} \f$.

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
   SummaryParser sp( getenv("NMHDIR") + "/data/ORCA_MC_summary_all_10Apr2018.root");        // init summary parser
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
   n_reco_err = TMath::Sqrt(n_reco_err);
   ```

   Note that `GetDetectedTrue()` is pseudo-code - this function needs to be implemented in the code that uses the `DetResponse`. It returns the number of 'detected' events (interacted per Mton * effective mass) for a given (flavor, is_cc, is_nubar, true energy, cos-theta and bjorken-y) combination. For example, such info is stored in the output histograms of `FluxChain.C` in directory `detflux`.

   NB! The binning of the true detected events needs to match the binning used for the `DetResponse`!

 */
class DetResponse : public EventFilter {

 public:
  DetResponse(reco reco_type, TString resp_name="",
  	      Int_t ebins  = 40, Double_t emin  =  1., Double_t emax  = 100.,
  	      Int_t ctbins = 40, Double_t ctmin = -1., Double_t ctmax = 1.,
  	      Int_t bybins =  1, Double_t bymin =  0., Double_t bymax = 1.);
  DetResponse(const DetResponse &detresp);
  ~DetResponse();

  std::vector<TrueB>& GetBinWeights(Double_t E_reco, Double_t ct_reco, Double_t by_reco);
  std::vector<TrueB>& GetBinWeights(SummaryEvent *evt);
  void                Fill(SummaryEvent *evt);
  void                WriteToFile(TString filename);
  void                ReadFromFile(TString filename);
  TCanvas*            DisplayResponse(Double_t e_reco, Double_t ct_reco);
  /// Return pointer to the 3D histogram with selected reco events
  TH3D*               GetHist3D() { return fHResp; }
  /// Get response name
  TString             Get_RespName() { return fRespName; }

 private:

  void InitResponse(Int_t ebins, Int_t ctbins, Int_t bybins);
  void CleanResponse();
  void Normalise();
  
  TString fRespName;   //!< name to identify the response
  Bool_t  fNormalised; //!< boolean to check that normalisation is called
  
  Int_t fEbins;        //!< number of reco energy bins
  Int_t fCtbins;       //!< number of reco cos-theta bins
  Int_t fBybins;       //!< number of reco bjorken-y bins

  /// map to create histogram names, flavors
  std::map<Int_t, TString> fFlavs = { {0, "elec"}, {1, "muon"}, {2, "tau"}, {3, "atmmu"}, {4, "noise"} };
  /// map to create histogram names, interaction types
  std::map<Int_t, TString> fInts  = { {0, "nc"}, {1, "cc"}  };
  /// map to create histogram names, neutrino/anti-neutrino
  std::map<Int_t, TString> fPols  = { {0, "nu"}, {1, "nub"} };
  /// map from gSeaGen MC_type (Geant4 particle code) to flavor
  std::map<Int_t, Int_t>   fType_to_Flav = { {12, 0}, {14, 1}, {16, 2}, {13, 3}, {0, 4} };
  
  TH3D    *fhSim[5][2][2];      //!< total numbers of simulated events [flavor][nc/cc][nu/nub]
  TH3D    *fHResp;              //!< 3D histogram to help with binning functionality; stores the events with reco observables that pass cuts
  std::vector<TrueB> ***fResp;  //!< Response structure in 3D in [Ereco][CTreco][BYreco] = vector<TrueB> {contributing true bins}

  TList fHeap;
  
};

#endif
