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
   class to store E_true, ct_true, by_true bins that contribute to E_reco, ct_reco, by_reco bin
 */
struct TrueB : public TObject {
  
  Int_t    fFlav;        //!< flavor index
  Int_t    fIsCC;        //!< is cc
  Int_t    fIsNB;        //!< is anti-neutrino
  Int_t    fE_true_bin;  //!< true energy bin index
  Int_t    fCt_true_bin; //!< true costheta bin index
  Int_t    fBy_true_bin; //!< true bjorken-y bin index
  Double_t fN;           //!< total number of MC events from this bin that passed the selection
  Double_t fFracTrue;    //!< fraction of events from true bin that contribute to RecoB
  Double_t fFracReco;    //!< relative signal contribution of true bin to RecoB

  /** Default constructor. */
  TrueB(): fFlav(0), fIsCC(0), fIsNB(0), fE_true_bin(0), fCt_true_bin(0), fBy_true_bin(0), fN(0), fFracTrue(0), fFracReco(0) {};
  
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

    //counters are initialized to 1
    fN        = 1;
    fFracTrue = 1;
    fFracReco = 1;
    
  };
  
  /** Destructor. */
  ~TrueB() {};

  /** Function to increment member counters */
  void Increment() {
    fN++;
    fFracTrue++;
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
    output << "Flavor, is-cc, is-nb; true e, ct, by bin; N, frac true, frac reco: "
	   << tb.fFlav << ' ' << tb.fIsCC << ' ' << tb.fIsNB << ' '
	   << tb.fE_true_bin << ' ' << tb.fCt_true_bin << ' ' << tb.fBy_true_bin << ' '
	   << tb.fN << ' ' << tb.fFracTrue << ' ' << tb.fFracReco;
    return output;
  }

  ClassDef(TrueB, 1)
};

//==========================================================================================

/**
   Class description here.
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
  void Fill(SummaryEvent *evt);
  void Normalise();
  void WriteToFile(TString filename);
  void ReadFromFile(TString filename);
  TCanvas* DisplayResponse(Double_t e_reco, Double_t ct_reco);
  TH3D* GetHist3D() { return fRespH; } //!< Return pointer to the 3D histogram with selected reco events
 private:

  TString fRespName; //!< name to identify the response
  
  Int_t fEbins;      //!< number of reco energy bins
  Int_t fCtbins;     //!< number of reco cos-theta bins
  Int_t fBybins;     //!< number of reco bjorken-y bins

  /// map to create histogram names, flavors
  std::map<Int_t, TString> fFlavs = { {0, "elec"}, {1, "muon"}, {2, "tau"}, {3, "atmmu"}, {4, "noise"} };
  /// map to create histogram names, interaction types
  std::map<Int_t, TString> fInts  = { {0, "nc"}, {1, "cc"}  };
  /// map to create histogram names, neutrino/anti-neutrino
  std::map<Int_t, TString> fPols  = { {0, "nu"}, {1, "nub"} };
  /// map from gSeaGen MC_type (Geant4 particle code) to flavor
  std::map<Int_t, Int_t>   fType_to_Flav = { {12, 0}, {14, 1}, {16, 2}, {13, 3}, {0, 4} };
  
  TH3D    *fhSim[5][2][2];      //!< total numbers of simulated events [flavor][nc/cc][nu/nub]
  TH3D    *fRespH;              //!< 3D histogram to help with binning functionality
  std::vector<TrueB> ***fResp;  //!< Response structure in 3D in [Ereco][CTreco][BYreco] = vector<TrueB> {contributing bins}

  TList fHeap;
  
};

#endif
