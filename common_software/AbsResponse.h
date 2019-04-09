#ifndef AbsResponse_h
#define AbsResponse_h

#include "EventFilter.h"
#include "TH3.h"
#include <map>

/** 
    This is an abstract base class for different detector responses implemented in MONA common software.

    It inherits from `EventFilter` and thus has the functionality for filtering events.

    It was created to allow for several detector responses to be implemented that have a common base class. The latter is essential for fitting purposes. The classes and functions that enable the use of `RooFit` for NMO analyses require well-defined function interfaces. Having a virtual base class for the responses means that the `RooFit` wrapper class `fitter_software/FitPDF.h/C` only needs to store a pointer to the abstract base class and does not need to be modified with the addition of responses. The use of the modified responses is handled by the user in `fitter_software/FitUtil::RecoEvts`, either by modification of the `FitUtil` class or by inheriting from it and overwriting some of the virtual functions in `FitUtil`.

    This class has two purely virtual functions `Fill(SummaryEvent *evt)` and `GetResponseType`. Both of these functions need to be overloaded when a class inherits from this class. The first function is a method to be written to fill the response and it takes `SummaryEvent`'s as input, which is the data format in `MONA`. The second function returns an identifier code of the new response class. For each new class the enum `resp` with the codes should be extended and the new class should return a code from that enum. This is required in fitter software. For example `FitUtil::RecoEvts` will receive a pointer to a response when it is called from inside of `FitPDF::evaluate()`. However, different responses will require a different calculation to be performed. The return code of `GetResponseType` can be used to identify the response and decide on further action.

*/
class AbsResponse : public EventFilter {

 public:

  //-------------------------------------------------------------------------------------
  // constructors & destructors
  //-------------------------------------------------------------------------------------
  
  AbsResponse(reco reco_type, TString resp_name,
	      Int_t ebins  = 40, Double_t emin  =  1., Double_t emax  = 100.,
  	      Int_t ctbins = 40, Double_t ctmin = -1., Double_t ctmax = 1.,
  	      Int_t bybins =  1, Double_t bymin =  0., Double_t bymax = 1.);

  AbsResponse(reco reco_type, TString resp_name, TH3* h_bins);

  AbsResponse(reco reco_type, TString resp_name,
	      Int_t t_ebins , Double_t t_emin , Double_t t_emax ,
  	      Int_t t_ctbins, Double_t t_ctmin, Double_t t_ctmax,
  	      Int_t t_bybins, Double_t t_bymin, Double_t t_bymax,
	      Int_t r_ebins , Double_t r_emin , Double_t r_emax ,
  	      Int_t r_ctbins, Double_t r_ctmin, Double_t r_ctmax,
  	      Int_t r_bybins, Double_t r_bymin, Double_t r_bymax );

  AbsResponse(reco reco_type, TString resp_name, TH3* h_bins_true, TH3* h_bins_reco);

  AbsResponse(TString name, const AbsResponse &other);
  
  virtual ~AbsResponse();

  //-------------------------------------------------------------------------------------
  // virtual methods to be over-loaded
  //-------------------------------------------------------------------------------------

  /// enumerator that defines various response types that are implemented in inheriting classes, which need to return a code from this enumerator by implementing the function `GetResponseType`. For each new response class added to common_software, an element should be added to this enumerator that could be returned by `GetResponseType` in the inheriting class. This is necessary for deciding what to do in `FitUtil::RecoEvts`, as different response types require different calculation implementation.
  enum resp { BinnedResponse, EvtResponse };
  
  /** Purely virtual method to fill the response that needs to be re-implemented in inheriting classes */
  virtual void Fill(SummaryEvent* evt) = 0;

  /** Purely virtual method to return the response type from the enumerator `AbsResponse::resp` in the inheriting class.
      This is required to be able to be to have a singular interface to RooFit functionality, but decide in the function `FitUtil::RecoEvts` how to proceed with the calculation.
   */
  virtual resp GetResponseType() = 0;

  //-------------------------------------------------------------------------------------
  // setters/getters
  //-------------------------------------------------------------------------------------
  
  /** Function that returns a pointer to the TH3 histogram with reco bin settings
      \return Pointer to a TH3 with reco binning setting
   */
  TH3D*   GetHist3DReco() { return fhBinsReco; }
  
  /** Function that returns a pointer to the TH3 histogram with true bin settings
      \return Pointer to a TH3 with true binning setting
   */  
  TH3D*   GetHist3DTrue() { return fhBinsTrue; }

  /** function that returns the response name
      \return Response name
   */
  TString GetRespName()   { return fRespName;  }

 protected:

  //-------------------------------------------------------------------------------------
  // protected members
  //-------------------------------------------------------------------------------------
  
  TH3D   *fhBinsTrue; //!< 3D histogram that defines the binning in true space
  TH3D   *fhBinsReco; //!< 3D histogram that defines the binning in reco space
  TString fRespName;  //!< response name
  
  // enumerator with recognised particle types
  enum supported_particles {ELEC=0, MUON, TAU, ATMMU, NOISE};

  /// map to create histogram names and flavors, for internal use in inheriting classes
  std::map<UInt_t, TString> fFlavs = { {ELEC, "elec"}, {MUON, "muon"}, {TAU, "tau"} };
  /// map to create histogram names, interaction types, for internal use in inheriting classes
  std::map<UInt_t, TString> fInts  = { {0, "nc"}, {1, "cc"}  };
  /// map to create histogram names, neutrino/anti-neutrino
  std::map<UInt_t, TString> fPols  = { {0, "nu"}, {1, "nub"} };
  /// map from gSeaGen MC_type (Geant4 particle code) to flavor map used internally. Exception if PGG code of the event not in this map
  std::map<UInt_t, UInt_t>   fType_to_Supported = { {12, ELEC}, {14, MUON}, {16, TAU}, {13, ATMMU}, {0, NOISE} };

 private:

  //-------------------------------------------------------------------------------------
  // private
  //-------------------------------------------------------------------------------------
  void InitHists(TString resp_name, TH3* h_bins_true, TH3* h_bins_reco);
  
};

#endif
