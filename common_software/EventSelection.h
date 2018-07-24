#ifndef EventSelection_h
#define EventSelection_h

#include "SummaryEvent.h"
#include "EventFilter.h"
#include "TH2.h"
#include "TTree.h"
#include "TFile.h"

/**
   This class inherits from `EventFilter` and can be used to define event selections when looping over `SummaryEvent`'s.


   At the init of the class the user selects the reconstruction type that is used for filling member histograms. Then, the function `EventFilter::AddCut` can be used to define selection cuts for the events of this selection. The method `EventSelection::Fill` stores a `SummaryEvent` for this selection if the event passes the selection cuts.

   The class always fills member histograms when `Fill` is called. Additionally, if a pointer to a `TTree` where the analyzed `SummaryEvent`'s are read from is provided, the class operates in 'Tree more' and copies over the events to the member `fTree`.
 */
class EventSelection : public EventFilter {

 public:
  
  /// enum to determine which reco type is used in filling hists
  enum reco {mc_truth, track, shower};

  EventSelection(reco reco_type, TString selection_name="", TTree *t = NULL,
		 Int_t ebins  = 40, Double_t emin  =  1., Double_t emax  = 100.,
		 Int_t ctbins = 40, Double_t ctmin = -1., Double_t ctmax = 1.);
  EventSelection(const EventSelection &evsel);
  ~EventSelection();

  void Fill(SummaryEvent *evt, Double_t w=1.);
  void SetTreeMode(TTree *t);
  void WriteToFile(TString fname);

  TString Get_SelName()   { return fSelName;   } //!< Get the selection name
  TH2D*   Get_h_E_costh() { return fh_E_costh; } //!< Get the pointer to E_costh histo

 private:
  reco           fRecoType;       //!< variable to determine which reco type is used
  TString        fSelName;        //!< name of the selection
  TFile         *fTreeFile;       //!< the temp file the tree will be associated with
  TTree         *fTree;           //!< tree with selected events
  TH2D          *fh_E_costh;      //!< Energy vs CosTheta histogram

};

#endif
