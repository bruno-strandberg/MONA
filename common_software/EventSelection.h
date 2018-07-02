#ifndef EventSelection_h
#define EventSelection_h

#include "SummaryEvent.h"
#include "EventFilter.h"
#include "TH2.h"
#include "TTree.h"
#include "TFile.h"

class EventSelection : public EventFilter {

 public:
  
  enum reco {mc_truth, track, shower}; //!< enum to determine which reco type is used in filling hists

  EventSelection(reco reco_type, TString selection_name="", TTree *t = NULL,
		 Int_t ebins  = 40, Double_t emin  =  1., Double_t emax  = 100.,
		 Int_t ctbins = 40, Double_t ctmin = -1., Double_t ctmax = 1.);
  EventSelection(const EventSelection &evsel);
  ~EventSelection();

  void Fill(SummaryEvent *evt, Double_t w=1.);
  void SetTreeMode(TTree *t);
  void WriteToFile(TString fname);

 private:
  reco           fRecoType;       //!< variable to determine which reco type is used
  TString        fSelName;        //!< name of the selection
  TFile         *fTreeFile;       //!< the temp file the tree will be associated with
  TTree         *fTree;           //!< tree with selected events
  TH2D          *fh_E_costh;      //!< Energy vs CosTheta histogram

};

#endif
