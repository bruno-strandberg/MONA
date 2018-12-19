#ifndef EventSelection_h
#define EventSelection_h

#include "SummaryEvent.h"
#include "EventFilter.h"
#include "TH2.h"
#include "TH3.h"
#include "TTree.h"
#include "TFile.h"

/**
   This class inherits from `EventFilter` and can be used to define event selections when looping over `SummaryEvent`'s.


   At the init of the class the user selects the reconstruction type that is used for filling member histograms. Then, the function `EventFilter::AddCut` can be used to define selection cuts for the events of this selection. The method `EventSelection::Fill` stores a `SummaryEvent` for this selection if the event passes the selection cuts.

   The class always fills member histograms when `Fill` is called. TO-BE-IMPLEMENTED: Additionally, when `TreeMode` flag is specified, the code saves the reco variables to a TTree, such that `fitter_software` can be used to perform an un-binned fit to data.
 */
class EventSelection : public EventFilter {

 public:
  
  EventSelection(reco reco_type, TString selection_name="", Bool_t TreeMode = kFALSE,
		 Int_t ebins  = 40, Double_t emin  =  1., Double_t emax  = 100.,
		 Int_t ctbins = 40, Double_t ctmin = -1., Double_t ctmax = 1.,
		 Int_t bybins =  1, Double_t bymin =  0., Double_t bymax = 1.);
  EventSelection(const EventSelection &evsel);
  ~EventSelection();

  void Fill(SummaryEvent *evt, Double_t w=1.);
  void WriteToFile(TString fname);

  TString Get_SelName()      { return fSelName;      } //!< Get the selection name
  TH2D*   Get_h_E_costh()    { return fh_E_costh;    } //!< Get the pointer to E_costh histo
  TH3D*   Get_h_E_costh_by() { return fh_E_costh_by; } //!< Get the pointer to E_costh_by histo

 private:
  TString        fSelName;        //!< name of the selection
  TH2D          *fh_E_costh;      //!< Energy vs CosTheta histogram
  TH3D          *fh_E_costh_by;   //!< Energy vs CosTheta vs by histogram

};

#endif
