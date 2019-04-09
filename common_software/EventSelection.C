#include "EventSelection.h"
#include "NMHUtils.h"

#include "TRandom3.h"
#include "TVector3.h"

#include <vector>
#include <iostream>
#include <stdexcept>

using namespace std;

/** 
    Constructor of an EventSelection.
    
    \param reco_type       Reco type used for filling histograms, e.g. `EventSelection::track`
    \param selection_name  String identifier for this selection
    \param TreeMode        Save the events to a `TTree` to perform un-binned fits to data (currently not implemented)
    \param ebins           Number of energy bins for histogram(s)
    \param emin            Minimum energy for histogram(s)
    \param emax            Maximum energy for histogram(s)
    \param ctbins          Number of \f$cos\theta\f$ bins for histogram(s)
    \param ctmin           Minimum \f$cos\theta\f$ for histogram(s)
    \param ctmax           Maximum \f$cos\theta\f$ for histogram(s)
    \param bybins          Number of bjorken-y bins for histogram(s)
    \param bymin           Minimum bjorken-y for histogram(s)
    \param bymax           Maximum bjorken-y for histogram(s)
 */
EventSelection::EventSelection(reco reco_type, TString selection_name, Bool_t TreeMode,
			       Int_t ebins, Double_t emin, Double_t emax,
			       Int_t ctbins, Double_t ctmin, Double_t ctmax,
			       Int_t bybins, Double_t bymin, Double_t bymax) : EventFilter(reco_type) {

  if (TreeMode) {
    throw std::invalid_argument( "ERROR! EventSelection::EventSelection() TreeMode is currently not supported." );
  }
  
  fSelName  = selection_name;

  //----------------------------------------------
  //init 2D hist
  //----------------------------------------------

  std::vector<Double_t> e_edges = NMHUtils::GetLogBins(ebins, emin, emax);

  TString hname = fSelName + "_E_vs_costh";
  fh_E_costh = new TH2D(hname, hname, ebins, &e_edges[0], ctbins, ctmin, ctmax); 
  fh_E_costh->SetDirectory(0);

  //----------------------------------------------
  //for 3D hist constructor all axes need to be defined the same way
  //----------------------------------------------

  vector<Double_t> ct_edges = NMHUtils::GetBins(ctbins, ctmin, ctmax);
  vector<Double_t> by_edges = NMHUtils::GetBins(bybins, bymin, bymax);
 
  fh_E_costh_by = new TH3D(hname, hname, ebins, &e_edges[0], ctbins, &ct_edges[0], bybins, &by_edges[0]);
  fh_E_costh_by->SetDirectory(0);
}

//***********************************************************************************

/**
   Copy constructor.

   The cuts are copied by the default copy constructor of the `EventFilter`, here the ROOT members (histograms, TTree) are taken care of.

   \param evsel   Another instance of EventSelection
 */
EventSelection::EventSelection(const EventSelection &evsel) : EventFilter(evsel) {

  fSelName  = evsel.fSelName;

  fh_E_costh = (TH2D*)evsel.fh_E_costh->Clone();
  fh_E_costh->SetDirectory(0);

  fh_E_costh_by = (TH3D*)evsel.fh_E_costh_by->Clone();
  fh_E_costh_by->SetDirectory(0);
}

//***********************************************************************************

/**
   Destructor.
 */
EventSelection::~EventSelection() {

  if (fh_E_costh) delete fh_E_costh;
  if (fh_E_costh_by) delete fh_E_costh_by;

}

//***********************************************************************************

/** 
    Fill function that fills a `SummaryEvent` to member histograms if the event passes the cuts set for this event selection.

    \param evt  Pointer to the `SummaryEvent`.
    \param w    Event weight when filled to histograms.

 */
void EventSelection::Fill(SummaryEvent *evt, Double_t w) {

  if ( !PassesCuts(evt) ) return;

  SetObservables(evt); //implemented in EventFilter.C

  // fill all member histograms
  fh_E_costh->Fill( fEnergy , -fDir.z(), w );
  fh_E_costh_by->Fill( fEnergy, -fDir.z(), fBy, w );

  // if tree mode is implemented, a tree should be filled here with reco variable names that
  // can be recognized by `fitter_software/FitUtil` - E_reco, ct_reco, by_reco.
  
}

//***********************************************************************************

/**
   Write the histograms to file.
   
   \param fname    Name of the root file data is written to.
 */
void EventSelection::WriteToFile(TString fname) {

  TFile fout(fname, "RECREATE");
  fh_E_costh->Write();
  fh_E_costh_by->Write();
  fout.Close();
  
}
