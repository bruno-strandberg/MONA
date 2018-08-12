#include "EventSelection.h"
#include "NMHUtils.h"
#include "TRandom3.h"
#include "TVector3.h"
#include<vector>

using namespace std;

/** 
    Constructor of an EventSelection.
    
    \param reco_type       Reco type used for filling histograms, e.g. `EventSelection::track`
    \param selection_name  String identifier for this selection
    \param t               Pointer to a `TTree/TChain` where events are selected from; if provided the selection is in "Tree mode" and events that pass the cut are cloned to the `fTree` of this class; if `t = NULL` the instance is in "Histogram mode" and events are not copied, only histograms of this instance are filled.
    \param ebins           Number of energy bins for histogram(s)
    \param emin            Minimum energy for histogram(s)
    \param emax            Maximum energy for histogram(s)
    \param ctbins          Number of \f$cos\theta\f$ bins for histogram(s)
    \param ctmin           Minimum \f$cos\theta\f$ for histogram(s)
    \param ctmax           Maximum \f$cos\theta\f$ for histogram(s)
 */
EventSelection::EventSelection(reco reco_type, TString selection_name, TTree *t,
			       Int_t ebins, Double_t emin, Double_t emax,
			       Int_t ctbins, Double_t ctmin, Double_t ctmax) : EventFilter(reco_type) {

  fSelName  = selection_name;

  if (t != NULL) SetTreeMode(t);
  else { fTree = NULL; fTreeFile = NULL; }

  std::vector<Double_t> e_edges = NMHUtils::GetLogBins(ebins, emin, emax);

  TString hname = fSelName + "_E_vs_costh";
  fh_E_costh = new TH2D(hname, hname, ebins, &e_edges[0], ctbins, ctmin, ctmax); 
  fh_E_costh->SetDirectory(0);

}

//***********************************************************************************

/**
   Copy constructor.

   The cuts are copied by the default copy constructor of the `EventFilter`, here the ROOT members (histograms, TTree) are taken care of.

   \param evsel   Another instance of EventSelection
 */
EventSelection::EventSelection(const EventSelection &evsel) : EventFilter(evsel) {

  fSelName  = evsel.fSelName;

  if (evsel.fTree != NULL) { SetTreeMode(evsel.fTree); }
  else { fTree = NULL; fTreeFile = NULL; }

  fh_E_costh = (TH2D*)evsel.fh_E_costh->Clone();
  fh_E_costh->SetDirectory(0);

}

//***********************************************************************************

/**
   Destructor.
 */
EventSelection::~EventSelection() {

  if (fh_E_costh) delete fh_E_costh;

  if (fTreeFile) {
    fTreeFile->Close();
    delete fTreeFile;
  }

}

//***********************************************************************************

/** 
    Fill function that fills a `SummaryEvent` to member histograms (and stores in `fTree` if in Tree mode) if the event passes the cuts set for this event selection.

    \param evt  Pointer to the `SummaryEvent`.
    \param w    Event weight when filled to histograms.

 */
void EventSelection::Fill(SummaryEvent *evt, Double_t w) {

  if ( !PassesCuts(evt) ) return;

  SetObservables(evt);

  // fill all member histograms
  fh_E_costh->Fill ( fEnergy , -fDir.z(), w );

  // if tree is set, fill this event to the tree
  if (fTree != NULL) fTree->Fill();

}

//***********************************************************************************

/**
   Function to initialize 'Tree mode'.

   When in tree mode, the summary event that passes the cut in `EventSelection::Fill` is cloned to the TTree of this event selection.

   \param t    Pointer to a `TTree` with `SummaryEvent`s.

 */
void EventSelection::SetTreeMode(TTree *t) {

  // generate a random number string
  TRandom3 rand(0);
  TString fname = (TString)getenv("NMHDIR") + (TString)"/tmp_" + 
    (TString)to_string( (int)rand.Integer(10E6) ) + ".root";

  // open a TFile and clone an empty tree. The cloned tree remains attached to the input tree t,
  // such that if one does t->GetEntry(i); fTree->Fill() the entry i is copied to fTree.
  fTreeFile = new TFile(fname, "RECREATE");
  fTree = t->CloneTree(0);
  
}

//***********************************************************************************

/**
   Write the histograms and TTree (if in tree mode) to file.
   
   \param fname    Name of the root file data is written to.
 */
void EventSelection::WriteToFile(TString fname) {

  TString tmp_fname = fTreeFile->GetName();
  
  fTree->Write();
  fh_E_costh->Write();

  fTreeFile->Close();
  delete fTreeFile;
  
  system((TString)"mv " + tmp_fname + " " + fname);

}
