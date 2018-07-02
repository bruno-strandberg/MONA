#include "EventSelection.h"
#include "NMHUtils.h"
#include "TRandom3.h"
#include "TVector3.h"
#include<vector>

using namespace std;

EventSelection::EventSelection(reco reco_type, TString selection_name, TTree *t,
			       Int_t ebins, Double_t emin, Double_t emax,
			       Int_t ctbins, Double_t ctmin, Double_t ctmax) {

  fRecoType = reco_type;
  fSelName  = selection_name;

  if (t != NULL) SetTreeMode(t);
  else { fTree = NULL; fTreeFile = NULL; }

  std::vector<Double_t> e_edges = NMHUtils::GetLogBins(ebins, emin, emax);

  TString hname = fSelName + "_E_vs_costh";
  fh_E_costh = new TH2D(hname, hname, ebins, &e_edges[0], ctbins, ctmin, ctmax); 

}

//***********************************************************************************

EventSelection::~EventSelection() {

  if (fh_E_costh) delete fh_E_costh;

  if (fTreeFile) {
    fTreeFile->Close();
    delete fTreeFile;
  }

}

//***********************************************************************************

void EventSelection::Fill(SummaryEvent *evt, Double_t w) {

  if ( !PassesCuts(evt) ) return;

  // here the correct variables from a reco need to be selected for filling histograms
  Double_t energy = 0;
  TVector3 dir, pos;

  switch (fRecoType) {

  case mc_truth:
    energy    = evt->Get_MC_energy();
    //bjorken_y = evt->Get_MC_bjorkeny();
    dir       = evt->Get_MC_dir();
    pos       = evt->Get_MC_pos();
    break;
  case track:
    energy    = evt->Get_track_energy();
    //bjorken_y = evt->Get_track_bjorkeny();
    dir       = evt->Get_track_dir();
    pos       = evt->Get_track_pos();
    break;
  case shower:
    energy    = evt->Get_shower_energy();
    //bjorken_y = evt->Get_shower_bjorkeny();
    dir       = evt->Get_shower_dir();
    pos       = evt->Get_shower_pos();
    break;
  }

  // fill all member histograms
  fh_E_costh->Fill ( energy , -dir.z(), w );

  // if tree is set, fill this event to the tree
  if (fTree != NULL) fTree->Fill();

}

//***********************************************************************************

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

void EventSelection::WriteToFile(TString fname) {

  TString tmp_fname = fTreeFile->GetName();
  
  fTree->Write();
  fh_E_costh->Write();

  fTreeFile->Close();
  delete fTreeFile;
  
  system((TString)"mv " + tmp_fname + " " + fname);

}
