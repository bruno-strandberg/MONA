#ifndef Asymmetry_h
#define Asymmetry_h

#include <map>
#include <vector>

#include "TH2.h"
#include "TH3.h"

#include "EventSelection.h"
#include "DetResponse.h"
#include "SummaryParser.h"

/**
   This namespace holds the variables for the 'Asymmetry.C' macro.
 */
namespace ASYM {

  SummaryParser *fSp; //!< summary parser instance

  /// enumerator for flavors
  enum Flavor {elec=0, muon, tau, nc};
  /// map for flavor identification
  std::map<Int_t, Int_t>   fFlav_pdg = { {elec, 12}, {muon, 14}, {tau, 16}, {nc, 12} };
  /// map for flavor strings
  std::map<Int_t, TString> fFlav_str = { {elec, "elec"}, {muon, "muon"}, {tau, "tau"}, {nc, "nc"} };
  /// enumerator for data storage in fFluxHists
  enum WHs {NH_nu = 0, NH_nub, IH_nu, IH_nub};

  /// flux hists with structure fFluxHists[NH_IH_pair_index][flavor]<NH_nu><NH_nub><IH_nu><IH_nub>
  std::vector< std::vector< std::tuple< TH3D*, TH3D*, TH3D*, TH3D* > > > fFluxHists;
  /// weight hists with structure fSimHists[NH_IH_pair_index][flavor], first is h_nu, second h_nub
  std::vector< std::vector< std::pair<TH3D*,TH3D*> > > fSimHists;
  /// event selections with structure fEvtSels[NH_IH_pair_index].first[selection_index], first is NH, second is IH
  std::vector< std::pair< vector<EventSelection*>, vector<EventSelection*> > > fEvtSels;
  /// detector responses with structure fResponse[NH_IH_pair_index][selection_index]
  std::vector< std::vector<DetResponse*> > fResponse;
  /// detected events histograms with structure fDetHists[NH_IH_pair_index].first[selection_index], first is NH, second is IH. These histograms serve to cross-check that the detector response and the simpler weighting method of this macro yield the same results
  std::vector< std::pair< vector<TH3D*>, vector<TH3D*> > > fDetHists;
  
  void InitFluxHists( vector< std::pair<TString, TString> > &NH_IH_pairs );
  void InitVars();
  void FillSimHists();
  void InitEvtSels_1();
  void FillSelections();
  void FillFromResponse();
  void GetAsymmetries(vector< std::pair<TString, TString> > NH_IH_pairs, TString idstr);

  void CleanUp();

};

#endif
