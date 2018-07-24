#ifndef Asymmetry_h
#define Asymmetry_h

#include <map>
#include <vector>

#include "TH2.h"

#include "EventSelection.h"
#include "SummaryParser.h"

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
  std::vector< std::vector< std::tuple< TH2D*, TH2D*, TH2D*, TH2D* > > > fFluxHists;
  /// weight hists with structure fSimHists[flavor], first in h_nu, second h_nub
  std::vector< std::pair<TH2D*,TH2D*> > fSimHists;
  /// event selections per each NH/IH pair
  std::vector< std::pair< vector<EventSelection*>, vector<EventSelection*> > > fEvtSels;
  
  void InitFluxHists( vector< std::pair<TString, TString> > &NH_IH_pairs );
  void InitVars(TH2D *h_template, TString summary_files);
  void FillSimHists();
  void InitEvtSels_1(Int_t NH_IH_pcount);
  void FillSelections();
  void GetAsymmetries(vector< std::pair<TString, TString> > NH_IH_pairs, TString idstr);

  void CleanUp();

};

#endif
