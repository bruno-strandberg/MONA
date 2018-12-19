#include "TSystem.h"
#include "TROOT.h"
#include "TDatime.h"
#include "TFile.h"
#include "TStopwatch.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"

#include "DetResponse.h"
#include "FitFunction.h"
#include "NMHUtils.h"
#include "SummaryParser.h"
#include "SummaryEvent.h"

#include <iostream>
using namespace std;

void reconstruction_counts(TString sum_file="../data/ORCA_MC_summary_all_10Apr2018.root") { 
  SummaryParser sp(sum_file);
  
  bool plot = true;
  int n_bins = 40;
  std::vector<Double_t> e_edges  = NMHUtils::GetLogBins(n_bins, 1, 100);
  std::vector<Double_t> ct_edges = NMHUtils::GetBins(n_bins, -1, 1);

  // First tracks
  Double_t q;
  std::vector<Int_t> good_track_and_shower(10);
  std::vector<Int_t> good_track(10);
  std::vector<Int_t> good_shower(10);
  std::vector<Int_t> bad_track(10);
  std::vector<Int_t> bad_shower(10);
  std::vector<Int_t> bad_track_and_shower(10);
  std::vector<Int_t> rejected(10);
  std::vector<Int_t> total_events(10);
  for (Int_t i = 0; i < sp.GetTree()->GetEntries(); i++) {

    // Filter shyte
    bool good_tr = false;
    bool good_sh = false;
    sp.GetTree()->GetEntry(i);

    // We order the histograms by quality: 0.0-0.1 --> 0, 0.1-0.2 --> 1, etc.
    q = sp.GetEvt()->Get_RDF_track_score();
    int index = (int)(TMath::Floor(q * 10));
    if (index == 10) { index = 9; } // A perfect track will get an index 10, which does not exist

    total_events[index]++;
    if ((sp.GetEvt()->Get_RDF_muon_score() > 0.05) or (sp.GetEvt()->Get_RDF_noise_score() > 0.18)) { rejected[index]++; continue; } // Filters for the events
    if ((sp.GetEvt()->Get_track_ql0() > 0.5) and (sp.GetEvt()->Get_track_ql1() > 0.5)) { good_tr = true; }
    if ((sp.GetEvt()->Get_shower_ql0() > 0.5) and (sp.GetEvt()->Get_shower_ql1() > 0.5)) { good_sh = true; }
    if ((good_tr) and (good_sh)) { good_track_and_shower[index]++; }
    if (good_tr) { good_track[index]++; } 
    if (good_sh) { good_shower[index]++; }
    if (not good_tr) { bad_track[index]++; }
    if (not good_sh) { bad_shower[index]++; }
    if ((not good_tr) and (not good_sh)) { bad_track_and_shower[index]++; } 

    // Fill events
  }
  cout << "Detector response for tracks" << endl;
  for (Int_t i = 0; i < 10; i++) {
    cout << "Quality     : " << i << endl;
    cout << "Total events: " << total_events[i] << endl;
    cout << "Good events : " << good_track_and_shower[i] << endl;
    cout << "Good tracks : " << good_track[i] << endl;
    cout << "Good showers: " << good_shower[i] << endl;
    cout << "Bad tracks  : " << bad_track[i] << endl;
    cout << "Bad showers : " << bad_shower[i] << endl;
    cout << "Bad both    : " << bad_track_and_shower[i] << endl;
    cout << "Rejected evt: " << rejected[i] << endl;
  }

  // Then showers
  // I do not trust the clear method, so for now I just overwrite all elements with 0.
  for (Int_t i = 0; i < 10; i++) {
    good_track_and_shower[i] = 0;
    good_track[i] = 0;
    good_shower[i] = 0;
    bad_track[i] = 0;
    bad_shower[i] = 0;
    bad_track_and_shower[i] = 0;
    rejected[i] = 0;
    total_events[i] = 0;
  }

  for (Int_t i = 0; i < sp.GetTree()->GetEntries(); i++) {

    // Filter shyte
    bool good_tr = false;
    bool good_sh = false;
    sp.GetTree()->GetEntry(i);

    // We order the histograms by quality: 0.0-0.1 --> 0, 0.1-0.2 --> 1, etc.
    q = sp.GetEvt()->Get_RDF_track_score();
    int index = (int)(TMath::Floor(q * 10));
    if (index == 10) { index = 9; } // A perfect track will get an index 10, which does not exist

    total_events[index]++;
    if ((sp.GetEvt()->Get_RDF_muon_score() > 0.05) or (sp.GetEvt()->Get_RDF_noise_score() > 0.50)) { rejected[index]++; continue; } // Filters for the events
    if ((sp.GetEvt()->Get_track_ql0() > 0.5) and (sp.GetEvt()->Get_track_ql1() > 0.5)) { good_tr = true; }
    if ((sp.GetEvt()->Get_shower_ql0() > 0.5) and (sp.GetEvt()->Get_shower_ql1() > 0.5)) { good_sh = true; }
    if ((good_tr) and (good_sh)) { good_track_and_shower[index]++; }
    if (good_tr) { good_track[index]++; } 
    if (good_sh) { good_shower[index]++; }
    if (not good_tr) { bad_track[index]++; }
    if (not good_sh) { bad_shower[index]++; }
    if ((not good_tr) and (not good_sh)) { bad_track_and_shower[index]++; } 

    // Fill events
  }
  cout << "Detector response for showers" << endl;
  for (Int_t i = 0; i < 10; i++) {
    cout << "Quality     : " << i << endl;
    cout << "Total events: " << total_events[i] << endl;
    cout << "Good events : " << good_track_and_shower[i] << endl;
    cout << "Good tracks : " << good_track[i] << endl;
    cout << "Good showers: " << good_shower[i] << endl;
    cout << "Bad tracks  : " << bad_track[i] << endl;
    cout << "Bad showers : " << bad_shower[i] << endl;
    cout << "Bad both    : " << bad_track_and_shower[i] << endl;
    cout << "Rejected evt: " << rejected[i] << endl;
  }
}
