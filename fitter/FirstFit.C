#include "FitFunction.h"

#include "TROOT.h"
#include "TSystem.h"
#include "TF3.h"
#include "TDatime.h"

#include "SummaryParser.h"
#include "EventSelection.h"
#include "DetResponse.h"
#include "SummaryEvent.h"
#include "NMHUtils.h"

#include <iostream>
using namespace std;

void FirstFit(TString experiment_file = "../evt_sampler/output/Experiments/Experiment_oscpars1_sample_0_NH.root", 
	      TString summary_data = "../data/ORCA_MC_summary_all_10Apr2018.root") {

  TDatime a;
  a.Print();

  gROOT->ProcessLine(".L FitFunction.C+");
  gSystem->Load("$OSCPROBDIR/libOscProb.so");

  //---------------------------------------------------------------
  // initialise the detector response and the event selection
  //---------------------------------------------------------------
  DetResponse dr(DetResponse::track, "track_response", 40, 1, 100, 40, -1, 1, 1, 0, 1);
  EventSelection t(EventSelection::track, "track", NULL, 40, 1, 100, 40, -1, 1, 1, 0, 1);

  vector<EventFilter*> P = {&dr, &t};

  for (auto &p: P) {

    p->AddCut( &SummaryEvent::Get_track_ql0       , std::greater<double>()   ,  0.5, true );
    p->AddCut( &SummaryEvent::Get_track_ql1       , std::greater<double>()   ,  0.5, true );
    p->AddCut( &SummaryEvent::Get_RDF_track_score , std::greater<double>()   ,  0.6, true );
    p->AddCut( &SummaryEvent::Get_RDF_muon_score  , std::less_equal<double>(), 0.05, true );
    p->AddCut( &SummaryEvent::Get_RDF_noise_score , std::less_equal<double>(), 0.18, true );

  }

  //---------------------------------------------------------------
  // fill the response or read from file
  //---------------------------------------------------------------

  TString response_name = "track_response.root";

  if ( NMHUtils::FileExists(response_name) ) {

    cout << "NOTICE Reading response from file" << endl;
    dr.ReadFromFile(response_name);

  }
  else {

    cout << "NOTICE Filling response from summary data" << endl;

    SummaryParser all_sdata(summary_data);

    for (Int_t i = 0; i < all_sdata.GetTree()->GetEntries(); i++) {
    
      all_sdata.GetTree()->GetEntry(i);
      SummaryEvent *evt = all_sdata.GetEvt();

      // skip muons and noise for now
      if ( evt->Get_MC_is_neutrino() < 0.5 ) continue;

      dr.Fill(evt);

    }

    dr.WriteToFile(response_name);

  }

  cout << "NOTICE Response finished" << endl;

  //---------------------------------------------------------------
  // fill the event selection
  //---------------------------------------------------------------

  SummaryParser exp(experiment_file);

  for (Int_t i = 0; i < exp.GetTree()->GetEntries(); i++) {
    exp.GetTree()->GetEntry(i);
    t.Fill( exp.GetEvt() );
  }

  cout << "NOTICE Selection finished" << endl;

  //---------------------------------------------------------------
  // init the fit function and fit
  //---------------------------------------------------------------
  FitFunction fitf(&dr, 3, "../data/eff_mass/EffMhists_elec_CC.root", "../data/eff_mass/EffMhists_muon_CC.root",
		   "../data/eff_mass/EffMhists_tau_CC.root", "../data/eff_mass/EffMhists_elec_NC.root");

  // note that fit range is defined in the function declaration
  TF3 *func = new TF3("fitfunc", fitf, 2, 80, -1, 0, 0, 1, 6);
  func->SetParameters(0.297, 0.0215, 0.425, 1.38, 7.37e-5, 2.56e-3);
  func->SetParNames("sinsq12","sinsq13","sinsq23","dcp","dm21","dm31");

  func->SetParLimits(0, 0.25, 0.354);
  func->SetParLimits(1, 0.019, 0.0242);
  func->SetParLimits(2, 0.381, 0.636);
  func->SetParLimits(3, 1, 1.9);
  func->SetParLimits(4, 6.93E-5, 7.96E-5);
  func->SetParLimits(5, -2.7E-3, 2.7E-3);
  
  t.Get_h_E_costh_by()->Fit("fitfunc", "VR");

  TDatime b;
  b.Print();

}
