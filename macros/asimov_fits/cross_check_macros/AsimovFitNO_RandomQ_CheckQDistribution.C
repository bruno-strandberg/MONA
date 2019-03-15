#include "HelperFunctions.h"

#include "TSystem.h"
#include "TROOT.h"
#include "TDatime.h"
#include "TH3.h"
#include "TFile.h"
#include "TStopwatch.h"
#include "TString.h"

#include "DetResponse.h"
#include "EventSelection.h"
#include "NMHUtils.h"
#include "SummaryParser.h"
#include "SummaryEvent.h"

#include "FitUtil.h"
#include "FitPDF.h"

#include "RooDataHist.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"
#include "RooFitResult.h"
#include "RooMinimizer.h"
#include "RooRealVar.h"

#include <iostream>
using namespace std;
using namespace RooFit;

void AsimovFitNO_RandomQ_CheckQDistribution() {

  Bool_t use_random_q = kTRUE;
  TString filefolder = "./detector_responses/pid_binning_2/cross_check_random_q/";

  // DetRes input values
  Int_t EBins = 24;
  Int_t EMin = 1;
  Int_t EMax = 100;
  Int_t ctBins = 40;
  Int_t ctMin = -1;
  Int_t ctMax = 1;
  Int_t byBins = 1;
  Int_t byMin = 0;
  Int_t byMax = 1;

  // Fitter ranges
  Int_t fitEMin  = 2;
  Int_t fitEMax  = 80;
  Int_t fitctMin = -1;
  Int_t fitctMax = 0;

  //----------------------------------------------------------
  // detector response for tracks and showers
  //----------------------------------------------------------
  DetResponse track_response(DetResponse::track, "track_response", EBins, EMin, EMax, ctBins, ctMin, ctMax, byBins, byMin, byMax);
  track_response.AddCut( &SummaryEvent::Get_track_ql0       , std::greater<double>()   ,  0.5, true );
  track_response.AddCut( &SummaryEvent::Get_track_ql1       , std::greater<double>()   ,  0.5, true );
  track_response.AddCut( &SummaryEvent::Get_RDF_track_score , std::greater<double>()   ,  0.6, true );
  track_response.AddCut( &SummaryEvent::Get_RDF_muon_score  , std::less_equal<double>(), 0.05, true );
  track_response.AddCut( &SummaryEvent::Get_RDF_noise_score , std::less_equal<double>(), 0.18, true );

  DetResponse shower_response(DetResponse::shower, "shower_response", EBins, EMin, EMax, ctBins, ctMin, ctMax, byBins, byMin, byMax);
  shower_response.AddCut( &SummaryEvent::Get_shower_ql0     , std::greater<double>()   ,  0.5, true );
  shower_response.AddCut( &SummaryEvent::Get_shower_ql1     , std::greater<double>()   ,  0.5, true );
  shower_response.AddCut( &SummaryEvent::Get_RDF_track_score, std::less_equal<double>(),  0.6, true );
  shower_response.AddCut( &SummaryEvent::Get_RDF_muon_score , std::less_equal<double>(), 0.05, true );
  shower_response.AddCut( &SummaryEvent::Get_RDF_noise_score, std::less_equal<double>(),  0.5, true );

  //-----------------------------------------------------
  // fill the detector response and event selection
  //-----------------------------------------------------

  TString track_file = "track_response.root";
  TString shower_file = "shower_response.root";

  //if ( !NMHUtils::FileExists(filefolder + track_file) or !NMHUtils::FileExists(filefolder + shower_file) ) {

  TFile* f_q_dist = TFile::Open("./detector_responses/pid_quality_distribution.root", "READ");
  TH1D* q_dist = (TH1D*)f_q_dist->Get("h_quality");

  TH1D* sh_q_dist = new TH1D("h_q_sh", "h_q_sh", 100, 0, 1);
  TH1D* tr_q_dist = new TH1D("h_q_tr", "h_q_tr", 100, 0, 1);


  auto summary_file = (TString)getenv("NMHDIR") + "/data/ORCA_MC_summary_all_10Apr2018.root";
  SummaryParser sp(summary_file);
  for (Int_t i = 0; i < sp.GetTree()->GetEntries(); i++) {
    if (i % (Int_t)1e6 == 0) cout << "Event: " << i << endl;
    SummaryEvent *evt = sp.GetEvt(i);

    // Randomize the q value of the event.
    if (use_random_q) {
      // Get the track score and throw a random from the dist.
      // If the event is a track (shower) and the random falls
      // in the area of shower (track), reroll until it doesnt.
      Double_t track_score = evt->Get_RDF_track_score();
      Double_t ran = q_dist->GetRandom();
      if (track_score <= 0.6) {
        while (ran > 0.6) ran = q_dist->GetRandom();
      }
      else {
        while (ran <= 0.6) ran = q_dist->GetRandom();
      }

      evt->Set_RDF_track_score(ran);

      if (track_score <= 0.6) { 
        sh_q_dist->Fill( evt->Get_RDF_track_score() );
      }
      else {
        tr_q_dist->Fill( evt->Get_RDF_track_score() );
      }
    }

    track_response.Fill(evt);
    shower_response.Fill(evt);
  }

  cout << "NOTICE: Finished filling response" << endl;

  TFile fout(filefolder + "q_distribution_after_randomizing.root", "RECREATE");
  sh_q_dist->Write();
  tr_q_dist->Write();
  fout.Write();
  fout.Close();

}
