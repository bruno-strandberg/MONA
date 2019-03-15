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

void AsimovFitNBinsNO_RandomQ_CheckQDistribution() {

  Bool_t use_random_q = kTRUE;
  const int N_PID_CLASSES = 3;

  std::map<Int_t, Double_t> pid_map = SetPIDCase(N_PID_CLASSES);

  TString filefolder = TString::Format("./detector_responses/pid_binning_%i/cross_check_random_q/", N_PID_CLASSES);

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
  std::vector<DetResponse> track_response_vector;
  std::vector<DetResponse> shower_response_vector;

  for (Int_t i = 0; i < N_PID_CLASSES; i++) {
    std::function<bool(double, double)> comparison_operator;
    if (i == 0) { comparison_operator = std::greater_equal<double>(); // The first bin needs to include the lower limit.
    } else { comparison_operator = std::greater<double>(); }

    DetResponse track_response(DetResponse::track, Form("track_response_%.2f", pid_map[i]), 
                               EBins, EMin, EMax, ctBins, ctMin, ctMax, byBins, byMin, byMax);
    track_response.AddCut( &SummaryEvent::Get_track_ql0       , std::greater<double>()   , 0.5         , true );
    track_response.AddCut( &SummaryEvent::Get_track_ql1       , std::greater<double>()   , 0.5         , true );
    track_response.AddCut( &SummaryEvent::Get_RDF_track_score , comparison_operator      , pid_map[i]  , true );
    track_response.AddCut( &SummaryEvent::Get_RDF_track_score , std::less_equal<double>(), pid_map[i+1], true );
    track_response.AddCut( &SummaryEvent::Get_RDF_muon_score  , std::less_equal<double>(), 0.05        , true );
    track_response.AddCut( &SummaryEvent::Get_RDF_noise_score , std::less_equal<double>(), 0.18        , true );
    track_response_vector.push_back(track_response);

    DetResponse shower_response(DetResponse::shower, Form("shower_response_%.2f", pid_map[i]), 
                                EBins, EMin, EMax, ctBins, ctMin, ctMax, byBins, byMin, byMax);
    shower_response.AddCut( &SummaryEvent::Get_shower_ql0     , std::greater<double>()   , 0.5         , true );
    shower_response.AddCut( &SummaryEvent::Get_shower_ql1     , std::greater<double>()   , 0.5         , true );
    shower_response.AddCut( &SummaryEvent::Get_RDF_track_score, comparison_operator      , pid_map[i]  , true );
    shower_response.AddCut( &SummaryEvent::Get_RDF_track_score, std::less_equal<double>(), pid_map[i+1], true );
    shower_response.AddCut( &SummaryEvent::Get_RDF_muon_score , std::less_equal<double>(), 0.05        , true );
    shower_response.AddCut( &SummaryEvent::Get_RDF_noise_score, std::less_equal<double>(), 0.5         , true );
    shower_response_vector.push_back(shower_response);
  }

  //-----------------------------------------------------
  // fill the detector response and event selection
  //-----------------------------------------------------

  //if (not files_exist) {

  TFile* f_q_dist = TFile::Open("./detector_responses/pid_quality_distribution.root", "READ");
  TH1D* q_dist = (TH1D*)f_q_dist->Get("h_quality");

  std::vector<TH1D*> q_dist_vector;
  for (Int_t i = 0; i < N_PID_CLASSES; i++) {
    TH1D* q_dist = new TH1D( Form("h_q_%i", i), Form("h_q_%i", i), 100, 0, 1);
    q_dist_vector.push_back(q_dist);
  }

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

      if (N_PID_CLASSES == 3) {
        if (track_score < 0.4) {
          q_dist_vector[0]->Fill( evt->Get_RDF_track_score() );
        }
        else if ((track_score >= 0.4) and (track_score < 0.6)) {
          q_dist_vector[1]->Fill( evt->Get_RDF_track_score() );
        }
        else if (track_score > 0.6) {
          q_dist_vector[2]->Fill( evt->Get_RDF_track_score() );
        }
        else cout << "ERROR: Undefined situation for track score." << endl;
        return 0;
      }
      else {
        Int_t j = track_score * N_PID_CLASSES;
        if (j == N_PID_CLASSES) j--; // If score == 1, we go out of range for the vector.
        q_dist_vector[j]->Fill( evt->Get_RDF_track_score() );
      }

    }
  }

  cout << "NOTICE: Finished filling response" << endl;

  TFile fout(filefolder + "q_distribution_after_randomizing.root", "RECREATE");
  for (Int_t i = 0; i < N_PID_CLASSES; i++) {
    q_dist_vector[i]->Write();
  }
  fout.Write();
  fout.Close();

}
