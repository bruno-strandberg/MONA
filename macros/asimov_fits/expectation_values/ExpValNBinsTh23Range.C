#include "../HelperFunctions.h"

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

/* Script to calculate the asimov sensitivity at the PDG central values under the assumption
 * that Nature is NO. The script uses any bin number between 3 and 10 PID bins, where the bins 
 * are evenly spaced the following way in the shower and track channels:
 */

void ExpValNBinsTh23Range() {

  const int N_PID_CLASSES = 2;
  const Double_t PID_CUT = 0.6;

  std::map<Int_t, Double_t> pid_map = SetPIDCase(N_PID_CLASSES);

  TString filefolder = DetectorResponseFolder(N_PID_CLASSES);

  TString s_rootfile = Form("ExpVal%iBinsTh23Range.root", N_PID_CLASSES);

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
  std::vector<DetResponse*> track_response_vector;
  std::vector<DetResponse*> shower_response_vector;

  for (Int_t i = 0; i < N_PID_CLASSES; i++) {
    std::function<bool(double, double)> comparison_operator;
    if (i == 0) { comparison_operator = std::greater_equal<double>(); // The first bin needs to include the lower limit.
    } else { comparison_operator = std::greater<double>(); }

    DetResponse* track_response = new DetResponse(DetResponse::track, Form("track_response_%.2f", pid_map[i]), 
                               EBins, EMin, EMax, ctBins, ctMin, ctMax, byBins, byMin, byMax);
    track_response->AddCut( &SummaryEvent::Get_track_ql0       , std::greater<double>()   , 0.5         , true );
    track_response->AddCut( &SummaryEvent::Get_track_ql1       , std::greater<double>()   , 0.5         , true );
    track_response->AddCut( &SummaryEvent::Get_RDF_track_score , comparison_operator      , pid_map[i]  , true );
    track_response->AddCut( &SummaryEvent::Get_RDF_track_score , std::less_equal<double>(), pid_map[i+1], true );
    track_response->AddCut( &SummaryEvent::Get_RDF_muon_score  , std::less_equal<double>(), 0.05        , true );
    track_response->AddCut( &SummaryEvent::Get_RDF_noise_score , std::less_equal<double>(), 0.18        , true );
    track_response_vector.push_back(track_response);

    DetResponse* shower_response = new DetResponse(DetResponse::shower, Form("shower_response_%.2f", pid_map[i]), 
                                EBins, EMin, EMax, ctBins, ctMin, ctMax, byBins, byMin, byMax);
    shower_response->AddCut( &SummaryEvent::Get_shower_ql0     , std::greater<double>()   , 0.5         , true );
    shower_response->AddCut( &SummaryEvent::Get_shower_ql1     , std::greater<double>()   , 0.5         , true );
    shower_response->AddCut( &SummaryEvent::Get_RDF_track_score, comparison_operator      , pid_map[i]  , true );
    shower_response->AddCut( &SummaryEvent::Get_RDF_track_score, std::less_equal<double>(), pid_map[i+1], true );
    shower_response->AddCut( &SummaryEvent::Get_RDF_muon_score , std::less_equal<double>(), 0.05        , true );
    shower_response->AddCut( &SummaryEvent::Get_RDF_noise_score, std::less_equal<double>(), 0.5         , true );
    shower_response_vector.push_back(shower_response);
  }

  //-----------------------------------------------------
  // fill the detector response and event selection
  //-----------------------------------------------------

  auto summary_file = (TString)getenv("MONADIR") + "/data/ORCA_MCsummary_SEv2_ORCA115_23x9m_ECAP180401.root";
  SummaryParser sp(summary_file);
  for (Int_t i = 0; i < sp.GetTree()->GetEntries(); i++) {
    if (i % (Int_t)1e6 == 0) cout << "Event: " << i << endl;
    SummaryEvent *evt = sp.GetEvt(i);
    for (Int_t i = 0; i < N_PID_CLASSES; i++) {
      track_response_vector[i]->Fill(evt);
      shower_response_vector[i]->Fill(evt);
    }
  }
  
  
  cout << "NOTICE: Finished filling response" << endl;

  //----------------------------------------------------------
  // Create and save the expectation value histograms
  //----------------------------------------------------------

  auto meff_file = (TString)getenv("MONADIR") + "/data/eff_mass/EffMass_ORCA115_23x9m_ECAP180401.root";

  // Open root file to save histograms
  TFile* fout = new TFile(s_rootfile, "RECREATE");
  TDirectory* bins = fout->mkdir("test");


  for (Int_t i = 0; i < 11; i++) {
    FitUtil *fitutil = new FitUtil(3, track_response_vector[0]->GetHist3D(), fitEMin, fitEMax, fitctMin, fitctMax, 0, 1, meff_file);

    Int_t th23 = 40 + i;
    Double_t sinSqTh23_true = TMath::Power(TMath::Sin(th23 * TMath::Pi()/180.), 2);

    for (int i = 0; i < N_PID_CLASSES; i++) {
      FitPDF pdf_tracks(  Form("pdf_tracks_%.2f", pid_map[i]),  "pdf_tracks",  fitutil, track_response_vector[i]);
      FitPDF pdf_showers( Form("pdf_showers_%.2f", pid_map[i]), "pdf_showers", fitutil, shower_response_vector[i]);

      // Set NO values and make expectation histograms
      fitutil->SetNOcentvals();
      fitutil->GetVar("SinsqTh23")->setVal( sinSqTh23_true );

      auto track_no_3 = (TH3D*)pdf_tracks.GetExpValHist();
      auto shower_no_3 = (TH3D*)pdf_showers.GetExpValHist();
      auto track_no = track_no_3->Project3D("yx");
      auto shower_no = shower_no_3->Project3D("yx");
      TString track_name_no  = Form("tracks_expval_no_%i_%i_of_%i", th23, i+1, N_PID_CLASSES); // Start at 1 for a change, see if this works out...
      TString shower_name_no = Form("showers_expval_no_%i_%i_of_%i", th23, i+1, N_PID_CLASSES);
      track_no->SetName(track_name_no);
      shower_no->SetName(shower_name_no);

      // Set IO values and make expectation histograms
      fitutil->SetIOcentvals();
      fitutil->GetVar("SinsqTh23")->setVal( sinSqTh23_true );

      auto track_io_3 = (TH3D*)pdf_tracks.GetExpValHist();
      auto shower_io_3 = (TH3D*)pdf_showers.GetExpValHist();
      auto track_io = track_io_3->Project3D("yx");
      auto shower_io = shower_io_3->Project3D("yx");
      TString track_name_io  = Form("tracks_expval_io_%i_%i_of_%i", th23, i+1, N_PID_CLASSES); 
      TString shower_name_io = Form("showers_expval_io_%i_%i_of_%i", th23, i+1, N_PID_CLASSES);
      track_io->SetName(track_name_io);
      shower_io->SetName(shower_name_io);

      // Asymmetry calculation
      TString asym_track_name  = Form("asym_track_%i_%i_of_%i", th23, i+1, N_PID_CLASSES);
      auto asym_track = NMHUtils::Asymmetry(track_no, track_io, asym_track_name, 
            fitEMin, fitEMax, fitctMin, fitctMax, 0, 1);
      auto h_asym_track = std::get<0>(asym_track);

      TString asym_shower_name  = Form("asym_shower_%i_%i_of_%i", th23, i+1, N_PID_CLASSES);
      auto asym_shower = NMHUtils::Asymmetry(shower_no, shower_io, asym_shower_name, 
            fitEMin, fitEMax, fitctMin, fitctMax, 0, 1);
      auto h_asym_shower = std::get<0>(asym_shower);

      if (pid_map[i] < PID_CUT) {
        bins->cd();
        shower_no->Write();
        shower_io->Write();
        h_asym_shower->Write();
      }
      else {
        bins->cd();
        track_no->Write();
        track_io->Write();
        h_asym_track->Write();
      }
    }
  }
  fout->Close();
}
