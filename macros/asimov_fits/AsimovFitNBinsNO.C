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

void AsimovFitNBinsNO() {

  const int N_PID_CLASSES = 10;
  const Double_t PID_STEP = 1 / float(N_PID_CLASSES);
  const Double_t PID_EDGE = 0.6 * N_PID_CLASSES;
  TString filefolder = TString::Format("./pid_detres/pid_binning_%i/", N_PID_CLASSES);

  // DetRes input values
  Int_t EBins = 40;
  Int_t EMin = 1;
  Int_t EMax = 100;
  Int_t ctBins = 40;
  Int_t ctMin = -1;
  Int_t ctMax = 1;
  Int_t byBins = 1;
  Int_t byMin = 0;
  Int_t byMax = 1;

  
  //----------------------------------------------------------
  // detector response for tracks and showers
  //----------------------------------------------------------
  std::vector<DetResponse> track_response_vector;
  std::vector<DetResponse> shower_response_vector;

  for (Int_t i = 0; i < N_PID_CLASSES; i++) {
    std::function<bool(double, double)> comparison_operator;
    if (i == 0) { comparison_operator = std::greater_equal<double>(); // The first bin needs to include the lower limit.
    } else { comparison_operator = std::greater<double>(); }

    DetResponse track_response(DetResponse::track, Form("track_response_%.2f", PID_STEP * i), 
                               EBins, EMin, EMax, ctBins, ctMin, ctMax, byBins, byMin, byMax);
    track_response.AddCut( &SummaryEvent::Get_track_ql0       , std::greater<double>()   , 0.5             , true );
    track_response.AddCut( &SummaryEvent::Get_track_ql1       , std::greater<double>()   , 0.5             , true );
    track_response.AddCut( &SummaryEvent::Get_RDF_track_score , comparison_operator      , PID_STEP * i    , true );
    track_response.AddCut( &SummaryEvent::Get_RDF_track_score , std::less_equal<double>(), PID_STEP * (i+1), true );
    track_response.AddCut( &SummaryEvent::Get_RDF_muon_score  , std::less_equal<double>(), 0.05            , true );
    track_response.AddCut( &SummaryEvent::Get_RDF_noise_score , std::less_equal<double>(), 0.18            , true );
    track_response_vector.push_back(track_response);

    DetResponse shower_response(DetResponse::shower, Form("shower_response_%.2f", PID_STEP * i), 
                                EBins, EMin, EMax, ctBins, ctMin, ctMax, byBins, byMin, byMax);
    shower_response.AddCut( &SummaryEvent::Get_shower_ql0     , std::greater<double>()   , 0.5             , true );
    shower_response.AddCut( &SummaryEvent::Get_shower_ql1     , std::greater<double>()   , 0.5             , true );
    shower_response.AddCut( &SummaryEvent::Get_RDF_track_score, comparison_operator      , PID_STEP * i    , true );
    shower_response.AddCut( &SummaryEvent::Get_RDF_track_score, std::less_equal<double>(), PID_STEP * (i+1), true );
    shower_response.AddCut( &SummaryEvent::Get_RDF_muon_score , std::less_equal<double>(), 0.05            , true );
    shower_response.AddCut( &SummaryEvent::Get_RDF_noise_score, std::less_equal<double>(), 0.5             , true );
    shower_response_vector.push_back(shower_response);
  }

  //-----------------------------------------------------
  // fill the detector response and event selection
  //-----------------------------------------------------

  // Start at true and if any one of the files is missing, become false
  Bool_t files_exist = kTRUE; 
  for (Int_t i = 0; i < N_PID_CLASSES; i++) {
    TString track_file  = Form("track_response_%.2f.root" , PID_STEP * i);
    TString shower_file = Form("shower_response_%.2f.root", PID_STEP * i);
    Bool_t track_exists  = NMHUtils::FileExists(filefolder + track_file);
    Bool_t shower_exists = NMHUtils::FileExists(filefolder + shower_file);
    files_exist = ((files_exist and track_exists) and shower_exists);
  }

  if (not files_exist) {

    auto summary_file = (TString)getenv("NMHDIR") + "/data/ORCA_MC_summary_all_10Apr2018.root";
    SummaryParser sp(summary_file);
    for (Int_t i = 0; i < sp.GetTree()->GetEntries(); i++) {
      if (i % (Int_t)1e6 == 0) cout << "Event: " << i << endl;
      SummaryEvent *evt = sp.GetEvt(i);
      for (Int_t i = 0; i < N_PID_CLASSES; i++) {
        track_response_vector[i].Fill(evt);
        shower_response_vector[i].Fill(evt);
      }
    }

    for (Int_t i = 0; i < N_PID_CLASSES; i++) {
      track_response_vector[i].WriteToFile(  filefolder +  Form("track_response_%.2f.root" , PID_STEP * i) );
      shower_response_vector[i].WriteToFile( filefolder + Form("shower_response_%.2f.root", PID_STEP * i) );
    }
    cout << "NOTICE: Finished filling response" << endl;
  }
  else {
    cout << "NOTICE: Reading responses from disk" << endl;
    for (Int_t i = 0; i < N_PID_CLASSES; i++) {
      track_response_vector[i].ReadFromFile(  filefolder + Form("track_response_%.2f.root" , PID_STEP * i) );
      shower_response_vector[i].ReadFromFile( filefolder + Form("shower_response_%.2f.root", PID_STEP * i) );
    }
  }
  
  cout << "NOTICE: Finished filling response" << endl;

  //----------------------------------------------------------
  // set up the PDFs and static oscillation parameters
  //----------------------------------------------------------

  auto meff_file = (TString)getenv("NMHDIR") + "/data/eff_mass/EffMass_ORCA115_23x9m_ECAP0418.root";
  FitUtil *fitutil = new FitUtil(3, track_response_vector[0].GetHist3D(), 2, 80, -1, 0, 0, 1, meff_file);

  std::vector<TH3D*> track_vector_io;
  std::vector<TH3D*> shower_vector_io;

  std::vector<FitPDF> pdf_tracks_vector;
  std::vector<FitPDF> pdf_showers_vector;

  for (int i = 0; i < N_PID_CLASSES; i++) {
    FitPDF pdf_tracks(  Form("pdf_tracks_%.2f", i*PID_STEP),  "pdf_tracks",  fitutil, &track_response_vector[i]);
    FitPDF pdf_showers( Form("pdf_showers_%.2f", i*PID_STEP), "pdf_showers", fitutil, &shower_response_vector[i]);
    pdf_tracks_vector.push_back( pdf_tracks );
    pdf_showers_vector.push_back( pdf_showers );

    // Set IO values and make expectation histograms
    fitutil->SetIOlims();
    fitutil->SetIOcentvals();

    track_vector_io.push_back(  (TH3D*)pdf_tracks.GetExpValHist() );
    shower_vector_io.push_back( (TH3D*)pdf_showers.GetExpValHist() );
    TString track_name_io  = Form("tracks_expval_io_%.2f", i*PID_STEP);
    TString shower_name_io = Form("showers_expval_io_%.2f", i*PID_STEP);
    track_vector_io[i] ->SetName(track_name_io);
    shower_vector_io[i]->SetName(shower_name_io);

    // Set parameters to constant
    fitutil->GetVar("SinsqTh12")->setConstant(kTRUE);
    fitutil->GetVar("SinsqTh13")->setConstant(kTRUE);
    fitutil->GetVar("dcp")->setConstant(kTRUE);
    fitutil->GetVar("Dm21")->setConstant(kTRUE);
  }

  //----------------------------------------------------------
  // set up data for simultaneous fitting and fit
  //----------------------------------------------------------
  cout << "NOTICE Fitter started fitting" << endl;
  
  TStopwatch timer;

  // Fit under NO model, IO data
  fitutil->SetNOlims();
  fitutil->SetNOcentvals();

  std::map<string, TH1*> hist_map_io;
  RooCategory cats_io("categories_io", "data categories");
  for (Int_t i = 0; i < N_PID_CLASSES; i++) {
    if (i < PID_EDGE) {
      hist_map_io.insert( {(string)shower_vector_io[i]->GetName(), shower_vector_io[i] } );
      cats_io.defineType( shower_vector_io[i]->GetName() );
    }
    else {
      cats_io.defineType( track_vector_io[i]->GetName() );
      hist_map_io.insert( {(string)track_vector_io[i] ->GetName(),  track_vector_io[i] } );
    }
  }

  RooSimultaneous simPdf_io("simPdf_io", "simultaneous Pdf for IO", cats_io);
  for (Int_t i = 0; i < N_PID_CLASSES; i++) {
    if (i < PID_EDGE) {
      simPdf_io.addPdf(pdf_showers_vector[i], shower_vector_io[i]->GetName() );
    }
    else {
      simPdf_io.addPdf(pdf_tracks_vector[i],  track_vector_io[i]->GetName() );
    }
  }

  RooDataHist data_hists_io("data_hists_io", "track and shower data", fitutil->GetObs(), cats_io, hist_map_io);
  RooFitResult *fitres_io = simPdf_io.fitTo( data_hists_io, Save(kTRUE) );
  cout << "NOTICE Fitter finished fitting, time duration [s]: " << (Double_t)timer.RealTime() << endl;

  RooArgSet result_io ( fitres_io->floatParsFinal() );
  cout << "*********Fit result comparison****************************" << endl;
  cout << "dm31       fitted: " << ((RooRealVar*)result_io.find("Dm31"))->getVal() << endl;
  cout << "sinsq_th23 fitted: " << ((RooRealVar*)result_io.find("SinsqTh23"))->getVal() << endl;
  cout << "*********Fit result comparison****************************" << endl;
  //----------------------------------------------------------
  // set hierarchy to fitted values
  //----------------------------------------------------------

  Double_t dm31      = ((RooRealVar*)result_io.find("Dm31"))->getVal();
  Double_t sinSqTh23 = ((RooRealVar*)result_io.find("SinsqTh23"))->getVal();
  fitutil->GetVar("Dm31")->setVal( dm31 );
  fitutil->GetVar("SinsqTh23")->setVal( sinSqTh23 );

  std::vector<TH3D*> fitted_io;
  for (Int_t i = 0; i < N_PID_CLASSES; i++) {
    if (i < PID_EDGE) {
      TH3D *showers_fitted_io = (TH3D*)pdf_showers_vector[i].GetExpValHist();
      showers_fitted_io->SetName( Form("showers_fitted_io_%.2f", i*PID_STEP) );
      fitted_io.push_back( showers_fitted_io );
    }
    else {
      TH3D *tracks_fitted_io = (TH3D*)pdf_tracks_vector[i].GetExpValHist();
      tracks_fitted_io->SetName( Form("tracks_fitted_io_%.2f", i*PID_STEP) );
      fitted_io.push_back( tracks_fitted_io );
    }
  }

  std::vector<Double_t> chi2_io;
  for (Int_t i = 0; i < N_PID_CLASSES; i++) {
    if (i < PID_EDGE) chi2_io.push_back( HistoChi2Test( shower_vector_io[i], fitted_io[i], 2, 80, -1, 0) );
    else              chi2_io.push_back( HistoChi2Test( track_vector_io[i],  fitted_io[i], 2, 80, -1, 0) );
  }

  Double_t chi2_tot_io = 0;
  for (Int_t i = 0; i < N_PID_CLASSES; i++) {
    cout << "NMHUtils: Chi2 between events IO and events fitted on NO is: " << chi2_io[i] << endl;
    chi2_tot_io += chi2_io[i] * chi2_io[i];
  }
  cout << "Squared sum is : " << std::sqrt( chi2_tot_io ) << endl;

}
