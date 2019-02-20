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

void AsimovFit3BinsIO() {

  const int N_PID_CLASSES = 3;
  const Double_t PID_CUT = 0.6;
  const Double_t PID_STEP = 1 / float(N_PID_CLASSES);
  const Double_t PID_EDGE = PID_CUT * N_PID_CLASSES;


  std::map<Int_t, Double_t> pid_map;
  pid_map.insert(std::make_pair(0, 0.0)); // shower
  pid_map.insert(std::make_pair(1, 0.4)); // middle group: shower
  pid_map.insert(std::make_pair(2, 0.6)); // track
  pid_map.insert(std::make_pair(3, 1.0)); // upper limit 

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

    DetResponse track_response(DetResponse::track, Form("track_response_%.2f", PID_STEP * i), 
                               EBins, EMin, EMax, ctBins, ctMin, ctMax, byBins, byMin, byMax);
    track_response.AddCut( &SummaryEvent::Get_track_ql0       , std::greater<double>()   , 0.5         , true );
    track_response.AddCut( &SummaryEvent::Get_track_ql1       , std::greater<double>()   , 0.5         , true );
    track_response.AddCut( &SummaryEvent::Get_RDF_track_score , comparison_operator      , pid_map[i]  , true );
    track_response.AddCut( &SummaryEvent::Get_RDF_track_score , std::less_equal<double>(), pid_map[i+1], true );
    track_response.AddCut( &SummaryEvent::Get_RDF_muon_score  , std::less_equal<double>(), 0.05        , true );
    track_response.AddCut( &SummaryEvent::Get_RDF_noise_score , std::less_equal<double>(), 0.18        , true );
    track_response_vector.push_back(track_response);

    DetResponse shower_response(DetResponse::shower, Form("shower_response_%.2f", PID_STEP * i), 
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
  FitUtil *fitutil = new FitUtil(3, track_response_vector[0].GetHist3D(), fitEMin, fitEMax, fitctMin, fitctMax, 0, 1, meff_file);

  std::vector<TH3D*> track_vector_true;
  std::vector<TH3D*> shower_vector_true;

  std::vector<FitPDF> pdf_tracks_vector;
  std::vector<FitPDF> pdf_showers_vector;

  for (int i = 0; i < N_PID_CLASSES; i++) {
    FitPDF pdf_tracks(  Form("pdf_tracks_%.2f", pid_map[i]),  "pdf_tracks",  fitutil, &track_response_vector[i]);
    FitPDF pdf_showers( Form("pdf_showers_%.2f", pid_map[i]), "pdf_showers", fitutil, &shower_response_vector[i]);
    pdf_tracks_vector.push_back( pdf_tracks );
    pdf_showers_vector.push_back( pdf_showers );

    // Set NO values and make expectation histograms
    fitutil->SetNOcentvals();

    track_vector_true.push_back(  (TH3D*)pdf_tracks.GetExpValHist() );
    shower_vector_true.push_back( (TH3D*)pdf_showers.GetExpValHist() );
    TString track_name_true  = Form("tracks_expval_true_%.2f", pid_map[i]);
    TString shower_name_true = Form("showers_expval_true_%.2f", pid_map[i]);
    track_vector_true[i] ->SetName(track_name_true);
    shower_vector_true[i]->SetName(shower_name_true);

    // Set parameters to constant
    fitutil->GetVar("SinsqTh12")->setConstant(kTRUE);
    fitutil->GetVar("SinsqTh13")->setConstant(kTRUE);
    fitutil->GetVar("dcp")->setConstant(kTRUE);
    fitutil->GetVar("Dm21")->setConstant(kTRUE);
  }

  // create quantiles for theta23
  fitutil->GetVar("SinsqTh23")->setRange("firstq" , 0., 0.5);
  fitutil->GetVar("SinsqTh23")->setRange("secondq", 0.5, 1.);
  
  //----------------------------------------------------------
  // set up data for simultaneous fitting and fit
  //----------------------------------------------------------
  cout << "NOTICE Fitter started fitting" << endl;
  
  TStopwatch timer;

  // Fit under IO model, NO data
  SetIOlimsChi2Fit(fitutil);
  fitutil->SetIOcentvals();

  std::map<string, TH1*> hist_map;
  RooCategory cats("categories", "data categories");
  for (Int_t i = 0; i < N_PID_CLASSES; i++) {
    if (pid_map[i] < PID_CUT) {
      hist_map.insert( {(string)shower_vector_true[i]->GetName(), shower_vector_true[i] } );
      cats.defineType( shower_vector_true[i]->GetName() );
      cout << "NOTICE: Added hist and cat to shower" << endl;
    }
    else {
      hist_map.insert( {(string)track_vector_true[i]->GetName(), track_vector_true[i] } );
      cats.defineType( track_vector_true[i]->GetName() );
      cout << "NOTICE: Added hist and cat to track" << endl;
    }
  }

  RooSimultaneous simPdf("simPdf", "simultaneous Pdf for NO", cats);
  for (Int_t i = 0; i < N_PID_CLASSES; i++) {
    if (pid_map[i] < PID_CUT) {
      simPdf.addPdf(pdf_showers_vector[i], shower_vector_true[i]->GetName() );
      cout << "NOTICE: Added simpdf to shower" << endl;
    }
    else {
      simPdf.addPdf(pdf_tracks_vector[i],  track_vector_true[i]->GetName() );
      cout << "NOTICE: Added simpdf to track" << endl;
    }
  }

  RooDataHist data_hists("data_hists", "track and shower data", fitutil->GetObs(), cats, hist_map);

  // Fit in both quadrants to find the real minimum of Th23.
  fitutil->SetIOcentvals();
  fitutil->GetVar("SinsqTh23")->setVal(0.4);
  RooFitResult *fitres_1q = simPdf.chi2FitTo( data_hists, Save(), Range("firstq"), DataError(RooAbsData::Poisson) );
  RooArgSet result_1q ( fitres_1q->floatParsFinal() );

  fitutil->SetIOcentvals();
  fitutil->GetVar("SinsqTh23")->setVal(0.6);
  RooFitResult *fitres_2q = simPdf.chi2FitTo( data_hists, Save(), Range("secondq"), DataError(RooAbsData::Poisson) );
  RooArgSet result_2q ( fitres_2q->floatParsFinal() );

  RooArgSet *result;
  Double_t fitChi2_1q = fitres_1q->minNll();
  Double_t fitChi2_2q = fitres_2q->minNll();
  cout << "first q " << TMath::Sqrt( fitChi2_1q ) << endl;
  cout << "second q" << TMath::Sqrt( fitChi2_2q ) << endl;
  if (fitChi2_1q == fitChi2_2q) cout << "NOTICE: Minimizer found same minimum for both quadrants." << endl;
  if (fitChi2_1q < fitChi2_2q) result = &result_1q;
  else                         result = &result_2q;

  cout << "NOTICE Fitter finished fitting, time duration [s]: " << (Double_t)timer.RealTime() << endl;

  cout << "*********Fit result comparison****************************" << endl;
  cout << "dm31       fitted: " << ((RooRealVar*)result->find("Dm31"))->getVal() << endl;
  cout << "sinsq_th23 fitted: " << ((RooRealVar*)result->find("SinsqTh23"))->getVal() << endl;
  cout << "*********Fit result comparison****************************" << endl;

  //----------------------------------------------------------
  // set hierarchy to fitted values
  //----------------------------------------------------------

  Double_t dm31      = ((RooRealVar*)result->find("Dm31"))->getVal();
  Double_t sinSqTh23 = ((RooRealVar*)result->find("SinsqTh23"))->getVal();
  fitutil->GetVar("Dm31")->setVal( dm31 );
  fitutil->GetVar("SinsqTh23")->setVal( sinSqTh23 );

  std::vector<TH3D*> fitted;
  for (Int_t i = 0; i < N_PID_CLASSES; i++) {
    if (pid_map[i] < PID_CUT) {
      TH3D *showers_fitted = (TH3D*)pdf_showers_vector[i].GetExpValHist();
      showers_fitted->SetName( Form("showers_fitted_no_%.2f", pid_map[i]) );
      fitted.push_back( showers_fitted );
    }
    else {
      TH3D *tracks_fitted = (TH3D*)pdf_tracks_vector[i].GetExpValHist();
      tracks_fitted->SetName( Form("tracks_fitted_no_%.2f", pid_map[i]) );
      fitted.push_back( tracks_fitted );
    }
  }

  std::vector< std::tuple<TH1*, Double_t, Double_t> > chi2;
  for (Int_t i = 0; i < N_PID_CLASSES; i++) {
    if (pid_map[i] < PID_CUT) chi2.push_back( NMHUtils::Asymmetry( shower_vector_true[i], fitted[i], Form("sensitivity_shower_%i", i),
                                                 fitEMin, fitEMax, fitctMin, fitctMax) );
    else                      chi2.push_back( NMHUtils::Asymmetry( track_vector_true[i],  fitted[i], Form("sensitivity_track_%i", i),
                                                 fitEMin, fitEMax, fitctMin, fitctMax) );
  }

  Double_t chi2_tot = 0;
  for (Int_t i = 0; i < N_PID_CLASSES; i++) {
    Double_t chi2_i = std::get<1>(chi2[i]);

    cout << "NMHUtils: Chi2 between events NO and events fitted on IO is: " << chi2_i << endl;
    chi2_tot += chi2_i * chi2_i;
  }
  cout << "Squared sum is : " << std::sqrt( chi2_tot ) << endl;

}
